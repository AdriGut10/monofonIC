// This file is part of monofonIC (MUSIC2)
// A software package to generate ICs for cosmological simulations
// Copyright (C) 2020 by Oliver Hahn
//
// monofonIC is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// monofonIC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <transfer_function_plugin.hh>
#include <math/interpolate.hh>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

class transfer_CLASS_file_plugin : public TransferFunction_plugin
{
private:
  enum class dark_matter_component
  {
    cdm,
    ncdm0
  };

  using column_map = std::map<std::string, std::size_t>;

  struct class_native_transfer_table
  {
    column_map columns;
    std::vector<std::vector<double>> rows;
    double redshift = std::numeric_limits<double>::quiet_NaN();
    std::size_t declared_rows = 0;
    bool has_initial_curvature_normalisation = false;
    bool has_hmpc_units = false;
  };

  using signed_interpolator = interpolated_function_1d<true, false, false>;

  signed_interpolator delta_b_, delta_dm_, delta_m_;
  signed_interpolator theta_b_, theta_dm_, theta_m_;

  dark_matter_component dm_component_ = dark_matter_component::cdm;
  cosmology::total_type_t total_type_ = cosmology::MATTER_;
  double f_b_ = 0.0;
  double f_dm_ = 0.0;
  double normalisation_ = 1.0;
  double velocity_normalisation_ = std::numeric_limits<double>::quiet_NaN();
  double kmin_ = 0.0;
  double kmax_ = 0.0;
  bool has_velocities_ = false;

  static std::string lower_case(std::string value)
  {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
  }

  static std::string trim(const std::string &value)
  {
    const auto begin = value.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos)
      return {};
    const auto end = value.find_last_not_of(" \t\r\n");
    return value.substr(begin, end - begin + 1);
  }

  dark_matter_component parse_dark_matter_component() const
  {
    const auto value = lower_case(
        pcf_->get_value_safe<std::string>("cosmology", "DarkMatterComponent", "cdm"));

    if (value == "cdm")
      return dark_matter_component::cdm;
    if (value == "ncdm[0]")
      return dark_matter_component::ncdm0;

    throw std::runtime_error(
        "Invalid cosmology/DarkMatterComponent '" + value +
        "' for CLASS_file; supported values are 'cdm' and 'ncdm[0]'");
  }

  std::string dm_density_column() const
  {
    return dm_component_ == dark_matter_component::cdm ? "d_cdm" : "d_ncdm[0]";
  }

  std::string dm_velocity_column() const
  {
    return dm_component_ == dark_matter_component::cdm ? "t_cdm" : "t_ncdm[0]";
  }

  static column_map parse_numbered_header(const std::string &line)
  {
    column_map columns;
    std::stringstream stream(line);
    std::string token;

    while (stream >> token)
    {
      const auto colon = token.find(':');
      if (colon == std::string::npos || colon == 0 || colon + 1 == token.size())
        continue;

      const auto number_text = token.substr(0, colon);
      if (!std::all_of(number_text.begin(), number_text.end(),
                       [](unsigned char c) { return std::isdigit(c) != 0; }))
        continue;

      std::size_t parsed = 0;
      const auto one_based = std::stoul(number_text, &parsed);
      if (parsed != number_text.size() || one_based == 0)
        throw std::runtime_error("Invalid numbered column in CLASS transfer header: '" + token + "'");

      const auto name = token.substr(colon + 1);
      if (!columns.emplace(name, one_based - 1).second)
        throw std::runtime_error("Duplicate column name '" + name + "' in CLASS transfer header");
    }

    return columns;
  }

  static bool parse_header_redshift(const std::string &line, double &redshift)
  {
    const auto lower = lower_case(line);
    constexpr const char marker[] = "at redshift z=";
    const auto position = lower.find(marker);
    if (position == std::string::npos)
      return false;

    const auto value = trim(line.substr(position + sizeof(marker) - 1));
    std::size_t parsed = 0;
    redshift = std::stod(value, &parsed);
    if (parsed == 0 || !std::isfinite(redshift))
      throw std::runtime_error("Invalid redshift in CLASS transfer header");
    return true;
  }

  static bool parse_declared_rows(const std::string &line, std::size_t &count)
  {
    const auto lower = lower_case(line);
    constexpr const char marker[] = "number of wavenumbers equal to";
    const auto position = lower.find(marker);
    if (position == std::string::npos)
      return false;

    const auto value = trim(line.substr(position + sizeof(marker) - 1));
    std::size_t parsed = 0;
    count = std::stoul(value, &parsed);
    if (parsed == 0)
      throw std::runtime_error("Invalid row count in CLASS transfer header");
    return true;
  }

  static class_native_transfer_table read_native_class_table(const std::string &filename)
  {
    std::ifstream input(filename);
    if (!input)
      throw std::runtime_error("Could not open CLASS transfer function file '" + filename + "'");

    class_native_transfer_table table;
    std::string line;
    std::size_t line_number = 0;

    while (std::getline(input, line))
    {
      ++line_number;
      const auto stripped = trim(line);
      if (stripped.empty())
        continue;

      if (stripped.front() == '#')
      {
        const auto lower = lower_case(stripped);
        table.has_initial_curvature_normalisation =
            table.has_initial_curvature_normalisation ||
            lower.find("normalized to initial curvature=1") != std::string::npos;

        double redshift = 0.0;
        if (parse_header_redshift(stripped, redshift))
          table.redshift = redshift;

        std::size_t declared_rows = 0;
        if (parse_declared_rows(stripped, declared_rows))
          table.declared_rows = declared_rows;

        if (lower.find("1:k (h/mpc)") != std::string::npos)
        {
          table.has_hmpc_units = true;
          table.columns = parse_numbered_header(stripped);
        }
        continue;
      }

      std::stringstream stream(stripped);
      std::vector<double> row;
      double value = 0.0;
      while (stream >> value)
      {
        if (!std::isfinite(value))
          throw std::runtime_error("Non-finite value in CLASS transfer table at line " +
                                   std::to_string(line_number));
        row.push_back(value);
      }
      stream.clear();
      stream >> std::ws;
      if (!stream.eof())
        throw std::runtime_error("Non-numeric value in CLASS transfer table at line " +
                                 std::to_string(line_number));
      if (!row.empty())
        table.rows.push_back(std::move(row));
    }

    if (input.bad())
      throw std::runtime_error("I/O error while reading CLASS transfer function file '" +
                               filename + "'");

    return table;
  }

  static bool has_column(const class_native_transfer_table &table, const std::string &name)
  {
    return table.columns.find(name) != table.columns.end();
  }

  static std::size_t column(const class_native_transfer_table &table, const std::string &name)
  {
    const auto found = table.columns.find(name);
    if (found == table.columns.end())
      throw std::runtime_error("CLASS transfer table is missing required column '" + name + "'");
    return found->second;
  }

  void validate_required_columns(const class_native_transfer_table &table) const
  {
    if (!table.has_initial_curvature_normalisation)
      throw std::runtime_error(
          "CLASS_file requires a native CLASS table normalized to initial curvature=1");
    if (!table.has_hmpc_units)
      throw std::runtime_error("CLASS_file requires a native CLASS table with k in h/Mpc");
    if (table.columns.empty() || !has_column(table, "k"))
      throw std::runtime_error("Could not find the numbered native CLASS column header");

    const std::size_t column_count = table.columns.size();
    std::vector<bool> numbered(column_count, false);
    for (const auto &entry : table.columns)
    {
      if (entry.second >= column_count || numbered[entry.second])
        throw std::runtime_error("CLASS transfer header columns must be numbered consecutively from 1");
      numbered[entry.second] = true;
    }
    if (!std::all_of(numbered.begin(), numbered.end(), [](bool value) { return value; }))
      throw std::runtime_error("CLASS transfer header columns must be numbered consecutively from 1");

    column(table, "d_b");
    column(table, dm_density_column());

    if (total_type_ == cosmology::TOTAL_)
      column(table, "d_tot");

    const bool has_tb = has_column(table, "t_b");
    const bool has_tdm = has_column(table, dm_velocity_column());
    const bool has_velocity_pair = has_tb && has_tdm;
    if (has_velocity_pair && total_type_ == cosmology::TOTAL_ && !has_column(table, "t_tot"))
      throw std::runtime_error(
          "CLASS transfer table requires 't_tot' for TransferComponent=total when velocities are present");

    const bool wants_relative_velocity =
        pcf_->get_value_safe<bool>("setup", "DoBaryonVrel", false);
    if (wants_relative_velocity && !has_velocity_pair)
      throw std::runtime_error(
          "DoBaryonVrel=yes requires 't_b' and '" + dm_velocity_column() +
          "' in the CLASS transfer table (generate it with output=dTk,vTk)");
  }

  void initialise_interpolators(const class_native_transfer_table &table)
  {
    if (table.rows.size() < 6)
      throw std::runtime_error("CLASS transfer table must contain at least six data rows");
    if (table.declared_rows != 0 && table.declared_rows != table.rows.size())
      throw std::runtime_error("CLASS transfer table row count does not match its header");

    const auto ik = column(table, "k");
    const auto idb = column(table, "d_b");
    const auto iddm = column(table, dm_density_column());
    const bool has_native_matter = has_column(table, "d_m");
    const auto idaggregate = total_type_ == cosmology::TOTAL_
                                 ? column(table, "d_tot")
                                 : (has_native_matter ? column(table, "d_m") : 0);

    has_velocities_ = has_column(table, "t_b") &&
                      has_column(table, dm_velocity_column());
    const auto itb = has_velocities_ ? column(table, "t_b") : 0;
    const auto itdm = has_velocities_ ? column(table, dm_velocity_column()) : 0;
    const auto ittot = has_velocities_ && total_type_ == cosmology::TOTAL_
                           ? column(table, "t_tot")
                           : 0;

    const auto expected_columns = table.columns.size();
    std::vector<double> k, db, ddm, dm, tb, tdm, tm;
    k.reserve(table.rows.size());
    db.reserve(table.rows.size());
    ddm.reserve(table.rows.size());
    dm.reserve(table.rows.size());
    if (has_velocities_)
    {
      tb.reserve(table.rows.size());
      tdm.reserve(table.rows.size());
      tm.reserve(table.rows.size());
    }

    for (const auto &row : table.rows)
    {
      if (row.size() != expected_columns)
        throw std::runtime_error("CLASS transfer table row width does not match its numbered header");

      const double kh = row[ik];
      if (!std::isfinite(kh) || kh <= 0.0 || (!k.empty() && kh <= k.back()))
        throw std::runtime_error("CLASS transfer table k values must be finite, positive, and increasing");

      const double density_b = row[idb];
      const double density_dm = row[iddm];
      const double inverse_k2_normalisation = -normalisation_ / (kh * kh);

      double density_aggregate = 0.0;
      switch (total_type_)
      {
      case cosmology::TOTAL_:
        density_aggregate = row[idaggregate];
        break;
      case cosmology::MATTER_:
        density_aggregate = has_native_matter ? row[idaggregate]
                                              : f_b_ * density_b + f_dm_ * density_dm;
        break;
      case cosmology::BPLUSC_:
        density_aggregate = f_b_ * density_b + f_dm_ * density_dm;
        break;
      default:
        throw std::runtime_error("Invalid TransferComponent in CLASS_file");
      }

      k.push_back(kh);
      db.push_back(density_b * inverse_k2_normalisation);
      ddm.push_back(density_dm * inverse_k2_normalisation);
      dm.push_back(density_aggregate * inverse_k2_normalisation);

      if (has_velocities_)
      {
        const double velocity_b = row[itb];
        const double velocity_dm = row[itdm];
        const double velocity_aggregate = total_type_ == cosmology::TOTAL_
                                              ? row[ittot]
                                              : f_b_ * velocity_b + f_dm_ * velocity_dm;
        tb.push_back(velocity_b * inverse_k2_normalisation);
        tdm.push_back(velocity_dm * inverse_k2_normalisation);
        tm.push_back(velocity_aggregate * inverse_k2_normalisation);
      }
    }

    delta_b_.set_data(k, db);
    delta_dm_.set_data(k, ddm);
    delta_m_.set_data(k, dm);
    if (has_velocities_)
    {
      theta_b_.set_data(k, tb);
      theta_dm_.set_data(k, tdm);
      theta_m_.set_data(k, tm);
    }

    // Avoid endpoints where cubic-spline interpolation becomes lower order.
    kmin_ = k[1];
    kmax_ = k[k.size() - 2];
  }

  double density_value(double k, tf_type type) const
  {
    switch (type)
    {
    case delta_matter:
    case delta_matter0:
      return delta_m_(k);
    case delta_cdm:
    case delta_cdm0:
      return delta_dm_(k);
    case delta_baryon:
    case delta_baryon0:
      return delta_b_(k);
    case delta_bc:
      return delta_b_(k) - delta_dm_(k);
    default:
      throw std::runtime_error("Invalid density transfer type requested from CLASS_file");
    }
  }

  double velocity_value(double k, tf_type type) const
  {
    if (!has_velocities_)
      return 0.0;
    if (!std::isfinite(velocity_normalisation_) || velocity_normalisation_ <= 0.0)
      throw std::runtime_error("CLASS_file velocity normalisation was not initialized");

    double value = 0.0;
    switch (type)
    {
    case theta_matter:
    case theta_matter0:
      value = theta_m_(k);
      break;
    case theta_cdm:
    case theta_cdm0:
      value = theta_dm_(k);
      break;
    case theta_baryon:
    case theta_baryon0:
      value = theta_b_(k);
      break;
    case theta_bc:
      value = theta_b_(k) - theta_dm_(k);
      break;
    default:
      throw std::runtime_error("Invalid velocity transfer type requested from CLASS_file");
    }
    return value / velocity_normalisation_;
  }

public:
  explicit transfer_CLASS_file_plugin(config_file &cf, const cosmology::parameters &cosmo_params)
      : TransferFunction_plugin(cf, cosmo_params)
  {
    dm_component_ = parse_dark_matter_component();
    total_type_ = cosmo_params_.get_total_type();
    f_b_ = cosmo_params_["f_b"];
    f_dm_ = cosmo_params_["f_c"];

    const double amplitude = cosmo_params_["A_s"];
    if (amplitude > 0.0)
    {
      const double pivot_hmpc = cosmo_params_["k_p"] / cosmo_params_["h"];
      if (!std::isfinite(pivot_hmpc) || pivot_hmpc <= 0.0)
        throw std::runtime_error("CLASS_file requires a positive cosmology/k_p");

      const double pi = std::acos(-1.0);
      normalisation_ = std::sqrt(2.0 * pi * pi * amplitude *
                                 std::pow(1.0 / pivot_hmpc, cosmo_params_["n_s"] - 1.0) /
                                 std::pow(2.0 * pi, 3.0));
      if (!std::isfinite(normalisation_) || normalisation_ <= 0.0)
        throw std::runtime_error("Could not compute CLASS_file primordial normalization from A_s");
      tf_isnormalised_ = true;
      music::ilog << "CLASS_file: using A_s=" << colors::CONFIG_VALUE << amplitude
                  << colors::RESET << " for transfer normalization." << std::endl;
    }
    else
    {
      const double sigma8 = cosmo_params_["sigma_8"];
      if (!std::isfinite(sigma8) || sigma8 <= 0.0)
        throw std::runtime_error("CLASS_file requires either a positive A_s or sigma_8");
      normalisation_ = 1.0;
      tf_isnormalised_ = false;
      music::ilog << "CLASS_file: using sigma_8=" << colors::CONFIG_VALUE << sigma8
                  << colors::RESET << " for numerical transfer normalization." << std::endl;
    }

    const auto filename = pcf_->get_value<std::string>("cosmology", "transfer_file");
    music::ilog << "Reading native CLASS transfer function data from:" << std::endl
                << "  '" << filename << "'" << std::endl;

    const auto table = read_native_class_table(filename);
    validate_required_columns(table);

    const double target_redshift =
        pcf_->get_value_safe<double>("cosmology", "ztarget", 0.0);
    if (!std::isfinite(table.redshift))
      throw std::runtime_error("CLASS transfer table header does not specify its redshift");
    const double redshift_tolerance =
        1.0e-8 * std::max(1.0, std::max(std::abs(table.redshift), std::abs(target_redshift)));
    if (std::abs(table.redshift - target_redshift) > redshift_tolerance)
      throw std::runtime_error("CLASS transfer table redshift z=" +
                               std::to_string(table.redshift) +
                               " does not match cosmology/ztarget=" +
                               std::to_string(target_redshift));

    initialise_interpolators(table);

    tf_distinct_ = true;
    tf_withvel_ = has_velocities_;
    tf_withtotal0_ = false;
    tf_velunits_ = false;

    const char *component =
        dm_component_ == dark_matter_component::cdm ? "cdm" : "ncdm[0]";
    music::ilog << "Read native CLASS table with " << table.rows.size()
                << " rows; logical dark matter is " << colors::CONFIG_VALUE << component
                << colors::RESET << ", k=" << colors::CONFIG_VALUE << kmin_ << colors::RESET
                << " to " << colors::CONFIG_VALUE << kmax_ << colors::RESET << " h/Mpc."
                << std::endl;
  }

  void set_velocity_normalisation(double fHa_target_Mpc_inv) override
  {
    if (!std::isfinite(fHa_target_Mpc_inv) || fHa_target_Mpc_inv <= 0.0)
      throw std::runtime_error("CLASS_file received an invalid f*a*H velocity normalization");
    velocity_normalisation_ = fHa_target_Mpc_inv;
    music::ilog.Print("CLASS_file: f*a*H/c = %.17g 1/Mpc", velocity_normalisation_);
  }

  double compute(double k, tf_type type) const override
  {
    if (!std::isfinite(k) || k <= 0.0 || k > kmax_)
      return 0.0;
    k = std::max(k, kmin_);

    switch (type)
    {
    case delta_matter:
    case delta_cdm:
    case delta_baryon:
    case delta_bc:
    case delta_matter0:
    case delta_cdm0:
    case delta_baryon0:
      return density_value(k, type);

    case theta_matter:
    case theta_cdm:
    case theta_baryon:
    case theta_bc:
    case theta_matter0:
    case theta_cdm0:
    case theta_baryon0:
      return velocity_value(k, type);

    default:
      throw std::runtime_error("Invalid transfer type requested from CLASS_file");
    }
  }

  double get_kmin(void) const override { return kmin_; }
  double get_kmax(void) const override { return kmax_; }
};

namespace
{
TransferFunction_plugin_creator_concrete<transfer_CLASS_file_plugin> creator("CLASS_file");
}
