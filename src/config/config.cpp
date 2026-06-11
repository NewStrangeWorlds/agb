
#include "config.h"

#include <iostream>
#include <string>
#include <cstdint>

#include <toml++/toml.hpp>

#include "../additional/physical_const.h"

namespace agb {


namespace {

//Lenient getters so an integer literal also satisfies a float key and vice versa
//(TOML / toml++ are type-strict; value_or would silently return the default on a
//type mismatch, e.g. resolution = 1000 for a double field). Missing keys/tables
//also fall through to the default.
double getDouble(toml::node_view<toml::node> n, double def)
{
  if (auto v = n.value<double>())  return *v;
  if (auto v = n.value<int64_t>()) return static_cast<double>(*v);
  return def;
}

int64_t getInt(toml::node_view<toml::node> n, int64_t def)
{
  if (auto v = n.value<int64_t>()) return *v;
  if (auto v = n.value<double>())  return static_cast<int64_t>(*v);
  return def;
}

bool getBool(toml::node_view<toml::node> n, bool def)
{
  return n.value<bool>().value_or(def);
}

std::string getStr(toml::node_view<toml::node> n, const std::string& def)
{
  return n.value<std::string>().value_or(def);
}

}



ModelConfig::ModelConfig(const std::string model_folder)
{
  loadConfigFile(model_folder);
}


bool ModelConfig::loadConfigFile(const std::string folder_path)
{
  model_folder = folder_path;
  if (model_folder.back() != '/')
    model_folder.append("/");

  const std::string file_path = model_folder + "config.toml";

  toml::table t;
  try
  {
    t = toml::parse_file(file_path);
  }
  catch (const toml::parse_error& err)
  {
    std::cout << "Couldn't parse config file " << file_path << ":\n" << err << "\n";
    return false;
  }

  std::cout << "\nParameters found in " << file_path << ":\n";

  //--- [star] (required; given in solar units, stored in cgs) ---
  stellar_radius         = getDouble(t["star"]["radius"], 0.0);
  stellar_mass           = getDouble(t["star"]["mass"], 0.0);
  stellar_luminosity     = getDouble(t["star"]["luminosity"], 0.0);
  stellar_mass_loss_rate = getDouble(t["star"]["mass_loss_rate"], 0.0);
  c_o_ratio              = getDouble(t["star"]["c_o_ratio"], 0.0);

  std::cout << "- star: R*=" << stellar_radius << " R_sun, M*=" << stellar_mass
            << " M_sun, L*=" << stellar_luminosity << " L_sun, Mdot="
            << stellar_mass_loss_rate << " M_sun/yr, C/O=" << c_o_ratio << "\n";

  if (stellar_radius <= 0. || stellar_mass <= 0. || stellar_luminosity <= 0.)
    std::cout << "  WARNING: a [star] parameter is missing or non-positive.\n";

  stellar_radius         *= constants::radius_sun;
  stellar_mass           *= constants::mass_sun;
  stellar_luminosity     *= constants::luminosity_sun;
  stellar_mass_loss_rate *= constants::mass_sun / constants::year;

  //--- [model] ---
  starting_model_path     = getStr(t["model"]["starting_model"], "");
  fastchem_parameter_file = getStr(t["model"]["fastchem_parameter_file"], "");
  std::cout << "- starting model: " << starting_model_path
            << "  | FastChem: " << fastchem_parameter_file << "\n";

  //--- [spectral_grid] ---
  min_wavelength      = getDouble(t["spectral_grid"]["min_wavelength"], min_wavelength);
  max_wavelength      = getDouble(t["spectral_grid"]["max_wavelength"], max_wavelength);
  spectral_resolution = getDouble(t["spectral_grid"]["resolution"], spectral_resolution);
  nb_core_impact_param = static_cast<size_t>(
    getInt(t["spectral_grid"]["nb_core_impact_param"],
           static_cast<int64_t>(nb_core_impact_param)));
  std::cout << "- spectral grid: " << min_wavelength << "-" << max_wavelength
            << " micron, R=" << spectral_resolution << "\n";

  //--- [radiative_transfer] ---
  nb_radiative_transfer_iter     = static_cast<unsigned int>(
    getInt(t["radiative_transfer"]["nb_iterations"], nb_radiative_transfer_iter));
  radiative_transfer_convergence = getDouble(t["radiative_transfer"]["convergence"], radiative_transfer_convergence);
  use_spline_discretisation      = getBool(t["radiative_transfer"]["use_spline_discretisation"], use_spline_discretisation);
  flux_from_divergence           = getBool(t["radiative_transfer"]["flux_from_divergence"], flux_from_divergence);
  std::cout << "- RT: " << nb_radiative_transfer_iter << " iter, conv "
            << radiative_transfer_convergence << ", spline=" << use_spline_discretisation
            << ", flux_from_divergence=" << flux_from_divergence << "\n";

  //--- [temperature] (incl. the corrector knobs promoted from code) ---
  nb_temperature_iter        = static_cast<unsigned int>(
    getInt(t["temperature"]["nb_iterations"], nb_temperature_iter));
  temperature_convergence    = getDouble(t["temperature"]["convergence"], temperature_convergence);
  temperature_max_change     = getDouble(t["temperature"]["max_relative_change"], temperature_max_change);
  smooth_temperature_profile = getBool(t["temperature"]["smooth_profile"], smooth_temperature_profile);
  force_monotonic_temperature= getBool(t["temperature"]["force_monotonic"], force_monotonic_temperature);

  use_linearisation                = getBool(t["temperature"]["use_linearisation"], use_linearisation);
  linearisation_relaxation         = getDouble(t["temperature"]["linearisation_relaxation"], linearisation_relaxation);
  linearisation_start_unsoeld_lucy = getBool(t["temperature"]["linearisation_start_unsoeld_lucy"], linearisation_start_unsoeld_lucy);
  linearisation_switch_dt_fraction = getDouble(t["temperature"]["linearisation_switch_dt_fraction"], linearisation_switch_dt_fraction);
  linearisation_switch_re_residual = getDouble(t["temperature"]["linearisation_switch_re_residual"], linearisation_switch_re_residual);
  linearisation_switch_count       = static_cast<unsigned int>(
    getInt(t["temperature"]["linearisation_switch_count"], linearisation_switch_count));
  linearisation_flux_constraint    = getBool(t["temperature"]["linearisation_flux_constraint"], linearisation_flux_constraint);
  linearisation_xi                 = getDouble(t["temperature"]["linearisation_xi"], linearisation_xi);
  linearisation_zeta_tau_scale     = getDouble(t["temperature"]["linearisation_zeta_tau_scale"], linearisation_zeta_tau_scale);

  temperature_relaxation_init = getDouble(t["temperature"]["relaxation_init"], temperature_relaxation_init);
  temperature_relaxation_min  = getDouble(t["temperature"]["relaxation_min"], temperature_relaxation_min);
  temperature_relaxation_max  = getDouble(t["temperature"]["relaxation_max"], temperature_relaxation_max);
  temperature_relaxation_down = getDouble(t["temperature"]["relaxation_down"], temperature_relaxation_down);
  temperature_relaxation_up   = getDouble(t["temperature"]["relaxation_up"], temperature_relaxation_up);
  use_anderson_acceleration   = getBool(t["temperature"]["use_anderson"], use_anderson_acceleration);
  anderson_window             = static_cast<unsigned int>(
    getInt(t["temperature"]["anderson_window"], anderson_window));
  anderson_start_iter         = static_cast<unsigned int>(
    getInt(t["temperature"]["anderson_start_iter"], anderson_start_iter));
  anderson_activation_fraction= getDouble(t["temperature"]["anderson_activation_fraction"], anderson_activation_fraction);
  std::cout << "- temperature: " << nb_temperature_iter << " iter, conv "
            << temperature_convergence << ", max change " << temperature_max_change
            << ", use_linearisation=" << use_linearisation
            << " (relax " << linearisation_relaxation
            << ", start_UL=" << linearisation_start_unsoeld_lucy << ")\n";

  //--- [hydrodynamics] ---
  nb_hydrodynamics_iter     = static_cast<unsigned int>(
    getInt(t["hydrodynamics"]["nb_iterations"], nb_hydrodynamics_iter));
  hydrodynamics_convergence = getDouble(t["hydrodynamics"]["convergence"], hydrodynamics_convergence);
  use_henyey_solver         = getBool(t["hydrodynamics"]["use_henyey_solver"], use_henyey_solver);
  std::cout << "- hydrodynamics: " << nb_hydrodynamics_iter << " iter, conv "
            << hydrodynamics_convergence << ", henyey=" << use_henyey_solver << "\n";

  //--- [dust] ---
  refractive_index_file = getStr(t["dust"]["refractive_index_file"], "");
  std::cout << "- dust refractive indices: " << refractive_index_file << "\n";

  //--- [movable_grid] (default off) ---
  use_movable_grid           = getBool(t["movable_grid"]["enabled"], use_movable_grid);
  regrid_frequency           = static_cast<unsigned int>(
    getInt(t["movable_grid"]["regrid_frequency"], regrid_frequency));
  grid_relaxation            = getDouble(t["movable_grid"]["grid_relaxation"], grid_relaxation);
  monitor_weight_opacity     = getDouble(t["movable_grid"]["monitor_weight_opacity"], monitor_weight_opacity);
  monitor_weight_temperature = getDouble(t["movable_grid"]["monitor_weight_temperature"], monitor_weight_temperature);
  monitor_weight_nucleation  = getDouble(t["movable_grid"]["monitor_weight_nucleation"], monitor_weight_nucleation);
  monitor_weight_velocity    = getDouble(t["movable_grid"]["monitor_weight_velocity"], monitor_weight_velocity);
  monitor_smoothing_passes   = static_cast<int>(
    getInt(t["movable_grid"]["monitor_smoothing_passes"], monitor_smoothing_passes));
  monitor_max                = getDouble(t["movable_grid"]["monitor_max"], monitor_max);
  monitor_rel_floor          = getDouble(t["movable_grid"]["monitor_rel_floor"], monitor_rel_floor);
  if (use_movable_grid)
    std::cout << "- movable grid: ON (regrid every " << regrid_frequency << ")\n";

  //--- [opacity] ---
  opacity_path = getStr(t["opacity"]["path"], "");
  if (!opacity_path.empty() && opacity_path.back() != '/')
    opacity_path.append("/");
  wavenumber_file_path = opacity_path + "wavenumber_full.dat";

  if (auto species = t["opacity"]["species"].as_array())
  {
    for (auto&& element : *species)
    {
      auto entry = element.as_table();
      if (!entry) continue;

      const std::string symbol = getStr((*entry)["symbol"], "");
      const std::string folder = getStr((*entry)["folder"], "");
      if (!symbol.empty() && !folder.empty())
      {
        opacity_species_symbol.push_back(symbol);
        opacity_species_folder.push_back(folder);
      }
    }
  }
  std::cout << "- opacity: " << opacity_path << " ("
            << opacity_species_symbol.size() << " species)\n";

  //--- [output] (paths relative to the model folder; "none"/absent -> disabled) ---
  auto outputPath = [&](const char* key) -> std::string
  {
    const std::string p = getStr(t["output"][key], "");
    if (p.empty() || p == "none" || p == "None") return "";
    return model_folder + p;
  };
  output_spectrum_path   = outputPath("spectrum");
  output_atmosphere_path = outputPath("atmosphere");
  output_dust_path       = outputPath("dust");
  output_hydro_path      = outputPath("hydro");

  std::cout << "\n";

  return true;
}


}
