 
#ifndef _config_h
#define _config_h

#include <string>
#include <vector>

namespace agb {

struct ModelConfig {
  ModelConfig(const std::string model_folder);
  bool loadConfigFile(const std::string folder);
  bool loadOutputConfigFile(const std::string folder);
  void readOpacityConfig(std::fstream& file);
  
  std::string model_folder = "";

  double stellar_radius = 0;
  double stellar_mass = 0;
  double stellar_luminosity = 0;
  double stellar_mass_loss_rate = 0;
  double c_o_ratio = 0;

  std::string fastchem_parameter_file = "";

  std::string starting_model_path = "";
  
  std::string opacity_path = "";
  std::string wavenumber_file_path = "";
  double spectral_resolution = 0;
  double min_wavelength = 0;
  double max_wavelength = 0;

  size_t nb_core_impact_param = 20;

  std::vector<std::string> opacity_species_symbol;
  std::vector<std::string> opacity_species_folder;

  std::string refractive_index_file = "";

  unsigned int nb_radiative_transfer_iter = 30;
  double radiative_transfer_convergence = 1e-4;
  bool use_spline_discretisation = true;
  
  unsigned int nb_temperature_iter = 200;
  double temperature_convergence = 1e-2;
  double temperature_max_change = 0.005;
  bool smooth_temperature_profile = true;

  //Hard-enforce a strictly outward-decreasing temperature profile after each
  //correction (clip T[i] -> 0.99 T[i-1]). Not read from model.config. Default OFF:
  //this constraint fights the Unsoeld-Lucy correction wherever radiative
  //equilibrium wants a flatter profile or a (e.g. dust-front) inversion, clipping
  //the just-applied heating back down every step and driving a +/-cap limit cycle
  //that never settles. Enable only if a run genuinely needs it.
  bool force_monotonic_temperature = false;

  //Robustness controls for the radiative-equilibrium temperature iteration.
  //These are not read from model.config; adjust here if needed. The iteration
  //uses an exact T^4 (Planck) update with an adaptive per-layer under-relaxation
  //omega (halved on sign flips to kill oscillations, grown back when the
  //correction is stable), optionally accelerated by Anderson mixing. The per-step
  //size is still capped by temperature_max_change so the change stays small enough
  //for the interleaved hydrodynamics-dust cycle to re-converge.
  double temperature_relaxation_init = 1.0;     //initial per-layer omega
  double temperature_relaxation_min  = 0.05;    //floor for omega when oscillating
  double temperature_relaxation_max  = 1.0;     //ceiling for omega
  double temperature_relaxation_down = 0.5;     //factor on omega when delta flips sign
  double temperature_relaxation_up   = 1.3;     //factor on omega when delta keeps its sign
  bool   use_anderson_acceleration = true;      //Anderson mixing on the T profile
  unsigned int anderson_window = 4;             //history depth m
  unsigned int anderson_start_iter = 2;         //plain (damped) iterations before accelerating
  //Only accelerate once the iteration has entered the settling regime, i.e. when
  //the max relative T change has dropped below this fraction of the per-step cap.
  //This keeps the early large-change transient (which the hydro cycle is sensitive
  //to) on plain damped steps and applies Anderson where the slow creep / r^2 H
  //oscillation actually occurs.
  double anderson_activation_fraction = 0.5;

  unsigned int nb_hydrodynamics_iter = 200;
  double hydrodynamics_convergence = 1e-2;

  //Structure solver for the stationary wind. Default: the Henyey-type global
  //Newton-Raphson eigenvalue solve (Setup A: dust decoupled / alpha frozen in the
  //Newton, alpha damped in the outer loop). Set false to fall back to the legacy
  //Melia-Phi shooting method (kept for cross-checks). Not read from model.config;
  //flip here if needed.
  bool use_henyey_solver = false;

  std::string output_spectrum_path = "";
  std::string output_atmosphere_path = "";
  std::string output_dust_path = "";
  std::string output_hydro_path = "";
};



}


#endif