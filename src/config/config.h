 
#ifndef _config_h
#define _config_h

#include <string>
#include <vector>

namespace agb {

struct ModelConfig {
  ModelConfig(const std::string model_folder);
  bool loadConfigFile(const std::string folder);   //reads <folder>/config.toml

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

  //Make the frequency-INTEGRATED flux eddington_flux_int conservative via the flux-
  //divergence form (thesis eq. 2.59): r^2 H_int = L/(16 pi^2) + int r^2 [int kappa_abs(B-J)
  //dnu] dr, computed in conservativeFluxIntegral() and overwriting the (non-conservative,
  //eq. 2.58) integral. This is consistent with the moment equation (2.61) the RT solves, so
  //r^2 H_int is conserved at radiative equilibrium (lets Unsoeld-Lucy's flux term and the
  //flux-convergence check work). The per-frequency eddington_flux stays on the eq. 2.58
  //transport form (calcFlux) so the monochromatic flux is physical (>=0) for the spectrum
  //and the flux-mean extinction. Set false for the pure legacy eq. 2.58 integral.
  //config.toml: [radiative_transfer] flux_from_divergence.
  bool flux_from_divergence = true;
  
  unsigned int nb_temperature_iter = 200;
  double temperature_convergence = 1e-2;
  double temperature_max_change = 0.005;
  bool smooth_temperature_profile = true;

  //Full-linearisation temperature corrector (thesis sec. 3.2.3 / App. B.1) instead
  //of the approximate Unsoeld-Lucy correction. Treats the converged RT operator
  //(Eddington/sphericality factors, opacities, geometry) as frozen and performs one
  //Newton step on (T_gas, T_dust) so the radiation field's response dJ/dT is the
  //EXACT moment-system response M_n^{-1}, not the tau-weighted Unsoeld-Lucy heuristic.
  //The per-frequency mean intensities are eliminated by a Rybicki-type reduction to a
  //dense 2D x 2D temperature system (D = nb_grid_points). Requires the Taylor moment
  //discretisation (use_spline_discretisation = false). Now a config.toml key.
  //Default OFF -> exact current Unsoeld-Lucy behaviour.
  bool   use_linearisation = true;
  double linearisation_relaxation = 0.5;   //global under-relaxation of the Newton step
  //Two-phase corrector: start with Unsoeld-Lucy and switch to the (Newton) linearisation
  //only once the profile is close enough that the frozen-opacity linearisation is valid.
  //Far from convergence the Newton step overshoots (the dust nucleation rate ~ exp(-/T) is
  //extremely temperature-sensitive, so a frozen-opacity step badly mispredicts the dust
  //response); Unsoeld-Lucy's damped fixed-point iteration is robust there. Switch once, for
  //linearisation_switch_count consecutive steps, BOTH the per-step max|dT|/T has settled
  //below linearisation_switch_dt_fraction*temperature_max_change AND the max radiative-
  //equilibrium residual (energy-balance deviation) is below linearisation_switch_re_residual.
  //Then latch (one-way). Set linearisation_start_unsoeld_lucy=false to linearise immediately
  //(e.g. a frozen-structure test starting near convergence).
  bool   linearisation_start_unsoeld_lucy = true;
  double linearisation_switch_dt_fraction = 0.5;
  double linearisation_switch_re_residual = 3e-2;
  unsigned int linearisation_switch_count = 3;
  //Composite radiative-equilibrium constraint (thesis eq. 3.64/3.65): the local
  //radiative-equilibrium ratio (weight xi, reliable in the optically thin region) plus
  //the flux-constancy ratio  [int d(f q r^2 J)/dX dnu]/(r*^2 H*) - 1  (weight zeta(r),
  //needed in the optically thick / diffusion region where local RE degenerates and the
  //temperature is otherwise left unpinned even though the flux is not conserved). With
  //the flux term off, this reduces to the plain full linearisation (local RE only).
  //Default OFF: enforcing flux constancy via the temperature on a fixed structure
  //diverges here, because the discrete local RE (from the moment-system J) and the
  //discrete flux constancy (from the separately discretised calcFlux H) are mutually
  //inconsistent by ~2.3% on the converged structure - no temperature satisfies both.
  //That residual is an RT J-vs-H discretisation floor, not a temperature-correction
  //error. Leave the flux term off; the plain (local-RE) full linearisation converges.
  bool   linearisation_flux_constraint = true;
  double linearisation_xi  = 1.0;          //weight of the local-RE term (constant, ~1)
  //zeta(r) = tau / (tau + tau_scale): ~1 deep (flux constancy), ->0 at the thin outer
  //edge (where the flux derivative is numerically noisy). tau is a grey radial optical
  //depth from the outer boundary, built from the flux-mean extinction.
  double linearisation_zeta_tau_scale = 1.0;

  //Hard-enforce a strictly outward-decreasing temperature profile after each
  //correction (clip T[i] -> 0.99 T[i-1]). Now a config.toml key. Default OFF:
  //this constraint fights the Unsoeld-Lucy correction wherever radiative
  //equilibrium wants a flatter profile or a (e.g. dust-front) inversion, clipping
  //the just-applied heating back down every step and driving a +/-cap limit cycle
  //that never settles. Enable only if a run genuinely needs it.
  bool force_monotonic_temperature = false;

  //Robustness controls for the radiative-equilibrium temperature iteration.
  //These are now config.toml keys; adjust here if needed. The iteration
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

  //Movable (adaptive) radial grid via equidistribution of a monitor function.
  //Now a config.toml key; adjust here. The grid moves as a step decoupled from
  //the solvers: every regrid_frequency outer iterations the nodes are redistributed
  //so that  w(r) = 1 + sum a_k |d ln q_k/d ln r|  (over flux-mean opacity, gas
  //temperature, nucleation rate) is equal-integral per cell, then under-relaxed,
  //state is interpolated, and the RT ray geometry is rebuilt. Default OFF.
  bool   use_movable_grid = false;
  unsigned int regrid_frequency = 5;          //regrid every N outer iterations
  double grid_relaxation = 0.1;               //omega_grid: node-motion under-relaxation
  double monitor_weight_opacity = 1.0;        //a_chi  (flux-mean extinction)
  double monitor_weight_temperature = 0.0;    //a_T    (gas temperature)
  double monitor_weight_nucleation = 1.0;     //a_J    (dust nucleation rate)
  double monitor_weight_velocity = 1.0;       //a_v    (wind velocity -> clusters at sonic point)
  int    monitor_smoothing_passes = 4;
  double monitor_max = 20.0;                   //cap on the monitor (max cell-size ratio)
  double monitor_rel_floor = 1e-8;            //per-quantity floor relative to its peak

  //Structure solver for the stationary wind. Default: the Henyey-type global
  //Newton-Raphson eigenvalue solve (Setup A: dust decoupled / alpha frozen in the
  //Newton, alpha damped in the outer loop). Set false to fall back to the legacy
  //Melia-Phi shooting method (kept for cross-checks). Now a config.toml key;
  //flip here if needed.
  bool use_henyey_solver = false;

  std::string output_spectrum_path = "";
  std::string output_atmosphere_path = "";
  std::string output_dust_path = "";
  std::string output_hydro_path = "";
};



}


#endif