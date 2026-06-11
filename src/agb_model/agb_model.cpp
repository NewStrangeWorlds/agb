 
#include "agb_model.h"

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "../additional/aux_functions.h"
#include "../additional/physical_const.h"
#include "../additional/movable_grid.h"


namespace agb{


namespace {

//Dense Gaussian elimination with partial pivoting for the small (m x m) Anderson
//least-squares system. Returns an empty vector if the matrix is singular.
std::vector<double> solveLinearSystem(
  std::vector<std::vector<double>> a,
  std::vector<double> b)
{
  const size_t m = b.size();

  for (size_t col=0; col<m; ++col)
  {
    size_t pivot = col;
    for (size_t row=col+1; row<m; ++row)
      if (std::abs(a[row][col]) > std::abs(a[pivot][col])) pivot = row;

    if (std::abs(a[pivot][col]) < 1e-300)
      return std::vector<double>{};

    std::swap(a[col], a[pivot]);
    std::swap(b[col], b[pivot]);

    for (size_t row=col+1; row<m; ++row)
    {
      const double factor = a[row][col] / a[col][col];
      for (size_t k=col; k<m; ++k) a[row][k] -= factor * a[col][k];
      b[row] -= factor * b[col];
    }
  }

  std::vector<double> x(m, 0.);
  for (size_t row=m; row-- > 0; )
  {
    double sum = b[row];
    for (size_t k=row+1; k<m; ++k) sum -= a[row][k] * x[k];
    x[row] = sum / a[row][row];
  }

  return x;
}

}

AGBStarModel::AGBStarModel(const std::string folder)
 : config(folder)
 , spectral_grid(&config)
 , atmosphere(&config, spectral_grid.nbSpectralPoints())
 , chemistry(config.model_folder+config.fastchem_parameter_file, 1.0, config.c_o_ratio)
 //, dust_species(new AnalyticDust(&config, &spectral_grid, &atmosphere, 0.05))
 , dust_species(new GailSedlmayrDust(&config, &spectral_grid, &atmosphere))
 , transport_coeff(&config, &spectral_grid, config.opacity_species_symbol, config.opacity_species_folder)
 , radiative_transfer(&config, &spectral_grid, &atmosphere)
 , temperature_correction(&config, &spectral_grid, radiative_transfer.radiation_field)
 , hydrodynamics(&config, &spectral_grid, &atmosphere, radiative_transfer.radiation_field)
{

}


void AGBStarModel::calcModel()
{
  std::vector<double> temperature_gas_old = atmosphere.temperature_gas;
  std::vector<double> temperature_dust_old = atmosphere.temperature_dust;

  //reset the persistent temperature-iteration state for this run
  temp_iter_count = 0;
  prev_max_rel_change = 1e30;
  anderson_was_active = false;
  relaxation_gas.clear();    relaxation_dust.clear();
  prev_delta_b_gas.clear();  prev_delta_b_dust.clear();
  anderson_x_gas.clear();    anderson_f_gas.clear();
  anderson_x_dust.clear();   anderson_f_dust.clear();

  std::vector<double> degree_of_condensation(atmosphere.nb_grid_points, 0);

  chemistry.calcChemicalComposition(
    std::vector<double>{}, 
    atmosphere.temperature_gas, 
    atmosphere.pressure_bar,
    degree_of_condensation,
    atmosphere.number_densities, 
    atmosphere.mean_molecuar_weight,
    atmosphere.total_element_density,
    atmosphere.total_h_density);

  atmosphere.equationOfState();

  dust_species->calcDistribution(
      chemistry.element_abundances[chemistry.fastchem_element_indices[_C]]
    - chemistry.element_abundances[chemistry.fastchem_element_indices[_O]]);

  chemistryDustIteration();
  chemistryHydroIteration();

  //two-phase corrector: unless told to linearise immediately, begin with Unsoeld-Lucy and
  //let temperatureIteration() latch on the linearisation once the profile is close enough
  linearisation_active = config.use_linearisation && !config.linearisation_start_unsoeld_lucy;

  for (unsigned int it=0; it<config.nb_temperature_iter; ++it)
  {
    //chemistryDustIteration();
    chemistryHydroIteration();

    bool temperature_converged = temperatureIteration();
    std::cout << "Temperature iteration converged: " << temperature_converged << "\n\n";

    if (temperature_converged == true) break;

    std::vector<double> temperature_change_gas(atmosphere.nb_grid_points, 0);
    std::vector<double> temperature_change_dust(atmosphere.nb_grid_points, 0);

    std::pair<double, size_t> temperature_convergence_gas = checkTemperatureConvergence(
      atmosphere.temperature_gas, 
      temperature_gas_old, 
      temperature_change_gas);

    std::pair<double, size_t> temperature_convergence_dust = checkTemperatureConvergence(
      atmosphere.temperature_dust, 
      temperature_dust_old, 
      temperature_change_dust);
    
    std::cout << "Global iteration " << it << "\n";
    std::cout << "Max T change "
              << temperature_convergence_gas.first << "  " << temperature_convergence_gas.second << "\t"
              << temperature_convergence_dust.first << "  " << temperature_convergence_dust.second << "\n\n";

    //track the per-iteration change (previously these were never updated, so the
    //reported change was cumulative drift from the initial guess)
    temperature_gas_old = atmosphere.temperature_gas;
    temperature_dust_old = atmosphere.temperature_dust;

    //movable grid: redistribute nodes to follow the steep (dust-front) features,
    //then refresh the change references so the next iteration is measured on the
    //new grid
    if (config.use_movable_grid && (it % config.regrid_frequency == 0))
    {
      applyMovableGrid();
      temperature_gas_old = atmosphere.temperature_gas;
      temperature_dust_old = atmosphere.temperature_dust;
    }

    //if (config.output_atmosphere_path != "")
      //atmosphere.writeStructure(config.output_atmosphere_path);

    // if (std::abs(temperature_convergence_gas.first) < 1e-3
    //     && std::abs(temperature_convergence_dust.first) < 1e-3)
    //   break;
  }


  if (config.output_spectrum_path != "")
    radiative_transfer.saveSpectrum(config.output_spectrum_path);

  if (config.output_atmosphere_path != "")
    atmosphere.writeStructure(config.output_atmosphere_path);

  if (config.output_dust_path != "")
    dust_species->saveOutput(config.output_dust_path);

  if (config.output_hydro_path != "")
    hydrodynamics.saveOutput(config.output_hydro_path);
}



bool AGBStarModel::chemistryDustIteration()
{
  //Gas-phase chemistry at the full (undepleted) element abundances. The carbon
  //consumed by dust formation is NOT removed here: the Gail & Sedlmayr moment sweep
  //(GailSedlmayrDust::calcDistribution) depletes the growth species DIFFERENTIALLY
  //along the outward integration, throttling growth/nucleation by (1 - fc) with fc the
  //running degree of condensation. That is the literature-standard, single-pass
  //treatment (GS book Sec. 14.3; Winters Eq. 5.4) and avoids the fragile outer
  //chemistry<->dust fixed-point iteration (which can oscillate, since the moment method
  //recomputes the dust from scratch with no memory of what already condensed).
  std::vector<double> degree_of_condensation(atmosphere.nb_grid_points, 0);

  chemistry.calcChemicalComposition(
    std::vector<double>{},
    atmosphere.temperature_gas,
    atmosphere.pressure_bar,
    degree_of_condensation,
    atmosphere.number_densities,
    atmosphere.mean_molecuar_weight,
    atmosphere.total_element_density,
    atmosphere.total_h_density);

  atmosphere.equationOfState();

  //condensable carbon abundance = eps_C - eps_O (all oxygen assumed locked in CO)
  const double condensable_carbon_abundance =
      chemistry.element_abundances[chemistry.fastchem_element_indices[_C]]
    - chemistry.element_abundances[chemistry.fastchem_element_indices[_O]];

  dust_species->calcDistribution(condensable_carbon_abundance);

  if (config.output_dust_path != "")
    dust_species->saveOutput(config.output_dust_path);

  return true;
}



bool AGBStarModel::chemistryHydroIteration()
{
  std::vector<double> alpha_old = hydrodynamics.alpha;

  for (unsigned int iter=0; iter<config.nb_hydrodynamics_iter; ++iter)
  {
    chemistryDustIteration();
    radiativeTransfer();

    //hand the (frozen) dust kernels to the hydrodynamics so the Henyey solver can
    //couple the dust-moment equations (no-op for the shooting path)
    hydrodynamics.setDustState(
      dust_species->nucleationRate(),
      dust_species->growthTimescale(),
      dust_species->secondMoment());

    hydrodynamics.calcWindVelocity();

    atmosphere.equationOfState();

    std::pair<double, size_t> convergence = checkConvergence(alpha_old, hydrodynamics.alpha);

    std::cout << "Chemistry-Hydro iteration: " << iter << "  Max alpha change " 
              << convergence.first << "  " << convergence.second << "  " 
              << alpha_old[convergence.second] << "  " << hydrodynamics.alpha[convergence.second] << "\n";

    alpha_old = hydrodynamics.alpha;

    //A rejected wind solve (no interior critical point / non-finite Mdot) leaves the
    //structure frozen, so alpha does not change -- that is NOT convergence. Require a
    //genuine, accepted solve before declaring the chemistry-hydro loop converged.
    const bool converged = std::abs(convergence.first) < config.hydrodynamics_convergence
                           && !hydrodynamics.last_solve_rejected;

    if (converged)
      std::cout << "Chemistry-hydrodynamics converged!\n\n";
    else if (hydrodynamics.last_solve_rejected)
      std::cout << "Chemistry-hydrodynamics not converged (wind solve rejected)!\n\n";
    else
      std::cout << "Chemistry-hydrodynamics not converged!\n\n";

    // std::string file_name = "hydro_test_" + std::to_string(iter) + ".dat";
    // hydrodynamics.saveOutput(file_name);

    if (converged)
      break;
  }
  
  //exit(0);

  // for (size_t i=0; i<atmosphere.nb_grid_points; ++i)
  //   std::cout << i << "\t" << degree_of_condensation[i] << "\n";
  // std::cout << chemistry.element_abundances[chemistry.fastchem_element_indices[_C]] << "\t" << chemistry.element_abundances[chemistry.fastchem_element_indices[_O]] << "\n";
  // exit(0);

  return true;
}



void AGBStarModel::radiativeTransfer()
{
  std::cout << "Calculating dust opacities\n";
  for (size_t i=0; i<atmosphere.nb_grid_points; ++i)
    dust_species->calcTransportCoefficients(
      i, 
      atmosphere.absorption_coeff_dust[i], 
      atmosphere.scattering_coeff_dust[i]);
  
  std::cout << "Calculating gas opacities\n\n";
  for (size_t i=0; i<atmosphere.nb_grid_points; ++i)
    transport_coeff.calculate(
      atmosphere.temperature_gas[i], 
      atmosphere.pressure_bar[i], 
      atmosphere.number_densities[i], 
      atmosphere.absorption_coeff_gas[i], 
      atmosphere.scattering_coeff_gas[i]);
  
  for (size_t i=0; i<atmosphere.nb_grid_points; ++i)
  {
    for (size_t j=0; j<spectral_grid.nbSpectralPoints(); ++j)
    {
      atmosphere.absorption_coeff[i][j] = atmosphere.absorption_coeff_gas[i][j] + atmosphere.absorption_coeff_dust[i][j];
      atmosphere.scattering_coeff[i][j] = atmosphere.scattering_coeff_gas[i][j] + atmosphere.scattering_coeff_dust[i][j];

      atmosphere.extinction_coeff[i][j] = atmosphere.absorption_coeff[i][j] + atmosphere.scattering_coeff[i][j];
      atmosphere.extinction_coeff_gas[i][j] = atmosphere.absorption_coeff_gas[i][j] + atmosphere.scattering_coeff_gas[i][j];
      atmosphere.extinction_coeff_dust[i][j] = atmosphere.absorption_coeff_dust[i][j] + atmosphere.scattering_coeff_dust[i][j];
    }
  }
  
  radiative_transfer.solveRadiativeTransfer();
}



void AGBStarModel::applyMovableGrid()
{
  const size_t np = atmosphere.nb_grid_points;

  //monitor quantities on the current grid: flux-mean extinction, gas temperature,
  //and dust nucleation rate
  std::vector<double> chi_h(np, 0.);
  for (size_t i=0; i<np; ++i)
    chi_h[i] = radiative_transfer.radiation_field[i].fluxWeightedExtinction(
                 atmosphere.extinction_coeff[i]);

  std::vector<double> jstar = dust_species->nucleationRate();
  if (jstar.size() != np) jstar.assign(np, 1.0);

  const std::vector<std::vector<double>> quantities =
    { chi_h, atmosphere.temperature_gas, jstar, atmosphere.velocity };
  const std::vector<double> weights =
    { config.monitor_weight_opacity,
      config.monitor_weight_temperature,
      config.monitor_weight_nucleation,
      config.monitor_weight_velocity };

  const std::vector<double> target = aux::equidistributedGrid(
    atmosphere.radius, quantities, weights,
    config.monitor_smoothing_passes, config.monitor_max, config.monitor_rel_floor);

  //under-relax the node motion; pin the boundaries and keep the grid strictly monotone
  std::vector<double> new_radius(np);
  for (size_t i=0; i<np; ++i)
    new_radius[i] = atmosphere.radius[i]
                  + config.grid_relaxation * (target[i] - atmosphere.radius[i]);
  new_radius[0]    = atmosphere.radius[0];
  new_radius[np-1] = atmosphere.radius[np-1];
  for (size_t i=1; i<np; ++i)
    if (new_radius[i] <= new_radius[i-1])
      new_radius[i] = new_radius[i-1] * (1.0 + 1e-6);

  //remap the structure and rebuild the grid-dependent geometry / warm starts
  atmosphere.remapToGrid(new_radius);
  radiative_transfer.rebuildGeometry();
  hydrodynamics.resetWarmStart();

  //reset per-node temperature-iteration state (re-initialises lazily on the new grid)
  relaxation_gas.clear();    relaxation_dust.clear();
  prev_delta_b_gas.clear();  prev_delta_b_dust.clear();
  anderson_x_gas.clear();    anderson_f_gas.clear();
  anderson_x_dust.clear();   anderson_f_dust.clear();
  anderson_was_active = false;
  prev_max_rel_change = 1e30;

  std::cout << "Movable grid: regridded; inner/outer r/R* = "
            << atmosphere.radius_grid.front() << " / " << atmosphere.radius_grid.back() << "\n\n";
}


bool AGBStarModel::temperatureIteration()
{
  bool converged = false;

  //dust opacities do not depend on temperature
  //so, we only calculate them once
  std::cout << "Calculating dust opacities\n";
  for (size_t i=0; i<atmosphere.nb_grid_points; ++i)
    dust_species->calcTransportCoefficients(
      i, 
      atmosphere.absorption_coeff_dust[i], 
      atmosphere.scattering_coeff_dust[i]);


  //for (unsigned int it=0; it<config.nb_temperature_iter; ++it)
  for (unsigned int it=0; it<1; ++it)
  {
    std::cout << "Calculating gas opacities\n\n";
    for (size_t i=0; i<atmosphere.nb_grid_points; ++i)
      transport_coeff.calculate(
        atmosphere.temperature_gas[i], 
        atmosphere.pressure_bar[i], 
        atmosphere.number_densities[i], 
        atmosphere.absorption_coeff_gas[i], 
        atmosphere.scattering_coeff_gas[i]);

    for (size_t i=0; i<atmosphere.nb_grid_points; ++i)
    {
      for (size_t j=0; j<spectral_grid.nbSpectralPoints(); ++j)
      {
        //atmosphere.scattering_coeff_gas[i][j] = 0;
        //atmosphere.scattering_coeff_dust[i][j] = 0;
        atmosphere.absorption_coeff[i][j] = atmosphere.absorption_coeff_gas[i][j] + atmosphere.absorption_coeff_dust[i][j];
        atmosphere.scattering_coeff[i][j] = atmosphere.scattering_coeff_gas[i][j] + atmosphere.scattering_coeff_dust[i][j];

        atmosphere.extinction_coeff[i][j] = atmosphere.absorption_coeff[i][j] + atmosphere.scattering_coeff[i][j];
        atmosphere.extinction_coeff_gas[i][j] = atmosphere.absorption_coeff_gas[i][j] + atmosphere.scattering_coeff_gas[i][j];
        atmosphere.extinction_coeff_dust[i][j] = atmosphere.absorption_coeff_dust[i][j] + atmosphere.scattering_coeff_dust[i][j];
      }
    }

    radiative_transfer.solveRadiativeTransfer();

    //current profiles (x_k of the fixed-point map)
    std::vector<double> temperature_gas_prev  = atmosphere.temperature_gas;
    std::vector<double> temperature_dust_prev = atmosphere.temperature_dust;

    std::vector<double> delta_temperature_gas;
    std::vector<double> delta_temperature_dust;

    if (linearisation_active)
    {
      //full-linearisation Newton step: the exact dJ/dT response of the converged RT
      //operator drives both local radiative-equilibrium residuals to zero, removing
      //the Unsoeld-Lucy accuracy floor. Anderson is bypassed (Newton already mixes).
      radiative_transfer.linearisedTemperatureCorrection(
        delta_temperature_gas, delta_temperature_dust);

      //Per-layer adaptive under-relaxation of the Newton step (same mechanism as the
      //Unsoeld-Lucy path): in the COUPLED loop the structure (dust formation, hydro)
      //updates between steps, and a full Newton step can over-correct, driving a
      //period-2 +/-cap limit cycle (seen most strongly just beyond the dust front
      //where nucleation is extremely T-sensitive). Halve omega on a sign flip, grow it
      //back when the step keeps its sign; apply omega AFTER capping the raw step so the
      //damping is visible even while the cap binds (a step pinned at +/-cap with a
      //flipping sign is exactly the oscillation we must damp).
      const size_t ng = atmosphere.nb_grid_points;
      if (relaxation_gas.size()    != ng) relaxation_gas.assign(ng, config.temperature_relaxation_init);
      if (relaxation_dust.size()   != ng) relaxation_dust.assign(ng, config.temperature_relaxation_init);
      if (prev_delta_b_gas.size()  != ng) prev_delta_b_gas.assign(ng, 0.);
      if (prev_delta_b_dust.size() != ng) prev_delta_b_dust.assign(ng, 0.);

      auto damp = [&](std::vector<double>& delta, std::vector<double>& relax,
                      std::vector<double>& prev_delta, const std::vector<double>& temperature)
      {
        for (size_t i=0; i<ng; ++i)
        {
          if (delta[i] * prev_delta[i] < 0.) relax[i] *= config.temperature_relaxation_down;
          else                               relax[i] *= config.temperature_relaxation_up;
          relax[i] = std::min(std::max(relax[i], config.temperature_relaxation_min),
                              config.temperature_relaxation_max);
          prev_delta[i] = delta[i];

          const double cap = config.temperature_max_change * temperature[i];
          if (std::abs(delta[i]) > cap) delta[i] = std::copysign(cap, delta[i]);
          delta[i] *= relax[i];   //omega applied AFTER the cap
        }
      };
      damp(delta_temperature_gas,  relaxation_gas,  prev_delta_b_gas,  temperature_gas_prev);
      damp(delta_temperature_dust, relaxation_dust, prev_delta_b_dust, temperature_dust_prev);
    }
    else
    {
      //relaxed Unsoeld-Lucy correction; f_k = G(x_k) - x_k = delta_temperature.
      //calculate() does the exact T^4 update, per-layer adaptive damping and
      //correction smoothing, and carries the relaxation/prev-correction state.
      delta_temperature_gas = temperature_correction.calculate(
        atmosphere.temperature_gas,
        atmosphere.radius,
        atmosphere.extinction_coeff,
        atmosphere.absorption_coeff_gas,
        relaxation_gas,
        prev_delta_b_gas);

      delta_temperature_dust = temperature_correction.calculate(
        atmosphere.temperature_dust,
        atmosphere.radius,
        atmosphere.extinction_coeff,
        atmosphere.absorption_coeff_dust,
        relaxation_dust,
        prev_delta_b_dust);
    }

    //Anderson acceleration only once we are in the settling regime (per-step change
    //already well below the cap), so the early large-change transient that the
    //hydro-dust cycle is sensitive to stays on plain damped steps. Restart the
    //history fresh on activation so it is not seeded by the transient.
    const bool accelerate =
      config.use_anderson_acceleration
      && !linearisation_active
      && temp_iter_count >= config.anderson_start_iter
      && prev_max_rel_change < config.anderson_activation_fraction * config.temperature_max_change;

    if (accelerate && !anderson_was_active)
    {
      anderson_x_gas.clear();  anderson_f_gas.clear();
      anderson_x_dust.clear(); anderson_f_dust.clear();
    }
    anderson_was_active = accelerate;

    std::vector<double> temperature_gas_acc, temperature_dust_acc;
    if (accelerate)
    {
      temperature_gas_acc = andersonStep(
        anderson_x_gas, anderson_f_gas, temperature_gas_prev, delta_temperature_gas);
      temperature_dust_acc = andersonStep(
        anderson_x_dust, anderson_f_dust, temperature_dust_prev, delta_temperature_dust);
    }

    for (size_t i=0; i<atmosphere.nb_grid_points; ++i)
    {
      atmosphere.temperature_gas[i] = accelerate ?
        temperature_gas_acc[i]  : temperature_gas_prev[i]  + delta_temperature_gas[i];
      atmosphere.temperature_dust[i] = accelerate ?
        temperature_dust_acc[i] : temperature_dust_prev[i] + delta_temperature_dust[i];
    }

    //Optionally enforce monotonicity before the cap (default off: it fights the
    //correction and drives a +/-cap limit cycle, see config). When enabled, doing
    //it here lets the cap below rein the clip back to <= the per-step change so a
    //non-monotonic feature is flattened gradually rather than all at once.
    if (config.force_monotonic_temperature)
    {
      forceMonotonicProfile(atmosphere.temperature_gas);
      forceMonotonicProfile(atmosphere.temperature_dust);
    }

    //Cap the actually applied step LAST, so nothing (Anderson overshoot or the
    //monotonicity clip) can produce a change larger than the per-step bound the
    //hydro-dust cycle needs. Early on this binds; near convergence the change is
    //below the cap, so it no longer interferes.
    for (size_t i=0; i<atmosphere.nb_grid_points; ++i)
    {
      const double max_gas = config.temperature_max_change * temperature_gas_prev[i];
      const double diff_gas = atmosphere.temperature_gas[i] - temperature_gas_prev[i];
      if (std::abs(diff_gas) > max_gas)
        atmosphere.temperature_gas[i] = temperature_gas_prev[i] + std::copysign(max_gas, diff_gas);
      if (atmosphere.temperature_gas[i] <= 0.)
        atmosphere.temperature_gas[i] = 0.5 * temperature_gas_prev[i];

      const double max_dust = config.temperature_max_change * temperature_dust_prev[i];
      const double diff_dust = atmosphere.temperature_dust[i] - temperature_dust_prev[i];
      if (std::abs(diff_dust) > max_dust)
        atmosphere.temperature_dust[i] = temperature_dust_prev[i] + std::copysign(max_dust, diff_dust);
      if (atmosphere.temperature_dust[i] <= 0.)
        atmosphere.temperature_dust[i] = 0.5 * temperature_dust_prev[i];
    }

    ++temp_iter_count;

    //largest relative change actually applied this step (a settling indicator used
    //in the convergence test so we never declare convergence mid-wobble)
    double max_rel_change = 0.;
    for (size_t i=1; i<atmosphere.nb_grid_points; ++i)
    {
      max_rel_change = std::max(max_rel_change,
        std::abs((atmosphere.temperature_gas[i]  - temperature_gas_prev[i])  / temperature_gas_prev[i]));
      max_rel_change = std::max(max_rel_change,
        std::abs((atmosphere.temperature_dust[i] - temperature_dust_prev[i]) / temperature_dust_prev[i]));
    }
    prev_max_rel_change = max_rel_change;

    auto flux_convergence = checkFluxConvergence();
    
    std::vector<double> energy_balance_dust;
    std::vector<double> energy_balance_gas;

    auto max_energy_balance_gas = checkEnergyBalance(
      atmosphere.temperature_gas, 
      atmosphere.absorption_coeff_gas,
      energy_balance_gas);
    auto max_energy_balance_dust = checkEnergyBalance(
      atmosphere.temperature_dust,
      atmosphere.absorption_coeff_dust,
      energy_balance_dust);

    //Two-phase corrector switch (Unsoeld-Lucy -> linearisation): latch on the Newton
    //linearisation only once the profile has settled (small per-step change) AND is near
    //radiative equilibrium (small energy-balance residual) for a few consecutive steps.
    //Far from convergence the frozen-opacity Newton overshoots (esp. via the very
    //temperature-sensitive dust nucleation), so UL is used to smooth the coarse errors first.
    if (config.use_linearisation && config.linearisation_start_unsoeld_lucy && !linearisation_active)
    {
      const double re_residual = std::max(std::abs(max_energy_balance_gas.first),
                                          std::abs(max_energy_balance_dust.first));
      const bool settling = max_rel_change
                          < config.linearisation_switch_dt_fraction * config.temperature_max_change;
      const bool near_re  = re_residual < config.linearisation_switch_re_residual;

      if (settling && near_re) ++linearisation_ready_count;
      else                     linearisation_ready_count = 0;

      if (linearisation_ready_count >= config.linearisation_switch_count)
      {
        linearisation_active = true;
        std::cout << "[corrector] Unsoeld-Lucy -> linearisation at temperature iteration "
                  << temp_iter_count << " (max|dT|/T=" << max_rel_change
                  << ", RE residual=" << re_residual << ")\n";
      }
    }

    // auto max_energy_balance_gas = checkEnergyBalance(
    //   atmosphere.temperature_gas, 
    //   atmosphere.absorption_coeff,
    //   energy_balance_gas);
    // auto max_energy_balance_dust = checkEnergyBalance(
    //   atmosphere.temperature_dust,
    //   atmosphere.absorption_coeff,
    //   energy_balance_dust);

    //RT_FLUX_CHECK: verify the RT is internally flux-consistent, i.e. that its flux
    //and its mean intensity satisfy the zeroth moment equation
    //  d(r^2 H)/dr = r^2 * integral kappa_abs (B - J) dnu
    //(gas emits at T_gas, dust at T_dust; scattering cancels out). A mismatch means
    //the moment-system J and the calcFlux H disagree across e.g. the steep dust
    //opacity gradient -> flux non-conservation that no temperature correction can fix.
    if (std::getenv("RT_FLUX_CHECK"))
    {
      const size_t np = atmosphere.nb_grid_points;
      const size_t nl = spectral_grid.nbSpectralPoints();
      std::vector<double> r2h(np, 0.), rhs(np, 0.);

      for (size_t i=0; i<np; ++i)
      {
        r2h[i] = radiative_transfer.radiation_field[i].eddington_flux_int_conservative
               * atmosphere.radius[i]*atmosphere.radius[i];

        std::vector<double> y(nl, 0.);
        for (size_t j=0; j<nl; ++j)
        {
          const double Jnu = radiative_transfer.radiation_field[i].mean_intensity[j];
          const double Bg  = aux::planckFunctionWavelength(atmosphere.temperature_gas[i],  spectral_grid.wavelength_list[j]);
          const double Bd  = aux::planckFunctionWavelength(atmosphere.temperature_dust[i], spectral_grid.wavelength_list[j]);
          y[j] = atmosphere.absorption_coeff_gas[i][j]  * (Bg - Jnu)
               + atmosphere.absorption_coeff_dust[i][j] * (Bd - Jnu);
        }
        rhs[i] = atmosphere.radius[i]*atmosphere.radius[i]
               * radiative_transfer.radiation_field[i].wavelengthIntegration(y);
      }

      std::cout << "\n[RTFC] i\tr2H\td(r2H)/dr\tr2*kabs(B-J)\tratio\n";
      for (size_t i=1; i+1<np; ++i)
      {
        const double dflux = (r2h[i+1]-r2h[i-1])/(atmosphere.radius[i+1]-atmosphere.radius[i-1]);
        std::cout << "[RTFC] " << i << "\t" << r2h[i] << "\t" << dflux << "\t" << rhs[i]
                  << "\t" << (rhs[i]!=0. ? dflux/rhs[i] : 0.) << "\n";
      }
      std::cout << "\n";
    }

    for (size_t i=0; i<atmosphere.nb_grid_points; ++i)
      std::cout << i << "\t" << atmosphere.temperature_gas[i] << "\t"
                     << atmosphere.temperature_dust[i] << "\t"
                     << delta_temperature_gas[i] << "\t"
                     << delta_temperature_dust[i] << "\t"
                     << radiative_transfer.radiation_field[i].eddington_flux_int_conservative*atmosphere.radius[i]*atmosphere.radius[i] << "\t"
                     << energy_balance_gas[i] << "\t"
                     << energy_balance_dust[i] << "\t"
                     << "\n";
    
    std::cout << "\nTemperature iteration: " << it << "\n";
    std::cout << "Max relative T change: " << max_rel_change << "\n";
    std::cout << "Flux convergence (vs L target): "
              << flux_convergence.first << "\t" << flux_convergence.second << "\t"
              << radiative_transfer.radiation_field[flux_convergence.second].eddington_flux_int
                *atmosphere.radius[flux_convergence.second]
                *atmosphere.radius[flux_convergence.second] << "\n";

    std::cout << "Energy balance gas: "
              << max_energy_balance_gas.first << "\t" << max_energy_balance_gas.second << "\t"
              << atmosphere.temperature_gas[max_energy_balance_gas.second] << "\t"
              << delta_temperature_gas[max_energy_balance_gas.second] << "\n";
    std::cout << "Energy balance dust: "
              << max_energy_balance_dust.first << "\t" <<max_energy_balance_dust.second << "\t"
              << atmosphere.temperature_dust[max_energy_balance_dust.second] << "\t"
              << delta_temperature_dust[max_energy_balance_dust.second] << "\n";

    std::cout << "\n";

    //combined criterion: the model must be flux-constant at the target luminosity,
    //in local energy balance, AND have stopped changing (no residual wobble).
    if (std::abs(flux_convergence.first) < config.temperature_convergence
        && std::abs(max_energy_balance_gas.first) < config.temperature_convergence
        && std::abs(max_energy_balance_dust.first) < config.temperature_convergence
        && max_rel_change < config.temperature_convergence)
    {
      converged = true;
      break;
    }
  }

  return converged;
}



std::pair<double, size_t> AGBStarModel::checkFluxConvergence()
{
  std::pair<double, size_t> max_difference{0.0, 0};

  //Measure deviation against the target luminosity r^2 H = L / (16 pi^2), the same
  //target the Unsoeld-Lucy correction drives toward. (Referencing r^2 H at the
  //inner boundary instead lets that point's own wobble pollute the yardstick.)
  const double target_flux = config.stellar_luminosity
                           / (16. * constants::pi * constants::pi);

  for (size_t i=0; i<atmosphere.nb_grid_points; ++i)
  {
    const double flux = radiative_transfer.radiation_field[i].eddington_flux_int_conservative
                        * atmosphere.radius[i]*atmosphere.radius[i];
    const double rel_difference = (flux - target_flux)/target_flux;

    if (std::abs(rel_difference) > std::abs(max_difference.first))
    {
      max_difference.first = rel_difference;
      max_difference.second = i;
    }
  }

  return max_difference;
}



std::vector<double> AGBStarModel::andersonStep(
  std::vector<std::vector<double>>& x_history,
  std::vector<std::vector<double>>& f_history,
  const std::vector<double>& x_k,
  const std::vector<double>& f_k)
{
  const size_t n = x_k.size();

  //record the current iterate/residual and trim to the window (m differences need
  //m+1 stored pairs)
  x_history.push_back(x_k);
  f_history.push_back(f_k);

  while (x_history.size() > config.anderson_window + 1)
  {
    x_history.erase(x_history.begin());
    f_history.erase(f_history.begin());
  }

  const size_t stored = x_history.size();

  //plain (damped) step G(x_k) = x_k + f_k until we have at least one difference
  std::vector<double> result(n);
  for (size_t i=0; i<n; ++i)
    result[i] = x_k[i] + f_k[i];

  if (stored < 2)
    return result;

  const size_t m = stored - 1; //number of difference columns

  //difference matrices: dF[j] = f_{j+1} - f_j, dX[j] = x_{j+1} - x_j
  std::vector<std::vector<double>> dF(m, std::vector<double>(n));
  std::vector<std::vector<double>> dX(m, std::vector<double>(n));

  for (size_t j=0; j<m; ++j)
    for (size_t i=0; i<n; ++i)
    {
      dF[j][i] = f_history[j+1][i] - f_history[j][i];
      dX[j][i] = x_history[j+1][i] - x_history[j][i];
    }

  //normal equations (dF^T dF) gamma = dF^T f_k  with light Tikhonov regularisation
  std::vector<std::vector<double>> a(m, std::vector<double>(m, 0.));
  std::vector<double> b(m, 0.);

  for (size_t j=0; j<m; ++j)
  {
    for (size_t k=0; k<m; ++k)
    {
      double sum = 0.;
      for (size_t i=0; i<n; ++i) sum += dF[j][i] * dF[k][i];
      a[j][k] = sum;
    }
    double sum = 0.;
    for (size_t i=0; i<n; ++i) sum += dF[j][i] * f_k[i];
    b[j] = sum;
  }

  double trace = 0.;
  for (size_t j=0; j<m; ++j) trace += a[j][j];
  const double reg = 1e-10 * (trace/m + 1e-300);
  for (size_t j=0; j<m; ++j) a[j][j] += reg;

  std::vector<double> gamma = solveLinearSystem(a, b);

  //if the solve failed (singular), fall back to the plain damped step
  if (gamma.empty())
    return result;

  //x_{k+1} = (x_k + f_k) - sum_j gamma_j (dX_j + dF_j)
  for (size_t j=0; j<m; ++j)
    for (size_t i=0; i<n; ++i)
      result[i] -= gamma[j] * (dX[j][i] + dF[j][i]);

  return result;
}



std::pair<double, size_t> AGBStarModel::checkEnergyBalance(
  std::vector<double>& temperature,
  std::vector<std::vector<double>>& absorption__coeff,
  std::vector<double>& deviation)
{
  std::pair<double, size_t> max_deviation{0.0, 0};
  deviation.assign(temperature. size(), 0.);

  for (size_t i=0; i<atmosphere.nb_grid_points; ++i)
  {
    std::vector<double> y1(spectral_grid.nbSpectralPoints(), 0);
    std::vector<double> y2(spectral_grid.nbSpectralPoints(), 0);

    for (size_t j=0; j<spectral_grid.nbSpectralPoints(); ++j)
    {
      y1[j] = absorption__coeff[i][j] * radiative_transfer.radiation_field[i].mean_intensity[j];
      y2[j] = absorption__coeff[i][j] * aux::planckFunctionWavelength(temperature[i], spectral_grid.wavelength_list[j]);
    }

    deviation[i] = radiative_transfer.radiation_field[i].wavelengthIntegration(y1) 
                           / radiative_transfer.radiation_field[i].wavelengthIntegration(y2) 
                           - 1.;

    if (std::abs(deviation[i]) > std::abs(max_deviation.first))
    {
      max_deviation.first = deviation[i];
      max_deviation.second = i;
    }
  }

  return max_deviation;
}



std::pair<double, size_t> AGBStarModel::checkTemperatureConvergence(
  const std::vector<double>& temperature,
  const std::vector<double>& temperature_old,
  std::vector<double>& change)
{
  std::pair<double, size_t> max_change{0.0, 0};
  change.assign(temperature. size(), 0.);

  for (size_t i=1; i<atmosphere.nb_grid_points; ++i)
  {
    const double rel_difference = (temperature[i] - temperature_old[i])/temperature_old[i];
    change[i] = rel_difference;

    if (std::abs(rel_difference) > std::abs(max_change.first))
    {
      max_change.first = rel_difference;
      max_change.second = i;
    }
  }

  return max_change;
}



std::pair<double, size_t> AGBStarModel::checkConvergence(
  const std::vector<double>& old_data,
  const std::vector<double>& new_data)
{
  std::pair<double, size_t> max_change{0.0, 0};

  for (size_t i=1; i<atmosphere.nb_grid_points; ++i)
  {
    const double rel_difference = (new_data[i] - old_data[i])/old_data[i];

    if (std::abs(rel_difference) > std::abs(max_change.first))
    {
      max_change.first = rel_difference;
      max_change.second = i;
    }
  }

  return max_change;
}




void AGBStarModel::forceMonotonicProfile(std::vector<double>& data)
{
  for (size_t i=1; i<data.size(); ++i)
    if (data[i] > data[i-1])
      data[i] = 0.99*data[i-1];
}


void AGBStarModel::smoothProfile(std::vector<double>& data)
{

  for (size_t i=1; i<data.size()-1; ++i)
    data[i] = 0.5 * data[i] + 0.25*(data[i+1] + data[i-1]);

}


}