 
#include "agb_model.h"

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

#include "../additional/aux_functions.h"


namespace agb{

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

    std::vector<double> delta_temperature_gas = temperature_correction.calculate(
      atmosphere.temperature_gas,
      atmosphere.radius,
      atmosphere.extinction_coeff,
      atmosphere.absorption_coeff_gas);

    std::vector<double> delta_temperature_dust = temperature_correction.calculate(
      atmosphere.temperature_dust,
      atmosphere.radius,
      atmosphere.extinction_coeff,
      atmosphere.absorption_coeff_dust);

    // std::vector<double> delta_temperature_gas = temperature_correction.calculate(
    //   atmosphere.temperature_gas,
    //   atmosphere.radius,
    //   atmosphere.extinction_coeff_gas,
    //   atmosphere.absorption_coeff_gas);

    // std::vector<double> delta_temperature_dust = delta_temperature_gas;

    std::vector<double> temperature_gas_new = atmosphere.temperature_gas;
    std::vector<double> temperature_dust_new = atmosphere.temperature_dust;
    
    for (size_t i=0; i<atmosphere.nb_grid_points; ++i)
    {
      temperature_gas_new[i] += delta_temperature_gas[i];
      temperature_dust_new[i] += delta_temperature_dust[i];
    }

    for (size_t i=0; i<atmosphere.nb_grid_points; ++i)
    {
      atmosphere.temperature_gas[i] = std::sqrt(atmosphere.temperature_gas[i] * temperature_gas_new[i]);
      atmosphere.temperature_dust[i] = std::sqrt(atmosphere.temperature_dust[i] * temperature_dust_new[i]);
    }


    // for (size_t i=0; i<atmosphere.nb_grid_points; ++i)
    // {
    //   atmosphere.temperature_gas[i] += delta_temperature_gas[i];
    //   atmosphere.temperature_dust[i] += delta_temperature_dust[i];
    // }


    if (config.smooth_temperature_profile)
    {
      smoothProfile(atmosphere.temperature_gas);
      smoothProfile(atmosphere.temperature_dust);
    }

    forceMonotonicProfile(atmosphere.temperature_gas);
    forceMonotonicProfile(atmosphere.temperature_dust);

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

    // auto max_energy_balance_gas = checkEnergyBalance(
    //   atmosphere.temperature_gas, 
    //   atmosphere.absorption_coeff,
    //   energy_balance_gas);
    // auto max_energy_balance_dust = checkEnergyBalance(
    //   atmosphere.temperature_dust,
    //   atmosphere.absorption_coeff,
    //   energy_balance_dust);

    for (size_t i=0; i<atmosphere.nb_grid_points; ++i)
      std::cout << i << "\t" << atmosphere.temperature_gas[i] << "\t" 
                     << atmosphere.temperature_dust[i] << "\t" 
                     << delta_temperature_gas[i] << "\t" 
                     << delta_temperature_dust[i] << "\t" 
                     << radiative_transfer.radiation_field[i].eddington_flux_int*atmosphere.radius[i]*atmosphere.radius[i] << "\t" 
                     << energy_balance_gas[i] << "\t" 
                     << energy_balance_dust[i] << "\t" 
                     << "\n";
    
    std::cout << "\nTemperature iteration: " << it << "\n";
    std::cout << "Flux convergence: " 
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

    // if (std::abs(flux_convergence.first) < config.temperature_convergence
    //     && std::abs(max_energy_balance_gas.first) < config.temperature_convergence
    //     && std::abs(max_energy_balance_dust.first) < config.temperature_convergence)
    // {
    //   converged = true;
    //   break;
    // }

    if (std::abs(flux_convergence.first) < config.temperature_convergence)
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

  const double constant_flux = radiative_transfer.radiation_field[0].eddington_flux_int 
                             * atmosphere.radius[0]*atmosphere.radius[0];

  for (size_t i=1; i<atmosphere.nb_grid_points; ++i)
  {
    const double flux = radiative_transfer.radiation_field[i].eddington_flux_int 
                        * atmosphere.radius[i]*atmosphere.radius[i];
    const double rel_difference = (constant_flux - flux)/constant_flux;

    if (std::abs(rel_difference) > std::abs(max_difference.first))
    {
      max_difference.first = rel_difference;
      max_difference.second = i;
    }
  }

  return max_difference;
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