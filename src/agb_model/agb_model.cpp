 
#include "agb_model.h"

#include <vector>
#include <string>

#include "../additional/aux_functions.h"


namespace agb{

AGBStarModel::AGBStarModel(const std::string folder)
 : config(folder)
 , spectral_grid(&config)
 , atmosphere(&config, spectral_grid.nbSpectralPoints())
 , chemistry(config.model_folder+config.fastchem_parameter_file)
 , dust_species(new AnalyticDust(&config, &spectral_grid, &atmosphere, 0.1))
 , transport_coeff(&config, &spectral_grid, config.opacity_species_symbol, config.opacity_species_folder)
 , radiative_transfer(&config, &spectral_grid, &atmosphere)
 , temperature_correction(&config, &spectral_grid, radiative_transfer.radiation_field)
{

}


void AGBStarModel::calcModel()
{
  chemistry.calcChemicalComposition(
    std::vector<double>{1.0, config.c_o_ratio}, 
    atmosphere.temperature_gas, 
    atmosphere.pressure_bar, 
    atmosphere.number_densities, 
    atmosphere.mean_molecuar_weight,
    atmosphere.total_element_density,
    atmosphere.total_h_density);

  atmosphere.equationOfState();

  dust_species->calcDistribution();

  bool temperature_converged = temperatureIteration();
  std::cout << "Temperature iteration converged: " << temperature_converged << "\n\n";


  if (config.output_spectrum_path != "")
    radiative_transfer.saveSpectrum(config.output_spectrum_path);

  if (config.output_atmosphere_path != "")
    atmosphere.writeStructure(config.output_atmosphere_path);
}



bool AGBStarModel::temperatureIteration()
{
  bool converged = false;

  //dust opacities do not depend on temperature
  //so, we only calculate them once
  for (size_t i=0; i<atmosphere.nb_grid_points; ++i)
    dust_species->calcTransportCoefficients(
      i, 
      atmosphere.absorption_coeff_dust[i], 
      atmosphere.scattering_coeff_dust[i]);


  for (unsigned int it=0; it<config.nb_temperature_iter; ++it)
  {
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

    //std::vector<double> delta_temperature_dust = delta_temperature_gas;

    for (size_t i=0; i<atmosphere.nb_grid_points; ++i)
    {
      atmosphere.temperature_gas[i] += delta_temperature_gas[i];
      atmosphere.temperature_dust[i] += delta_temperature_dust[i];
    }


    if (config.smooth_temperature_profile)
    {
      smoothProfile(atmosphere.temperature_gas);
      smoothProfile(atmosphere.temperature_dust);
    }

    //forceMonotonicProfile(atmosphere.temperature_gas);
    //forceMonotonicProfile(atmosphere.temperature_dust);

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

    if (std::abs(flux_convergence.first) < config.temperature_convergence 
        && std::abs(max_energy_balance_gas.first) < config.temperature_convergence 
        && std::abs(max_energy_balance_dust.first) < config.temperature_convergence)
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



void AGBStarModel::forceMonotonicProfile(std::vector<double>& data)
{
  /*for (size_t i=1; i<data.size(); ++i)
    if (data[i] > data[i-1])
      data[i] = 0.99*data[i-1];*/
  if (data[1] > data[0]) data[0] = 1.01 * data[1];
}


void AGBStarModel::smoothProfile(std::vector<double>& data)
{

  for (size_t i=1; i<data.size()-1; ++i)
    data[i] = 0.5 * data[i] + 0.25*(data[i+1] + data[i-1]);

}


}