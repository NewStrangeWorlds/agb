 
#include "agb_model.h"

#include <vector>
#include <string>



namespace agb{

AGBStarModel::AGBStarModel(const std::string folder)
 : config(folder)
 , spectral_grid(&config)
 , atmosphere(&config)
 , chemistry(config.model_folder+config.fastchem_parameter_file)
 , dust_species(&config, &spectral_grid, &atmosphere)
 , transport_coeff(&config, &spectral_grid, config.opacity_species_symbol, config.opacity_species_folder)
 , radiative_transfer(&config, &spectral_grid, &atmosphere)
 , temperature_correction(&config, &spectral_grid, &radiative_transfer.radiation_field)
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
  
  dust_species.calcFormation();

  for (size_t i=0; i<atmosphere.nb_grid_points; ++i)
    dust_species.calcTransportCoefficients(
      i, 
      atmosphere.absorption_coeff_dust[i], 
      atmosphere.scattering_coeff_dust[i]);

  for (size_t i=0; i<atmosphere.nb_grid_points; ++i)
    transport_coeff.calculate(
      atmosphere.temperature_gas[i], 
      atmosphere.pressure_bar[i], 
      atmosphere.number_densities[i], 
      atmosphere.absorption_coeff_gas[i], 
      atmosphere.scattering_coeff_gas[i]);

  atmosphere.absorption_coeff.assign(
    atmosphere.nb_grid_points, 
    std::vector<double>(spectral_grid.nbSpectralPoints(), 0.));

  atmosphere.scattering_coeff.assign(
    atmosphere.nb_grid_points, 
    std::vector<double>(spectral_grid.nbSpectralPoints(), 0.));

  for (size_t i=0; i<atmosphere.nb_grid_points; ++i)
  {
    for (size_t j=0; j<spectral_grid.nbSpectralPoints(); ++j)
    {
      atmosphere.absorption_coeff[i][j] = atmosphere.absorption_coeff_gas[i][j] + atmosphere.absorption_coeff_dust[i][j];
      atmosphere.scattering_coeff[i][j] = atmosphere.scattering_coeff_gas[i][j] + atmosphere.scattering_coeff_dust[i][j];
    }
  }


  radiative_transfer.solveRadiativeTransfer();

  if (config.output_spectrum_path != "")
    radiative_transfer.saveSpectrum(config.output_spectrum_path);
}



}