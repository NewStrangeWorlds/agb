 
#include "agb_model.h"

#include <vector>
#include <string>



namespace agb{

AGBStarModel::AGBStarModel(const std::string folder)
 : config(folder)
 , spectral_grid(&config)
 , atmosphere(&config)
 , chemistry(config.model_folder+config.fastchem_parameter_file)
 , transport_coeff(&config, &spectral_grid, config.opacity_species_symbol, config.opacity_species_folder)
{
  // std::cout << "spectral points " << spectral_grid.nbSpectralPoints() << "\n";

  // for (size_t i=0; i<spectral_grid.nbSpectralPoints(); ++i)
  // {
  //   std::cout << spectral_grid.wavelength_list[i] << "\t" << spectral_grid.wavenumber_list[i] << "\n";
  // }
  
  chemistry.calcChemicalComposition(
    std::vector<double>{1.0, config.c_o_ratio}, 
    atmosphere.temperature_gas, 
    atmosphere.pressure_bar, 
    atmosphere.number_densities, 
    atmosphere.mean_molecuar_weight);
  
  std::vector<double> pressure_old = atmosphere.pressure;

  atmosphere.equationOfState();
  
  for (size_t i=0; i<atmosphere.nb_grid_points; ++i)
    transport_coeff.calculate(
      atmosphere.temperature_gas[i], 
      atmosphere.pressure_bar[i], 
      atmosphere.number_densities[i], 
      atmosphere.absorption_coeff[i], 
      atmosphere.scattering_coeff[i]);
}



}