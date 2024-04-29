 
#include "atmosphere.h"


#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <omp.h>
#include <stdlib.h>

#include "../config/config.h"
#include "../additional/physical_const.h"
#include "../additional/exceptions.h"

namespace agb {

Atmosphere::Atmosphere(ModelConfig* config_)
 : config(config_)
{ 
  std::string file_path = config->model_folder + config->starting_model_path;
  readStructure(file_path);
  
  std::cout << "\nRead-in atmosphere structure:\n";
  for (size_t i=0; i<nb_grid_points; ++i)
  {
    std::cout << std::setprecision(5) << std::scientific  
              << radius_grid[i]
              << "\t" << radius[i] 
              << "\t" << mass_density[i] 
              << "\t" << pressure[i] 
              << "\t" << temperature_gas[i] 
              << "\t" << temperature_dust[i] << "\t" 
              << velocity[i] << "\n";
  }

  std::cout << "\n";

  nb_grid_points = radius_grid.size();

  number_densities.assign(nb_grid_points, std::vector<double>(nb_chemistry_species, 0.));
  mean_molecuar_weight.assign(nb_grid_points, 0.);
  total_element_density.assign(nb_grid_points, 0.);
  total_h_density.assign(nb_grid_points, 0.);

  absorption_coeff.resize(nb_grid_points);
  scattering_coeff.resize(nb_grid_points);

  absorption_coeff_gas.resize(nb_grid_points);
  scattering_coeff_gas.resize(nb_grid_points);
  absorption_coeff_dust.resize(nb_grid_points);
  scattering_coeff_dust.resize(nb_grid_points);
}



std::vector<double> Atmosphere::cgsToBar(const std::vector<double> pressure_data)
{
  std::vector<double> pressure_data_bar = pressure_data;

  for (auto & p : pressure_data_bar)
    p *= 1e-6;

  return pressure_data_bar;
}


void Atmosphere::equationOfState()
{
  for (size_t i=0; i<nb_grid_points; ++i)
    pressure[i] = mass_density[i] / (mean_molecuar_weight[i] * constants::mass_proton) * constants::boltzmann_k * temperature_gas[i];

  pressure_bar = cgsToBar(pressure);
}


}