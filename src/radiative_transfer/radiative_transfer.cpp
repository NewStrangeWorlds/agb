
#include "radiative_transfer.h" 

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <omp.h>
#include <sstream>
#include <algorithm> 
#include <assert.h> 
#include <cmath>

#include "../additional/aux_functions.h"
#include "../spectral_grid/spectral_grid.h"
#include "../config/config.h"
#include "../atmosphere/atmosphere.h"

namespace agb{


RadiativeTransfer::RadiativeTransfer(
  ModelConfig* config_,
  SpectralGrid* spectral_grid_,
  Atmosphere* atmosphere_)
  : config(config_)
  , spectral_grid(spectral_grid_)
  , atmosphere(atmosphere_)
  , nb_core_impact_param(config->nb_core_impact_param)
  , nb_impact_param(nb_core_impact_param+atmosphere->nb_grid_points)
  , nb_grid_points(atmosphere->nb_grid_points)
  , nb_spectral_points(spectral_grid->nbSpectralPoints())
  , radiation_field(nb_grid_points, RadiationField(nb_spectral_points))
{
  
  eddington_factors.assign(nb_grid_points, std::vector<double>(nb_spectral_points, 0.0));

  extinction_coeff.assign(nb_spectral_points, std::vector<double>(nb_grid_points, 0.0));

  source_function.assign(nb_spectral_points, std::vector<double>(nb_grid_points, 0.0));

  createImpactParameterGrid();

  for (size_t i=0; i<nb_grid_points; ++i)
    radiation_field[i].createAngleGrid(impact_parameter_grid, atmosphere->radius[i], i);
}


void RadiativeTransfer::calcSourceFunction()
{
  for (size_t i=0; i<nb_grid_points; ++i)
    for (size_t j=0; j<nb_spectral_points; ++j)
      source_function[j][i] = atmosphere->absorption_coeff[i][j] 
        * aux::planckFunctionWavelength(
            atmosphere->temperature_gas[i], 
            spectral_grid->wavelength_list[j])
        + atmosphere->scattering_coeff[i][j]/extinction_coeff[j][i] 
        * radiation_field[i].mean_intensity[j];
}



void RadiativeTransfer::solveRadiativeTransfer()
{
  for (size_t i=0; i<nb_grid_points; ++i)
    for (size_t j=0; j<nb_spectral_points; ++j)
      extinction_coeff[j][i] = atmosphere->absorption_coeff[i][j] + atmosphere->scattering_coeff[i][j];

  calcSourceFunction();

  impact_parameter_grid[0].solveRadiativeTransfer(200, extinction_coeff[200], source_function[200]);

  exit(0);
  
  for (size_t i=0; i<nb_spectral_points; ++i)
    for (auto & ip : impact_parameter_grid)
      ip.solveRadiativeTransfer(
        i, 
        extinction_coeff[i],
        source_function[i]);
}



}