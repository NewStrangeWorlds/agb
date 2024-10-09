
#include "analytic_dust.h" 

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
#include "../additional/exceptions.h"
#include "../atmosphere/atmosphere.h"
#include "../additional/quadrature.h"
#include "../additional/physical_const.h"
#include "../../_deps/lx_mie-src/mie/mie.h"



namespace agb{


AnalyticDust::AnalyticDust(
  ModelConfig* config_,
  SpectralGrid* spectral_grid_,
  Atmosphere* atmosphere_,
  const double particle_radius)
  : DustSpecies(config_, spectral_grid_, atmosphere_)
  , const_radius(particle_radius)
{


}



void AnalyticDust::calcDistribution()
{
  const double condensation_temperature = 1100.;
  const double max_number_density = 5.0e-13;
  
  unsigned int condensation_radius_idx = 0;
  double condensation_radius = 0;

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    if (atmosphere->temperature_gas[i] < condensation_temperature)
    {
      condensation_radius_idx = i;
      condensation_radius = atmosphere->radius_grid[i];

      break;
    }
  }

  number_density.assign(nb_grid_points, 0.);
  
  //analytical fit to Winters standard model
  for (size_t i=0; i<nb_grid_points; ++i)
  {
    number_density[i] = max_number_density 
                       / (std::exp((condensation_radius - atmosphere->radius_grid[i])/0.06) + 1.);
    number_density[i] *= atmosphere->total_h_density[i];
  }

  //don't let the dust density increase towards smaller r/R*
  for (int i=condensation_radius_idx; i>-1; i--)
  {
    if (number_density[i] > number_density[i+1])
      number_density[i] = number_density[i+1];
  }

  //fixed, mono-dispersed size distribution throughout the wind
  size_distribution.assign(nb_grid_points, std::vector<double>(1, 1.0));
  particle_radius.assign(nb_grid_points, std::vector<double>(1, const_radius * 1e-4)); //in cm
}


void AnalyticDust::saveOutput(const std::string file_path)
{

  //nothing to do here

}


std::vector<double> AnalyticDust::degreeOfCondensation(
  const double element_abundance)
{
  return std::vector<double>(atmosphere->nb_grid_points, 0);
}


}