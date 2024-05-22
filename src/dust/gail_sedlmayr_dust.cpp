
#include "gail_sedlmayr_dust.h" 

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



namespace agb{


GailSedlmayrDust::GailSedlmayrDust(
  ModelConfig* config_,
  SpectralGrid* spectral_grid_,
  Atmosphere* atmosphere_)
  : DustSpecies(config_, spectral_grid_, atmosphere_)
{


}


//log(!) of saturation vapour pressure of graphite
//p_vap is in units of dyn cm-2
double GailSedlmayrDust::saturationVapourPressure(const double temperature)
{

  return (18.6516 - 85906.1/temperature) + log(1e6);

}



//returns the logarithm of the saturation ratio of graphite
double GailSedlmayrDust::saturationRatio(
  const double temperature,
  const double number_density_carbon)
{
  const double p_vap = saturationVapourPressure(temperature);

  const double partial_pressure = number_density_carbon 
                                * constants::boltzmann_k * temperature;

  double saturation_ratio = partial_pressure / std::exp(p_vap);

  if (saturation_ratio > 1)
    return std::log(saturation_ratio);
  else
    return -9999.;
}



void GailSedlmayrDust::calcDistribution()
{
  

  //fixed, mono-dispersed size distribution throughout the wind
  //size_distribution.assign(nb_grid_points, std::vector<double>(1, 1.0));
  //particle_radius.assign(nb_grid_points, std::vector<double>(1, const_radius * 1e-4)); //in cm
}


}