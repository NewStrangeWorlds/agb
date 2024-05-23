
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

  theta_infinity = surface_tension * 4. * constants::pi 
                 * monomer_radius*monomer_radius / constants::boltzmann_k;

}



double GailSedlmayrDust::nucleationRate(
  const double temperature,
  const double number_density_c,
  const double number_density_c2,
  const double number_density_c2h,
  const double number_density_c2h2)
{
  
  const double ln_saturation_ratio = saturationRatio(
    temperature,
    number_density_c);

  if (ln_saturation_ratio <= 0)
    return 1e-100;


  const double n_star = criticalClusterSize(
    temperature,
    ln_saturation_ratio);

  const double delta_f = freeEnergyOfFormation(
    temperature,
    ln_saturation_ratio,
    n_star);

  const double beta = monomerGrowthRate(
    temperature,
    number_density_c,
    number_density_c2,
    number_density_c2h,
    number_density_c2h2);

  const double z = zeldovichFactor(
    temperature,
    n_star,
    ln_saturation_ratio);

  const double c0 = equilibriumClusterDistribution(
    temperature,
    delta_f,
    number_density_c,
    ln_saturation_ratio);

  const double a_star = monomer_surface_area * std::pow(n_star, 2./3.);
  const double nucleation_rate = beta * a_star * z * c0;

  return nucleation_rate;
}



double GailSedlmayrDust::growthRate(
  const double temperature,
  const double number_density_c,
  const double number_density_c2,
  const double number_density_c2h,
  const double number_density_c2h2)
{
  double tau = monomer_surface_area * std::sqrt(constants::boltzmann_k * temperature)
              * ( sticking_coeff[0]/std::sqrt(mass_c) * number_density_c
              + 2 * sticking_coeff[1]/std::sqrt(mass_c2) * number_density_c2
              + 2 * sticking_coeff[1]/std::sqrt(mass_c2h) * number_density_c2h
              + 2 * sticking_coeff[1]/std::sqrt(mass_c2h2) * number_density_c2h2);

  return 1./tau;
}


void GailSedlmayrDust::calcDistribution()
{
  for (size_t i=0; i<nb_grid_points; ++i)
  {
    
  }

  //fixed, mono-dispersed size distribution throughout the wind
  //size_distribution.assign(nb_grid_points, std::vector<double>(1, 1.0));
  //particle_radius.assign(nb_grid_points, std::vector<double>(1, const_radius * 1e-4)); //in cm
}


}