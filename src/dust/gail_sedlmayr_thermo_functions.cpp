
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

//ln(!) of saturation vapour pressure of graphite
//p_vap is in units of dyn cm-2
double GailSedlmayrDust::saturationVapourPressure(const double temperature)
{

  return (18.6516 - 85906.1/temperature) + std::log(1e6);
  // const double theta = 5040./temperature;

  // return  3.24595e+1 - 1.68624e+1 * theta - 5.17205e-2 * theta*theta 
  //         + 3.99686e-3 * std::pow(theta, 3)
  //         -1.00638e-4 * std::pow(theta, 4);
}



//returns the ln of the saturation ratio of graphite
double GailSedlmayrDust::saturationRatio(
  const double temperature,
  const double number_density_carbon)
{
  const double ln_p_vap = saturationVapourPressure(temperature);

  const double partial_pressure = number_density_carbon 
                                * constants::boltzmann_k * temperature;

  double saturation_ratio = partial_pressure / std::exp(ln_p_vap);

  if (saturation_ratio >= 1)
    return std::log(saturation_ratio);
  else
    return -9999.;
}



double GailSedlmayrDust::criticalClusterSize(
  const double temperature,
  const double ln_saturation_ratio)
{
  if (ln_saturation_ratio <= 0)
    return 0.;

  const double n_star_inf = std::pow(
    (2. * theta_infinity) / (3. * temperature * ln_saturation_ratio), 
    3);

  const double part1 = std::sqrt(1 + 2.* std::pow(n_l/n_star_inf, 1./3.));
  const double part2 = 2.*std::pow(n_l/n_star_inf, 1./3.);

  const double n_star = 1 + n_star_inf/8. * std::pow(1 + part1 - part2, 3);

  if (n_star - 1 < 1.e-10)
    return n_star + 1.e-10;

  return n_star;
}



double GailSedlmayrDust::freeEnergyOfFormation(
  const double temperature,
  const double ln_saturation_ratio,
  const double critical_cluster_size)
{
  if (ln_saturation_ratio <= 0)
    return 0.;

  const double theta_n = theta_infinity 
                         / (1. + std::pow(n_l/(critical_cluster_size - 1.), 1./3.));

  const double delta_f = - constants::boltzmann_k * temperature * (critical_cluster_size - 1.)
                         * ln_saturation_ratio
                         + constants::boltzmann_k * theta_n 
                         * std::pow(critical_cluster_size - 1., 2./3.);

  return delta_f;
}


//beta factor from G&S
double GailSedlmayrDust::monomerGrowthRate(
  const double temperature,
  const double number_density_c,
  const double number_density_c2,
  const double number_density_c2h,
  const double number_density_c2h2)
{
  
  double beta = 1./std::sqrt(mass_c) * sticking_coeff[0] * number_density_c
              + 2./std::sqrt(mass_c2) * std::pow(2., 2./3.) * sticking_coeff[1] * number_density_c2
              + 2./std::sqrt(mass_c2h) * std::pow(2., 2./3.) * sticking_coeff[1] * number_density_c2h
              + 2./std::sqrt(mass_c2h2) * std::pow(2., 2./3.) * sticking_coeff[1] * number_density_c2h2;

  beta *= std::sqrt(constants::boltzmann_k * temperature/(2. * constants::pi));
  
  return beta;
}



double GailSedlmayrDust::zeldovichFactor(
  const double temperature,
  const double critical_cluster_size,
  const double ln_saturation_ratio)
{
  if (ln_saturation_ratio <= 0)
    return 0.;

  double z = 2./9. * theta_infinity / temperature * std::pow(critical_cluster_size - 1, 2./3.)
            * (1. / std::pow( std::pow(n_l, 1./3.) + std::pow(critical_cluster_size - 1, 1./3.), 2)
            + std::pow(n_l, 1./3.)
            * 1./std::pow(std::pow(n_l, 1./3.) + std::pow(critical_cluster_size - 1, 1./3.), 3));
  
  z /= constants::pi * 2;

  z = std::sqrt(z);

  return z;
}


//size distribution of the critical cluster in equilibrium
double GailSedlmayrDust::equilibriumClusterDistribution(
  const double temperature,
  const double free_energy,
  const double monomer_number_density,
  const double ln_saturation_ratio)
{
  if (ln_saturation_ratio <= 0)
    return 0;

  double c = monomer_number_density 
           * std::exp(-free_energy / (constants::boltzmann_k*temperature));
  
  return c;
}


}