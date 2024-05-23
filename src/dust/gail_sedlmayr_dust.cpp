
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
#include "../chemistry/chem_species.h"



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
  double nucleation_rate = beta * a_star * z * c0;

  if (nucleation_rate < 1e-100)
    nucleation_rate = 1e-100;

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


std::vector<double> GailSedlmayrDust::dustMomentZero()
{
  std::vector<double> K0(nb_grid_points, 0);

  auto linearInterpolation = [&](
    const double x_l, 
    const double x_r, 
    const double y_l, 
    const double y_r, 
    const double x) {
      return y_l + (y_r - y_l)/(x_r - x_l) * (x - x_l);};

  for (size_t i=1; i<nb_grid_points; ++i)
  {
    const double v_l = atmosphere->velocity[i-1];
    const double v_r = atmosphere->velocity[i];

    const double n_l = nucleation_rate[i-1] / atmosphere->total_h_density[i-1];
    const double n_r = nucleation_rate[i] / atmosphere->total_h_density[i];

    const double h = atmosphere->radius[i] - atmosphere->radius[i-1];

    const double j = std::pow(
      10, 
      linearInterpolation(
        atmosphere->radius[i-1], 
        atmosphere->radius[i], 
        std::log10(n_l), 
        std::log10(n_r), 
        atmosphere->radius[i-1] + h/2.));

    const double v = linearInterpolation(
      atmosphere->radius[i-1], 
      atmosphere->radius[i], 
      v_l, 
      v_r, 
      atmosphere->radius[i-1] + h/2.);
    
    const double rk1 = h * n_l / v_l;
    const double rk2 = h * j / v;
    const double rk3 = rk2;
    const double rk4 = h * n_r / v_r;

    K0[i] = K0[i-1] + 1./6.*rk1 + 1./3.*rk2 + 1./3. *rk3 + 1./6.*rk4;
  }

  //for (auto & k : K0)
    //if (k < 1e-45) k = 1e-45;

  return K0;
}



std::vector<double> GailSedlmayrDust::dustMoment(
  const int order,
  std::vector<double>& prev_moment,
  const int n_lower)
{
  std::vector<double> K(nb_grid_points, 0);

  auto linearInterpolation = [&](
    const double x_l, 
    const double x_r, 
    const double y_l, 
    const double y_r, 
    const double x) {
      return y_l + (y_r - y_l)/(x_r - x_l) * (x - x_l);};

  K[0] = 0;

  for (size_t i=1; i<nb_grid_points; ++i)
  {
    const double v_l = atmosphere->velocity[i-1];
    const double v_r = atmosphere->velocity[i];

    const double tau_l = growth_rate[i-1];
    const double tau_r = growth_rate[i];

    const double n_l = nucleation_rate[i-1]/atmosphere->total_h_density[i-1];
    const double n_r = nucleation_rate[i]/atmosphere->total_h_density[i];

    const double h = atmosphere->radius[i] - atmosphere->radius[i-1];

    const double bracket_l = 1./v_l * (order/3. / tau_l * prev_moment[i-1] + std::pow(n_lower,order/3.) * n_l);
    const double bracket_r = 1./v_r * (order/3. / tau_r * prev_moment[i] + std::pow(n_lower,order/3.) * n_r);

    const double bracket = linearInterpolation(
      atmosphere->radius[i-1], 
      atmosphere->radius[i], 
      bracket_l, 
      bracket_r, 
      atmosphere->radius[i-1] + h/2.);

    const double rk1 = h * bracket_l;
    const double rk2 = h * bracket;
    const double rk3 = rk2;
    const double rk4 = h * bracket_r;

    K[i] = K[i-1] + 1./6.*rk1 + 1./3.*rk2 + 1./3. *rk3 + 1./6.*rk4;
  }

  //for (auto & k : K)
    //if (k < 1e-45) k = 1e-45;

  return K;
}



void GailSedlmayrDust::calcDistribution()
{
  nucleation_rate.assign(nb_grid_points, 0);
  growth_rate.assign(nb_grid_points, 0);

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    const double number_density_c = atmosphere->number_densities[i][_C];
    const double number_density_c2 = atmosphere->number_densities[i][_C2];
    const double number_density_c2h = atmosphere->number_densities[i][_C2H];
    const double number_density_c2h2 = atmosphere->number_densities[i][_C2H2];
    const double temperature = atmosphere->temperature_gas[i];

    nucleation_rate[i] = nucleationRate(
      temperature,
      number_density_c,
      number_density_c2,
      number_density_c2h,
      number_density_c2h2);

    growth_rate[i] = growthRate(
      temperature,
      number_density_c,
      number_density_c2,
      number_density_c2h,
      number_density_c2h2);
  }


  dust_moments.assign(nb_moments, std::vector<double>(nb_grid_points, 0));

  dust_moments[0] = dustMomentZero();

  for (size_t i=1; i<nb_moments; ++i)
    dust_moments[i] = dustMoment(i, dust_moments[i-1], minimum_monomer_number);


  number_density.assign(nb_grid_points, 0);
  size_distribution.assign(nb_grid_points, std::vector<double>(1, 1.0));
  particle_radius.assign(nb_grid_points, std::vector<double>(1, 0)); //in cm

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    const double number_density_c = atmosphere->number_densities[i][_C];
    const double temperature = atmosphere->temperature_gas[i];

    const double ln_s = saturationRatio(temperature, number_density_c);

    if (ln_s < 0)
    {
      number_density[i] = 0.;
      particle_radius[i][0] = 0.;
    }
    else
    {
      number_density[i] = dust_moments[0][i] * atmosphere->total_h_density[i];
      particle_radius[i][0] = monomer_radius * dust_moments[1][i] / dust_moments[0][i];
    }


    std::cout << i << "\t" << ln_s << "\t" << number_density[i] << "\t" << particle_radius[i][0] 
              << "\t" << nucleation_rate[i] << "\t" << growth_rate[i] << "\t" << dust_moments[0][i] << "\t" << dust_moments[1][i] << "\t" << dust_moments[2][i] 
              //<< "\t" << constants::pi * monomer_radius*monomer_radius * dust_moments[2][i] * atmosphere->total_h_density[i] << "\t" << constants::pi * particle_radius[i][0]*particle_radius[i][0] * number_density[i]
              << "\n";
  }
}


}