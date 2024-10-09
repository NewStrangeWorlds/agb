
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
  double tau = monomer_surface_area * std::sqrt(constants::boltzmann_k * temperature/(2* constants::pi))
              * ( sticking_coeff[0]/std::sqrt(mass_c) * number_density_c
              + 2 * sticking_coeff[1]/std::sqrt(mass_c2) * number_density_c2
              + 2 * sticking_coeff[1]/std::sqrt(mass_c2h) * number_density_c2h
              + 2 * sticking_coeff[1]/std::sqrt(mass_c2h2) * number_density_c2h2);

  return 1./tau;
}


//calculation of K0 using a 4th-order Runge-Kutta method
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

    const double j_l = nucleation_rate[i-1] / atmosphere->total_h_density[i-1];
    const double j_r = nucleation_rate[i] / atmosphere->total_h_density[i];

    const double h = atmosphere->radius[i] - atmosphere->radius[i-1];

    const double j = std::pow(
      10, 
      linearInterpolation(
        atmosphere->radius[i-1], 
        atmosphere->radius[i], 
        std::log10(j_l), 
        std::log10(j_r), 
        atmosphere->radius[i-1] + h/2.));

    const double v = linearInterpolation(
      atmosphere->radius[i-1], 
      atmosphere->radius[i], 
      v_l, 
      v_r, 
      atmosphere->radius[i-1] + h/2.);
    
    const double rk1 = h * j_l / v_l;
    const double rk2 = h * j / v;
    const double rk3 = rk2;
    const double rk4 = h * j_r / v_r;

    K0[i] = K0[i-1] + 1./6.*rk1 + 1./3.*rk2 + 1./3. *rk3 + 1./6.*rk4;
  }

  return K0;
}


//calculation of dust moments using a 4th-order Runge-Kutta method
std::vector<double> GailSedlmayrDust::dustMoment(
  const int order,
  std::vector<double>& prev_moment,
  const int n_lower)
{
  //just in case someone tries to call it
  //with order 0 :-)
  if (order == 0)
    return dustMomentZero();

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

    const double j_l = nucleation_rate[i-1]/atmosphere->total_h_density[i-1];
    const double j_r = nucleation_rate[i]/atmosphere->total_h_density[i];

    const double h = atmosphere->radius[i] - atmosphere->radius[i-1];

    const double bracket_l = 1./v_l * (order/3. / tau_l * prev_moment[i-1] + std::pow(n_lower,order/3.) * j_l);
    const double bracket_r = 1./v_r * (order/3. / tau_r * prev_moment[i] + std::pow(n_lower,order/3.) * j_r);

    const double bracket = std::pow(
      10,
      linearInterpolation(
        atmosphere->radius[i-1], 
        atmosphere->radius[i], 
        std::log10(bracket_l), 
        std::log10(bracket_r), 
        atmosphere->radius[i-1] + h/2.));

    const double rk1 = h * bracket_l;
    const double rk2 = h * bracket;
    const double rk3 = rk2;
    const double rk4 = h * bracket_r;

    K[i] = K[i-1] + 1./6.*rk1 + 1./3.*rk2 + 1./3. *rk3 + 1./6.*rk4;
  }

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


  //the dust moments are normalised to n<H>
  //here, we scale them back to their normal values
  for (size_t i=0; i<nb_grid_points; ++i)
  {
    for (auto & k : dust_moments)
    {
      k[i] *= atmosphere->total_h_density[i];

      if (k[i] < 1e-45) k[i] = 1e-45;
    }
  }


  number_density.assign(nb_grid_points, 0);
  size_distribution.assign(nb_grid_points, std::vector<double>(1, 1.0));
  particle_radius.assign(nb_grid_points, std::vector<double>(1, 0)); //in cm


  std::cout << "Gail&Sedlmayr dust calculation\n";

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
      number_density[i] = dust_moments[0][i];
      
      //compute the particle radius from the mean volume
      particle_radius[i][0] = monomer_radius * std::pow(dust_moments[3][i] / dust_moments[0][i], 1./3.);
    }


    //set minumim number densities and particle radii
    //otherwise, radiative transfer and temperature calculations crash
    if (number_density[i] < 1e-30)
    {
      number_density[i] = 1e-30;
      particle_radius[i][0] = 1e-8;
    }


    // std::cout << i << "\t" 
    //           << atmosphere->velocity[i] << "\t"
    //           << ln_s << "\t" 
    //           << number_density[i] << "\t" 
    //           << number_density[i]/atmosphere->total_h_density[i] << "\t"
    //           << particle_radius[i][0]*1e4 << "\t" 
    //           << nucleation_rate[i] << "\t" 
    //           << growth_rate[i] << "\t" 
    //           << dust_moments[0][i] << "\t" 
    //           << dust_moments[1][i] << "\t" 
    //           << dust_moments[2][i] << "\t"
    //           << dust_moments[3][i] << "\t"
    //           << "\n";
  }

  std::cout << "\n";
}



std::vector<double> GailSedlmayrDust::degreeOfCondensation(
  const double carbon_abundance)
{
  std::vector<double> degree_of_condensation(atmosphere->nb_grid_points, 0);

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    //the element abundances from FastChem are given with respect to n<tot>, not n<H>
    degree_of_condensation[i] = 
      dust_moments[3][i] /(atmosphere->total_element_density[i] * (carbon_abundance));

    if (degree_of_condensation[i] >=1) degree_of_condensation[i] = 0.9999999;
  }

  return degree_of_condensation;
}


void GailSedlmayrDust::saveOutput(const std::string file_path)
{
  std::fstream file;
  file.open(file_path.c_str(), std::ios::out);

  if (file.fail()) 
  {
    std::cout << "Couldn't open dust output file " << file_path << "\n";
    return;
  }

  std::cout << "Saving dust output to " << file_path << "\n\n";

  file << std::setprecision(10) << std::scientific << "#r/R*\tn_<H>(cm-3)\tnumber_density(cm-3)\tradius(micron)\tnucleation_rate\tgrowth_rate\tK0\tK1\tK2\tK3\n";
  
  for (size_t i=0; i<nb_grid_points; ++i)
  {  
     const double number_density_c = atmosphere->number_densities[i][_C];
     const double number_density_c2 = atmosphere->number_densities[i][_C2];
     const double number_density_c2h = atmosphere->number_densities[i][_C2H];
     const double number_density_c2h2 = atmosphere->number_densities[i][_C2H2];
     const double temperature = atmosphere->temperature_gas[i];

     file << atmosphere->radius_grid[i] << "\t"
          << atmosphere->total_h_density[i] << "\t"
          << temperature << "\t"
          << number_density_c << "\t"
          << number_density_c2 << "\t"
          << number_density_c2h << "\t"
          << number_density_c2h2 << "\t"
          << number_density[i] << "\t"
          << particle_radius[i][0]*1e4 << "\t"
          << nucleation_rate[i] << "\t"
          << growth_rate[i] << "\t"
          << dust_moments[0][i] << "\t"
          << dust_moments[1][i] << "\t"
          << dust_moments[2][i] << "\t"
          << dust_moments[3][i] << "\t"
          << dust_moments[3][i]/dust_moments[0][i] * 4./3. * constants::pi * std::pow(monomer_radius, 3) << "\t"
          << "\n";
  }

  file.close();
}



}