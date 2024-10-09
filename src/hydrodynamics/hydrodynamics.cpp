
#include "hydrodynamics.h"

#include "../additional/aux_functions.h"
#include "../spectral_grid/spectral_grid.h"
#include "../config/config.h"
#include "../additional/exceptions.h"
#include "../atmosphere/atmosphere.h"
#include "../additional/quadrature.h"
#include "../additional/physical_const.h"
#include "../radiative_transfer/radiative_transfer.h"

#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h> 
#include <iostream>
#include <fstream>
#include <iomanip>


namespace agb {


Hydrodynamics::Hydrodynamics(
  ModelConfig* config_,
  SpectralGrid* spectral_grid_,
  Atmosphere* atmosphere_,
  std::vector<RadiationField>& radiation_field_)
  : config(config_)
  , spectral_grid(spectral_grid_)
  , atmosphere(atmosphere_)
  , radiation_field(radiation_field_)
  , nb_grid_points(atmosphere->nb_grid_points)
{
  isothermal_sound_speed.assign(nb_grid_points, 0);
  alpha.assign(nb_grid_points, 0);
  sound_speed_derivative.assign(nb_grid_points, 0);
  phi.assign(nb_grid_points, 0);
}


void Hydrodynamics::saveOutput(const std::string file_path)
{
  std::fstream file;
  file.open(file_path.c_str(), std::ios::out);

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    file << std::setprecision(10) << std::scientific 
         << atmosphere->radius_grid[i] << "\t"
         << atmosphere->velocity[i] << "\t" 
         << isothermal_sound_speed[i] << "\t" 
         << alpha[i] << "\t" 
         << phi[i] << "\n";
  }

  file.close();
}



void Hydrodynamics::calcWindVelocity()
{
  isothermal_sound_speed.assign(nb_grid_points, 0);
  alpha.assign(nb_grid_points, 0);
  sound_speed_derivative.assign(nb_grid_points, 0);
  phi.assign(nb_grid_points, 0);

  speedOfSound();

  std::vector<double> flux_weighted_extinction(nb_grid_points, 0);

  for (size_t i=0; i<flux_weighted_extinction.size(); ++i)
    flux_weighted_extinction[i] = radiation_field[i].fluxWeightedExtinction(atmosphere->extinction_coeff[i]);
  
  calcAlpha(flux_weighted_extinction);
  
  sound_speed_derivative = soundSpeedDerivative();
  int critical_point = findCriticalPoint();

  phi = calcPhi(critical_point);

  std::vector<double> wind_velocity = windVelocity(phi, critical_point);

  mass_loss_rate = calcMassLossRate(
    critical_point, 
    sound_speed_derivative[critical_point], 
    flux_weighted_extinction[critical_point]) / constants::mass_sun * constants::year;


  const double max_change = 1e-2;

  atmosphere->velocity = wind_velocity;

  std::vector<double> mass_density = massDensity(mass_loss_rate);

  // for (size_t i=0; i<mass_density.size(); ++i)
  // { 
  //   double delta_rho = atmosphere->mass_density[i] - mass_density[i];

  //   if (std::abs(delta_rho)/atmosphere->mass_density[i] > max_change) 
  //   {
  //     if (delta_rho < 0)
  //       mass_density[i] = atmosphere->mass_density[i] * (1 + max_change);
  //     else
  //       mass_density[i] = atmosphere->mass_density[i] * (1 - max_change);
  //   }
  // }


  //atmosphere->mass_density = mass_density;

  for (size_t i=0; i<nb_grid_points; ++i)
    atmosphere->mass_density[i] = std::sqrt(mass_density[i] * atmosphere->mass_density[i]);

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    std::cout << i << "  " 
              << atmosphere->velocity[i] << "\t" 
              << isothermal_sound_speed[i] << "\t" 
              << alpha[i] << "\t" 
              << phi[i] << "\t" 
              << wind_velocity[i] << "\t"
              << mass_density[i] << "\t"
              << atmosphere->mass_density[i] << "\n";
  }
  
  std::cout << "Critical point: " << critical_point << "\n";
  std::cout << "Mass loss rate: " << mass_loss_rate << "\n";
}



void Hydrodynamics::speedOfSound()
{
  //for (size_t i=0; i<nb_grid_points; ++i)
    //isothermal_sound_speed[i] = sqrt(atmosphere->pressure[i]/atmosphere->mass_density[i]);

  for (size_t i=0; i<nb_grid_points; ++i)
    isothermal_sound_speed[i] = std::sqrt(1.0 / (atmosphere->mean_molecuar_weight[i] * constants::mass_proton) 
      * constants::boltzmann_k * atmosphere->temperature_gas[i]);
}


void Hydrodynamics::calcAlpha(
  const std::vector<double>& flux_weighted_extinction)
{
  for (size_t i=0; i<nb_grid_points; ++i)
  {
    alpha[i] = flux_weighted_extinction[i]/atmosphere->mass_density[i] 
             * config->stellar_luminosity
             /(4*constants::pi * constants::light_c * constants::gravitation_const * config->stellar_mass);
  }
}



std::vector<double> Hydrodynamics::windVelocity(
  const std::vector<double>& phi,
  const int critical_point)
{
  std::vector<double> wind_velocity(nb_grid_points, 0);

  for (size_t i=0; i<critical_point; ++i)
    wind_velocity[i] = phi[i] - std::sqrt(phi[i]*phi[i] - isothermal_sound_speed[i] * isothermal_sound_speed[i]);
  
  for (size_t i=critical_point+1; i<nb_grid_points; ++i)
    wind_velocity[i] = phi[i] + std::sqrt(phi[i]*phi[i] - isothermal_sound_speed[i] * isothermal_sound_speed[i]);
  
  //to make the transition between the upward and downward branch smooth
  //we interpolate the wind velocity at the critical point between the two end points
  const double x_left = atmosphere->radius[critical_point-1];
  const double x_right = atmosphere->radius[critical_point+1];
  
  double y_left = wind_velocity[critical_point-1];
  double y_right = wind_velocity[critical_point+1];
  wind_velocity[critical_point] = y_left + (atmosphere->radius[critical_point] - x_left) * (y_right - y_left) / (x_right - x_left);

  return wind_velocity;
}


//integration of phi from the inner boundary to the critical point
//uses an explicit Euler scheme
void Hydrodynamics::integratePhiOutward(
  const int critical_point,
  const double start_velocity,
  std::vector<double>& phi)
{ 
  double v = start_velocity;
  phi[0] = 0.5 * (v + isothermal_sound_speed[0]*isothermal_sound_speed[0]/v);

  for (size_t i=1; i<critical_point+1; ++i)
  {
    const double r = atmosphere->radius[i-1];

    double phi_deriv = - 0.5/v 
             * constants::gravitation_const * config->stellar_mass 
             /r/r * (1 - alpha[i-1])
             + isothermal_sound_speed[i-1] * isothermal_sound_speed[i-1] 
             / v/r;

    phi[i] = phi[i-1] + (atmosphere->radius[i] - r) * phi_deriv;
    
    v = phi[i] - std::sqrt(phi[i]*phi[i] - isothermal_sound_speed[i]*isothermal_sound_speed[i]);
  }
 
}


//integration of phi from the outer boundary inwards to the critical point
//uses an explicit Euler scheme
void Hydrodynamics::integratePhiInward(
  const int critical_point,
  const double end_velocity,
  std::vector<double>& phi)
{
  double v = end_velocity;
  
  phi.back() = 0.5 * (v + isothermal_sound_speed.back()*isothermal_sound_speed.back() / v);

  for (int i=phi.size()-2; i>critical_point-1; --i)
  {
    const double r = atmosphere->radius[i+1];

    double phi_deriv = - 0.5/v 
             * constants::gravitation_const * config->stellar_mass 
             /r/r * (1 - alpha[i+1])
             + isothermal_sound_speed[i+1] * isothermal_sound_speed[i+1] 
             / v/r;

    phi[i] = phi[i+1] - (r - atmosphere->radius[i]) * phi_deriv;
    
    v = phi[i] + std::sqrt(phi[i]*phi[i] - isothermal_sound_speed[i]*isothermal_sound_speed[i]);
  }

}


std::vector<double> Hydrodynamics::calcPhi(const int critical_point)
{
  std::vector<double> phi(nb_grid_points, 0.0);

  //first we solve phi from the inner boundary to the critical point
  //we use a shooting method to find the correct velocity at inner boundary
  double v_0 = 1e-2;
  double v_1 = isothermal_sound_speed[critical_point];

  for (int it=0; it<1000; ++it)
  {
    double v_start = (v_1 - v_0) * 0.5 + v_0;
    
    integratePhiOutward(critical_point, v_start, phi);

    if (std::isnan(phi[critical_point]))
      v_1 = v_start;
    else
    {
      if (phi[critical_point] > isothermal_sound_speed[critical_point])
        v_0 = v_start;
      else 
        v_1 = v_start;
    }

    if (std::abs(v_1 - v_0)/v_1 < 1e-5 && !std::isnan(phi[critical_point]))
      break;
  }


  //now we solve phi from the outer boundary to the critical point
  //we use a shooting method to find the correct velocity at outer boundary
  v_0 = isothermal_sound_speed[critical_point];
  v_1 = isothermal_sound_speed[critical_point] * 100;

  for (int it=0; it<1000; ++it)
  {
    double v_start = (v_1 - v_0) * 0.5 + v_0;
    
    integratePhiInward(critical_point, v_start, phi);

    if (std::isnan(phi[critical_point]))
      v_0 = v_start;
    else
    {
       if (phi[critical_point] < isothermal_sound_speed[critical_point])
         v_0 = v_start;
       else
        v_1 = v_start;
    }

    if (std::abs(v_1 - v_0)/v_1 < 1e-5  && !std::isnan(phi[critical_point]))
      break;
  }
  
  //set phi at the critical point to its actual value
  phi[critical_point] = isothermal_sound_speed[critical_point];

  return phi;
}



int Hydrodynamics::findCriticalPoint()
{
  std::vector<double> phi_deriv(nb_grid_points, 0);

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    phi_deriv[i] = constants::gravitation_const * config->stellar_mass 
                                / (atmosphere->radius[i]*atmosphere->radius[i])
                                * (1. - alpha[i])
                                - 2 * isothermal_sound_speed[i] * isothermal_sound_speed[i]/atmosphere->radius[i];
  }

  int i = 0;

  while (i < nb_grid_points-1 && phi_deriv[i] > 0)
     ++i;

  return i;
}



std::vector<double> Hydrodynamics::soundSpeedDerivative()
{
  std::vector<double> sound_speed_derivative(nb_grid_points, 0);

  for (size_t i=0; i<nb_grid_points-1; ++i)
  {
    sound_speed_derivative[i] = (isothermal_sound_speed[i+1]*isothermal_sound_speed[i+1] 
                              - isothermal_sound_speed[i]*isothermal_sound_speed[i])
                              / (atmosphere->radius[i+1] - atmosphere->radius[i]);
  }

  sound_speed_derivative[nb_grid_points-1] = sound_speed_derivative[nb_grid_points-2];

  return sound_speed_derivative;
}


double Hydrodynamics::calcMassLossRate(
  const int critical_point,
  const double sound_speed_derivative,
  const double flux_weighted_extinction)
{
  const double radius = atmosphere->radius[critical_point];
  const double sound_speed = isothermal_sound_speed[critical_point];
 
  double mass_loss_rate = config->stellar_luminosity * flux_weighted_extinction * sound_speed / constants::light_c
                          * 1./(constants::gravitation_const * config->stellar_mass/radius/radius 
                                - 2*sound_speed*sound_speed/radius - sound_speed_derivative);

  return mass_loss_rate;
}


std::vector<double> Hydrodynamics::massDensity(
  const double mass_loss_rate1)
{
  std::vector<double> mass_density(nb_grid_points, 0);

  const double mass_loss_rate = config->stellar_mass_loss_rate;
  
  for (size_t i=0; i<nb_grid_points; ++i)
    mass_density[i] = mass_loss_rate 
      / (4*constants::pi*atmosphere->radius[i]*atmosphere->radius[i] * atmosphere->velocity[i]);

  return mass_density;
}

}