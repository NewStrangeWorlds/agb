
#include "temperature_correction.h" 

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
#include "../radiative_transfer/radiative_transfer.h"
#include "../config/config.h"
#include "../atmosphere/atmosphere.h"
#include "../additional/quadrature.h"
#include "../additional/physical_const.h"

namespace agb{


std::vector<double> TemperatureCorrection::calculate(
  std::vector<double>& temperature,
  std::vector<double>& radius,
  std::vector<std::vector<double>>& extinction_coeff,
  std::vector<std::vector<double>>& absorption_coeff)
{
  std::vector<double> fq_j(temperature.size(), 0.);
  std::vector<double> qx_h(temperature.size(), 0.);
  std::vector<double> kappa_j(temperature.size(), 0.);
  std::vector<double> kappa_b(temperature.size(), 0.);
  std::vector<double> kappa_h(temperature.size(), 0.);

  calcIntegratedQuantities(
    radius,
    temperature,
    extinction_coeff,
    absorption_coeff,
    fq_j,
    qx_h,
    kappa_j,
    kappa_b,
    kappa_h);


  std::vector<double> delta_temperature = unsoeldLucyCorrection(
    temperature,
    radius,
    fq_j,
    qx_h,
    kappa_j,
    kappa_b,
    kappa_h);

  // std::vector<double> delta_temperature2 = lambdaIteration(
  //   temperature,
  //   absorption_coeff);

  for (size_t i=0; i<radius.size(); ++i)
  { 
    if (std::abs(delta_temperature[i])/temperature[i] > config->temperature_max_change) 
    {
      if (delta_temperature[i] < 0)
        delta_temperature[i] = -config->temperature_max_change * temperature[i];
      else
        delta_temperature[i] = config->temperature_max_change * temperature[i];
    }
  }

  return delta_temperature;
}


void TemperatureCorrection::calcIntegratedQuantities(
  const std::vector<double>& radius,
  const std::vector<double>& temperature,
  const std::vector<std::vector<double>>& extinction_coeff,
  const std::vector<std::vector<double>>& absorption_coeff,
  std::vector<double>& fq_j,
  std::vector<double>& qx_h,
  std::vector<double>& kappa_j,
  std::vector<double>& kappa_b,
  std::vector<double>& kappa_h)
{

  for (size_t i=0; i<radiation_field.size(); ++i)
  {
    std::vector<double> y(spectral_grid->nbSpectralPoints(), 0.);

    for (size_t j=0; j<spectral_grid->nbSpectralPoints(); ++j)
      y[j] = radiation_field[i].eddington_factor[j] * radiation_field[i].sphericality_factor[j] * radiation_field[i].mean_intensity[j];

    fq_j[i] = radiation_field[i].wavelengthIntegration(y) / radiation_field[i].mean_intensity_int;

    for (size_t j=0; j<spectral_grid->nbSpectralPoints(); ++j)
      y[j] = radiation_field[i].sphericality_factor[j] * extinction_coeff[i][j] * radiation_field[i].eddington_flux[j];

    qx_h[i] = radiation_field[i].wavelengthIntegration(y) / radiation_field[i].eddington_flux_int;


    for (size_t j=0; j<spectral_grid->nbSpectralPoints(); ++j)
      y[j] = absorption_coeff[i][j] * radiation_field[i].mean_intensity[j];

    kappa_j[i] = radiation_field[i].wavelengthIntegration(y) / radiation_field[i].mean_intensity_int;

    for (size_t j=0; j<spectral_grid->nbSpectralPoints(); ++j)
      y[j] = absorption_coeff[i][j] * aux::planckFunctionWavelength(temperature[i], spectral_grid->wavelength_list[j]);

    kappa_b[i] = radiation_field[i].wavelengthIntegration(y) / aux::planckFunctionIntegrated(temperature[i]);

    for (size_t j=0; j<spectral_grid->nbSpectralPoints(); ++j)
      y[j] = extinction_coeff[i][j] * radiation_field[i].eddington_flux[j];

    kappa_h[i] = radiation_field[i].wavelengthIntegration(y) / radiation_field[i].eddington_flux_int;
  }

}



std::vector<double> TemperatureCorrection::unsoeldLucyCorrection(
  const std::vector<double>& temperature,
  const std::vector<double>& radius,
  const std::vector<double>& fq_j,
  const std::vector<double>& qx_h,
  const std::vector<double>& kappa_j,
  const std::vector<double>& kappa_b,
  const std::vector<double>& kappa_h)
{
  const size_t nb_grid_points = radius.size();

  std::vector<double> tau_h(nb_grid_points, 0);

  for (int i=nb_grid_points-2; i>-1; --i)
    tau_h[i] = tau_h[i+1] + 0.5*(radius[i+1] - radius[i])*(qx_h[i] + qx_h[i+1]);
  

  std::vector<double> delta_h(nb_grid_points, 0.);

  for (size_t i=0; i<nb_grid_points; ++i)
    delta_h[i] = (radiation_field[i].eddington_flux_int * radius[i]*radius[i] - config->stellar_luminosity / (16. * constants::pi * constants::pi ));

  std::vector<double> integral(nb_grid_points, 0.);

  integral.back() = fq_j.back() * delta_h.back() /(radiation_field.back().eddington_flux_int/radiation_field.back().mean_intensity_int);

  for (int i=nb_grid_points-2; i>-1; --i)
  {
    double qi = 0.5 * (tau_h[i] - tau_h[i+1]) 
              * (delta_h[i] + delta_h[i+1]);
    integral[i] = integral[i+1] + qi;
  }


  std::vector<double> delta_temperature(nb_grid_points, 0.);

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    double delta_b = (kappa_j[i]/kappa_b[i] * radiation_field[i].mean_intensity_int 
                    - aux::planckFunctionIntegrated(temperature[i]))
                    - kappa_j[i]/kappa_b[i]/(fq_j[i] * radius[i]*radius[i]) * integral[i];

    delta_temperature[i] = constants::pi / (4*constants::stefan_boltzmann * std::pow(temperature[i], 3)) * delta_b;
    //delta_temperature[i] = std::pow(std::pow(temperature[i], 4) + constants::pi/constants::stefan_boltzmann * delta_b, 0.25) - temperature[i];

    //std::cout << "delta " << i << "\t" << delta_h[i] << "\t" << tau_b[i] << "\t" << integral[i] << "\t" << kappa_j[i]/kappa_b[i] * radiation_field[i].mean_intensity_int - aux::planckFunctionIntegrated(temperature[i]) << "\t" << delta_b << "\t" << delta_temperature[i] << "\n";
  }

  return delta_temperature;
}




std::vector<double> TemperatureCorrection::lambdaIteration(
  const std::vector<double>& temperature,
  const std::vector<std::vector<double>>& absorption_coeff)
{
  std::vector<double> delta_temperature(temperature.size(), 0.);
  
  for (size_t i=0; i<temperature.size(); ++i)
  {
    std::vector<double> y(spectral_grid->nbSpectralPoints(), 0.);

    for (size_t j=0; j<spectral_grid->nbSpectralPoints(); ++j)
      y[j] = absorption_coeff[i][j] * (radiation_field[i].mean_intensity[j] - aux::planckFunctionWavelength(temperature[i], spectral_grid->wavelength_list[i]));

    delta_temperature[i] = radiation_field[i].wavelengthIntegration(y);

    for (size_t j=0; j<spectral_grid->nbSpectralPoints(); ++j)
      y[j] = absorption_coeff[i][j] * aux::planckFunctionDerivWavelength(temperature[i], spectral_grid->wavelength_list[i]);

    delta_temperature[i] /= radiation_field[i].wavelengthIntegration(y);
  }
  
  return delta_temperature;
}



}