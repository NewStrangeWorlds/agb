
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
  std::vector<std::vector<double>>& absorption_coeff,
  std::vector<double>& relaxation,
  std::vector<double>& prev_delta_b)
{
  const size_t nb_grid_points = temperature.size();

  //(re)initialise the persistent per-layer state on first use or grid changes
  if (relaxation.size() != nb_grid_points)
    relaxation.assign(nb_grid_points, config->temperature_relaxation_init);
  if (prev_delta_b.size() != nb_grid_points)
    prev_delta_b.assign(nb_grid_points, 0.);

  std::vector<double> fq_j(nb_grid_points, 0.);
  std::vector<double> qx_h(nb_grid_points, 0.);
  std::vector<double> kappa_j(nb_grid_points, 0.);
  std::vector<double> kappa_b(nb_grid_points, 0.);
  std::vector<double> planck_function_int(nb_grid_points, 0.);

  calcIntegratedQuantities(
    radius,
    temperature,
    extinction_coeff,
    absorption_coeff,
    fq_j,
    qx_h,
    kappa_j,
    kappa_b,
    planck_function_int);

  //Planck-integrated correction delta_B per layer
  std::vector<double> delta_b = unsoeldLucyCorrection(
    temperature,
    radius,
    fq_j,
    qx_h,
    kappa_j,
    kappa_b,
    planck_function_int);

  //smooth the correction (not the profile): grid-scale noise is removed but the
  //filter vanishes as delta_B -> 0, so it cannot diffuse the converged solution
  if (config->smooth_temperature_profile)
    smoothCorrection(delta_b);

  std::vector<double> delta_temperature(nb_grid_points, 0.);

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    //adapt the per-layer under-relaxation: shrink omega when the correction
    //flips sign (oscillation), grow it back when the sign is stable (creep)
    if (delta_b[i] * prev_delta_b[i] < 0.)
      relaxation[i] *= config->temperature_relaxation_down;
    else
      relaxation[i] *= config->temperature_relaxation_up;

    if (relaxation[i] < config->temperature_relaxation_min)
      relaxation[i] = config->temperature_relaxation_min;
    if (relaxation[i] > config->temperature_relaxation_max)
      relaxation[i] = config->temperature_relaxation_max;

    prev_delta_b[i] = delta_b[i];

    //exact Planck (T^4) update with the FULL correction (omega is applied below,
    //after the cap): T_new^4 = T^4 + (pi/sigma) * delta_B. This is self-limiting
    //for large corrections, unlike the linearised pi/(4 sigma T^3) form. Floor the
    //argument so T cannot collapse below T/2.
    const double t4 = std::pow(temperature[i], 4);
    double t4_new = t4 + constants::pi / constants::stefan_boltzmann * delta_b[i];

    const double t4_floor = 0.0625 * t4; //(0.5 T)^4
    if (t4_new < t4_floor) t4_new = t4_floor;

    double step = std::pow(t4_new, 0.25) - temperature[i];

    //cap the per-step change so the interleaved hydrodynamics-dust cycle can
    //re-converge from the perturbed profile (the magnitude bound the hydro coupling
    //needs)
    const double max_change = config->temperature_max_change * temperature[i];
    if (std::abs(step) > max_change)
      step = std::copysign(max_change, step);

    //Apply the adaptive under-relaxation AFTER the cap. This is essential: when the
    //cap binds (model far from radiative equilibrium) applying omega before the cap
    //leaves the step pinned at +/-cap and the damping is invisible, so a layer whose
    //correction keeps flipping sign oscillates forever at the cap. Damping the
    //already-capped step lets omega actually shrink an oscillating layer's motion
    //while a steadily-converging layer (stable sign, omega -> 1) still moves at the cap.
    delta_temperature[i] = relaxation[i] * step;
  }

  return delta_temperature;
}


void TemperatureCorrection::smoothCorrection(std::vector<double>& correction)
{
  if (correction.size() < 3) return;

  std::vector<double> smoothed = correction;

  for (size_t i=1; i<correction.size()-1; ++i)
    smoothed[i] = 0.5 * correction[i] + 0.25 * (correction[i+1] + correction[i-1]);

  correction = smoothed;
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
  std::vector<double>& planck_function_int)
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
      y[j] = aux::planckFunctionWavelength(temperature[i], spectral_grid->wavelength_list[j]);

    planck_function_int[i] = radiation_field[i].wavelengthIntegration(y);

    for (size_t j=0; j<spectral_grid->nbSpectralPoints(); ++j)
      y[j] = absorption_coeff[i][j] * aux::planckFunctionWavelength(temperature[i], spectral_grid->wavelength_list[j]);

    kappa_b[i] = radiation_field[i].wavelengthIntegration(y) / planck_function_int[i];
  }

}



std::vector<double> TemperatureCorrection::unsoeldLucyCorrection(
  const std::vector<double>& temperature,
  const std::vector<double>& radius,
  const std::vector<double>& fq_j,
  const std::vector<double>& qx_h,
  const std::vector<double>& kappa_j,
  const std::vector<double>& kappa_b,
  const std::vector<double>& planck_function_int)
{
  std::cout << "Performing Unsöld-Lucy correction.\n\n";
  const size_t nb_grid_points = radius.size();

  std::vector<double> tau_h(nb_grid_points, 0);

  for (int i=nb_grid_points-2; i>-1; --i)
    tau_h[i] = tau_h[i+1] + 0.5*(radius[i+1] - radius[i])*(qx_h[i] + qx_h[i+1]);


  std::vector<double> delta_h(nb_grid_points, 0.);

  //flux deviation from the target. Use the SAME (eq. 2.58) integral as the rest of the
  //Unsoeld-Lucy factors (qx_h, the h-factor): mixing the conservative (eq. 2.59) integral
  //here with the 2.58 weighting made the correction internally inconsistent and oscillate
  //in the coupled loop. UL is only the robust phase-1 starter (its ~2.3% flux floor is
  //fine); the linearisation phase reaches true flux conservation.
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


  //Planck-integrated correction delta_B per layer. The conversion to a
  //temperature change (exact T^4 update + relaxation) is done by the caller.
  std::vector<double> delta_b(nb_grid_points, 0.);

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    delta_b[i] = (kappa_j[i]/kappa_b[i] * radiation_field[i].mean_intensity_int
                - planck_function_int[i])
                - kappa_j[i]/kappa_b[i]/(fq_j[i] * radius[i]*radius[i]) * integral[i];
  }

  return delta_b;
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