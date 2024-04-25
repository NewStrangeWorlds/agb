
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
#include "../additional/quadrature.h"
#include "../additional/physical_const.h"

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
  eddington_factor_f.assign(nb_spectral_points, std::vector<double>(nb_grid_points, 0.0));
  boundary_eddington_factor_h.assign(nb_spectral_points, 0.);
  extinction_coeff.assign(nb_spectral_points, std::vector<double>(nb_grid_points, 0.0));
  scattering_coeff.assign(nb_spectral_points, std::vector<double>(nb_grid_points, 0.0));
  source_function.assign(nb_spectral_points, std::vector<double>(nb_grid_points, 0.0));
  sphericality_factor.assign(nb_spectral_points, std::vector<double>(nb_grid_points, 0.0));

  createImpactParameterGrid();

  for (size_t i=0; i<nb_grid_points; ++i)
    radiation_field[i].createAngleGrid(impact_parameter_grid, atmosphere->radius[i], i);
}


void RadiativeTransfer::calcSourceFunction()
{
  for (size_t i=0; i<nb_grid_points; ++i)
    for (size_t j=0; j<nb_spectral_points; ++j)
      source_function[j][i] = atmosphere->absorption_coeff[i][j] / extinction_coeff[j][i]
        * aux::planckFunctionWavelength(
            atmosphere->temperature_gas[i], 
            spectral_grid->wavelength_list[j])
        + atmosphere->scattering_coeff[i][j]/extinction_coeff[j][i] 
        * radiation_field[i].mean_intensity[j];
}


double RadiativeTransfer::boundaryFluxCorrection()
{
  std::vector<double> y(nb_spectral_points, 0);

  for (size_t i=0; i<nb_spectral_points; ++i)
    y[i] = aux::planckFunctionDerivWavelength(
      atmosphere->temperature_gas[0], 
      spectral_grid->wavelength_list[i])
      / extinction_coeff[i][0];
  
  const double integral = -aux::quadratureTrapezoidal(spectral_grid->wavelength_list_cm, y);
  const double boundary_flux = config->stellar_luminosity 
                               / (16. * constants::pi * constants::pi * std::pow(atmosphere->radius[0],2.));

  // for (size_t i=0; i<nb_spectral_points; ++i)
  //   y[i] = aux::planckFunctionDerivWavenumber(
  //     atmosphere->temperature_gas[0], 
  //     spectral_grid->wavenumber_list[i])
  //     / extinction_coeff[i][0];

  // const double integral2 = aux::quadratureTrapezoidal(spectral_grid->wavenumber_list, y);

  return boundary_flux/integral;
}



void RadiativeTransfer::calcEddingtonFactors()
{
  std::vector<double> x = radiation_field.back().angles;
  std::vector<double> y(x.size(), 0);


  for (size_t i=0; i<nb_spectral_points; ++i)
  {
    for (size_t j=0; j<x.size(); ++j)
      y[j] = radiation_field.back().angle_grid[j].u[i]*x[j];

    const double boundary_flux = - aux::quadratureTrapezoidal(x, y);

    boundary_eddington_factor_h[i] = boundary_flux/radiation_field.back().mean_intensity_impact[i];

    for (size_t j=0; j<nb_grid_points; ++j)
      eddington_factor_f[i][j] = radiation_field[j].eddington_k_impact[i]/radiation_field[j].mean_intensity_impact[i];
  }

}


void RadiativeTransfer::calcSphericalityFactor()
{

  for (size_t i=0; i<nb_spectral_points; ++i)
  { 
    std::vector<double> y(nb_grid_points, 0);

    for (size_t j=0; j<nb_grid_points; ++j)
      y[j] = (3*eddington_factor_f[i][j] - 1)/(atmosphere->radius[j]*eddington_factor_f[i][j]);

    sphericality_factor[i][0] = 1.0;
    
    double integral = 0;

    for (size_t j=1; j<nb_grid_points; ++j)
    {
      double delta = (atmosphere->radius[j] - atmosphere->radius[j-1])
                       * (y[j] + y[j+1]) / 2.0;
      integral += delta;

      sphericality_factor[i][j] = std::exp(integral + std::log(atmosphere->radius[0]*atmosphere->radius[0]))
                                  /(atmosphere->radius[j]*atmosphere->radius[j]);
    }
  }

}


void RadiativeTransfer::solveRadiativeTransfer()
{
  for (size_t i=0; i<nb_grid_points; ++i)
    for (size_t j=0; j<nb_spectral_points; ++j)
    {
      extinction_coeff[j][i] = atmosphere->absorption_coeff[i][j] + atmosphere->scattering_coeff[i][j];
      scattering_coeff[j][i] = atmosphere->scattering_coeff[i][j];
    }

  double boundary_flux_correction = boundaryFluxCorrection();
  std::vector<double> boundary_planck_derivative(nb_spectral_points, 0.);

  for (size_t i=0; i<nb_spectral_points; ++i)
    boundary_planck_derivative[i] = aux::planckFunctionDerivWavelength(
      atmosphere->temperature_gas[0], 
      spectral_grid->wavelength_list[i]);

  std::vector<double> radius(nb_grid_points, 0);
  std::vector<double> radius2(nb_grid_points, 0);

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    radius[i] = atmosphere->radius[i];
    radius2[i] = radius[i]*radius[i];
  }


  //intial values
  eddington_factor_f.assign(nb_spectral_points, std::vector<double>(nb_grid_points, 0.33333));
  boundary_eddington_factor_h.assign(nb_spectral_points, 0.5);

  calcSphericalityFactor();


  std::vector<std::vector<double>> eddington_factor_prev = eddington_factor_f;

  for (size_t it=0; it<100; ++it)
  { 
    #pragma omp parallel for
    for (size_t i=0; i<nb_spectral_points; ++i)
      solveMomentSystem(
        i,
        radius,
        radius2,
        boundary_planck_derivative[i], 
        boundary_flux_correction);

    calcSourceFunction();
    
    #pragma omp parallel for
    for (size_t i=0; i<nb_spectral_points; ++i)
    {
      for (auto & ip : impact_parameter_grid)
        ip.solveRadiativeTransfer(
          i,
          boundary_planck_derivative[i],
          boundary_flux_correction,
          extinction_coeff[i],
          source_function[i]);
    }

    for (size_t i=0; i<nb_grid_points; ++i)
      radiation_field[i].angularIntegration();
  
    calcEddingtonFactors();
    calcSphericalityFactor();

    double convergence = checkConvergence(eddington_factor_prev, eddington_factor_f);
    
    std::cout << it << "  " << convergence << "\n";
    
    if (convergence < 1e-4) break;

    eddington_factor_prev = eddington_factor_f;
  }

  #pragma omp parallel for
  for (size_t i=0; i<nb_spectral_points; ++i)
    calcFlux(
      i,
      radius,
      radius2,
      source_function[i]);

  for (size_t i=0; i<nb_grid_points; ++i)
    std::cout << i << "\t" << radiation_field[i].mean_intensity_impact[200] << "\t" << radiation_field[i].mean_intensity[200] << "\t" << eddington_factor_f[200][i] << "\t" << radiation_field[i].eddington_flux[200] << "\n";

  exit(0);
}


double RadiativeTransfer::checkConvergence(
  std::vector<std::vector<double>>& old_values,
  std::vector<std::vector<double>>& new_values)
{
  double max_rel_difference = 0;

  for (size_t i=0; i<old_values.size(); ++i)
    for (size_t j=0; j<old_values[i].size(); ++j)
    {
      double rel_difference = std::abs((old_values[i][j] - new_values[i][j])/old_values[i][j]);

      if (rel_difference > max_rel_difference)
        max_rel_difference = rel_difference;
    }
  
  return max_rel_difference;
}



}