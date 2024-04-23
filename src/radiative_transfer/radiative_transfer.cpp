
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
  eddington_factor_h.assign(nb_spectral_points, std::vector<double>(nb_grid_points, 0.0));
  extinction_coeff.assign(nb_spectral_points, std::vector<double>(nb_grid_points, 0.0));
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
      source_function[j][i] = atmosphere->absorption_coeff[i][j] 
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
  
  return boundary_flux/integral;
}



void RadiativeTransfer::calcEddingtonFactors()
{
  for (size_t i=0; i<nb_spectral_points; ++i)
  {
    for (size_t j=0; j<nb_grid_points; ++j)
    {
      eddington_factor_f[i][j] = radiation_field[j].eddington_k_impact[i]/radiation_field[j].mean_intensity_impact[i];
      eddington_factor_h[i][j] = radiation_field[j].eddington_flux_impact[i]/radiation_field[j].mean_intensity_impact[i];
    }
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
      double delta_q = (atmosphere->radius[j] - atmosphere->radius[j-1])
                       * (y[j] + y[j+1]) / 2.0;
      integral += delta_q;

      sphericality_factor[i][j] = std::exp(integral 
                                 + std::log(atmosphere->radius[0]*atmosphere->radius[0]
                                          /(atmosphere->radius[j]*atmosphere->radius[j])));
    }
  }

}


void RadiativeTransfer::solveRadiativeTransfer()
{
  for (size_t i=0; i<nb_grid_points; ++i)
    for (size_t j=0; j<nb_spectral_points; ++j)
      extinction_coeff[j][i] = atmosphere->absorption_coeff[i][j] + atmosphere->scattering_coeff[i][j];

  calcSourceFunction();

  double boundary_flux_correction = boundaryFluxCorrection();
  std::vector<double> boundary_planck_derivative(nb_spectral_points, 0.);

  for (size_t i=0; i<nb_spectral_points; ++i)
    boundary_planck_derivative[i] = aux::planckFunctionDerivWavelength(
      atmosphere->temperature_gas[0], 
      spectral_grid->wavelength_list_cm[i]);

  // impact_parameter_grid[30].solveRadiativeTransfer(
  //   100,
  //   boundary_planck_derivative[100],
  //   boundary_flux_correction,
  //   extinction_coeff[100],
  //   source_function[100]);

  // exit(0);

  for (size_t i=0; i<nb_spectral_points; ++i)
    for (auto & ip : impact_parameter_grid)
      ip.solveRadiativeTransfer(
        i,
        boundary_planck_derivative[i],
        boundary_flux_correction,
        extinction_coeff[i],
        source_function[i]);
  
  for (size_t i=0; i<nb_grid_points; ++i)
    radiation_field[i].angularIntegration();

  calcEddingtonFactors();
  calcSphericalityFactor();

  for (size_t i=0; i<nb_grid_points; ++i)
    std::cout << i << "\t" << radiation_field[i].mean_intensity_impact[100] << "\t" << radiation_field[i].eddington_flux_impact[100] << "\t" << radiation_field[i].eddington_k_impact[100] << "\t" << radiation_field[i].eddington_k_impact[100]/radiation_field[i].mean_intensity_impact[100] << "\t" << sphericality_factor[100][i] << "\n";

  exit(0);
}



}