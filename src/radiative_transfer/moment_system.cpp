 
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

#include "../spectral_grid/spectral_grid.h"
#include "../config/config.h"
#include "../atmosphere/atmosphere.h"
#include "../additional/tri_diagonal_matrix.h"
#include "../additional/aux_functions.h"


namespace agb{


void RadiativeTransfer::solveMomentSystem(
  const size_t nu,
  const std::vector<double>& radius,
  const std::vector<double>& radius2,
  const double boundary_planck_derivative,
  const double boundary_flux_correction)
{
  const std::vector<double> x_grid = generateXGrid(
    extinction_coeff[nu], 
    radius,
    sphericality_factor[nu]);

  std::vector<double> emission_coeff(nb_grid_points, 0);

  for (size_t i=0; i<nb_grid_points; ++i)
    emission_coeff[i] = aux::planckFunctionWavelength(
                          atmosphere->temperature_gas[i], 
                          spectral_grid->wavelength_list[nu]) 
                      * atmosphere->absorption_coeff[i][nu];

  aux::TriDiagonalMatrix m(nb_grid_points);
  std::vector<double> rhs(nb_grid_points, 0.);

  assembleMomentSystem(
    x_grid,
    radius,
    radius2,
    emission_coeff,
    extinction_coeff[nu],
    scattering_coeff[nu],
    eddington_factor_f[nu],
    boundary_eddington_factor_h[nu],
    sphericality_factor[nu],
    boundary_planck_derivative,
    boundary_flux_correction,
    m,
    rhs);

  std::vector<double> result = m.solve(rhs);

  for (size_t i=0; i<nb_grid_points; ++i)
    radiation_field[i].mean_intensity[nu] = result[i];
}


std::vector<double> RadiativeTransfer::generateXGrid(
  const std::vector<double>& extinction_coeff,
  const std::vector<double>& radius,
  const std::vector<double>& sphericality_factor)
{
  std::vector<double> x_grid(nb_grid_points, 0.);

  for (int i=nb_grid_points-2; i>-1; --i)
  {
    double integral = - (radius[i] - radius[i+1]) * (extinction_coeff[i]*sphericality_factor[i] + extinction_coeff[i+1]*sphericality_factor[i+1]) / 2.0;

    x_grid[i] = x_grid[i+1] + integral;
  }

  return x_grid;
}


//Spline discretisation
void RadiativeTransfer::assembleMomentSystem(
  const std::vector<double>& x_grid,
  const std::vector<double>& radius,
  const std::vector<double>& radius2,
  const std::vector<double>& emission_coeff,
  const std::vector<double>& extinction_coeff,
  const std::vector<double>& scattering_coeff,
  const std::vector<double>& eddington_factor,
  const double boundary_eddington_h,
  const std::vector<double>& sphericality_factor,
  const double boundary_planck_derivative,
  const double boundary_flux_correction,
  aux::TriDiagonalMatrix& m,
  std::vector<double>& rhs)
{
  std::vector<double> a(nb_grid_points, 0.);

  for (size_t i=0; i<nb_grid_points; ++i)
    a[i] = eddington_factor[i] * sphericality_factor[i] * radius2[i];

  auto b = [&](size_t i) {
    return radius2[i]/sphericality_factor[i] * (1. - scattering_coeff[i]/extinction_coeff[i]);};

  auto c = [&](size_t i) {
    return radius2[i]/sphericality_factor[i] / extinction_coeff[i];};

  const double hr = x_grid[1] - x_grid.front();

  m.b[0] = hr/3. * b(0) + 1/hr * a[0];
  m.c[0] = hr/6. * b(1) - 1/hr * a[1];

  rhs[0] = hr/3. * c(0) * emission_coeff[0]
         + hr/6. * c(1) * emission_coeff[1]
         - radius2[0]/extinction_coeff[0] * boundary_flux_correction * boundary_planck_derivative;
  
  for (size_t i=1 ; i<nb_grid_points-1; ++i)
  {
    const double hr = x_grid[i+1] - x_grid[i];
	  const double hl = x_grid[i] - x_grid[i-1];

    m.a[i] =  hl/6. * b(i-1) - 1./hl * a[i-1];
	  m.b[i] = (hr+hl)/3. * b(i) + 1./hr * a[i] + 1./hl * a[i];
	  m.c[i] = hr/6. * b(i+1) - 1./hr * a[i+1];

	  rhs[i] = hl/6. * c(i-1) * emission_coeff[i-1]
	         + (hr+hl)/3. * c(i) * emission_coeff[i]
		       + hr/6. * c(i+1) * emission_coeff[i+1];
  }

  const double hl = x_grid.back() - x_grid[nb_grid_points-2];

  m.a.back() = hl/6. * b(nb_grid_points-2) - 1./hl * a[nb_grid_points-2];
  m.b.back() = hl/3. * b(nb_grid_points-1) + 1./hl * a[nb_grid_points-1] - radius2.back() * boundary_eddington_h;
  rhs.back() = hl/6. * c(nb_grid_points-2) * emission_coeff[nb_grid_points-2] + hl/3. * c(nb_grid_points-1) * emission_coeff[nb_grid_points-1];
}



void RadiativeTransfer::calcFlux(
  const size_t nu,
  const std::vector<double>& radius,
  const std::vector<double>& radius2,
  const std::vector<double>& source_function)
{
  const std::vector<double> x_grid = generateXGrid(
    extinction_coeff[nu], 
    radius,
    sphericality_factor[nu]);

  std::vector<double> mean_intensity(nb_grid_points, 0);

  for (size_t i=0; i<nb_grid_points; ++i)
    mean_intensity[i] = radiation_field[i].mean_intensity[nu];

  aux::TriDiagonalMatrix m(nb_grid_points);
  std::vector<double> rhs(nb_grid_points, 0.);

  assembleMomentSystemFlux(
    x_grid,
    radius2,
    source_function,
    mean_intensity,
    eddington_factor_f[nu],
    sphericality_factor[nu],
    m,
    rhs);

  std::vector<double> result = m.solve(rhs);

  for (size_t i=0; i<nb_grid_points; ++i)
    radiation_field[i].eddington_flux[nu] = result[i]/radius2[i];
}


//spline discretisation
void RadiativeTransfer::assembleMomentSystemFlux(
  const std::vector<double>& x_grid,
  const std::vector<double>& radius2,
  const std::vector<double>& source_function,
  const std::vector<double>& mean_intensity,
  const std::vector<double>& eddington_factor,
  const std::vector<double>& sphericality_factor,
  aux::TriDiagonalMatrix& m,
  std::vector<double>& rhs)
{
  std::vector<double> a(nb_grid_points, 0.);

  for (size_t i=0; i<nb_grid_points; ++i)
    a[i] = eddington_factor[i] * sphericality_factor[i] * radius2[i] * mean_intensity[i];
  
  const double hr = x_grid[1] - x_grid.front();
  
  m.c[0] = 2./hr;
  m.b[0] = 4./hr;
  rhs[0] = 6./hr/hr * (a[1] - a[0])
           - radius2[0]/sphericality_factor[0] * (mean_intensity[0] - source_function[0]); 

  for (size_t i=1; i<nb_grid_points-1; ++i)
  {
    const double hr = x_grid[i+1] - x_grid[i];
	  const double hl = x_grid[i] - x_grid[i-1];

    m.a[i] = 1./hl;
	  m.b[i] = 2. * (1./hl + 1./hr);
    m.c[i] = 1./hr;
    rhs[i] = 3./(hl*hl) * (a[i] - a[i-1]) + 3./(hr*hr) * (a[i+1] - a[i]);
  }

  const double hl = x_grid.back() - x_grid[nb_grid_points-2];

  m.a.back() = 2./hl;
  m.b.back() = 4./hl;
  rhs.back() = 6./(hl*hl) * (a.back() - a[nb_grid_points-2])
               + radius2.back()/sphericality_factor.back()
               * (mean_intensity.back() - source_function.back());
}



//Taylor discretisation
/*void RadiativeTransfer::assembleMomentSystem(
  const std::vector<double>& x_grid,
  const std::vector<double>& radius,
  const std::vector<double>& radius2,
  const std::vector<double>& emission_coeff,
  const std::vector<double>& extinction_coeff,
  const std::vector<double>& scattering_coeff,
  const std::vector<double>& eddington_factor,
  const double boundary_eddington_h,
  const std::vector<double>& sphericality_factor,
  const double boundary_planck_derivative,
  const double boundary_flux_correction,
  aux::TriDiagonalMatrix& m,
  std::vector<double>& rhs)
{
  const double hr = x_grid[1] - x_grid.front();

  std::vector<double> a(nb_grid_points, 0.);

  for (size_t i=0; i<nb_grid_points; ++i)
    a[i] = eddington_factor[i] * sphericality_factor[i] * radius2[i];

  auto b = [&](size_t i) {
    return radius2[i]/sphericality_factor[i] * (1 - scattering_coeff[i]/extinction_coeff[i]);};

  auto c = [&](size_t i) {
    return radius2[i]/sphericality_factor[i] / extinction_coeff[i];};

  m.b[0] = -1/hr * a[0] - hr/2 * b(0);
  m.c[0] = 1/hr * a[1];

  rhs[0] = radius2[0]/extinction_coeff[0] * boundary_flux_correction * boundary_planck_derivative
           - hr/2. * c(0) * emission_coeff[0];
  
  for (size_t i=1 ; i<nb_grid_points-1; ++i)
  {
    const double hr = x_grid[i+1] - x_grid[i];
	  const double hl = x_grid[i] - x_grid[i-1];

    m.a[i] = 2./(hl*(hl+hr)) * a[i-1];
	  m.b[i] = - 2./(hl*hr) * a[i] - b(i);
	  m.c[i] = 2./(hr*(hr+hl)) * a[i+1];

	  rhs[i] = - c(i) * emission_coeff[i];
  }

  const double hl = x_grid.back() - x_grid[nb_grid_points-2];

  m.a.back() = 1./hl * a[nb_grid_points-2];
  m.b.back() = -1./hl * a.back() - hl/2 * b(nb_grid_points-1) + radius2.back() * boundary_eddington_h;

  rhs.back() = -hl/2. * c(nb_grid_points-1) * emission_coeff.back();
}*/


/*void RadiativeTransfer::assembleMomentSystemFlux(
  const std::vector<double>& x_grid,
  const std::vector<double>& radius2,
  const std::vector<double>& source_function,
  const std::vector<double>& mean_intensity,
  const std::vector<double>& eddington_factor,
  const std::vector<double>& sphericality_factor,
  aux::TriDiagonalMatrix& m,
  std::vector<double>& rhs)
{
  std::vector<double> a(nb_grid_points, 0.);

  for (size_t i=0; i<nb_grid_points; ++i)
    a[i] = eddington_factor[i] * sphericality_factor[i] * radius2[i];
  
  const double hr = x_grid[1] - x_grid.front();
  
  m.c[0] = 0.;
  m.b[0] = 1.;
  
  rhs[0] = 1./hr * (a[1]*mean_intensity[1] - a[0]*mean_intensity[0])
           - hr/2. * radius2[0]/sphericality_factor[0] * (mean_intensity[0] - source_function[0]); 

  for (size_t i=1; i<nb_grid_points-1; ++i)
  {
    const double hr = x_grid[i+1] - x_grid[i];
	  const double hl = x_grid[i] - x_grid[i-1];

    m.a[i] = 0;
	  m.b[i] = 1;
    m.c[i] = 0;
	  
    rhs[i] = - hr/(hl*(hl+hr)) * a[i-1] * mean_intensity[i-1]
	           + (hr-hl)/(hr*hl) * a[i] * mean_intensity[i]
			       + hl/(hr*(hr+hl)) * a[i+1] * mean_intensity[i+1];
  }

  const double hl = x_grid.back() - x_grid[nb_grid_points-2];

  m.a.back() = 0;
  m.b.back() = 1;
  rhs.back() = 1/hl*(a.back() * mean_intensity.back() - a[nb_grid_points-2] * mean_intensity[nb_grid_points-2])
	           + hl/2 * radius2.back()/sphericality_factor.back()
             * (mean_intensity.back() - source_function.back());
}*/


}