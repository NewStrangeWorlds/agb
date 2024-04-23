
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


void RadiativeTransfer::createImpactParameterGrid()
{
  impact_parameter_grid.resize(nb_impact_param);

  impact_parameter_grid[0].p = 0;
  impact_parameter_grid[0].nb_z_points = nb_grid_points;


  for (size_t i=1; i<nb_core_impact_param; ++i)
  {
    double a = 1.0 - 1.0/nb_core_impact_param*i;

    impact_parameter_grid[i].p = std::sqrt(
      atmosphere->radius[0] * atmosphere->radius[0]
      - atmosphere->radius[0] * atmosphere->radius[0]*a*a);
    
    impact_parameter_grid[i].nb_z_points = nb_grid_points;
  }


  for (size_t i=0; i<nb_grid_points; i++)
  {
    impact_parameter_grid[i+nb_core_impact_param].p = atmosphere->radius[i];
    impact_parameter_grid[i+nb_core_impact_param].nb_z_points = nb_grid_points - i;
  }

  createZGrids();
}



void RadiativeTransfer::createZGrids()
{
  for (auto & ip : impact_parameter_grid)
    ip.z_grid.resize(ip.nb_z_points);

  for (size_t k=0; k<nb_core_impact_param; k++)
  {
    for (size_t i = 0; i<impact_parameter_grid[k].nb_z_points; i++)
    {
      size_t j = nb_grid_points - impact_parameter_grid[k].nb_z_points + i;
      double d = atmosphere->radius[j] * atmosphere->radius[j] - impact_parameter_grid[k].p * impact_parameter_grid[k].p;
      
      impact_parameter_grid[k].z_grid[i].z = std::sqrt(d);
      impact_parameter_grid[k].z_grid[i].radius_index = j;
      impact_parameter_grid[k].z_grid[i].radius = atmosphere->radius[j];
    }
  }

  for (size_t k=nb_core_impact_param; k<nb_impact_param; k++)
    for (size_t i = 0; i<impact_parameter_grid[k].nb_z_points; i++)
    {
      size_t j = nb_grid_points - impact_parameter_grid[k].nb_z_points + i;
      double d = atmosphere->radius[j] * atmosphere->radius[j] - impact_parameter_grid[k].p * impact_parameter_grid[k].p;

      impact_parameter_grid[k].z_grid[i].z = std::sqrt(d);
      impact_parameter_grid[k].z_grid[i].radius_index = j;
      impact_parameter_grid[k].z_grid[i].radius = atmosphere->radius[j];
    }
}


//calculate the optical depth along an impact parameter from 
//the outside towards the inside
//uses a trapezoidal rule
std::vector<double> ImpactParam::opticalDepth(
  const std::vector<double>& extinction_coeff)
{
  std::vector<double> optical_depth(nb_z_points, 0.0);
  
  for (int i=nb_z_points-2; i>-1; --i)
  {
    double q = -(z_grid[i].radius - z_grid[i+1].radius) 
               * (extinction_coeff[z_grid[i].radius_index] 
                + extinction_coeff[z_grid[i+1].radius_index]) / 2.0;
    optical_depth[i] = optical_depth[i+1] + q;
  }

  return optical_depth;
}


void ImpactParam::assembleSystem(
  const std::vector<double>& optical_depth,
  const std::vector<double>& source_function,
  const double boundary_planck_derivative,
  const double boundary_flux_correction,
  const double boundary_exctinction_coeff,
  aux::TriDiagonalMatrix& M,
  std::vector<double>& rhs)
{
  double hr = optical_depth[1] - optical_depth[0];

  M.b[0] = - 1/hr - hr/2;
  M.c[0] = 1/hr;

  if (z_grid[0].z == 0)
  {
    rhs[0] = - hr/2 * source_function[0];
  }
  else
  {
    rhs[0] = - hr/2 * source_function[0] + z_grid[0].angle_point->angle
		        * 3 * boundary_planck_derivative * boundary_flux_correction / boundary_exctinction_coeff;
  }
  
  for (size_t i=1; i<nb_z_points-1; ++i)
  {
    double hr = optical_depth[i+1] - optical_depth[i];
    double hl = optical_depth[i] - optical_depth[i-1];

    M.a[i] = 2/(hl*(hl+hr));
    M.b[i] = - 2/(hr*hl) - 1;
    M.c[i] = 2/(hr*(hl+hr));
    rhs[i] = - source_function[i];
  }

  double hl = optical_depth.back() - optical_depth[nb_z_points-2];

  M.a.back() = 1/hl;
  M.b.back() = -1/hl - hl/2 + 1;
  rhs.back() = -hl/2 *  source_function.back();
}


//solves the radiative transfer equation along an impact parameter
//for a given spectral index nu
void ImpactParam::solveRadiativeTransfer(
  const size_t nu,
  const double boundary_planck_derivative,
  const double boundary_flux_correction,
  const std::vector<double>& extinction_coeff,
  const std::vector<double>& source_function)
{ 
  //the uppermost impact parameter
  if (nb_z_points == 1)
  {
    z_grid[0].angle_point->u[nu] = 0;

    return;
  }


  std::vector<double> optical_depth = opticalDepth(
    extinction_coeff);
   
  std::vector<double> source_function_z;
  source_function_z.reserve(nb_z_points);

  for (auto & z : z_grid)
    source_function_z.push_back(source_function[z.radius_index]);

  aux::TriDiagonalMatrix M(nb_z_points);
  std::vector<double> rhs(nb_z_points, 0.);

  assembleSystem(
    optical_depth,
    source_function_z,
    boundary_planck_derivative,
    boundary_flux_correction,
    extinction_coeff[0],
    M,
    rhs);

  std::vector<double> u = M.solve(rhs);

  // for (size_t i=0; i<nb_z_points; ++i)
  //   std::cout << i << "\t" << optical_depth[i] << "\t" << M.a[i] << "\t" << M.b[i] << "\t" << M.c[i] << "\t" << rhs[i] << "\t" << u[i] << "\n";

  for (size_t i=0; i<nb_z_points; ++i)
  {
    //if (u[i] < 1e-45) u[i] = 1e-45;

    z_grid[i].angle_point->u[nu] = u[i];
  }

}


} 

