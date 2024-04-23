
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
    }
}


//calculate the optical depth along an impact parameter from 
//the outside towards the inside
//uses a trapezoidal rule
std::vector<double> ImpactParam::opticalDepth(
  const size_t nu,
  const std::vector<double>& extinction_coeff)
{
  std::vector<double> optical_depth(nb_z_points, 0.0);
  
  for (int i=nb_z_points-2; i>-1; --i)
  {
    double q = -(z_grid[i].radius - z_grid[i+1].radius) 
               * (extinction_coeff[z_grid[i].radius_index] 
                + extinction_coeff[z_grid[i+1].radius_index]);
    optical_depth[i] = optical_depth[i+1] + q;
  }

  return optical_depth;
}



void ImpactParam::solveRadiativeTransfer(
  const size_t nu,
  const std::vector<double>& extinction_coeff,
  const std::vector<double>& source_function)
{
  
  std::vector<double> optical_depth = opticalDepth(
    nu,
    extinction_coeff);
 
  for (size_t i=0; i<nb_z_points; ++i)
    std::cout << i << "\t" << optical_depth[i] << "\t" << extinction_coeff[z_grid[i].radius_index] << "\n";
}


/*static void SolveImpactDGL(int imp, int freq)
{
  int r,a,j;
  double* TauMesh = (double* ) calloc(ImpactMesh[imp].zPointNumber,sizeof(double));
  double* b = (double* ) calloc(ImpactMesh[imp].zPointNumber,sizeof(double));
  double* x = (double* ) calloc(ImpactMesh[imp].zPointNumber,sizeof(double));
  double* SourceFunction = (double* ) calloc(ImpactMesh[imp].zPointNumber,sizeof(double));
  struct TriDiagMatrix A;

  InitTriDiagMatrix(&A,ImpactMesh[imp].zPointNumber);

  if (Config.ImpactSystem.MeshCreation == RFSpline)
  {
    if (ImpactMesh[imp].zPointNumber > 2)  
	  TauMeshSplines(imp,freq,TauMesh);
    else  
	  TauMeshTrapez(imp,freq,TauMesh);
  }

  if (Config.ImpactSystem.MeshCreation == RFTrapezoid)  TauMeshTrapez(imp,freq,TauMesh);


  for (j=0; j<ImpactMesh[imp].zPointNumber; j++)
	SourceFunction[j] = GetSourceFunc(ImpactMesh[imp].zMesh[j].RadiusIndex,freq);

  
  if (Config.ImpactSystem.DiscType == RFCubicSplines)  ImpactCoeffSpline(A,b,TauMesh,SourceFunction,imp,freq);
  if (Config.ImpactSystem.DiscType == RFHermiteFormula)  ImpactCoeffHermite(A,b,TauMesh,SourceFunction,imp,freq);
  if (Config.ImpactSystem.DiscType == RFTaylorDiff)  ImpactCoeffDiff(A,b,TauMesh,SourceFunction,imp,freq);

  
  SolveTriDiagEq(A,b,x,ImpactMesh[imp].zPointNumber);

  
  for (j=0; j<ImpactMesh[imp].zPointNumber; j++)
  {
    r = ImpactMesh[imp].zMesh[j].RadiusIndex;
    a = ImpactMesh[imp].zMesh[j].AngleIndex;

    if (x[j] < 1e-45) x[j] = 1e-45;
	RadField[r].AngleMesh[a].U[freq] = x[j];
  }

  free(TauMesh);
  free(x);
  free(b);
  free(SourceFunction);
  FreeTriDiagMatrix(&A);
}*/




} 

