
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
#include "../additional/quadrature.h"


namespace agb{


RadiationField::RadiationField(const size_t nb_spectral_points)
{
  eddington_factor.assign(nb_spectral_points, 0.0);
  eddington_flux.assign(nb_spectral_points, 0.0);
  mean_intensity.assign(nb_spectral_points, 0.0);

  mean_intensity_impact.assign(nb_spectral_points, 0.);
  eddington_flux_impact.assign(nb_spectral_points, 0.);
  eddington_k_impact.assign(nb_spectral_points, 0.);
}


void RadiationField::createAngleGrid(
  std::vector<ImpactParam>& impact_parameter_grid,
  const double radius,
  const size_t radius_index_)
{
  radius_index = radius_index_;

  //count the number of angles per radius shell
  //there is one for every impact parameter crossing it
  nb_angles = 0;

  for (auto & ip : impact_parameter_grid)
  { 
    for (size_t k=0; k<ip.nb_z_points; ++k)
    {
      if (ip.z_grid[k].radius_index == radius_index)
        nb_angles++;
    }
      
  }

  const size_t nb_spectral_points = eddington_flux.size();
  angle_grid.assign(nb_angles, AnglePoint(nb_spectral_points));


  auto angle_it = angle_grid.begin();

  for (auto & ip : impact_parameter_grid)
  {
    for (size_t j=0; j<ip.nb_z_points; j++)
    {
      if (ip.z_grid[j].radius_index == radius_index)
      {
        double angle = ip.z_grid[j].z / radius;
        
        (*angle_it).angle = angle;

        ip.z_grid[j].angle_point = &(*angle_it);
        angle_it++;

        break;
      }
    }

    if (angle_it == angle_grid.end())
      break;
  }

  angles.assign(nb_angles, 0);

  for (size_t i=0; i<nb_angles; ++i)
    angles[i] = angle_grid[i].angle;
}


void RadiationField::angularIntegration()
{
  std::vector<double> y(nb_angles, 0);
  
  for (size_t i=0; i<mean_intensity.size(); ++i)
  {
    for (size_t j=0; j<nb_angles; ++j)
      y[j] = angle_grid[j].u[i];
    
    mean_intensity_impact[i] = - 0.5 * aux::quadratureTrapezoidal(angles, y);

    for (size_t j=0; j<nb_angles; ++j)
      y[j] *= angles[j];

    eddington_flux_impact[i] = - 0.5 * aux::quadratureTrapezoidal(angles, y);

    for (size_t j=0; j<nb_angles; ++j)
      y[j] *= angles[j];

    eddington_k_impact[i] = - 0.5 * aux::quadratureTrapezoidal(angles, y);
  }

}


}
