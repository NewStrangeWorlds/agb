#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "spectral_grid.h"

#include "../config/config.h"
#include "../additional/exceptions.h"


namespace agb{


SpectralGrid::SpectralGrid(ModelConfig* config_)
{
  config = config_;

  loadWavenumberList();

  wavelength_list_full = wavenumberToWavelength(wavenumber_list_full);

  createHighResGrid();
}



void SpectralGrid::loadWavenumberList()
{
  std::string file_name = config->wavenumber_file_path;


  std::fstream file;

  file.open(file_name.c_str(), std::ios::in);


  if (file.fail())
    throw FileNotFound(std::string ("SpectralGrid::loadWavenumberList"), file_name);


  size_t nb_wavenumbers;
  file >> nb_wavenumbers;

  wavenumber_list_full.resize(nb_wavenumbers);


  for (std::vector<double>::iterator it = wavenumber_list_full.begin(); it != wavenumber_list_full.end(); ++it)
    file >> *it;


  file.close();
}


std::vector<double> SpectralGrid::createConstantResolutionGrid(
  const double min_wavelength,
  const double max_wavelength,
  const double spectral_resolution)
{
  std::vector<double> grid;
  grid.reserve(100000);
  
  grid.push_back(min_wavelength);

  while (grid.back() < max_wavelength)
  {
    if (grid.size() == grid.capacity() - 1)
      grid.reserve(grid.size() + 100000);

    const double delta_mu = grid.back() / spectral_resolution;
    grid.push_back(grid.back() + delta_mu); 
  }
  
  grid.shrink_to_fit();
  
  //since internally we work in wavenumbers, we reverse the wavelength grid
  std::reverse(grid.begin(), grid.end());

  return grid;
}



void SpectralGrid::createHighResGrid()
{
  std::vector<double> wavelength_grid = createConstantResolutionGrid(
    config->min_wavelength,
    config->max_wavelength,
    config->spectral_resolution);

  std::vector<double> wavenumber_grid = wavelengthToWavenumber(wavelength_grid);

  
  //find the closest points in the tabulated wavenumber list
  std::vector<int> included_points(wavenumber_list_full.size(), 0);
  auto it_start = wavenumber_list_full.begin();

  for (auto & nu : wavenumber_grid)
  {
    const size_t idx = findClosestIndex(nu, wavenumber_list_full, it_start);

    included_points[idx] = 1;
  }


  //and create the final grids
  index_list.resize(0);
  index_list.reserve(wavelength_list_full.size());

  for (size_t i=0; i<included_points.size(); ++i)
  {
    if (included_points[i] == 1)
      index_list.push_back(i);
  }


  index_list.shrink_to_fit();

  wavenumber_list.assign(index_list.size(), 0);
  wavelength_list.assign(index_list.size(), 0);

  for (size_t i=0; i<index_list.size(); ++i)
  {
    wavenumber_list[i] = wavenumber_list_full[index_list[i]];
    wavelength_list[i] = wavelength_list_full[index_list[i]];
  }
  
  nb_spectral_points = wavelength_list.size();
}




size_t SpectralGrid::findClosestIndexAsc(
  const double x,
  std::vector<double>& data,
  std::vector<double>::iterator start)
{
  auto iter_geq = std::lower_bound(
    start, 
    data.end(),
    x);

  if (iter_geq == start)
    return start - data.begin();

  double a = *(iter_geq - 1);
  double b = *(iter_geq);

  if (std::fabs(x - a) < fabs(x - b)) 
    return iter_geq - data.begin() - 1;

  return iter_geq - data.begin();
}


size_t SpectralGrid::findClosestIndexDesc(
  const double x,
  std::vector<double>& data,
  std::vector<double>::iterator start)
{
  auto iter_geq = std::lower_bound(
    start, 
    data.end(),
    x,
    std::greater<double>());

  if (iter_geq == start)
    return start - data.begin();

  double a = *(iter_geq - 1);
  double b = *(iter_geq);

  if (std::fabs(x - a) < fabs(x - b)) 
    return iter_geq - data.begin() - 1;

  return iter_geq - data.begin();
}


size_t SpectralGrid::findClosestIndex(
  const double x,
  std::vector<double>& data,
  std::vector<double>::iterator start)
{
  if (data.size() < 2)
    return 0;
  
  if (data[0] < data[1])
    return findClosestIndexAsc(x, data, start);

  if (data[0] > data[1])
    return findClosestIndexDesc(x, data, start);

  return 0;
}



/*std::vector<double> SpectralGrid::wavenumberList(const std::vector<size_t>& indices)
{
  std::vector<double> output(indices.size(), 0.0);

  for (size_t i=0; i<indices.size(); ++i)
    output[i] = wavenumber_list[indices[i]];

  return output;
}



std::vector<double> SpectralGrid::wavelengthList(const std::vector<size_t>& indices)
{
  std::vector<double> output(indices.size(), 0.0);

  for (size_t i=0; i<indices.size(); ++i)
    output[i] = wavelength_list[indices[i]];

  return output;
}*/



}






