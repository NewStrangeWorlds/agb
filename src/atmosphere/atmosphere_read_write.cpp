 
#include "atmosphere.h"


#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <omp.h>
#include <stdlib.h>

#include "../config/config.h"
#include "../additional/exceptions.h"

namespace agb {

void Atmosphere::readStructure(const std::string file_path)
{
  std::fstream file;
  file.open(file_path.c_str(), std::ios::in);

  if (file.fail()) 
  {
    if (file.fail())
    throw FileNotFound(std::string ("Atmospheric structure file "), file_path);
  }

  std::string line;
  std::string input;

  //header
  std::getline(file, line);
  std::getline(file, line);


  while (std::getline(file, line))
  {
    std::istringstream input(line);

    double radius_grid_, radius_, rho, p, t_gas, t_dust, v;

    input >> radius_grid_ >> radius_ >> rho >> p >> t_gas >> t_dust >> v;

    radius_grid.push_back(radius_grid_);
    radius.push_back(radius_);
    mass_density.push_back(rho);
    pressure.push_back(p);
    temperature_gas.push_back(t_gas);
    temperature_dust.push_back(t_dust);
    velocity.push_back(v);
  }
  
  nb_grid_points = radius_grid.size();

  pressure_bar = cgsToBar(pressure);
}


}