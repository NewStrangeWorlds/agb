 
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
#include "../additional/physical_const.h"

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



void Atmosphere::writeStructure(const std::string file_path)
{
  std::fstream file(file_path.c_str(), std::ios::out);
  
  if (file.fail())
  {
    std::cout << "Unable to open atmosphere output file " << file_path << "\n";
    return;
  }

  std::cout << "Saving atmospheric structure to " << file_path << "\n\n";

  file << std::setprecision(4) << std::fixed << "#Model: "
       << "R* = " << config->stellar_radius/constants::radius_sun << "   "
       << "M* = " << config->stellar_mass/constants::mass_sun << "  "
       << "L* = " << std::scientific << config->stellar_luminosity/constants::luminosity_sun << "  "
       << "C/O = " << std::fixed << config->c_o_ratio << "\n";

  file << std::setw(16) << std::left << "#r/R*" << "\t"
       << std::setw(16) << std::left << "r(cm)" << "\t"
       << std::setw(16) << std::left << "rho(g/cm3)" << "\t"
       << std::setw(16) << std::left << "p(dyn/cm2)" << "\t"
       << std::setw(16) << std::left << "T_gas(K)" << "\t"
       << std::setw(16) << std::left << "T_dust(K)" << "\t"
       << std::setw(16) << std::left << "v(cm/s)" << "\n";

  for (size_t i=0; i<radius.size(); i++)
  {
    file << std::setprecision(10) << std::scientific
         << radius_grid[i] << "\t"
         << radius[i] << "\t"
         << mass_density[i] << "\t"
         << pressure[i] << "\t"
         << temperature_gas[i] << "\t"
         << temperature_dust[i] << "\t"
         << velocity[i] << "\n";
  }

  file.close();
}


}