
#include "config.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <omp.h>
#include <stdlib.h>

#include "../additional/physical_const.h"

namespace agb {

ModelConfig::ModelConfig(const std::string model_folder)
{
  loadConfigFile(model_folder);
}


bool ModelConfig::loadConfigFile(const std::string folder_path)
{
  model_folder = folder_path;

  if (model_folder.back() != '/')
    model_folder.append("/");

  std::string file_path = model_folder;
  file_path.append("model.config");

  
  std::fstream file;
  file.open(file_path.c_str(), std::ios::in);

  if (file.fail()) 
  {
    std::cout << "Couldn't open model options file " << file_path << "\n";
    
    return false;
  }

  
  std::cout << "\nParameters found in model.config:\n";

  std::string line;
  std::string input;
  
  std::getline(file, line);

  std::cout << "General Model Parameters\n";

  file >> stellar_radius;
  std::cout << "- R*: " << stellar_radius << "\n";
  stellar_radius *= constants::radius_sun;

  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  file >> stellar_mass;

  std::cout << "- M*: " << stellar_mass << "\n";
  stellar_mass *= constants::mass_sun;

  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  
  file >> stellar_luminosity;
  std::cout << "- L*: " << stellar_luminosity << "\n";
  stellar_luminosity *= constants::luminosity_sun;
  
  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);

  file >> c_o_ratio;
  std::cout << "- C/O: " << c_o_ratio << "\n";


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  
  file >> opacity_path;

  if (opacity_path.back() != '/')
    opacity_path.append("/");

  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  file >> min_wavelength >> max_wavelength >> spectral_resolution;

  std::cout << "- opacity folder: " << opacity_path << "\n";
  std::cout << "- min wavelength: " << min_wavelength << "\n";
  std::cout << "- max wavelength: " << max_wavelength << "\n";
  std::cout << "- spectral resolution: " << spectral_resolution << "\n";

  wavenumber_file_path = opacity_path + "wavenumber_full.dat";

  
  // std::getline(file, line);
  // file >> spectral_resolution >> line;
  // std::cout << "- Spectral resolution: " << spectral_resolution << "\n";

  // std::getline(file, line);

  // file >> input >> line;
  // std::cout << "- Opacity data folder: " << input << "\n";
  // cross_section_file_path = input;
  
  // if (cross_section_file_path.back() != '/')
  //   cross_section_file_path.append("/");
  
  // wavenumber_file_path = cross_section_file_path + "wavenumber_full.dat";


  // std::getline(file, line);

  
  
  std::cout << "\n";


  return true;
}


} 
