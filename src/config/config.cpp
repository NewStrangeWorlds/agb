
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
  loadOutputConfigFile(model_folder);
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

  file >> fastchem_parameter_file;
  std::cout << "- FastChem parameter file: " << fastchem_parameter_file << "\n";


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);

  file >> starting_model_path;
  std::cout << "- starting model: " << starting_model_path << "\n";


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  file >> min_wavelength >> max_wavelength >> spectral_resolution;
  
  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  
  file >> refractive_index_file;


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  
  file >> opacity_path;

  if (opacity_path.back() != '/')
    opacity_path.append("/");
  
  std::cout << "- refractive indices: " << refractive_index_file << "\n";
  std::cout << "- min wavelength: " << min_wavelength << "\n";
  std::cout << "- max wavelength: " << max_wavelength << "\n";
  std::cout << "- spectral resolution: " << spectral_resolution << "\n";
  std::cout << "- opacity folder: " << opacity_path << "\n";

  wavenumber_file_path = opacity_path + "wavenumber_full.dat";
  
  std::getline(file, line);
  std::getline(file, line);
  readOpacityConfig(file);

  
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


void ModelConfig::readOpacityConfig(std::fstream& file)
{
  std::string line;
  std::getline(file, line);
  
  
  while(std::getline(file, line))
  {
    std::istringstream input(line);

    std::string species, folder;

    input >> species >> folder;
    
    if (species.length() > 0 && folder.length() > 0)
    {
      opacity_species_symbol.push_back(species);
      opacity_species_folder.push_back(folder);
    }
    
  }


  std::cout << "- Opacity species:\n";
  for (size_t i=0; i<opacity_species_symbol.size(); ++i)
    std::cout << "   species " << opacity_species_symbol[i] << "\t folder: " << opacity_species_folder[i] << "\n"; 
  
  
  std::cout << "\n";
}



bool ModelConfig::loadOutputConfigFile(const std::string folder_path)
{
  model_folder = folder_path;

  if (model_folder.back() != '/')
    model_folder.append("/");

  std::string file_path = model_folder;
  file_path.append("output.config");

  
  std::fstream file;
  file.open(file_path.c_str(), std::ios::in);

  if (file.fail()) 
  {
    std::cout << "Couldn't open model output options file " << file_path << "\n";
    
    return false;
  }

  
  std::cout << "\nParameters found in output.config:\n";

  std::string line;
  std::string input;
  
  std::getline(file, line);

  file >> output_spectrum_path;
  std::cout << "- output spectrum: " << output_spectrum_path << "\n";
  
  if (output_spectrum_path == "None" || output_spectrum_path == "none")
    output_spectrum_path = "";
  else
    output_spectrum_path = model_folder + output_spectrum_path;

  std::cout << "\n";

  return true;
} 

}
