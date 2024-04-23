 
#ifndef _config_h
#define _config_h

#include <string>

namespace agb {

struct ModelConfig {
  ModelConfig(const std::string model_folder);
  bool loadConfigFile(const std::string folder);
  
  std::string model_folder = "";

  double stellar_radius = 0;
  double stellar_mass = 0;
  double stellar_luminosity = 0;
  double c_o_ratio = 0;
  
  std::string opacity_path = "";
  std::string wavenumber_file_path = "";
  double spectral_resolution = 0;
  double min_wavelength = 0;
  double max_wavelength = 0;
};



}


#endif