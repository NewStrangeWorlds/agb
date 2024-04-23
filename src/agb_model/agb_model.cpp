 
#include "agb_model.h"

#include <vector>
#include <string>



namespace agb{

AGBStarModel::AGBStarModel(const std::string folder)
 : config(folder)
 , spectral_grid(&config)
{
  std::cout << "spectral points " << spectral_grid.nbSpectralPoints() << "\n";

  for (size_t i=0; i<spectral_grid.nbSpectralPoints(); ++i)
  {
    std::cout << spectral_grid.wavelength_list[i] << "\t" << spectral_grid.wavenumber_list[i] << "\n";
  }
}



}