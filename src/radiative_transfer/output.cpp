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

#include "../additional/aux_functions.h"
#include "../spectral_grid/spectral_grid.h"
#include "../config/config.h"
#include "../atmosphere/atmosphere.h"
#include "../additional/quadrature.h"
#include "../additional/physical_const.h"

namespace agb{

void RadiativeTransfer::saveSpectrum(
  const std::string file_path)
{
  std::fstream file;
  file.open(file_path.c_str(), std::ios::out);

  if (file.fail()) 
  {
    std::cout << "Couldn't open spectrum output file " << file_path << "\n";
    return;
  }

  file << std::setprecision(10) << std::scientific << "#Wavelength (micron)\tFlux (erg s-1 cm-2 micron-1)\t r2 F (cm2 erg s-1 cm-2 micron-1)\n";
  
  for (size_t i=0; i<nb_spectral_points; ++i)
  {  
    const double flux = radiation_field.back().eddington_flux[i] * 4 * constants::pi;
    const double radius2 = atmosphere->radius.back() * atmosphere->radius.back();

     file << spectral_grid->wavelength_list[i] << "\t" << flux << "\t" << flux*radius2 << "\n";
  }

  file.close();
}


}
