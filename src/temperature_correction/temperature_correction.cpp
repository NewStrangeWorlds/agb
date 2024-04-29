
#include "temperature_correction.h" 

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


std::vector<double> TemperatureCorrection::calculate(
  std::vector<double>& temperature,
  std::vector<std::vector<double>>& absorption_coeff,
  std::vector<std::vector<double>>& scattering_coeff)
{


}



}