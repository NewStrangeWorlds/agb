 
#ifndef _temperature_correction_h
#define _temperature_correction_h

#include <vector>
#include <string>
#include <iostream>

namespace agb {

//forward declaration
class ModelConfig;
class SpectralGrid;
class RadiationField;


class TemperatureCorrection{
  public:
    TemperatureCorrection(
      ModelConfig* config_,
      SpectralGrid* spectral_grid_,
      std::vector<RadiationField>* radiation_field_)
      : config(config_)
      , spectral_grid(spectral_grid_)
      , radiation_field(radiation_field_) {}
    ~TemperatureCorrection() {}

    std::vector<double> calculate(
      std::vector<double>& temperature,
      std::vector<std::vector<double>>& absorption_coeff,
      std::vector<std::vector<double>>& scattering_coeff);
  protected:
    ModelConfig* config;
    SpectralGrid* spectral_grid;
    std::vector<RadiationField>* radiation_field;
};



}

#endif