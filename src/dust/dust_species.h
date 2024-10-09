 
#ifndef _dust_species_h
#define _dust_species_h

#include <vector>
#include <string>
#include <iostream>
#include <complex>
#include <cmath>


namespace agb {

//forward declaration
class ModelConfig;
class SpectralGrid;
class Atmosphere;


class DustSpecies{
  public:
    DustSpecies(
      ModelConfig* config_,
      SpectralGrid* spectral_grid_,
      Atmosphere* atmosphere_);
    ~DustSpecies() {}
    
    std::vector<double> number_density;
    std::vector<std::vector<double>> size_distribution;
    std::vector<std::vector<double>> particle_radius;

    virtual void calcDistribution() = 0;
    virtual std::vector<double> degreeOfCondensation(
      const double element_abundance) = 0;
    void calcTransportCoefficients(
      const size_t radius_idx,
      std::vector<double>& absorption_coeff,
      std::vector<double>& scattering_coeff);
    virtual void saveOutput(const std::string file_path) = 0;
  protected:
    ModelConfig* config;
    SpectralGrid* spectral_grid;
    Atmosphere* atmosphere;
    
    const size_t nb_grid_points = 0;
    const size_t nb_spectral_points = 0;
    const double monomer_radius = 1.28e-8;

    std::vector<std::complex<double>> refractive_index; 

    void readRefractiveIndexFile(const std::string file_path);
    void opticalProperties(
      double radius, 
      std::vector<double>& absorption_efficiency, 
      std::vector<double>& scattering_efficiency, 
      std::vector<double>& asymmetry_parameter);
};



}

#endif