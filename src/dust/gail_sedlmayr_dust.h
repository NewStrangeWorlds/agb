 
#ifndef _gail_sedlmayr_dust_h
#define _gail_sedlmayr_dust__h

#include <vector>
#include <string>
#include <iostream>
#include <complex>
#include <cmath>

#include "dust_species.h"

#include "../additional/physical_const.h"


namespace agb {

//forward declaration
class ModelConfig;
class SpectralGrid;
class Atmosphere;


class GailSedlmayrDust: public DustSpecies{
  public:
    GailSedlmayrDust(
      ModelConfig* config_,
      SpectralGrid* spectral_grid_,
      Atmosphere* atmosphere_);
    ~GailSedlmayrDust() {}

    void calcDistribution();
    void calcTransportCoefficients(
      const size_t radius_idx,
      std::vector<double>& absorption_coeff,
      std::vector<double>& scattering_coeff);
    std::vector<double> degreeOfCondensation(
      const double carbon_abundance);
    void saveOutput(const std::string file_path);
  private:
    //fixed values for graphite from Gail&Sedlmayr (1984)
    const double monomer_radius = 1.28e-8;
    const double monomer_surface_area = 
      4.*constants::pi * monomer_radius*monomer_radius; //cm^2
    const int minimum_monomer_number = 1000;
    const double critical_saturation_ratio = 3.0;
    const double n_l = 5.;
    const double surface_tension = 1400; // erg/ cm^2
    const double theta_infinity = 
      surface_tension * monomer_surface_area / constants::boltzmann_k; //in K
    const std::vector<double> sticking_coeff{0.37, 0.34};

    const double mass_c = 12.01 * constants::amu;
    const double mass_c2 = 2*12.01 * constants::amu;
    const double mass_c2h = 25.029 * constants::amu;
    const double mass_c2h2 = 26.04 * constants::amu;
    
    const size_t nb_moments = 6;
    std::vector<std::vector<double>> dust_moments;

    std::vector<double> nucleation_rate;
    std::vector<double> growth_rate;

    double nucleationRate(
      const double temperature,
      const double number_density_c,
      const double number_density_c2,
      const double number_density_c2h,
      const double number_density_c2h2);
    double growthRate(
      const double temperature,
      const double number_density_c,
      const double number_density_c2,
      const double number_density_c2h,
      const double number_density_c2h2);

    std::vector<double> dustMomentZero();
    std::vector<double> dustMoment(
      const int order,
      std::vector<double>& moment_prev,
      const int n_lower);

    double saturationVapourPressure(
      const double temperature);
    double saturationRatio(
      const double temperature, 
      const double number_density_carbon);
    double criticalClusterSize(
      const double temperature,
      const double ln_saturation_ratio);
    double freeEnergyOfFormation(
      const double temperature,
      const double ln_saturation_ratio,
      const double critical_cluster_size);
    double monomerGrowthRate(
      const double temperature,
      const double number_density_c,
      const double number_density_c2,
      const double number_density_c2h,
      const double number_density_c2h2);
    double zeldovichFactor(
      const double temperature,
      const double critical_cluster_size,
      const double ln_saturation_ratio);
    double equilibriumClusterDistribution(
      const double temperature,
      const double free_energy,
      const double monomer_number_density,
      const double ln_saturation_ratio);
};



}

#endif