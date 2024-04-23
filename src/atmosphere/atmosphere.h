 
#ifndef _atmosphere_h
#define _atmosphere_h

#include <vector>
#include <string>
#include <iostream>

#include "../chemistry/chem_species.h"


namespace agb {

//forward declaration
class ModelConfig;


class Atmosphere{
  public:
    Atmosphere(ModelConfig* config_);
    ~Atmosphere() {}

    size_t nb_grid_points = 0;
    const size_t nb_chemistry_species = constants::species_data.size();

    std::vector<double> radius_grid;
    std::vector<double> pressure;
    std::vector<double> pressure_bar;
    std::vector<double> radius;
    std::vector<double> temperature_gas;
    std::vector<double> temperature_dust;
    std::vector<double> mass_density;
    std::vector<double> velocity;

    std::vector<std::vector<double>> number_densities;
    std::vector<double> mean_molecuar_weight;

    std::vector<std::vector<double>> absorption_coeff;
    std::vector<std::vector<double>> scattering_coeff;

    std::vector<std::vector<double>> eddington_flux;
    std::vector<std::vector<double>> flux;
    std::vector<std::vector<double>> mean_intensity;
    std::vector<std::vector<double>> eddington_k;

    void equationOfState();
  protected:
    ModelConfig* config;

    std::vector<double> cgsToBar(const std::vector<double> pressure_data);

    void readStructure(const std::string file_path);
    void writeStructure(const std::string file_path);
};



}

#endif