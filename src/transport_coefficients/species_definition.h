
#ifndef SPECIES_DEFINITION_H
#define SPECIES_DEFINITION_H

#include "opacity_species.h"

#include <vector>
#include <iostream>

#include "../chemistry/chem_species.h"
#include "../additional/physical_const.h"
#include "../spectral_grid/spectral_grid.h"


namespace agb{


class ModelConfig;
class SpectralGrid;


class GasGeneric : public OpacitySpecies {
  public:
    GasGeneric(
      ModelConfig* config_ptr, 
      SpectralGrid* spectral_grid_ptr, 
      const unsigned int index, 
      const std::string name, 
      const std::string folder) 
        : OpacitySpecies(index, name, folder) 
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          init();
        }
    GasGeneric(
      ModelConfig* config_ptr, 
      SpectralGrid* spectral_grid_ptr, 
      const unsigned int index, 
      const std::string name, 
      const std::string folder, 
      const size_t reference_species) 
        : OpacitySpecies(index, name, folder) 
        {
          config = config_ptr; spectral_grid = spectral_grid_ptr; 
          pressure_reference_species = reference_species; 
          init();
        }
    GasGeneric(
      ModelConfig* config_ptr, 
      SpectralGrid* spectral_grid_ptr, 
      const unsigned int index, 
      const std::string name, 
      const std::string folder, 
      const std::vector<size_t> cia_collision_species) 
        : OpacitySpecies(index, name, folder)
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          cia_collision_partner = cia_collision_species; 
          init();
        }
    virtual ~GasGeneric() {}
};



class GasH : public OpacitySpecies {
  public:
    GasH(ModelConfig* config_ptr, SpectralGrid* spectral_grid_ptr, const std::string folder) 
        : OpacitySpecies(_H, "H", folder)
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          init();
        }
    GasH(ModelConfig* config_ptr, SpectralGrid* spectral_grid_ptr) 
        : OpacitySpecies(_H, "H", "")
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          init();
        }
    virtual ~GasH() {}
  protected:
    virtual bool calcRalyleighCrossSections(std::vector<double>& cross_sections);
};



class GasHm : public OpacitySpecies {
  public:
    GasHm(ModelConfig* config_ptr, SpectralGrid* spectral_grid_ptr) 
        : OpacitySpecies(_Hm, "H-", "")
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          init();
        }
    virtual ~GasHm() {}
  protected:
    virtual bool calcContinuumAbsorption(const double temperature, const std::vector<double>& number_densities, std::vector<double>& absorption_coeff);
  private:
    std::vector<double> boundFreeAbsorption(const double temperature);
    std::vector<double> freeFreeAbsorption(const double temperature);
};


class GasH2m : public OpacitySpecies {
  public:
    GasH2m(ModelConfig* config_ptr, SpectralGrid* spectral_grid_ptr) 
         : OpacitySpecies(_H2, "H2-ff", "")
         {config = config_ptr; spectral_grid = spectral_grid_ptr; init();}
    virtual ~GasH2m() {}
  protected:
    virtual bool calcContinuumAbsorption(const double temperature, const std::vector<double>& number_densities, std::vector<double>& absorption_coeff);
  private:
    const std::vector<double> fit_coeff = {-0.9301, 2.031, 1.321, -0.03775, -0.1817, -1.098, 0.1173, 0.2116, -0.0425, 0.515, -0.26, 0.02659, -0.02431, 0.01914, -0.1201, 0.1463, -0.08034, 0.02281, -0.003203, -0.001883, 0.0109};
    const double min_temperature = 1400.11;
    const double max_theta =  5040.4/min_temperature;
    const double max_temperature = 10080.8;
    const double min_theta =  5040.4/max_temperature;
    const double min_lambda = 0.3505; //in micron
    const double max_lambda = 15.1883; //in micron

    std::vector<double> freeFreeAbsorption(const double temperature);
};



class GasH2 : public OpacitySpecies {
  public:
    GasH2(ModelConfig* config_ptr, SpectralGrid* spectral_grid_ptr, const std::string folder) 
        : OpacitySpecies(_H2, "H2", folder)
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          init();
        }
    GasH2(ModelConfig* config_ptr, SpectralGrid* spectral_grid_ptr) 
        : OpacitySpecies(_H2, "H2", "")
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          init();
        }
    virtual ~GasH2() {}
  protected:
    virtual bool calcRalyleighCrossSections(std::vector<double>& cross_sections);
};


class GasHe : public OpacitySpecies {
  public:
    GasHe(ModelConfig* config_ptr, SpectralGrid* spectral_grid_ptr, const std::string folder)
        : OpacitySpecies(_He, "He", folder)
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          init();
        }
    GasHe(ModelConfig* config_ptr, SpectralGrid* spectral_grid_ptr) 
        : OpacitySpecies(_He, "He", "")
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          init();
        }
    virtual ~GasHe() {}
  protected:
    virtual bool calcRalyleighCrossSections(std::vector<double>& cross_sections);
};


}

#endif
