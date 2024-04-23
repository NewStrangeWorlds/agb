 
#ifndef _atmosphere_h
#define _atmosphere_h

#include <vector>
#include <string>
#include <iostream>


namespace agb {


class Atmosphere{
  public:
    Atmosphere();
    ~Atmosphere() {}

    std::vector<double> pressure;
    std::vector<double> radius;
    std::vector<double> temperature_gas;
    std::vector<double> temperature_dust;
    std::vector<double> mass_density;
    std::vector<double> velocity;
  protected:
    void readStructure(const std::string file_path);
    void writeStructure(const std::string file_path);
};



}

#endif