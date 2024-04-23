#ifndef _spectral_grid_h
#define _spectral_grid_h

#include <vector>
#include <string>


namespace agb {

//forward declaration
class ModelConfig;


class SpectralGrid{
  public:
    SpectralGrid (ModelConfig* config);

    std::vector<double> wavenumber_list;      //wavenumber list used to calculate the high-res spectra
    std::vector<double> wavelength_list;      //wavelength list used to calculate the high-res spectra in units of microns
    std::vector<double> wavelength_list_cm;      //wavelength list used to calculate the high-res spectra in units of cm

    std::vector<double> wavelengthToWavenumber(
      const std::vector<double>& wavelengths);
    std::vector<double> wavenumberToWavelength(
      const std::vector<double>& wavenumbers);

    double wavelengthToWavenumber(const double wavelength)
     {return 1.0/wavelength * 1e4;}
    double wavenumberToWavelength(const double wavenumber)
     {return 1.0/wavenumber * 1e4;}

    size_t nbSpectralPointsFull() {
      return nb_spectral_points_full;}
    size_t nbSpectralPoints() {
      return nb_spectral_points;}

    std::vector<size_t> spectralIndexList() {
      return index_list;}

    size_t findClosestIndex(
      const double search_value,
      std::vector<double>& data,
      std::vector<double>::iterator it_start);

    std::vector<double> interpolateToWavenumberGrid(
      const std::vector<double>& data_x,
      const std::vector<double>& data_y,
      const bool log_interpolation);
    std::vector<double> interpolateToWavelengthGrid(
      const std::vector<double>& data_x,
      const std::vector<double>& data_y,
      const bool log_interpolation);
    std::vector<double> interpolateToWavelengthGrid(
      const std::vector<double>& data_x,
      const std::vector<double>& data_y,
      const std::vector<double>& new_x,
      const bool log_interpolation);
  private:
    ModelConfig* config;

    std::vector<double> wavenumber_list_full; //the full, global wavenumber list, the opacities have been calculated at
    std::vector<double> wavelength_list_full; //the full, global wavelength list, the opacities have been calculated at

    std::vector<size_t> index_list;

    size_t nb_spectral_points_full;           //number of points in the global wavenumber list
    size_t nb_spectral_points;                //number of points in the spectral grid

    void loadWavenumberList();

    std::vector<double> createConstantResolutionGrid(
      const double min_wavelength,
      const double max_wavelength,
      const double spectral_resolution);

    void createHighResGrid();

    size_t findClosestIndexDesc(
      const double search_value,
      std::vector<double>& data,
      std::vector<double>::iterator it_start);

    size_t findClosestIndexAsc(
      const double search_value,
      std::vector<double>& data,
      std::vector<double>::iterator it_start);
};



}


#endif
