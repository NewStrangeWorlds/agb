#ifndef _aux_functions_h
#define _aux_functions_h

#include <vector>


namespace agb{ namespace aux{


double planckFunctionWavenumber(const double temperature, const double wavenumber);
double planckFunctionWavelength(const double temperature, const double wavelength);
double planckFunctionDerivWavenumber(const double temperature, const double wavenumber);
double planckFunctionDerivWavelength(const double temperature, const double wavelength);
double linearInterpolation(const double x1, const double x2, const double y1, const double y2, const double x);
double voigtProfile(const double x, const double gaussian_width, const double lorentz_width);

double normalDistribution(const double mu, const double sigma, const double x);
double normalDistribution(const double sigma, const double x);

std::vector<double> interpolateToWavenumberGrid(const std::vector<double>& wavenumber_data, const std::vector<double>& data,
                                                const std::vector<double>& wavenumber_grid,
                                                const double outside_range_value,
                                                const bool interpolate_log);

}}



#endif
