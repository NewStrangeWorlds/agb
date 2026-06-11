
#ifndef _interpolation_h
#define _interpolation_h

#include <vector>
#include <cmath>
#include <algorithm>

namespace agb { namespace aux {


//Linear interpolation of (x, y) at a single query point xq.
//x must be sorted ascending; values are clamped at the ends.
inline double interpolateLinearPoint(
  const std::vector<double>& x,
  const std::vector<double>& y,
  const double xq)
{
  if (xq <= x.front()) return y.front();
  if (xq >= x.back())  return y.back();

  const auto it = std::upper_bound(x.begin(), x.end(), xq);
  const size_t i = static_cast<size_t>(it - x.begin()) - 1;   //x[i] <= xq < x[i+1]

  const double t = (xq - x[i]) / (x[i+1] - x[i]);
  return y[i] + t * (y[i+1] - y[i]);
}


//Linear interpolation of (x_old, y_old) onto x_new. With log_interp the
//interpolation is done on log10(y) (for strictly positive quantities such as
//density, pressure, opacity). x_old must be sorted ascending.
inline std::vector<double> interpolate(
  const std::vector<double>& x_old,
  const std::vector<double>& y_old,
  const std::vector<double>& x_new,
  const bool log_interp = false)
{
  std::vector<double> y_new(x_new.size(), 0.);

  if (log_interp)
  {
    std::vector<double> ly(y_old.size());
    for (size_t i=0; i<y_old.size(); ++i)
      ly[i] = std::log10(y_old[i] > 0. ? y_old[i] : 1e-300);

    for (size_t i=0; i<x_new.size(); ++i)
      y_new[i] = std::pow(10.0, interpolateLinearPoint(x_old, ly, x_new[i]));
  }
  else
  {
    for (size_t i=0; i<x_new.size(); ++i)
      y_new[i] = interpolateLinearPoint(x_old, y_old, x_new[i]);
  }

  return y_new;
}


}}

#endif
