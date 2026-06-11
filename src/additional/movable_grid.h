
#ifndef _movable_grid_h
#define _movable_grid_h

#include <vector>
#include <cmath>

#include "interpolation.h"

namespace agb { namespace aux {


//Equidistribution of an adaptive radial grid (de Boor; Dorfi & Drury).
//Builds a dimensionless monitor function
//   w(r) = 1 + sum_k weight_k * |d ln q_k / d ln r|
//(smoothed with a few 1-2-1 passes), then returns target node radii such that the
//integral of w dr is equal across every cell - i.e. nodes cluster where any of the
//monitor quantities q_k varies rapidly in a relative sense (orders of magnitude).
//The number of nodes is preserved and the boundaries radius.front()/radius.back()
//are pinned. The caller is expected to under-relax the motion towards this target.
inline std::vector<double> equidistributedGrid(
  const std::vector<double>& radius,
  const std::vector<std::vector<double>>& quantities,
  const std::vector<double>& weights,
  const int smoothing_passes = 2,
  const double max_monitor = 20.0,
  const double rel_floor = 1e-8)
{
  const size_t n = radius.size();
  const double tiny = 1e-300;

  //Per-quantity floor RELATIVE to that quantity's own peak. This bounds each
  //quantity's usable dynamic range to ~|log10(rel_floor)| decades, so an absolute
  //numerical floor (e.g. nucleation rate clamped to 1e-100) cannot create a huge
  //spurious log-gradient that would pull every node onto one spike and starve the
  //rest of the grid (incl. the sonic-point region).
  std::vector<double> qfloor(quantities.size(), tiny);
  for (size_t k=0; k<quantities.size(); ++k)
  {
    double qmax = tiny;
    for (const double q : quantities[k]) qmax = std::max(qmax, q);
    qfloor[k] = qmax * rel_floor;
  }

  //dimensionless monitor at every node
  std::vector<double> w(n, 1.0);

  for (size_t i=1; i+1<n; ++i)
  {
    const double dlnr = std::abs(std::log(radius[i+1]) - std::log(radius[i-1])) + tiny;

    double extra = 0.0;
    for (size_t k=0; k<quantities.size(); ++k)
    {
      const double q1 = std::log(std::max(quantities[k][i+1], qfloor[k]));
      const double q0 = std::log(std::max(quantities[k][i-1], qfloor[k]));
      extra += weights[k] * std::abs(q1 - q0) / dlnr;
    }
    w[i] = 1.0 + extra;
  }
  w[0]   = w[1];
  w[n-1] = w[n-2];

  //smooth the monitor so the mesh tracks features, not point noise
  for (int s=0; s<smoothing_passes; ++s)
  {
    std::vector<double> ws = w;
    for (size_t i=1; i+1<n; ++i)
      ws[i] = 0.5*w[i] + 0.25*(w[i-1] + w[i+1]);
    w = ws;
  }

  //cap the monitor so the densest cells are at most max_monitor times finer than the
  //baseline - bounds the cell-size ratio and guarantees the smooth/inner regions
  //(and the sonic point) keep enough nodes
  for (double & wi : w) wi = std::min(wi, max_monitor);

  //Work in x = ln(r): the baseline (w = 1) then reproduces a LOGARITHMIC grid
  //(dense towards the inner boundary, as the model grid is built), and the monitor
  //only ADDS clustering on top of that. Integrating in linear r instead would make
  //the baseline uniform in r and slowly drain nodes out of the inner atmosphere.
  std::vector<double> x(n);
  for (size_t i=0; i<n; ++i) x[i] = std::log(radius[i]);

  //cumulative W(x) = integral of w dx (monotone increasing)
  std::vector<double> W(n, 0.0);
  for (size_t i=1; i<n; ++i)
    W[i] = W[i-1] + 0.5*(w[i] + w[i-1])*(x[i] - x[i-1]);

  //place nodes at equal increments of W (invert W -> x); pin the boundaries
  std::vector<double> target(n);
  target[0]   = radius[0];
  target[n-1] = radius[n-1];

  const double w_total = W[n-1];
  for (size_t j=1; j+1<n; ++j)
  {
    const double w_target = w_total * static_cast<double>(j) / static_cast<double>(n-1);
    target[j] = std::exp(interpolateLinearPoint(W, x, w_target));
  }

  return target;
}


}}

#endif
