 
#include "agb_model.h"

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "../additional/aux_functions.h"
#include "../additional/physical_const.h"
#include "../additional/movable_grid.h"


namespace agb{


void AGBStarModel::applyMovableGrid()
{
  const size_t np = atmosphere.nb_grid_points;

  //monitor quantities on the current grid: flux-mean extinction, gas temperature,
  //and dust nucleation rate
  std::vector<double> chi_h(np, 0.);
  for (size_t i=0; i<np; ++i)
    chi_h[i] = radiative_transfer.radiation_field[i].fluxWeightedExtinction(
                 atmosphere.extinction_coeff[i]);

  std::vector<double> jstar = dust_species->nucleationRate();
  if (jstar.size() != np) jstar.assign(np, 1.0);

  const std::vector<std::vector<double>> quantities =
    { chi_h, atmosphere.temperature_gas, jstar, atmosphere.velocity };
  const std::vector<double> weights =
    { config.monitor_weight_opacity,
      config.monitor_weight_temperature,
      config.monitor_weight_nucleation,
      config.monitor_weight_velocity };

  const std::vector<double> target = aux::equidistributedGrid(
    atmosphere.radius, quantities, weights,
    config.monitor_smoothing_passes, config.monitor_max, config.monitor_rel_floor);

  //under-relax the node motion; pin the boundaries and keep the grid strictly monotone
  std::vector<double> new_radius(np);
  for (size_t i=0; i<np; ++i)
    new_radius[i] = atmosphere.radius[i]
                  + config.grid_relaxation * (target[i] - atmosphere.radius[i]);
  new_radius[0]    = atmosphere.radius[0];
  new_radius[np-1] = atmosphere.radius[np-1];
  for (size_t i=1; i<np; ++i)
    if (new_radius[i] <= new_radius[i-1])
      new_radius[i] = new_radius[i-1] * (1.0 + 1e-6);

  //remap the structure and rebuild the grid-dependent geometry / warm starts
  atmosphere.remapToGrid(new_radius);
  radiative_transfer.rebuildGeometry();
  hydrodynamics.resetWarmStart();

  //reset per-node temperature-iteration state (re-initialises lazily on the new grid)
  relaxation_gas.clear();    relaxation_dust.clear();
  prev_delta_b_gas.clear();  prev_delta_b_dust.clear();
  anderson_x_gas.clear();    anderson_f_gas.clear();
  anderson_x_dust.clear();   anderson_f_dust.clear();
  anderson_was_active = false;
  prev_max_rel_change = 1e30;

  std::cout << "Movable grid: regridded; inner/outer r/R* = "
            << atmosphere.radius_grid.front() << " / " << atmosphere.radius_grid.back() << "\n\n";
}



}