==================================
Dust formation (Gail & Sedlmayr)
==================================

*Module:* ``src/dust/`` — :cpp:class:`agb::GailSedlmayrDust` (derived from the
abstract :cpp:class:`agb::DustSpecies`), with the thermochemistry in
``gail_sedlmayr_thermo_functions.cpp``.  An :cpp:class:`agb::AnalyticDust`
prescription is also available for testing.

Carbon condenses out of the cooling, expanding gas into small amorphous-carbon
grains.  The code follows the classical **Gail & Sedlmayr moment method**:
homogeneous nucleation of seed clusters followed by grain growth, tracked through
the moments of the grain-size distribution.

The moment equations
====================

Let :math:`K_j(r)` be the :math:`j`-th moment of the grain-size distribution
(number of monomers per grain to the power :math:`j/3`), normalised per hydrogen
nucleus.  In the stationary, spherically symmetric wind they obey the coupled
system

.. math::
   :label: dust-moments

   v\,\frac{\mathrm{d}K_0}{\mathrm{d}r} &= \frac{J_\star}{n_{\langle H\rangle}}, \\
   v\,\frac{\mathrm{d}K_j}{\mathrm{d}r} &=
   \frac{j}{3}\,\frac{1}{\tau}\,K_{j-1}
   + N_\ell^{\,j/3}\,\frac{J_\star}{n_{\langle H\rangle}},
   \qquad j = 1,\dots,5 ,

where

* :math:`J_\star` is the **stationary nucleation rate** (seed clusters per unit
  volume and time),
* :math:`1/\tau` is the **net monomer growth rate** (so :math:`\tau` is the
  grain-growth timescale),
* :math:`N_\ell` (``minimum_monomer_number`` = 1000) is the lower cluster-size
  bound of the integration, and
* :math:`n_{\langle H\rangle}` is the hydrogen number density.

The code integrates **six** moments (:math:`K_0\dots K_5`).  Physically:
:math:`K_0` is the grain number density, and :math:`K_3` is the **condensed
volume per hydrogen nucleus** — i.e. the carbon locked into dust.  The mean grain
radius follows from :math:`\langle a\rangle = a_0\,(K_3/K_0)^{1/3}` with the
monomer radius :math:`a_0=1.28\times10^{-8}\,\mathrm{cm}`.

Nucleation rate
===============

The nucleation rate is built from classical nucleation theory
(:cpp:func:`nucleationRate`):

.. math::

   J_\star = c_0\,\beta\,A_\star\,Z ,

assembled from

* the **supersaturation ratio** :math:`S(T,n_C)` and its logarithm
  (``saturationRatio``, ``saturationVapourPressure``); nucleation is shut off
  (:math:`J_\star\to0`) where :math:`\ln S\le0`;
* the **critical cluster size** :math:`n_\star` (``criticalClusterSize``);
* the **free energy of cluster formation** :math:`\Delta F`
  (``freeEnergyOfFormation``), using the surface tension and the temperature
  scale :math:`\theta_\infty=\sigma A_1/k_B`;
* the **monomer growth rate** onto the critical cluster :math:`\beta`
  (``monomerGrowthRate``);
* the **Zeldovich factor** :math:`Z` (``zeldovichFactor``);
* the **equilibrium cluster distribution** :math:`c_0`
  (``equilibriumClusterDistribution``); and
* the critical-cluster surface area
  :math:`A_\star = A_1 n_\star^{2/3}`.

Growth is fed by collisions of the grains with the carbon-bearing species
:math:`\mathrm{C}, \mathrm{C_2}, \mathrm{C_2H}, \mathrm{C_2H_2}`, weighted by
their sticking coefficients and thermal velocities (:cpp:func:`growthRate`).  The
material constants (graphite/amorphous-carbon: surface tension
:math:`1400\,\mathrm{erg\,cm^{-2}}`, critical saturation ratio 3, sticking
coefficients 0.37/0.34) are the Gail & Sedlmayr (1984) values.

.. _dust-depletion:

Carbon depletion in a single causal sweep
=========================================

*The decisive physics:* the moment equations :eq:`dust-moments` are integrated
**outward in one coupled RK4 sweep** (:cpp:func:`agb::GailSedlmayrDust::calcDistribution`),
and at each radius the growth species that feed **both** :math:`1/\tau` and
:math:`J_\star` are throttled by the **running degree of condensation**

.. math::

   f_c = \frac{K_3}{(\varepsilon_C-\varepsilon_O)\,n_{\langle H\rangle}/n_{\langle H\rangle}}
        \;\equiv\; \frac{K_3}{c_\mathrm{cond}} ,
   \qquad
   \text{throttle factor } = 1 - f_c .

Both the (linear) growth rate and the (strongly nonlinear) nucleation rate are
evaluated from the **depleted** densities :math:`n\cdot(1-f_c)`, so condensation
**self-limits** as carbon is consumed.  Because the sweep marches outward and
:math:`K_3` only grows, :math:`f_c` rises monotonically, the throttle decreases
smoothly, and there is **no oscillation** and no outer fixed-point iteration.

.. important::

   This differential, in-sweep depletion (Gail & Sedlmayr book §14.3; Winters
   thesis Eq. 5.4) was the missing physics that earlier caused the dust opacity —
   and hence the radiative acceleration :math:`\alpha` — to grow without bound,
   collapsing the wind into a spurious "converged" state.  With it in place the
   self-consistent eigenvalue mass-loss rate for IRC+10216 settles at
   :math:`\dot{M}\approx7.8\times10^{-5}\,M_\odot\,\mathrm{yr}^{-1}` (a ~3 %
   match to the observed value), with a degree of condensation
   :math:`f_c\approx0.13`.

The integration detail: gas-phase densities, the growth timescale and the
hydrogen density are interpolated to the RK4 midpoints (logarithmically for
positive-definite quantities, linearly for temperature/velocity); the moments are
seeded to zero at the inner boundary (no dust there).  After the sweep the
moments are de-normalised by :math:`n_{\langle H\rangle}` and floored to avoid
underflow in the downstream RT.

Outputs and dust opacity
========================

For each layer the module stores the (depleted) nucleation rate, growth
timescale, grain number density, mean grain radius, and the degree of
condensation, written to the dust output file.  The grain optical properties are
then obtained from **Mie theory** (LX-MIE) using the material's complex
refractive index (:doc:`opacities`); :cpp:func:`calcTransportCoefficients`
returns the per-layer dust absorption and scattering coefficients that enter the
total opacity.

Coupling to the hydrodynamics
=============================

The frozen kernels needed by the Henyey structure solver are exposed through the
:cpp:class:`agb::DustSpecies` accessors :cpp:func:`nucleationRate`,
:cpp:func:`growthTimescale` and :cpp:func:`secondMoment`, handed to the
hydrodynamics via :cpp:func:`agb::Hydrodynamics::setDustState` (a no-op for the
default shooting path).
