========
Overview
========

Scientific context
===================

Carbon-rich asymptotic giant branch (AGB) stars lose mass through a slow,
massive, dust-driven wind.  In the cool, extended atmosphere, carbon that is not
locked into carbon monoxide condenses into small amorphous-carbon grains.  These
grains are extremely opaque to the stellar radiation field; the momentum they
absorb is shared with the gas through collisions, and the resulting radiative
acceleration, once it exceeds gravity, drives a stationary outflow.  The wind is
therefore an intrinsically **coupled** problem: the temperature structure sets
the chemistry and the dust nucleation; the dust sets the opacity; the opacity
sets both the temperature (through radiative transfer) and the radiative
acceleration (through the wind dynamics); and the wind density feeds back on the
chemistry and dust.

This code solves the **stationary, spherically symmetric** version of that
problem.  All time derivatives are dropped; the model seeks the steady-state
structure :math:`\{T_\mathrm{gas}(r),\,T_\mathrm{dust}(r),\,\rho(r),\,v(r)\}`
together with the mass-loss rate :math:`\dot{M}` that is consistent with the
radiation field, the chemistry, and the dust.

The scheme follows the coupled-model approach of Winters et al. and Gail &
Sedlmayr, with the singularity-free wind formulation of Melia.

What the code computes
======================

Given a small number of stellar parameters (radius, mass, luminosity, C/O ratio,
and either a prescribed or eigenvalue mass-loss rate) and a starting structure,
the model produces:

* the **radial temperature structure** of gas and dust, in radiative
  equilibrium with the frequency-dependent radiation field;
* the **gas-phase chemical composition** at every radius (FastChem);
* the **dust distribution** — nucleation rate, grain-growth timescale, grain
  number density, mean grain radius, and the degree of carbon condensation —
  from the Gail & Sedlmayr moment method;
* the **wind structure** — velocity :math:`v(r)`, density :math:`\rho(r)`, the
  location of the sonic (critical) point, and the mass-loss rate
  :math:`\dot{M}`;
* the **emergent spectrum** and the radial run of the radiation moments.

The solution strategy in one paragraph
======================================

The model is solved by **fixed-point (Picard) iteration** over two nested
loops, as sketched in Winters' coupled-model diagram:

#. **Outer loop** (*global iteration*): update the chemistry, dust and wind
   structure at the current temperature, then perform one radiative-equilibrium
   temperature correction.  Repeat until the temperature has converged.

#. **Inner loop** (*chemistry–hydrodynamics iteration*): at fixed temperature,
   alternate the chemistry+dust calculation, the radiative transfer, and the
   wind solve until the radiative acceleration parameter :math:`\alpha` is
   self-consistent.

Each physics block is a separate module with a narrow interface; see
:doc:`implementation/architecture`.  The numerical robustness machinery
(adaptive under-relaxation, Anderson acceleration, the two-phase
Unsöld–Lucy → linearisation corrector, the optional movable grid) is layered on
top of this loop and is described in :doc:`theory/numerics`.

Status and validation
======================

* The radiative-transfer module reproduces the thesis scheme term by term; the
  Taylor (finite-difference) moment discretisation is unconditionally stable,
  whereas the cubic-spline discretisation loses positivity for optically thick
  cells and is provided mainly for cross-checks.
* With self-consistent carbon depletion in the dust moment sweep, the
  eigenvalue mass-loss rate for IRC+10216 settles at
  :math:`\dot{M}\approx 7.8\times10^{-5}\,M_\odot\,\mathrm{yr}^{-1}`, a ~3 %
  match to the observed value.
* The default wind solver is the Melia ``Φ`` shooting method with an eigenvalue
  :math:`\dot{M}`; the Henyey global-Newton structure solver is implemented and
  available but is *not* required for a converged model (see
  :doc:`theory/hydrodynamics`).
