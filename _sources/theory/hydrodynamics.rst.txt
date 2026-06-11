=========================
Wind hydrodynamics
=========================

*Module:* ``src/hydrodynamics/`` — :cpp:class:`agb::Hydrodynamics`, with the
optional global-Newton :cpp:class:`agb::StructureSolver`
(``structure_solver.cpp``, ``henyey_solver.cpp``).

The hydrodynamics solves the **stationary, spherically symmetric momentum
equation** for the wind velocity, the density, the location of the sonic point,
and the mass-loss rate, given the flux-mean extinction from the radiative
transfer.

The radiative acceleration
==========================

The radiation force is parametrised by the ratio of radiative to gravitational
acceleration (:cpp:func:`agb::Hydrodynamics::calcAlpha`):

.. math::
   :label: alpha-hydro

   \alpha(r) = \frac{\chi_F(r)}{\rho(r)}\,
   \frac{L_\star}{4\pi c\,G M_\star} ,

with :math:`\chi_F` the **flux-mean extinction** coefficient supplied by the RT
module.  When :math:`\alpha>1` the radiation overcomes gravity and a wind is
driven.  In the cool dust-forming region :math:`\chi_F` is dominated by the grain
opacity, so :math:`\alpha` rises sharply across the dust front — this is the
engine of the outflow.

The Melia :math:`\Phi` transform
================================

The isothermal-sound-speed momentum equation has a **critical (sonic) point**
where the wind speed equals the isothermal sound speed
:math:`v=c_T`, at which the naïve equation is singular.  Following **Melia
(1988)**, the code works with the variable

.. math::
   :label: phi-def

   \Phi(r) = \tfrac12\!\left(v + \frac{c_T^2}{v}\right) ,

which is **regular** through the sonic point (its minimum is exactly
:math:`\Phi=c_T` at :math:`v=c_T`).  The velocity is recovered from :math:`\Phi`
by inverting :eq:`phi-def`,

.. math::

   v = \Phi \mp \sqrt{\Phi^2 - c_T^2} ,

taking the subsonic (``−``, :math:`v<c_T`) branch **inside** the critical point
and the supersonic (``+``) branch **outside** it
(:cpp:func:`agb::Hydrodynamics::windVelocity`).  :math:`\Phi` obeys the
first-order equation

.. math::
   :label: phi-ode

   \frac{\mathrm{d}\Phi}{\mathrm{d}r} =
   -\frac{1}{2v}\,\frac{G M_\star}{r^2}\,(1-\alpha)
   + \frac{c_T^2}{v\,r} ,

integrated with an explicit Euler scheme.

Critical point and the mass-loss eigenvalue
===========================================

At the sonic point the regularity condition (numerator and denominator of the
wind equation vanishing together) reads

.. math::
   :label: critical

   \frac{G M_\star}{r_c^2}\,(1-\alpha)
   - \frac{2 c_T^2}{r_c}
   + \frac{\mathrm{d}c_T^2}{\mathrm{d}r}\bigg|_{r_c} = 0 .

:cpp:func:`findCriticalPoint` locates the interior node where this condition is
met.  Solving the same regularity condition for the mass-loss rate (through
:math:`\alpha\propto\chi_F/\rho` and :math:`\rho=\dot M/4\pi r^2 v`) gives the
**eigenvalue mass-loss rate** (:cpp:func:`calcMassLossRate`):

.. math::
   :label: mdot-eigen

   \dot{M} = \frac{L_\star\,\chi_F\,c_T / c}
   {\dfrac{G M_\star}{r_c^2} - \dfrac{2 c_T^2}{r_c}
    + \dfrac{\mathrm{d}c_T^2}{\mathrm{d}r}\bigg|_{r_c}} .


Shooting solution (default)
===========================

The default solver (:cpp:func:`calcWindVelocity`, ``use_henyey_solver = false``)
proceeds:

#. compute :math:`\alpha` and the sound-speed derivative;
#. find the interior critical point :eq:`critical`;
#. integrate :math:`\Phi` :eq:`phi-ode` **outward** from the inner boundary to
   the critical point and **inward** from the outer boundary to the critical
   point, shooting on the boundary velocity so the two branches match
   (:cpp:func:`calcPhi`, ``integratePhiOutward`` / ``integratePhiInward``);
#. recover :math:`v(r)`, the eigenvalue :math:`\dot{M}` :eq:`mdot-eigen`, and the
   density :math:`\rho=\dot M/4\pi r^2 v`.

The velocity update is **under-relaxed in log space** with an adaptive factor,
and a solve with no interior critical point (or a non-finite :math:`\dot{M}`) is
**rejected**, leaving the structure frozen (``last_solve_rejected``).  A rejected
solve must not be read as convergence by the outer loop.

.. note::

   Shooting is known to be "slow but stable" (Dominik 1990).  With the carbon
   depletion of :doc:`dust` in place, the shooting solver plus the eigenvalue
   :math:`\dot{M}` converges the IRC+10216 coupling cleanly in ~12 chemistry–
   hydrodynamics iterations, with :math:`\alpha` settling to physical values
   (≈1.3–5) and reproducing the observed mass-loss rate.  **This is the
   recommended configuration.**

The grey (Lucy) starting model
==============================

When ``starting_model = "grey"`` (:cpp:func:`agb::Atmosphere::buildGreyStart`),
the model builds a logarithmic radial grid and a hydrostatic **Lucy spherical-
grey** structure as a seed:

.. math::

   T^4(r) = T_\mathrm{eff}^4\!\left[W(r) + \tfrac34\tau(r)\right],
   \qquad
   W(r) = \tfrac12\!\left(1 - \sqrt{1 - (R_\star/r)^2}\right) ,

with the grey optical depth :math:`\tau` integrated inward from a grey gas
opacity (``[grid] gas_opacity``).  This provides a structure inside the
convergence radius of the Newton-type solvers; it is also what the Henyey
bootstrap uses.

Optional: the Henyey global-Newton structure solver
===================================================

*Class:* :cpp:class:`agb::StructureSolver` (``use_henyey_solver = true``).

As an alternative to shooting, a **Henyey-type global Newton–Raphson** solver
relaxes the discretised boundary-value problem (the :math:`\Phi` equation, and
optionally the Gail & Sedlmayr dust moments :math:`K_0\dots K_5`) simultaneously
on all nodes, with the mass-loss rate as an eigenvalue closed by the
critical-point regularity row.  The Jacobian is obtained by **automatic
differentiation** (CppAD), validated against finite differences in the
``--selftest`` path.  Key design points:

* The residual is templated to evaluate in both ``double`` and
  ``CppAD::AD<double>`` (so the same code gives the value and the exact
  Jacobian).
* The radiative acceleration :math:`\alpha` is computed **inside** the residual
  from the frozen flux-mean extinction, so :math:`\dot{M}` is a true eigenvalue
  (``s = ln Mdot`` plus the regularity row).
* **Per-block residual/variable scaling** is essential: the unknowns span ~20
  orders of magnitude (:math:`\Phi\sim10^5` vs. :math:`K_j\sim10^{-13}`), so each
  block is scaled to :math:`O(1)` before the dense Eigen solve.  Without scaling
  the Newton fails.
* A **grey bootstrap** provides the seed; the relaxation-Picard intermediate
  structure is *not* a valid seed (its transient :math:`\chi_F` is garbage).
* Outer **log-space :math:`\chi_F` (α) damping** between solves stabilises the
  critical-point location.

.. important::

   The Henyey solver is **not required** for a converged model and is **off by
   default**.  Extensive experiments showed that the exact global Newton with a
   freely relocating critical point is ill-conditioned for the *marginal*
   dust-driven wind (:math:`\alpha\approx1` at the sonic point): coupling the
   dust opacity back into :math:`\alpha` inside the Newton (the "Stage 2"
   feedback) drives :math:`\alpha` super-Eddington on any moment excursion,
   losing the critical point.  The looser shooting solver with adaptive velocity
   relaxation "muddles through" and is more robust here.  Winters and Dominik
   reach the same conclusion in prescribing :math:`\dot{M}` and freezing
   :math:`\alpha` in the Newton, damping it only in the outer loop.

  
Prescribed vs. eigenvalue mass-loss rate
========================================

Two setups are possible:

* **Eigenvalue** :math:`\dot{M}` (default) — :math:`\dot{M}` follows from the
  critical-point regularity :eq:`mdot-eigen`; :math:`R_\star`/:math:`T_\star`
  float.  With carbon depletion this is stable and self-consistent for
  IRC+10216.
* **Prescribed** :math:`\dot{M}` (``[star] mass_loss_rate``) — Winters/Dominik
  fix :math:`\dot{M}` as a parameter.  A consistent prescribed value yields a
  regular transonic solution; an inconsistent one (e.g. far from the
  self-consistent eigenvalue) has no transonic solution and the solve is
  rejected.
