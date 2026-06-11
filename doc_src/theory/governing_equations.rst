=========================================
The coupled problem and the solution loop
=========================================

The model determines the stationary, spherically symmetric structure of a
carbon-star wind.  Four physics blocks are coupled:

#. **Radiative transfer** — given the temperature structure and the opacities,
   solve for the radiation field :math:`J_\nu(r)`, the flux :math:`H_\nu(r)`, and
   the variable Eddington factors.
#. **Radiative equilibrium** — given the radiation field and the opacities,
   correct the gas and dust temperatures so that each is in local radiative
   equilibrium (energy balance) and the integrated flux carries the stellar
   luminosity.
#. **Chemistry and dust** — given the temperature and density, compute the
   equilibrium gas composition (FastChem) and the carbon-dust nucleation,
   growth and condensation (Gail & Sedlmayr moment method), which set the dust
   opacity.
#. **Hydrodynamics** — given the (gas + dust) opacity and the radiation field,
   solve the stationary wind equation for the velocity, density, the sonic
   point, and the mass-loss rate.

The coupling closes a loop: temperature → chemistry/dust → opacity → radiation
field → temperature, and opacity + radiation → wind → density → chemistry/dust.

The two physical structure equations
====================================

**Stationary continuity** fixes the density once the velocity and mass-loss rate
are known:

.. math::

   \dot{M} = 4\pi r^2\,\rho(r)\,v(r)
   \quad\Longrightarrow\quad
   \rho(r) = \frac{\dot{M}}{4\pi r^2 v(r)} .

**Stationary momentum** (the wind equation) balances inertia, pressure gravity
and the radiative force.  Writing the radiative acceleration as a fraction
:math:`\alpha` of gravity,

.. math::
   :label: alpha-def

   \alpha(r) = \frac{\chi_F(r)\,L_\star}
   {4\pi c\,G M_\star\,\rho(r)} ,

where :math:`\chi_F` is the flux-mean extinction coefficient, the momentum
equation for an isothermal-sound-speed gas becomes a transonic problem with a
critical (sonic) point.  The code removes the sonic-point singularity with the
**Melia** :math:`\Phi` transform; see :doc:`hydrodynamics`.

Energy: radiative equilibrium
=============================

Because the wind is optically thin to thick and radiation dominates the energy
budget, the gas and dust temperatures are fixed by **radiative equilibrium**
rather than by an energy equation with advection.  The condition is local energy
balance for each absorbing component :math:`s\in\{\text{gas},\text{dust}\}`:

.. math::
   :label: re-balance

   \int \kappa_{\nu,s}\,\bigl(J_\nu - B_\nu(T_s)\bigr)\,\mathrm{d}\nu = 0 ,

together with global flux constancy
:math:`r^2 H_\mathrm{int}(r) = L_\star/(16\pi^2)`.  The temperature corrector
(:doc:`temperature_correction`) drives :eq:`re-balance` and the flux constraint
to zero.

The solution loop
=================

The driver is :cpp:func:`agb::AGBStarModel::calcModel`.  The structure is two
nested fixed-point iterations:

.. code-block:: text

   for it in 0 .. nb_temperature_iter:            # OUTER (global) iteration
       chemistryHydroIteration():                 #   INNER iteration
           repeat up to nb_hydrodynamics_iter:
               chemistryDustIteration()           #     FastChem + Gail-Sedlmayr dust
               radiativeTransfer()                #     opacities + RT solve
               hydrodynamics.calcWindVelocity()   #     wind + Mdot
               until max |Δα| < hydrodynamics_convergence
       temperatureIteration()                     #   one radiative-equilibrium step
       if converged: break
       (optionally) applyMovableGrid()            #   adaptive regrid

The **inner loop** holds the temperature fixed and iterates chemistry, dust,
radiative transfer and the wind until the radiative acceleration :math:`\alpha`
is self-consistent.  The **outer loop** then performs a single
radiative-equilibrium temperature correction and tests global convergence.  This
two-level Picard structure follows the coupled-model scheme of Winters et al.

.. note::

   A subtle but important design choice: the carbon consumed by dust formation
   is **not** removed in a separate outer chemistry↔dust fixed point.  Instead
   the Gail & Sedlmayr moment sweep depletes the growth species *differentially*
   along the outward integration (throttling growth and nucleation by
   :math:`1-f_c`, with :math:`f_c` the running degree of condensation).  This
   single-pass, causal treatment is robust where an outer Picard on the degree
   of condensation would oscillate.  See :doc:`dust`.

Convergence criteria
====================

The outer iteration is declared converged when **all** of the following are
simultaneously below ``temperature_convergence``:

* the flux deviation from the target :math:`L_\star/(16\pi^2)`,
* the gas energy-balance residual (ratio form of :eq:`re-balance`),
* the dust energy-balance residual,
* the maximum relative temperature change actually applied this step.

Requiring the last term (a *settling* indicator) prevents declaring convergence
mid-oscillation.  The inner loop additionally requires a genuinely *accepted*
wind solve — a rejected wind solve freezes the structure, so an unchanged
:math:`\alpha` must **not** be mistaken for convergence.
