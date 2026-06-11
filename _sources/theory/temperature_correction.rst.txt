============================================
Radiative-equilibrium temperature correction
============================================

*Module:* ``src/temperature_correction/`` —
:cpp:class:`agb::TemperatureCorrection`.

The temperature structure is fixed by **radiative equilibrium**, not by an
energy equation.  Two correctors are implemented and can be combined in a
two-phase scheme:

#. the **Unsöld–Lucy** correction — a robust, approximate fixed-point update;
#. the **full-linearisation** (Newton) correction — an exact-Jacobian step that
   removes the Unsöld–Lucy accuracy floor.

Both correct the **gas** and the **dust** temperatures independently, each toward
its own local radiative-equilibrium condition :eq:`re-balance`.

Why two temperatures
====================

Gas and dust exchange energy with the radiation field through *different*
opacities (molecular/atomic lines and continua for the gas; the amorphous-carbon
Mie opacity for the dust) and are not collisionally locked in the thin outer
wind.  Each component therefore has its own radiative-equilibrium temperature:
the energy-balance condition :eq:`re-balance` is solved separately for
:math:`T_\mathrm{gas}` (with :math:`\kappa_{\nu,\mathrm{gas}}`) and
:math:`T_\mathrm{dust}` (with :math:`\kappa_{\nu,\mathrm{dust}}`), sharing the
**same** radiation field :math:`J_\nu`.

Unsöld–Lucy correction
======================

*Method:* :cpp:func:`agb::TemperatureCorrection::calculate`.

The Unsöld–Lucy scheme produces a Planck-integrated correction
:math:`\Delta B` per layer from a combination of the **local** energy-balance
residual and the **non-local** flux-constancy residual, weighted by the
absorption-, Planck- and flux-mean opacities (``calcIntegratedQuantities``
assembles :math:`f q J`, :math:`q\chi H`, :math:`\kappa_J`, :math:`\kappa_B` and
the integrated Planck function).  The correction is converted to a temperature
change through an **exact** :math:`T^4` (Planck) update rather than a linearised
one.

The scheme is only first-order accurate and leaves a residual **accuracy
floor** (a few × :math:`10^{-3}` in the energy balance, ~2 % in flux
conservation for this problem).  It is, however, very robust far from
convergence, which is exactly where it is used.

.. _temp-ratio-form:

The decisive subtlety: the **ratio** form of the constraint
===========================================================

Both correctors enforce local radiative equilibrium in **ratio** form
(thesis 3.40/3.41), per species :math:`s`:

.. math::
   :label: re-ratio

   g_s = \frac{\displaystyle\sum_\nu w_\nu\,\kappa_{\nu,s}\,J_\nu}
              {\displaystyle\sum_\nu w_\nu\,\kappa_{\nu,s}\,B_\nu(T_s)} - 1 = 0 ,

**not** the difference form
:math:`\sum_\nu w_\nu\kappa_{\nu,s}(J_\nu-B_\nu)`.

.. warning::

   This choice is essential for the linearisation.  The difference-form integrals
   span many decades with depth, so the dense Newton matrix becomes
   catastrophically ill-scaled (deep rows huge, thin rows tiny) and produces
   garbage — uniform unresolved residuals, wild dust-front swings, and a flux
   that bulges in the interior.  Dividing each species row by its denominator
   :math:`\mathrm{den}_{s,i}=\sum_\nu w_\nu\kappa_{\nu,s}B_\nu` makes **every row
   :math:`O(1)`**, and the Newton solve is well conditioned.

Full-linearisation (Newton) correction
======================================

*Method:* :cpp:func:`agb::TemperatureCorrection::linearisedCorrection`
(``temperature_correction_linearisation.cpp``), with the RT half supplied by
:cpp:func:`agb::RadiativeTransfer::buildLinearisedMomentSystem`.

One Newton step is taken on :math:`(T_\mathrm{gas}, T_\mathrm{dust})` with the
converged RT operator **frozen** (Eddington/sphericality factors, opacities,
geometry).  The key point is that the intensity response is the **exact** moment-
system response, not the :math:`\tau`-weighted Unsöld–Lucy heuristic.

Per frequency, the moment system :math:`V_\nu J_\nu = \mathrm{rhs}_\nu(T)` is
re-assembled with the same Taylor stencil used in the RT solve, so the intensity
response to a temperature change is

.. math::

   \delta J_\nu = V_\nu^{-1}\bigl(K_\nu + S_g\,\delta T_g + S_d\,\delta T_d\bigr),
   \qquad
   S_{s,i} = \frac{\partial\,\mathrm{rhs}_i}{\partial T_{s,i}}
           = -c_i\,\kappa_{\mathrm{abs},s}\,\frac{\mathrm{d}B}{\mathrm{d}T},

with :math:`c_i=r_i^2/(q_i\chi_i)` and the residual
:math:`K_\nu=\mathrm{rhs}_\nu - V_\nu J_\nu`.  The diagonal (in the Taylor
discretisation) source derivatives make :math:`S_{s}` cheap.

The mean intensities are eliminated **Rybicki-style** — solving :math:`V_\nu`
against each unit node vector gives the inverse columns — collapsing the problem
to a dense :math:`2D\times2D` temperature system (:math:`D=` number of grid
points, factor 2 for gas+dust), solved by Gaussian elimination.

.. note::

   The frequency reduction into the dense Newton matrix uses a **deterministic**
   (static-scheduled, ascending-thread-order) reduction.  An earlier
   arrival-order reduction summed the ~5000 frequencies in a different order each
   run, perturbing the matrix at the :math:`10^{-15}` level; amplified by the
   stiff dust-nucleation coupling this changed the iteration count run-to-run
   (same solution).  The reduction is now bit-reproducible at a fixed thread
   count.

Local RE alone is not enough: the flux-constancy term
=====================================================

A frozen-structure test revealed that driving the **local** RE residual
:eq:`re-ratio` to zero does **not** conserve the flux: in the optically thick
interior local RE is degenerate (:math:`J\approx B` regardless of :math:`T`), so
the Newton converges without pinning the flux (:math:`r^2 H` bulged by ~13 %).

The fix is the **composite** constraint (thesis §3.2.4, eq. 3.64/3.65): per
species,

.. math::
   :label: composite-constraint

   \xi\,\underbrace{g_{\text{local},s}}_{\text{local-RE ratio}}
   \;+\;
   \zeta(r)\,\underbrace{g_{\text{flux},s}}_{\text{flux constancy}} = 0 ,

with the flux-constancy residual written in the **integrated (eq. 2.59)** form

.. math::

   g_{\text{flux},i} = \frac{1}{r_\star^2 H_\star}
   \int_{r_1}^{r_i} r'^2\,\Bigl[\textstyle\sum_\nu w_\nu\bigl(
   \kappa_g (J_\nu-B_g) + \kappa_d (J_\nu-B_d)\bigr)\Bigr]\,\mathrm{d}r' .

Because the divergence of this integral is exactly the local-balance source that
``solveMomentSystem`` discretises (eq. 2.61 = eq. 2.59 in :math:`X`-coordinates),
the flux term is **consistent** with the :math:`J`-solve and does not fight local
RE — and the radial integral with the :math:`1/(r_\star^2 H_\star)`
normalisation keeps it :math:`O(1)`-scaled.

The depth weight is :math:`\zeta(r)=\tau/(\tau+\tau_\mathrm{scale})`, built from a
grey radial optical depth (flux-mean extinction integrated inward): :math:`\zeta
\to1` deep where flux constancy must dominate, :math:`\zeta\to0` at the thin
outer edge where the flux derivative is numerically noisy (and forced to 0 at the
two boundary nodes).  The weights are ``linearisation_xi`` (:math:`\xi`, ≈1) and
``linearisation_zeta_tau_scale``.

.. important::

   Only the **integrated (eq. 2.59)** flux operator works.  Linearising the
   node-centred eq. 2.58 (transport) flux fails — its discrete divergence is
   inconsistent with the :math:`J`-solve stencil and the composite Jacobian
   cancels where the two constraints pull oppositely (determinant → 0,
   divergence).  The flux term is applied on the **gas row only** (dust keeps
   pure local RE) to avoid gas/dust collinearity.

Two-phase corrector and coupled-loop stability
==============================================

*Driver:* :cpp:func:`agb::AGBStarModel::temperatureIteration`.

The frozen-opacity Newton step over-corrects far from convergence — the dust
nucleation rate :math:`J_\star\sim\exp(-/T)` is so temperature-sensitive that a
frozen-opacity step badly mispredicts the dust response.  Two mechanisms keep the
coupled iteration stable:

#. **Two-phase start.**  Unless told otherwise
   (``linearisation_start_unsoeld_lucy``), the iteration begins with Unsöld–Lucy
   and **latches** onto the linearisation only once, for
   ``linearisation_switch_count`` consecutive steps, both the per-step
   ``max|dT|/T`` has settled below
   ``linearisation_switch_dt_fraction × max_relative_change`` **and** the maximum
   RE residual is below ``linearisation_switch_re_residual``.  The switch is
   one-way (no revert).

#. **Per-layer adaptive under-relaxation of the Newton step.**  In the coupled
   loop the structure (dust, hydro) updates between temperature steps, so a full
   Newton step can drive a period-2 ±cap limit cycle just beyond the dust front.
   The Newton :math:`\Delta T` is damped with the same adaptive scheme as
   Unsöld–Lucy: a per-layer factor :math:`\omega` halved on a sign flip and grown
   back when the step keeps its sign, **applied after the per-step cap** so the
   damping is visible even while the cap binds.  A global
   ``linearisation_relaxation`` (default 0.5) further under-relaxes the step.

The step order in :cpp:func:`temperatureIteration` is, deliberately:

#. compute the raw correction (UL or linearisation);
#. (linearisation) per-layer adaptive damping;
#. (UL only, and only in the settling regime) Anderson acceleration;
#. optional monotonicity clip;
#. **cap** :math:`|\Delta T|\le` ``max_relative_change`` :math:`\times T`
   **last**, and floor :math:`T>0`.

Capping last guarantees that nothing (an Anderson overshoot, a monotonicity clip)
produces a step larger than the bound the hydro–dust cycle tolerates.

.. note::

   Two historical oscillation bugs are worth recording, as they recur if the
   ordering is disturbed:

   * ``force_monotonic`` fought the correction (the clip cascaded the
     just-applied heating back down every step) → ±cap limit cycle.  Default off.
   * The under-relaxation :math:`\omega` was applied **before** the cap, so while
     the cap was binding the step stayed pinned at ±cap regardless of
     :math:`\omega` and an oscillating layer never damped.  Fix: cap first, then
     multiply by :math:`\omega`.

Anderson acceleration
=====================

For the Unsöld–Lucy phase, Anderson mixing
(:cpp:func:`agb::AGBStarModel::andersonStep`) accelerates the slow creep toward
radiative equilibrium.  It is gated to the **settling regime** only
(``prev_max_rel_change < anderson_activation_fraction × max_relative_change``),
so the early large-change transient — to which the hydro cycle is sensitive —
stays on plain damped steps, and the history is restarted fresh on activation so
it is not seeded by the transient.  Anderson is bypassed once the linearisation
is active (Newton already mixes).

Accuracy achieved
=================

With the integrated flux term, the linearisation corrector reaches flux
conservation at the :math:`\sim10^{-5}` level on a frozen structure (versus the
~2.3 % Unsöld–Lucy floor) and converges the coupled IRC+10216 model in ~20 outer
iterations from a warm start.
