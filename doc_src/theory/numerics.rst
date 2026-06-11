================================
Numerical methods and robustness
================================

This page collects the numerical machinery that wraps the physics blocks and
makes the stiff coupled iteration converge: the convergence tests, adaptive
under-relaxation, Anderson acceleration, and the optional movable grid.

Fixed-point structure and convergence tests
===========================================

The model is a two-level Picard iteration (:doc:`governing_equations`).  Three
residuals govern convergence, computed each outer step in
:cpp:func:`agb::AGBStarModel::temperatureIteration`:

* **Flux convergence** (:cpp:func:`checkFluxConvergence`) — the deviation of the
  conservative integrated flux :math:`r^2 H_\mathrm{int}` from the target
  :math:`L_\star/(16\pi^2)`, measured **against the target itself** (referencing
  the inner-boundary value would let that node's own wobble pollute the
  yardstick).
* **Energy balance** (:cpp:func:`checkEnergyBalance`) — the **ratio-form** local
  radiative-equilibrium residual :eq:`re-ratio`, separately for gas and dust.
* **Settling** — the maximum relative temperature change actually applied this
  step.

Convergence requires **all** of them below ``temperature_convergence``
*simultaneously*.  The settling term prevents declaring convergence mid-wobble;
the inner loop additionally requires a genuinely accepted wind solve.

Adaptive per-layer under-relaxation
===================================

Both temperature correctors apply a per-layer under-relaxation factor
:math:`\omega_i` that adapts to the local behaviour of the correction:

* **halve** :math:`\omega_i` (×``relaxation_down``) when the correction flips sign
  at layer :math:`i` (an oscillation), and
* **grow** it (×``relaxation_up``) when the correction keeps its sign,

clamped to ``[relaxation_min, relaxation_max]``.  Crucially, :math:`\omega` is
applied **after** the per-step cap (see below), so the damping remains effective
even while the cap is binding — a step pinned at ±cap with a flipping sign is
exactly the oscillation that must be damped.

The per-step cap
================

The actually applied temperature change is capped **last**, after any Anderson
acceleration or monotonicity clip:

.. math::

   |\Delta T_i| \le \texttt{max\_relative\_change}\times T_i ,

and :math:`T` is floored to stay positive.  This guarantees the change never
exceeds the bound the interleaved hydro–dust cycle can tolerate.  Early in the
run the cap binds; near convergence the change is below the cap and the cap no
longer interferes.

Anderson acceleration
=====================

*Method:* :cpp:func:`agb::AGBStarModel::andersonStep` (``anderson_step.cpp``).

Anderson (Pulay) mixing accelerates the Unsöld–Lucy fixed-point iteration by
extrapolating from a short history (depth ``anderson_window``) of profiles
:math:`x_k` and relaxed corrections :math:`f_k` (so that
:math:`G(x_k)=x_k+f_k`).  It is applied **only in the settling regime**
(``prev_max_rel_change < anderson_activation_fraction × max_relative_change``),
after ``anderson_start_iter`` plain damped steps, and the history is **restarted
on activation** so it is not seeded by the early transient.  Anderson is bypassed
once the Newton linearisation is active (Newton already provides the mixing).

Movable (adaptive) radial grid
==============================

*Method:* :cpp:func:`agb::AGBStarModel::applyMovableGrid`
(``movable_grid.cpp``).  **Default off** (``[movable_grid] enabled = false``).

To resolve the steep dust front, the radial nodes can be redistributed by
**equidistribution of a monitor function**.  Every ``regrid_frequency`` outer
iterations the monitor

.. math::
   :label: monitor

   w(r) = 1 + \sum_k a_k\,
   \left|\frac{\mathrm{d}\ln q_k}{\mathrm{d}\ln r}\right|

is formed from the flux-mean opacity, the gas temperature, the dust nucleation
rate, and the wind velocity (weights :math:`a_k = ` ``monitor_weight_*``; the
velocity term clusters nodes at the sonic point).  The monitor is smoothed
(``monitor_smoothing_passes``), floored (``monitor_rel_floor``) and capped
(``monitor_max``, the maximum cell-size ratio); the new nodes are placed so that
:math:`w` has equal integral per cell (``aux::equidistributedGrid``).

The node motion is then **under-relaxed** (``grid_relaxation``) with the
boundaries pinned and monotonicity preserved, the structure
(:math:`T_\mathrm{gas}, T_\mathrm{dust}, \rho, p, v`) is remapped onto the new
grid (:cpp:func:`agb::Atmosphere::remapToGrid`), the RT ray geometry is rebuilt
(:cpp:func:`agb::RadiativeTransfer::rebuildGeometry`), the hydrodynamics warm
start is discarded (:cpp:func:`agb::Hydrodynamics::resetWarmStart`), and the stale
per-node iteration state (relaxation factors, Anderson history) is reset.  The
chemistry, dust, opacities and radiation field are recomputed downstream and need
no remapping.

The grid is a step **decoupled** from the solvers — it moves between outer
iterations rather than inside any solve, which keeps each solver on a fixed grid.

Reproducibility
===============

For a fixed thread count the model is **bit-reproducible run-to-run**.  This
required care in the parallel reductions:

* The dense Newton matrix reduction over frequencies uses a static schedule with
  per-thread buffers summed in ascending thread order (an earlier arrival-order
  ``critical`` summed in a non-deterministic order, perturbing the matrix at
  :math:`10^{-15}` and changing the iteration count run-to-run).
* The Planck-power optimisation (``pow(x,5)`` → repeated multiply) was
  deliberately **not** applied because it can differ by 1 ULP and break
  bit-safety.

Results may still differ at the last digits when the thread count changes (the
floating-point summation order changes), but the converged physical solution is
the same.

Performance notes
=================

The radiative transfer is ~90 % of the runtime.  It has been optimised ~1.7×
(byte-identical output) by precomputing the iteration-invariant thermal emission
once per solve, reusing per-thread scratch buffers (a non-allocating tridiagonal
solve ``solveInto``), and parallelising the previously serial per-iteration
sections (angular integration, Eddington and sphericality factors).  The only
added build flag is ``-funroll-loops`` (see :doc:`../installation` for why
``-march=native`` must **not** be used).
