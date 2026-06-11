===================
Radiative transfer
===================

*Module:* ``src/radiative_transfer/`` —
:cpp:class:`agb::RadiativeTransfer`, :cpp:class:`agb::RadiationField`,
:cpp:class:`agb::ImpactParam`.

The radiative transfer is solved with the **Rybicki–Hummer variable-Eddington-
factor (VEF) moment method** for a spherically symmetric, static medium, as
described in the diploma thesis of D. Kitzmann.  The implementation matches that
scheme term by term.

Method overview
===============

The VEF method splits the angle-dependent transfer problem into two coupled
pieces that are iterated to consistency:

#. A **formal solution along tangent rays** (impact parameters) gives the
   angular distribution of the specific intensity at each radius, and from it the
   closure quantities — the **Eddington factor** :math:`f_\nu=K_\nu/J_\nu` and
   the **sphericality factor** :math:`q_\nu`.
#. A **moment system** for :math:`J_\nu(r)`, closed with the just-computed
   :math:`f_\nu` and :math:`q_\nu`, is a cheap tridiagonal problem that yields a
   new mean intensity and hence a new source function.

Iterating (1) ↔ (2) converges because the Eddington factors are only weakly
dependent on the source function, while the expensive angular detail is needed
only to update them.

Geometry: impact parameters and tangent rays
============================================

For a spherical atmosphere the natural coordinates are **impact parameters**
:math:`p` (tangent rays).  The set of rays comprises:

* ``nb_core_impact_param`` rays that strike the stellar core (:math:`p<R_\star`),
  carrying the inner boundary intensity, and
* one tangent ray per radial shell (:math:`p=r_i`), so that each shell is sampled.

Along each ray the geometric coordinate is :math:`z=\sqrt{r^2-p^2}`; the set of
``zPoint`` crossings of a ray with the radial shells defines the spatial grid for
the **Feautrier** variables

.. math::

   u_\nu = \tfrac12 (I_\nu^+ + I_\nu^-), \qquad
   v_\nu = \tfrac12 (I_\nu^+ - I_\nu^-),

the mean and difference of the outgoing/incoming intensities.  The ray geometry
is built once and rebuilt (:cpp:func:`agb::RadiativeTransfer::rebuildGeometry`)
only when the radial grid moves (movable grid).

Formal solution on a ray (Feautrier)
====================================

Along each impact parameter the second-order Feautrier equation for
:math:`u_\nu` is discretised on the ray's optical-depth grid and solved as a
tridiagonal system (:cpp:class:`agb::ImpactParam`).  Two discretisations are
available:

* **Taylor** finite differences (``assembleSystemTaylor``) — the default.
* **Cubic splines** (``assembleSystemSpline``).

The difference variable :math:`v_\nu` is then recovered (``calcV``), and angular
integration over all rays passing through a shell yields the moments
:math:`J_\nu, H_\nu, K_\nu` and thus the closure factors at that shell.

The moment system for :math:`J_\nu`
===================================

With :math:`f_\nu` and :math:`q_\nu` frozen from the formal solution, the
zeroth/first moment equations combine into a single second-order ODE for
:math:`J_\nu`.  In the optical-depth-like coordinate :math:`X` (built by
``generateXGrid`` from the extinction and the sphericality factor) this is

.. math::
   :label: moment-eq

   \frac{\partial}{\partial X}\!\left(\frac{\partial (f_\nu q_\nu r^2 J_\nu)}
   {\partial X}\right) = q_\nu r^2 (J_\nu - S_\nu) ,

with the source function

.. math::

   S_\nu = \frac{\kappa_\nu B_\nu(T) + \sigma_\nu J_\nu}{\chi_\nu}
   \quad\Longrightarrow\quad
   \chi_\nu(J_\nu - S_\nu) = \kappa_\nu (J_\nu - B_\nu) .

Here :math:`\kappa_\nu B_\nu` is the **thermal emission**; for a two-temperature
medium it is the sum over gas and dust,
:math:`\kappa_{\nu,\mathrm{gas}}B_\nu(T_\mathrm{gas})+
\kappa_{\nu,\mathrm{dust}}B_\nu(T_\mathrm{dust})`.  This emission is
iteration-invariant within one RT solve and is therefore precomputed once
(``precomputeEmission`` → ``planck_emission[nu][i]``).

:eq:`moment-eq` is assembled as a tridiagonal system
(``assembleMomentSystemTaylor`` / ``…Spline``) and solved with the Thomas
algorithm (:cpp:class:`agb::aux::TriDiagonalMatrix`).  Iterating with the
scattering term :math:`\sigma_\nu J_\nu` updated each pass converges the field;
``solveMomentSystem`` drives this for every frequency.

.. _rt-taylor-vs-spline:

Taylor vs. spline discretisation
================================

The two discretisations differ in their stability:

* **Taylor** finite differences keep the system matrix an **M-matrix**
  (positivity-preserving) unconditionally, so :math:`J_\nu\ge0` always.
* **Cubic splines** lose the M-matrix property once the optical depth per cell
  exceeds :math:`\sqrt 6\approx2.45`, producing *negative* mean intensities —
  exactly as the thesis (§A.1.3) warns.

.. warning::

   Use the **Taylor** discretisation (``use_spline_discretisation = false``) for
   production.  The spline option is retained for cross-checks only.  (Historic
   negative-intensity crashes were ultimately traced to a *density-collapse* bug
   upstream in the hydrodynamics, not to the RT scheme: a transparent atmosphere
   makes both the moment system and the Feautrier solve degenerate.)

Two flux definitions and why they coexist
=========================================

The frequency-*integrated* flux is needed in two roles that pull toward
different discretisations, so the code keeps **two** integrated-flux fields:

* **Transport form (eq. 2.58)**, ``eddington_flux`` / ``eddington_flux_int`` —
  the per-frequency flux from the gradient of :math:`f_\nu q_\nu r^2 J_\nu`
  (``calcFlux``).  This is physical (non-negative) per frequency, so it is used
  for the emergent **spectrum** and as the denominator of the **flux-mean
  extinction** and other ratios.
* **Divergence form (eq. 2.59)**, ``eddington_flux_int_conservative`` —
  obtained by integrating

  .. math::
     :label: flux-div

     \frac{\partial (r^2 H_\mathrm{int})}{\partial r}
     = r^2 \!\int \kappa_{\nu,\mathrm{abs}}\,(B_\nu - J_\nu)\,\mathrm{d}\nu ,

  outward from the diffusion inner boundary condition
  :math:`r^2 H_\mathrm{int}(r_1)=r_1^2\,\mathrm{bfc}\cdot
  (\mathrm{d}B/\mathrm{d}T)/\chi`.  Its divergence equals the local balance that
  the moment-equation :math:`J`-solve enforces, so :math:`r^2 H_\mathrm{int}` is
  **conserved at radiative equilibrium**.  This is the quantity compared against
  the target luminosity in the flux-convergence test and used by the
  temperature corrector's flux term.

.. note::

   Keeping these separate is essential.  ``eddington_flux_int`` is simultaneously
   a *ratio denominator* (which must stay consistent with its eq. 2.58
   numerator) and a *flux value* (which wants the conservative eq. 2.59 form).
   Conflating them corrupts the flux-mean extinction, which feeds the wind, and
   makes the coupled loop diverge.  The behaviour is controlled by
   ``flux_from_divergence`` (:doc:`../configuration`).

Outputs of the RT module
========================

For each radial shell, :cpp:class:`agb::RadiationField` stores the spectral
moments (``mean_intensity``, ``eddington_flux``, ``eddington_k``), their
wavelength integrals, the Eddington and sphericality factors, and the
angle-grid intensities.  It also provides the **flux-weighted (flux-mean)
extinction** :math:`\chi_F` (``fluxWeightedExtinction``), which is the single
most important quantity handed to the hydrodynamics — it sets the radiative
acceleration :eq:`alpha-def`.

The linearised moment operator
==============================

For the full-linearisation temperature corrector the RT module can also expose
its operator in linearised form (``buildLinearisedMomentSystem``): per frequency
it returns the frozen Taylor moment matrix :math:`V_\nu`, the diagonal source
derivative :math:`\mathrm{d}S/\mathrm{d}T`, and the moment-equation residual
:math:`K_\nu=\mathrm{rhs}_\nu - V_\nu J_\nu`.  The temperature corrector consumes
these; see :doc:`temperature_correction`.
