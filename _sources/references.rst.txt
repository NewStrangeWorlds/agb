==========
References
==========

The methods in this code draw on the following works.

Radiative transfer
==================

* **Kitzmann, D.**, *Diploma thesis* — the primary reference for the spherical
  variable-Eddington-factor moment scheme, the Taylor/spline discretisations, the
  flux forms (eq. 2.58/2.59/2.61), and the full-linearisation temperature
  correction (§3.2.3/3.2.4, App. B). Equation numbers in this documentation refer
  to it.
* **Rybicki, G. B. & Hummer, D. G.** — the variable-Eddington-factor /
  impact-parameter (Feautrier) formulation for spherical atmospheres.

Dust formation
==============

* **Gail, H.-P. & Sedlmayr, E. (1984)** and the textbook *Physics and Chemistry
  of Circumstellar Dust Shells* (Gail & Sedlmayr) — the moment method, classical
  nucleation theory, and the differential (in-sweep) carbon depletion (§14.3).
* **Winters, J. M.**, *PhD thesis* — the coupled stationary-wind model scheme
  (§4.3, Fig. 4.1), the differential carbon depletion (Eq. 5.4), and the
  prescribed-:math:`\dot{M}` setup (App. A).

Wind hydrodynamics
==================

* **Melia, F. (1988)** — the :math:`\Phi=\tfrac12(v+c_T^2/v)` transform that
  removes the sonic-point singularity.
* **Dominik, C. (1990)** and PhD thesis — the Henyey-type global relaxation of
  the wind + dust-moment boundary-value problem, and the observation that
  shooting is "unstable and slow."

Chemistry and opacities
=======================

* **FastChem** (Stock et al.) — equilibrium gas-phase chemistry.
  https://github.com/exoclime/fastchem
* **HELIOS-K** — the molecular/atomic opacity cross-section database used for the
  gas line opacities.
* **LX-MIE** (Kitzmann & Heng) — the Mie-theory code for the dust optical
  properties. https://github.com/daniel-kitzmann/LX-MIE

Reference object
================

* **IRC+10216 (CW Leonis)** — the prototypical carbon star and the reference
  application; observed :math:`\dot{M}\approx8\times10^{-5}\,M_\odot\,
  \mathrm{yr}^{-1}`.

Software libraries
==================

* **Eigen** (linear algebra), **CppAD** (automatic differentiation), and
  **toml++** (configuration parsing). See :doc:`installation`.

.. note::

   The diploma thesis PDF, and additional literature on the coupled model, the
   Melia transform and the Henyey relaxation, are collected in the project's
   ``Literature/`` directory alongside the source.
