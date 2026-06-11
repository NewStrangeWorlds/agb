============
Architecture
============

This page maps the physics onto the C++ code: the top-level class, the module
objects it owns, and the flow of control.

Entry point
===========

``src/model_main/model_main.cpp`` contains ``main``.  It takes the model folder
as the single argument, constructs an :cpp:class:`agb::AGBStarModel`, and calls
:cpp:func:`calcModel`.  The special argument ``--selftest`` instead runs the
Henyey structure-solver validation (Jacobian AD-vs-FD, Newton solve, eigenvalue
solve, grey bootstrap) and exits.

The top-level model object
==========================

:cpp:class:`agb::AGBStarModel` (``src/agb_model/``) owns one instance of each
physics module and orchestrates the iteration.  The members, in construction
order (which is also the dependency order), are:

.. list-table::
   :header-rows: 1
   :widths: 26 30 44

   * - Member
     - Class
     - Role
   * - ``config``
     - :cpp:class:`agb::ModelConfig`
     - Parsed ``config.toml``; owns all parameters.
   * - ``spectral_grid``
     - :cpp:class:`agb::SpectralGrid`
     - Wavenumber/wavelength grids and interpolation.
   * - ``atmosphere``
     - :cpp:class:`agb::Atmosphere`
     - The persistent structure (radius, :math:`T_\mathrm{gas}`,
       :math:`T_\mathrm{dust}`, :math:`\rho`, :math:`p`, :math:`v`), number
       densities, and all opacity arrays; the equation of state; structure I/O.
   * - ``chemistry``
     - :cpp:class:`agb::FastChemChemistry`
     - Equilibrium gas-phase composition (FastChem).
   * - ``dust_species``
     - :cpp:class:`agb::DustSpecies` ⟵ :cpp:class:`agb::GailSedlmayrDust`
     - Dust nucleation/growth (moment method) and dust opacity (Mie).
   * - ``transport_coeff``
     - :cpp:class:`agb::TransportCoefficients`
     - Gas opacities (lines, CIA, H⁻/H₂⁻, Rayleigh).
   * - ``radiative_transfer``
     - :cpp:class:`agb::RadiativeTransfer`
     - Spherical VEF moment radiative transfer; owns the
       :cpp:class:`agb::RadiationField` per shell.
   * - ``temperature_correction``
     - :cpp:class:`agb::TemperatureCorrection`
     - Unsöld–Lucy and full-linearisation correctors.
   * - ``hydrodynamics``
     - :cpp:class:`agb::Hydrodynamics`
     - Wind structure, sonic point, mass-loss rate; optional
       :cpp:class:`agb::StructureSolver` (Henyey).

The modules communicate almost entirely **through the shared**
:cpp:class:`agb::Atmosphere` **object** (temperatures, density, number densities,
opacities) and through the :cpp:class:`agb::RadiationField` (mean intensity, flux,
flux-mean extinction).  This keeps the interfaces narrow.

Control flow
============

.. code-block:: text

   main
   └─ AGBStarModel::calcModel                         (outer / global loop)
      ├─ chemistryHydroIteration                      (inner loop)
      │  ├─ chemistryDustIteration
      │  │  ├─ FastChemChemistry::calcChemicalComposition
      │  │  ├─ Atmosphere::equationOfState
      │  │  └─ GailSedlmayrDust::calcDistribution      (carbon-depleting RK4 sweep)
      │  ├─ radiativeTransfer
      │  │  ├─ dust  + gas opacities  (per shell)
      │  │  └─ RadiativeTransfer::solveRadiativeTransfer  (VEF moment method)
      │  ├─ Hydrodynamics::setDustState
      │  ├─ Hydrodynamics::calcWindVelocity            (Melia Φ / shooting / Henyey)
      │  └─ Atmosphere::equationOfState
      ├─ temperatureIteration                          (one RE correction step)
      │  ├─ RadiativeTransfer::solveRadiativeTransfer
      │  ├─ TemperatureCorrection::calculate           (Unsöld–Lucy)   ── or ──
      │  ├─ TemperatureCorrection::linearisedCorrection (full linearisation)
      │  ├─ adaptive damping / Anderson / cap
      │  └─ convergence tests (flux, energy balance, settling)
      └─ applyMovableGrid                              (optional, every N iters)

After the outer loop converges (or hits ``nb_temperature_iter``), the spectrum,
structure, dust and hydro outputs are written.

Persistent iteration state
==========================

:cpp:class:`agb::AGBStarModel` carries the state that must survive across outer
iterations: the temperature-iteration counter, the two-phase corrector latch
(``linearisation_active``, ``linearisation_ready_count``), the per-layer
relaxation factors and previous corrections (``relaxation_gas/dust``,
``prev_delta_b_gas/dust``), and the Anderson history
(``anderson_x_*``, ``anderson_f_*``).  All of this is reset at the start of
:cpp:func:`calcModel` and (the per-node parts) when the movable grid fires.

Third-party libraries
=====================

See :doc:`../installation` for versions.  FastChem (chemistry), LX-MIE (dust Mie
theory), Eigen + CppAD (Henyey linear algebra and autodiff), and toml++ (config
parsing) are fetched and built by CMake.
