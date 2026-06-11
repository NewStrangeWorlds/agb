=================
Source tree guide
=================

A directory-by-directory guide to ``src/``, for orientation when reading or
extending the code.  The physics behind each module is in :doc:`../theory/index`.

``model_main/``
===============

* ``model_main.cpp`` — ``main``; argument handling, the ``--selftest`` path, and
  construction of :cpp:class:`agb::AGBStarModel`.

``agb_model/``
==============

The orchestrator (:cpp:class:`agb::AGBStarModel`).

* ``agb_model.cpp`` — construction, ``calcModel`` (outer loop),
  ``chemistryHydroIteration`` (inner loop), ``chemistryDustIteration``,
  ``radiativeTransfer``, ``temperatureIteration``, and the convergence tests
  (``checkFluxConvergence``, ``checkEnergyBalance``, ``checkConvergence``).
* ``anderson_step.cpp`` — Anderson mixing of the temperature profile.
* ``movable_grid.cpp`` — the equidistribution regrid step.

``config/``
===========

* ``config.{h,cpp}`` — :cpp:class:`agb::ModelConfig`: every parameter, its
  default, and the toml++ parsing.  This is the authoritative list of
  configuration keys (:doc:`../configuration`).

``spectral_grid/``
==================

:cpp:class:`agb::SpectralGrid`.

* ``spectral_grid.cpp`` — load the global wavenumber list, build the
  constant-resolution working grid.
* ``convert.cpp`` — wavelength ↔ wavenumber conversions.
* ``spectral_grid_interpolate.cpp`` — interpolation onto the grids.

``atmosphere/``
===============

:cpp:class:`agb::Atmosphere` — the shared structure and opacity store.

* ``atmosphere.cpp`` — construction, equation of state, the grey (Lucy) starting
  model (``buildGreyStart``), and the movable-grid remap (``remapToGrid``).
* ``atmosphere_read_write.cpp`` — read/write the structure file
  (:ref:`format <running-structure-file>`).

``chemistry/``
==============

:cpp:class:`agb::FastChemChemistry` — the FastChem wrapper
(``fastchem_chemistry.cpp``); ``chem_species.h`` holds the species identifiers
and the physical-constants species table.

``transport_coefficients/``
===========================

The gas opacities (:cpp:class:`agb::TransportCoefficients`,
:cpp:class:`agb::OpacitySpecies`).

* ``transport_coeff.cpp`` — assemble all gas species' contributions.
* ``opacity_species.cpp`` — sampled line cross sections, interpolation,
  Rayleigh/continuum hooks.
* ``cross_section_file.cpp``, ``sampled_data.cpp`` — read and store the tabulated
  cross sections (:cpp:class:`agb::SampledData`).
* ``h_m.cpp``, ``h2_m.cpp`` — H⁻ and H₂⁻ continuum absorption.
* ``species_rayleigh_cross_sections.cpp`` — Rayleigh scattering.
* ``species_definition.h`` — per-species data.

``radiative_transfer/``
=======================

The VEF moment radiative transfer (:cpp:class:`agb::RadiativeTransfer`,
:cpp:class:`agb::RadiationField`, :cpp:class:`agb::ImpactParam`).

* ``radiative_transfer.cpp`` — driver ``solveRadiativeTransfer``, geometry setup,
  Eddington/sphericality factors, emission precompute.
* ``impact_parameter.cpp`` — the tangent-ray (Feautrier) formal solution.
* ``moment_system.cpp`` — the :math:`J_\nu` moment system (Taylor/spline
  assembly), ``calcFlux``, the conservative flux integral, and
  ``buildLinearisedMomentSystem`` for the corrector.
* ``radiation_field.cpp`` — moment storage, angular/wavelength integration, the
  flux-mean extinction.
* ``output.cpp`` — ``saveSpectrum``.

``temperature_correction/``
===========================

:cpp:class:`agb::TemperatureCorrection`.

* ``temperature_correction.cpp`` — Unsöld–Lucy (``calculate``), the integrated
  quantities, correction smoothing, Λ-iteration helper.
* ``temperature_correction_linearisation.cpp`` — the full-linearisation Newton
  corrector (RE constraints, Rybicki elimination, dense :math:`2D\times2D`
  solve).

``dust/``
=========

* ``dust_species.{h,cpp}`` — abstract :cpp:class:`agb::DustSpecies` base.
* ``gail_sedlmayr_dust.cpp`` — the moment method with carbon depletion
  (``calcDistribution``), nucleation/growth rates, the Mie opacity.
* ``gail_sedlmayr_thermo_functions.cpp`` — the nucleation thermochemistry
  (saturation, critical cluster, free energy, Zeldovich factor, …).
* ``analytic_dust.{h,cpp}`` — a simple analytic dust prescription for testing.

``hydrodynamics/``
==================

:cpp:class:`agb::Hydrodynamics` and the optional
:cpp:class:`agb::StructureSolver`.

* ``hydrodynamics.cpp`` — the Melia Φ shooting solver, critical-point search,
  mass-loss eigenvalue, alpha, density, output.
* ``structure_solver.cpp`` — the Henyey global-Newton solver (residuals,
  scaling, grey bootstrap, eigenvalue solve).
* ``henyey_solver.cpp`` — the CppAD Newton/Jacobian driver and self-tests.

``additional/``
===============

Shared utilities: ``aux_functions`` (Planck function, integration,
``equidistributedGrid``), ``physical_const.h`` (cgs constants and the species
table), ``tri_diagonal_matrix.h`` (Thomas solver with non-allocating
``solveInto``), ``interpolation.h``, ``quadrature.h``, ``movable_grid.h``,
``solve_linear_system.h``, ``exceptions.h``.
