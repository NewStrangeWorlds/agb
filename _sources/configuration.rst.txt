===========================
Configuration (config.toml)
===========================

The model is configured by a single `TOML <https://toml.io>`__ file named
:file:`config.toml` in the model folder.  It is parsed by
:cpp:class:`agb::ModelConfig` (``src/config/config.cpp``) using toml++.

Conventions
===========

* **Stellar parameters are given in solar units** and converted to cgs
  internally.
* The parser is *lenient about numeric types*: an integer literal also satisfies
  a floating-point key and vice versa (e.g. ``resolution = 1000`` is accepted for
  a ``double`` field).
* **Any missing key falls back to a built-in default** defined in
  ``src/config/config.h``.  Only the ``[star]`` block, the starting model, the
  spectral grid, and the opacity list are practically required.
* Output and data paths are interpreted **relative to the model folder** (except
  the opacity ``path`` and dust ``refractive_index_file``, which are full paths).

A complete annotated example is shown at the :ref:`bottom of this page
<configuration-example>`.

``[star]`` — stellar parameters
===============================

.. list-table::
   :header-rows: 1
   :widths: 26 14 60

   * - Key
     - Units
     - Description
   * - ``radius``
     - :math:`R_\odot`
     - Stellar (reference) radius :math:`R_\star`.
   * - ``mass``
     - :math:`M_\odot`
     - Stellar mass :math:`M_\star` (sets gravity).
   * - ``luminosity``
     - :math:`L_\odot`
     - Stellar luminosity :math:`L_\star` (the flux target the temperature
       correction drives toward and the driving luminosity of the wind).
   * - ``mass_loss_rate``
     - :math:`M_\odot\,\mathrm{yr}^{-1}`
     - Mass-loss rate.  Used as a *prescribed* value or as the seed for the
       eigenvalue solve, depending on the wind setup.
   * - ``c_o_ratio``
     - —
     - Carbon-to-oxygen ratio.  Must be > 1 for a carbon star; the condensable
       carbon is :math:`\varepsilon_C-\varepsilon_O` (all oxygen is assumed
       locked in CO).

``[model]``
===========

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Key
     - Description
   * - ``starting_model``
     - A path to a structure file (:ref:`format <running-structure-file>`), or
       the literal ``"grey"`` to build the Lucy grey bootstrap (see ``[grid]``).
   * - ``fastchem_parameter_file``
     - FastChem parameter / element-abundance file (relative to the model
       folder).

``[spectral_grid]``
===================

.. list-table::
   :header-rows: 1
   :widths: 30 14 56

   * - Key
     - Units
     - Description
   * - ``min_wavelength``
     - μm
     - Lower wavelength bound of the spectral grid.
   * - ``max_wavelength``
     - μm
     - Upper wavelength bound.
   * - ``resolution``
     - —
     - Constant spectral resolving power :math:`R=\lambda/\Delta\lambda` of the
       internal grid.
   * - ``nb_core_impact_param``
     - —
     - Number of impact parameters that intersect the stellar core, used by the
       ray (impact-parameter) radiative-transfer geometry. Default 20.

``[radiative_transfer]``
========================

.. list-table::
   :header-rows: 1
   :widths: 34 12 54

   * - Key
     - Default
     - Description
   * - ``nb_iterations``
     - 30
     - Maximum number of variable-Eddington-factor iterations per RT solve.
   * - ``convergence``
     - 1e-4
     - Relative convergence criterion of the RT iteration.
   * - ``use_spline_discretisation``
     - false
     - Use cubic-spline (``true``) instead of Taylor finite-difference
       (``false``) moment discretisation.  **Keep false**: splines lose
       positivity for optically thick cells (see :doc:`theory/radiative_transfer`).
   * - ``flux_from_divergence``
     - true
     - Compute the frequency-*integrated* flux from the flux-divergence form
       (thesis eq. 2.59), consistent with the moment-equation solve, so
       :math:`r^2 H` is conserved at radiative equilibrium.  The per-frequency
       flux stays on the transport form (eq. 2.58) so the monochromatic spectrum
       is non-negative.

``[temperature]``
=================

This block controls the radiative-equilibrium temperature determination
(:doc:`theory/temperature_correction`) and its robustness machinery
(:doc:`theory/numerics`).

Core
----

.. list-table::
   :header-rows: 1
   :widths: 34 12 54

   * - Key
     - Default
     - Description
   * - ``nb_iterations``
     - 200
     - Maximum number of global (outer) temperature iterations.
   * - ``convergence``
     - 1e-2
     - Convergence threshold applied jointly to the flux deviation, the gas and
       dust energy-balance residuals, and the maximum relative temperature
       change.
   * - ``max_relative_change``
     - 0.005
     - Per-step cap on :math:`|\Delta T|/T`.  Keeps each step small enough for
       the interleaved hydro–dust cycle to re-converge.
   * - ``smooth_profile``
     - true
     - Apply a 1–2–1 smoothing pass to the temperature profile.
   * - ``force_monotonic``
     - false
     - Hard-clip the profile to be strictly outward-decreasing
       (:math:`T_i\to0.99\,T_{i-1}`).  **Default off** — it fights the
       correction wherever radiative equilibrium wants a flatter profile or a
       dust-front inversion and drives a ±cap limit cycle.

Full-linearisation corrector
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 40 10 50

   * - Key
     - Default
     - Description
   * - ``use_linearisation``
     - true
     - Enable the full-linearisation (Newton) corrector in addition to / instead
       of Unsöld–Lucy.
   * - ``linearisation_relaxation``
     - 0.5
     - Global under-relaxation of the Newton step (dust-front stability).
   * - ``linearisation_start_unsoeld_lucy``
     - true
     - Two-phase corrector: start with the robust Unsöld–Lucy iteration and latch
       onto the linearisation only once it is safe.  Set ``false`` to linearise
       immediately (e.g. a warm restart near convergence).
   * - ``linearisation_switch_dt_fraction``
     - 0.5
     - Switch when ``max|dT|/T`` < this × ``max_relative_change`` …
   * - ``linearisation_switch_re_residual``
     - 3e-2
     - … **and** the max radiative-equilibrium residual is below this …
   * - ``linearisation_switch_count``
     - 3
     - … for this many consecutive steps (then latch one-way).
   * - ``linearisation_flux_constraint``
     - true
     - Add the composite flux-constancy term (thesis eq. 3.64/3.65) to the
       local-RE constraint.  See :doc:`theory/temperature_correction`.
   * - ``linearisation_xi``
     - 1.0
     - Weight :math:`\xi` of the local radiative-equilibrium term.
   * - ``linearisation_zeta_tau_scale``
     - 1.0
     - Optical-depth scale in :math:`\zeta(r)=\tau/(\tau+\tau_\mathrm{scale})`,
       the depth-dependent weight of the flux-constancy term.

Robustness (Unsöld–Lucy under-relaxation and Anderson)
------------------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 36 12 52

   * - Key
     - Default
     - Description
   * - ``relaxation_init``
     - 1.0
     - Initial per-layer under-relaxation factor :math:`\omega`.
   * - ``relaxation_min`` / ``relaxation_max``
     - 0.05 / 1.0
     - Floor / ceiling for :math:`\omega`.
   * - ``relaxation_down``
     - 0.5
     - Factor applied to :math:`\omega` when the correction flips sign
       (oscillation).
   * - ``relaxation_up``
     - 1.3
     - Factor applied to :math:`\omega` when the correction keeps its sign.
   * - ``use_anderson``
     - true
     - Enable Anderson mixing on the temperature profile.
   * - ``anderson_window``
     - 4
     - History depth :math:`m` of the Anderson mixing.
   * - ``anderson_start_iter``
     - 2
     - Number of plain damped iterations before acceleration starts.
   * - ``anderson_activation_fraction``
     - 0.5
     - Accelerate only once ``max|dT|/T`` < this × ``max_relative_change`` (the
       settling regime).

``[hydrodynamics]``
===================

.. list-table::
   :header-rows: 1
   :widths: 30 12 58

   * - Key
     - Default
     - Description
   * - ``nb_iterations``
     - 200
     - Maximum number of chemistry–hydrodynamics inner iterations.
   * - ``convergence``
     - 1e-2
     - Relative convergence threshold on the radiative acceleration
       :math:`\alpha`.
   * - ``use_henyey_solver``
     - false
     - Use the Henyey global-Newton structure solver instead of the Melia ``Φ``
       shooting method.  See :doc:`theory/hydrodynamics`.

``[grid]`` — only used for ``starting_model = "grey"``
======================================================

.. list-table::
   :header-rows: 1
   :widths: 26 12 62

   * - Key
     - Default
     - Description
   * - ``n_points``
     - 260
     - Number of logarithmically spaced radial cells.
   * - ``inner_radius``
     - 1.0
     - Inner boundary in units of :math:`R_\star`.
   * - ``outer_radius``
     - 25.0
     - Outer boundary in units of :math:`R_\star`.
   * - ``gas_opacity``
     - 1.0e-4
     - Grey gas mass-opacity (:math:`\mathrm{cm^2\,g^{-1}}`) that sets the
       optical-depth scale and hence the photospheric density of the seed.

``[dust]``
==========

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Key
     - Description
   * - ``refractive_index_file``
     - Path to the optical constants (complex refractive index vs. wavelength)
       of the dust material, e.g. amorphous carbon (Rouleau, Maron).  Fed to Mie
       theory to produce the grain absorption/scattering efficiencies.

``[movable_grid]`` — adaptive radial grid (default off)
=======================================================

See :doc:`theory/numerics` for the equidistribution scheme.

.. list-table::
   :header-rows: 1
   :widths: 36 12 52

   * - Key
     - Default
     - Description
   * - ``enabled``
     - false
     - Master switch for the movable grid.
   * - ``regrid_frequency``
     - 5
     - Redistribute nodes every N outer iterations.
   * - ``grid_relaxation``
     - 0.1
     - Under-relaxation :math:`\omega_\mathrm{grid}` of the node motion.
   * - ``monitor_weight_opacity``
     - 1.0
     - Weight of the flux-mean extinction in the monitor function.
   * - ``monitor_weight_temperature``
     - 0.0
     - Weight of the gas temperature.
   * - ``monitor_weight_nucleation``
     - 1.0
     - Weight of the dust nucleation rate.
   * - ``monitor_weight_velocity``
     - 1.0
     - Weight of the wind velocity (clusters nodes at the sonic point).
   * - ``monitor_smoothing_passes``
     - 4
     - Number of smoothing passes applied to the monitor.
   * - ``monitor_max``
     - 20.0
     - Cap on the monitor (maximum cell-size ratio).
   * - ``monitor_rel_floor``
     - 1e-8
     - Per-quantity floor relative to its peak.

``[opacity]``
=============

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - Key
     - Description
   * - ``path``
     - Absolute path to the opacity database (HELIOS-K cross-section folders).
       The wavenumber grid is read from ``<path>/wavenumber_full.dat``.
   * - ``species``
     - An array of ``{ symbol = "...", folder = "..." }`` tables.  ``symbol`` is
       the species identifier (e.g. ``"CO"``, ``"H2O"``, ``"CIA-H2-H2"``,
       ``"H-"``); ``folder`` is the sub-folder of the opacity database holding
       its cross sections, or ``"none"`` for species with built-in opacity
       (continuum / collision-induced, e.g. ``H-``, ``H2-``, Rayleigh).

``[output]``
============

Each value is a path relative to the model folder; an empty value or ``"none"``
disables that output.

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Key
     - Output
   * - ``atmosphere``
     - Converged structure (:ref:`format <running-structure-file>`).
   * - ``spectrum``
     - Emergent spectrum.
   * - ``dust``
     - Dust distribution and moments.
   * - ``hydro``
     - Wind structure.

.. _configuration-example:

Complete example (IRC+10216)
============================

.. code-block:: toml

   # AGB dust-driven wind model configuration (IRC+10216)

   [star]
   radius         = 1050.0     # R_sun
   mass           = 0.7        # M_sun
   luminosity     = 3.0e4      # L_sun
   mass_loss_rate = 8.0e-5     # M_sun / yr
   c_o_ratio      = 1.40

   [model]
   starting_model          = "atmosphere_converged.dat"   # a path, or "grey"
   fastchem_parameter_file = "fastchem_parameters.dat"

   [spectral_grid]
   min_wavelength = 0.3        # micron
   max_wavelength = 50.0       # micron
   resolution     = 1000

   [radiative_transfer]
   nb_iterations             = 40
   convergence               = 1.0e-4
   use_spline_discretisation = false
   flux_from_divergence      = true

   [temperature]
   nb_iterations       = 300
   convergence         = 1.0e-3
   max_relative_change = 0.01
   smooth_profile      = false
   use_linearisation                = true
   linearisation_relaxation         = 0.5
   linearisation_start_unsoeld_lucy = true
   linearisation_switch_dt_fraction = 0.5
   linearisation_switch_re_residual = 0.03
   linearisation_switch_count       = 3
   linearisation_flux_constraint    = true
   linearisation_xi                 = 1.0
   linearisation_zeta_tau_scale     = 1.0

   [hydrodynamics]
   nb_iterations     = 60
   convergence       = 1.0e-3
   use_henyey_solver = false

   [dust]
   refractive_index_file = "../opacity_data/Amorphous_C_Rouleau.dat"

   [movable_grid]
   enabled = false

   [opacity]
   path = "/media/data/opacity_data/helios-k/"
   species = [
     { symbol = "CIA-H2-H2", folder = "CIA/H2-H2" },
     { symbol = "CIA-H2-He", folder = "CIA/H2-He" },
     { symbol = "CIA-H-He",  folder = "CIA/H-He" },
     { symbol = "H-",        folder = "none" },
     { symbol = "H2",        folder = "none" },
     { symbol = "He",        folder = "none" },
     { symbol = "H",         folder = "none" },
     { symbol = "H2O",       folder = "Molecules/1H2-16O__POKAZATEL_e2b" },
     { symbol = "CO",        folder = "Molecules/12C-16O__Li2015_e2b" },
     { symbol = "CO2",       folder = "Molecules/12C-16O2__UCL-4000_e2b" },
     { symbol = "C2H2",      folder = "Molecules/12C2-1H2__aCeTY_e2b" },
   ]

   [output]
   spectrum   = "spectrum.dat"
   atmosphere = "atmosphere.dat"
   dust       = "dust.dat"
   hydro      = "hydro.dat"
