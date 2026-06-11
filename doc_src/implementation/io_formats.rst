============
File formats
============

All input/output tables are plain whitespace-separated ASCII with a comment
header.  Paths are relative to the model folder unless noted.

Input
=====

``config.toml``
   The model configuration; see :doc:`../configuration`.

Starting structure (``[model] starting_model``)
   A structure table (:ref:`format <running-structure-file>`), unless the value
   is ``"grey"`` (built-in Lucy bootstrap).

FastChem parameter file (``[model] fastchem_parameter_file``)
   FastChem's own element-abundance / parameter file.

Opacity database (``[opacity] path``)
   The HELIOS-K cross-section folders and ``wavenumber_full.dat``; external to
   the model folder.

Dust refractive indices (``[dust] refractive_index_file``)
   Complex refractive index (n, k) vs. wavelength of the dust material.

Output
======

Atmosphere structure (``[output] atmosphere``)
   Two header lines, then one row per grid point.  Columns:
   ``r/R*``, ``r(cm)``, ``rho(g/cm3)``, ``p(dyn/cm2)``, ``T_gas(K)``,
   ``T_dust(K)``, ``v(cm/s)``.  The first header line also records the stellar
   parameters and :math:`\dot{M}`.  This file can be fed back in as a starting
   model.

Spectrum (``[output] spectrum``)
   One header line, then columns:
   ``Wavelength (micron)``, ``Flux (erg s-1 cm-2 micron-1)``,
   ``r^2 F``.  The flux is the emergent (outermost-shell) astrophysical flux.

Dust (``[output] dust``)
   One header line, then per grid point:
   ``r/R*``, ``n<H>(cm-3)``, ``Tgas(K)``, ``f_cond`` (degree of condensation),
   the carbon-species number densities ``n_C, n_C2, n_C2H, n_C2H2``,
   ``n_d`` (grain number density), ``a_d(micron)`` (mean grain radius),
   ``J*(s-1 cm-3)`` (nucleation rate), ``tau_growth`` (growth timescale),
   the moments ``K0, K1, K2, K3``, and the mean condensed volume per grain
   ``V_d(cm3)``.

Hydrodynamics (``[output] hydro``)
   One header line, then per grid point:
   ``r/R*``, ``v(cm/s)``, ``c_T(cm/s)`` (isothermal sound speed),
   ``alpha`` (radiative/gravitational acceleration ratio),
   ``phi`` (the Melia :math:`\Phi` variable).

.. note::

   A ``hydro_debug.dat`` file may also appear in the model folder: it is an
   optional per-grid-point breakdown of the wind-equation terms (gravity,
   radiative :math:`\alpha`, sound-speed terms), appended whenever the critical
   point lands in the inner region or jumps, for diagnosing a spurious inner
   critical point.  It is not part of the standard output and can be deleted.
