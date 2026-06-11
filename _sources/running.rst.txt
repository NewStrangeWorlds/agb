===============
Running a model
===============

Invocation
==========

The executable takes exactly one argument, the **model folder**:

.. code-block:: console

   $ ./agb_wind <model_folder>

For example, with the reference IRC+10216 setup shipped in the repository:

.. code-block:: console

   $ ./agb_wind IRC10216/

The model folder must contain a :file:`config.toml` (see
:doc:`configuration`).  All other input and output paths in that file are
interpreted **relative to the model folder**.

Parallelism is controlled through OpenMP; set ``OMP_NUM_THREADS`` to choose the
number of threads.  The radiative transfer dominates the runtime (~90 %), and is
parallelised over spectral points and grid points.

.. note::

   Results are **bit-reproducible run-to-run for a fixed thread count**.  They
   may differ at the last few digits when the thread count changes, because the
   parallel floating-point reductions sum in a different order; the converged
   physical solution is the same.

Model folder layout
====================

A typical model folder (the IRC+10216 example) looks like:

.. code-block:: text

   IRC10216/
   ├── config.toml                 # the model configuration (required)
   ├── fastchem_parameters.dat     # FastChem element-abundance/parameter file
   ├── atmosphere_converged.dat    # starting structure (or use "grey")
   ├── atmosphere.dat              # output: converged structure
   ├── spectrum.dat                # output: emergent spectrum
   ├── dust.dat                    # output: dust distribution
   └── hydro.dat                   # output: wind structure

External data that is *not* inside the model folder:

* the **opacity database** (HELIOS-K cross sections), at the absolute path given
  by ``[opacity] path``;
* the **dust refractive indices** file, at ``[dust] refractive_index_file``.

Starting model
==============

The ``[model] starting_model`` key selects how the initial structure is built:

* **A path** (e.g. ``"atmosphere_converged.dat"``) — read an existing structure
  file.  This is the fastest route when restarting from a near-converged model.
* **The literal string** ``"grey"`` — build a logarithmic radial grid and a
  hydrostatic **Lucy grey** starting structure on the fly, with no input file.
  The grid spans ``[grid] inner_radius`` … ``outer_radius`` (in units of
  :math:`R_\star`) with ``[grid] n_points`` cells; see
  :doc:`theory/hydrodynamics`.

.. _running-structure-file:

Structure file format
=====================

The atmospheric structure file (both the input starting model and the
:file:`atmosphere.dat` output) is a whitespace-separated table with **two header
lines** followed by one row per radial grid point.  The columns are:

.. list-table::
   :header-rows: 1
   :widths: 8 18 74

   * - #
     - Column
     - Meaning
   * - 1
     - ``r/R*``
     - Radius in units of the stellar radius (the normalised radial grid).
   * - 2
     - ``r(cm)``
     - Radius in cm.
   * - 3
     - ``rho(g/cm3)``
     - Gas mass density.
   * - 4
     - ``p(dyn/cm2)``
     - Gas pressure (cgs).
   * - 5
     - ``T_gas(K)``
     - Gas temperature.
   * - 6
     - ``T_dust(K)``
     - Dust temperature.
   * - 7
     - ``v(cm/s)``
     - Wind velocity.

The grid is ordered from the inner boundary (index 0) outward.  The first header
line of the output additionally records the stellar parameters and the
mass-loss rate.

Output files
============

Each output is written only if the corresponding ``[output]`` key is set (a
missing key, an empty string, or ``"none"`` disables it).  Paths are relative to
the model folder.

* ``[output] atmosphere`` — the converged structure, same format as the input
  (:ref:`running-structure-file`).
* ``[output] spectrum`` — the emergent spectrum (wavelength vs. flux); see
  :doc:`implementation/io_formats`.
* ``[output] dust`` — the dust distribution: degree of condensation, the carbon
  species number densities, grain number density and radius, nucleation rate
  :math:`J_\star`, growth timescale, and the moments :math:`K_0\dots K_3`.
* ``[output] hydro`` — the wind structure (velocity, density, sound speed,
  radiative acceleration :math:`\alpha`, flux-mean extinction).

The outputs are written once, after the global iteration has converged (or after
the maximum number of iterations is reached).

Console diagnostics
===================

During a run the code prints, per global iteration:

* the chemistry–hydrodynamics sub-iteration progress and the maximum change of
  :math:`\alpha`;
* the corrector in use (``Unsöld–Lucy`` or the latched ``linearisation``);
* the maximum relative temperature change, the flux-convergence deviation from
  the target luminosity, and the gas/dust energy-balance (radiative-equilibrium)
  residuals.

These diagnostics are the primary way to monitor convergence; see
:doc:`theory/numerics` for how to read them.
