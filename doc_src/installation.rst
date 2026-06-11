=========================
Installation and building
=========================

Requirements
============

* A C++17 compiler (GCC or Clang) with **OpenMP** support.
* **CMake** ≥ 3.10.
* Network access at configure time: the build automatically downloads several
  dependencies through CMake ``FetchContent`` / ``ExternalProject``.

Bundled / fetched dependencies
==============================

The :file:`CMakeLists.txt` pulls in the following automatically:

.. list-table::
   :header-rows: 1
   :widths: 22 16 62

   * - Dependency
     - Version
     - Purpose
   * - `FastChem <https://github.com/exoclime/fastchem>`__
     - 4.0.2
     - Equilibrium gas-phase chemistry (built as ``fastchem_lib``).
   * - `LX-MIE <https://github.com/daniel-kitzmann/LX-MIE>`__
     - pinned commit
     - Mie theory for dust absorption/scattering efficiencies.
   * - `Eigen <https://gitlab.com/libeigen/eigen>`__
     - 3.4.0 (header-only)
     - Dense linear algebra for the Henyey structure solver.
   * - `toml++ <https://github.com/marzer/tomlplusplus>`__
     - 3.4.0 (header-only)
     - Parsing of the ``config.toml`` model configuration.
   * - `CppAD <https://github.com/coin-or/CppAD>`__
     - 20250000.3
     - Automatic differentiation for the Henyey Jacobian. Built and installed
       in isolation via ``ExternalProject`` into ``cppad-install/``.

FastChem internally uses Eigen, which leads to an important build constraint:

.. warning::

   **Do not add** ``-march=native`` **to the build flags.**  FastChem's bundled
   Eigen breaks (alignment / vectorisation faults) when compiled with it.  The
   Release configuration deliberately uses only ``-O3 -DNDEBUG -funroll-loops``,
   and ``-ffast-math`` is also omitted to keep the radiative-transfer results
   bit-reproducible.  If architecture-specific vectorisation is ever needed,
   scope it to individual radiative-transfer translation units — never apply it
   globally where FastChem/Eigen is compiled.

Building
========

The project follows the standard out-of-source CMake workflow:

.. code-block:: console

   $ mkdir build
   $ cd build
   $ cmake ..
   $ make -j

The first ``cmake ..`` will clone and build FastChem, LX-MIE, Eigen, toml++ and
CppAD; subsequent builds reuse them.  The executable is written **into the
source directory** as :file:`agb_wind` (the CMake ``RUNTIME_OUTPUT_DIRECTORY`` is
set to the source tree, and the runtime ``RPATH`` is pointed at
``cppad-install/lib`` so the loader finds ``libcppad_lib.so``).

The default and only configured build type is ``Release``.

Self-test
=========

A built-in self-test validates the Henyey structure solver's automatic-
differentiation Jacobian and Newton/eigenvalue machinery without needing a full
model:

.. code-block:: console

   $ ./agb_wind --selftest

This checks the analytic (CppAD) Jacobian against central finite differences for
the ``Φ`` equation, the coupled ``Φ``+dust-moment system and the mass-loss
eigenvalue row, runs a transonic eigenvalue solve, and exercises the grey
(Lucy) bootstrap.  It exits ``0`` on success.

Running a model
===============

Pass the path to a model folder (containing a :file:`config.toml`) as the single
command-line argument:

.. code-block:: console

   $ ./agb_wind IRC10216/

See :doc:`running` for the directory layout and :doc:`configuration` for the
full configuration reference.
