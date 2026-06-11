.. AGB Wind Model documentation master file

==========================================================
Atmospheric Model for Asymptotic Giant Branch Stars
==========================================================

This is the documentation of a self-consistent model for the **stationary,
spherically symmetric, dust-driven wind of a carbon-rich asymptotic giant branch
(AGB) star**.  The code couples

* a **spherical, frequency-dependent radiative transfer** solver using the
  Rybicki–Hummer variable-Eddington-factor moment method,
* a **radiative-equilibrium temperature** determination (Unsöld–Lucy correction
  and a full-linearisation Newton corrector),
* **equilibrium gas-phase chemistry** (FastChem),
* **carbon-dust nucleation and growth** following Gail & Sedlmayr, using the
  moment method with self-consistent carbon depletion, and
* a **stationary wind hydrodynamics** solver (Melia ``Φ`` transform with a
  mass-loss eigenvalue, plus an optional Henyey global-Newton structure solver),

iterated to a self-consistent stationary outflow.

.. note::

   The model is research software and under active development.  This
   documentation reflects the implementation on the ``main`` branch.  The
   underlying radiative-transfer scheme and notation follow the diploma thesis
   of D. Kitzmann; equation numbers of the form *(2.59)*, *(3.64)*, … below
   refer to that thesis.

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   overview
   installation
   running
   configuration

.. toctree::
   :maxdepth: 2
   :caption: Theory and modelling

   theory/index
   theory/governing_equations
   theory/radiative_transfer
   theory/temperature_correction
   theory/chemistry
   theory/dust
   theory/opacities
   theory/hydrodynamics
   theory/numerics

.. toctree::
   :maxdepth: 2
   :caption: Implementation

   implementation/architecture
   implementation/modules
   implementation/io_formats

.. toctree::
   :maxdepth: 1
   :caption: Appendix

   glossary
   references


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
