=====================
Theory and modelling
=====================

This part documents the physics and the numerical methods of each module, in
roughly the order in which they appear in the solution loop.

.. toctree::
   :maxdepth: 2

   governing_equations
   radiative_transfer
   temperature_correction
   chemistry
   dust
   opacities
   hydrodynamics
   numerics

Notation
========

Throughout, :math:`r` is the radial coordinate, :math:`R_\star` the stellar
reference radius, and a subscript :math:`\star` denotes a stellar quantity.  The
radiation field is described by its frequency-dependent moments — mean intensity
:math:`J_\nu`, Eddington (first) flux :math:`H_\nu`, and second moment
:math:`K_\nu` — and the variable Eddington factor :math:`f_\nu=K_\nu/J_\nu`.
The extinction (total), absorption and scattering coefficients are
:math:`\chi_\nu=\kappa_\nu+\sigma_\nu`.  The Planck function is
:math:`B_\nu(T)`.  Dust quantities use the Gail & Sedlmayr moments
:math:`K_j` (with :math:`K_3` the condensed volume per hydrogen nucleus) and the
nucleation rate :math:`J_\star`.

Equation numbers such as *(2.59)* refer to the diploma thesis of D. Kitzmann,
which is the primary reference for the radiative-transfer scheme.
