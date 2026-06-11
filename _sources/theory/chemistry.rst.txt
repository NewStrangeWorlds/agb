==============================
Gas-phase chemistry (FastChem)
==============================

*Module:* ``src/chemistry/`` — :cpp:class:`agb::FastChemChemistry`, wrapping
`FastChem <https://github.com/exoclime/fastchem>`__ 4.0.2.

The equilibrium gas-phase chemical composition is computed with **FastChem**, a
fast equilibrium chemistry solver.  At each radial grid point, given the gas
temperature and pressure (and the elemental abundances, with the prescribed C/O
ratio), FastChem returns the number densities of all gas-phase species, the mean
molecular weight, and the total element and hydrogen number densities.

Interface
=========

:cpp:func:`agb::FastChemChemistry::calcChemicalComposition` takes the
temperature and pressure profiles and a per-layer **degree of carbon
condensation** and fills:

* ``number_densities[i][species]`` — per-layer species number densities;
* ``mean_molecular_weight`` — used in the equation of state and the isothermal
  sound speed;
* ``total_element_density`` and ``total_h_density`` — reference densities used to
  normalise the dust moments and the condensable-carbon budget.

The species the rest of the model needs (the carbon molecules
:math:`\mathrm{C}, \mathrm{C_2}, \mathrm{C_2H}, \mathrm{C_2H_2}` that feed dust
formation, plus the opacity carriers) are indexed through
``fastchem_species_indices`` / ``fastchem_element_indices`` and the
``chem_species.h`` identifiers.

Carbon depletion coupling
=========================

The C/O ratio sets the **condensable carbon** abundance
:math:`\varepsilon_C-\varepsilon_O` (all oxygen is assumed locked into the very
stable CO molecule).  Only this excess carbon is available to form
amorphous-carbon dust.

In the current scheme the gas-phase chemistry is evaluated at the **full
(un-depleted)** element abundances; the carbon actually consumed by dust is *not*
fed back into a separate FastChem call.  Instead the Gail & Sedlmayr moment sweep
depletes the growth species differentially as it integrates outward (see
:doc:`dust`).  This avoids a fragile outer chemistry↔dust fixed point — the
moment method is memoryless, so an outer Picard on the degree of condensation can
oscillate.

.. note::

   The ``degree_of_condensation_c`` argument to
   :cpp:func:`calcChemicalComposition` exists so that a future, more detailed
   coupling (re-equilibrating the free C ↔ C₂ ↔ C₂H₂ subsystem at the depleted
   carbon abundance — "B2" in the development notes) can be enabled cheaply.  In
   the default path it is passed as zero.

Equation of state
=================

After the composition is known, :cpp:func:`agb::Atmosphere::equationOfState`
closes the thermodynamics (pressure ↔ density ↔ mean molecular weight), and the
isothermal sound speed used by the hydrodynamics follows from

.. math::

   c_T = \sqrt{\frac{k_B\,T_\mathrm{gas}}{\mu\,m_p}} ,

with :math:`\mu` the mean molecular weight per particle.
