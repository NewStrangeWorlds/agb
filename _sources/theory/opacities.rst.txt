==================================
Opacities (transport coefficients)
==================================

*Module:* ``src/transport_coefficients/`` —
:cpp:class:`agb::TransportCoefficients`, :cpp:class:`agb::OpacitySpecies`,
:cpp:class:`agb::SampledData`.

The radiative transfer needs the **absorption** and **scattering** coefficients
of the gas and the dust at every radius and every spectral point.  These are
assembled separately and summed:

.. math::

   \chi_\nu = \underbrace{\kappa_{\nu,\mathrm{gas}}+\sigma_{\nu,\mathrm{gas}}}_{\text{gas}}
            + \underbrace{\kappa_{\nu,\mathrm{dust}}+\sigma_{\nu,\mathrm{dust}}}_{\text{dust}} ,

with the gas/dust split kept explicitly (``absorption_coeff_gas`` /
``…_dust`` etc.) because the temperature corrector and the energy balance act on
each component's **absorption** separately (:doc:`temperature_correction`).

Gas opacities
=============

:cpp:func:`agb::TransportCoefficients::calculate` loops over the configured gas
species (one :cpp:class:`agb::OpacitySpecies` each) and accumulates their
contributions for the local temperature, pressure and number densities.  Three
kinds of opacity source are handled:

* **Molecular/atomic line absorption** — pre-computed cross sections sampled from
  a database (HELIOS-K), one folder per species (e.g. CO, H₂O, CO₂, C₂H₂).  The
  cross sections are tabulated on a pressure–temperature grid; for the local
  conditions the code locates the bracketing tabulated points
  (``findClosestDataPoints``) and interpolates (``calcAbsorptionCrossSections``).
  The :cpp:class:`agb::SampledData` objects hold the sampled cross sections on the
  global wavenumber grid.

* **Continuum / collision-induced absorption (CIA)** — e.g. H₂–H₂, H₂–He, H–He
  (``calcContinuumAbsorption``), and the negative ions **H⁻** and **H₂⁻**
  (``h_m.cpp``, ``h2_m.cpp``).  These have no line list ("``none``" folder) and
  are evaluated from built-in physics; their density depends on collision
  partners (``cia_collision_partner``).

* **Rayleigh scattering** — wavelength-dependent scattering cross sections per
  species (``species_rayleigh_cross_sections.cpp``,
  ``calcScatteringCrossSections`` / ``generalRayleighCrossSection`` with the
  King correction factor).

The species list, the database path, and the per-species folder are set in the
``[opacity]`` block (:doc:`../configuration`).  A folder of ``"none"`` marks a
species whose opacity is computed internally rather than read from a line list.

Spectral grid and sampling
==========================

*Module:* ``src/spectral_grid/`` — :cpp:class:`agb::SpectralGrid`.

The opacity database is tabulated on a fixed, high-resolution global **wavenumber
grid** (read from ``wavenumber_full.dat``).  The model builds a working
**constant-resolving-power** grid between ``min_wavelength`` and
``max_wavelength`` at the requested ``resolution``
(``createConstantResolutionGrid``), and selects the nearest database points
(``index_list``).  All radiative-transfer quantities live on this working grid;
helpers convert between wavelength (μm and cm) and wavenumber and interpolate data
onto either grid.

Dust opacities
==============

The dust absorption and scattering coefficients come from **Mie theory** applied
to the grain-size information from the dust module (:doc:`dust`).  Given the
complex refractive index of the dust material (``[dust] refractive_index_file``,
e.g. amorphous carbon) and the grain radius at each layer, the LX-MIE routine
returns the extinction/scattering efficiencies, which are scaled by the grain
geometric cross section and number density
(:cpp:func:`agb::GailSedlmayrDust::calcTransportCoefficients`).  Because the dust
opacity does not depend on temperature, it is computed **once** per temperature
iteration, not per Unsöld–Lucy sub-step.

.. note::

   A non-finite or non-positive mean grain radius would overflow the Mie series
   (a ``length_error`` was observed historically when a NaN ⟨a⟩ propagated from a
   collapsed density).  The dust module floors the grain radius and number
   density to guard against this.
