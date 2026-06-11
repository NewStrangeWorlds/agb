========
Glossary
========

.. glossary::

   AGB
      Asymptotic giant branch. A late evolutionary stage of low- and
      intermediate-mass stars, characterised by a cool, extended atmosphere and
      strong mass loss.

   Carbon star
      A star whose photospheric carbon-to-oxygen ratio exceeds unity
      (:term:`C/O ratio` > 1), so that carbon is available to form amorphous-carbon
      dust after oxygen is locked in CO.

   C/O ratio
      The carbon-to-oxygen elemental abundance ratio. Sets the condensable
      carbon :math:`\varepsilon_C-\varepsilon_O`.

   Critical point / sonic point
      The radius :math:`r_c` where the wind velocity equals the isothermal sound
      speed, :math:`v=c_T`. The wind equation is singular there; the Melia
      :math:`\Phi` transform removes the singularity.

   Degree of condensation
      :math:`f_c`, the fraction of condensable carbon locked into dust grains;
      :math:`f_c=K_3/c_\mathrm{cond}`.

   Eddington factor
      :math:`f_\nu=K_\nu/J_\nu`, the closure of the radiation moment system in
      the variable-Eddington-factor method.

   Equidistribution
      The principle behind the movable grid: place nodes so that a monitor
      function has equal integral in every cell.

   Feautrier method
      A second-order formulation of radiative transfer along a ray in terms of
      the mean (:math:`u`) and difference (:math:`v`) of the in/out intensities,
      solved as a tridiagonal system.

   Flux-mean extinction
      :math:`\chi_F`, the flux-weighted average of the extinction coefficient;
      the single radiation quantity that sets the radiative acceleration
      :math:`\alpha`.

   Henyey method
      A global Newton–Raphson relaxation of a discretised boundary-value problem
      on all grid nodes simultaneously. Here an optional wind structure solver.

   Impact parameter
      A tangent ray at perpendicular distance :math:`p` from the centre; the
      natural coordinate for spherical radiative transfer.

   Mass-loss rate
      :math:`\dot{M}=4\pi r^2\rho v`; either prescribed or obtained as an
      eigenvalue from the critical-point regularity condition.

   Moment method (dust)
      Tracking the moments :math:`K_j` of the grain-size distribution instead of
      the full distribution; due to Gail & Sedlmayr.

   Nucleation rate
      :math:`J_\star`, the rate of formation of seed clusters per unit volume and
      time (classical nucleation theory).

   Radiative equilibrium
      The condition that each absorbing component emits as much as it absorbs,
      :math:`\int\kappa_\nu(J_\nu-B_\nu)\,\mathrm{d}\nu=0`; fixes the temperature.

   Sphericality factor
      :math:`q_\nu`, the geometric factor in the spherical moment equation that
      accounts for the curvature of the radiation field.

   Unsöld–Lucy correction
      A robust, approximate temperature correction toward radiative equilibrium,
      combining the local energy balance and the flux constancy.

   Variable Eddington factor (VEF)
      The radiative-transfer method in which the moment system is closed with
      Eddington factors computed from a formal solution and iterated to
      consistency.
