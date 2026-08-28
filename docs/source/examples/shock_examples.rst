Shock Examples
==============

Here we describe the example shock problems.

Velocity and Mach conventions
-----------------------------

The Fortran, C, Python, MATLAB, and Excel shock interfaces report stations
``[initial, incident, reflected]``, numbered ``[1, 2, 5]`` in legacy CEA.
Speeds are positive magnitudes. The unshocked gas and reflecting wall are
stationary in the laboratory frame.

* ``velocity[0] = u1`` is the upstream gas speed in the incident-shock frame,
  also the incident wave speed in the laboratory frame.
* ``velocity[1] = u2`` is the downstream gas speed in the incident-shock frame.
  The laboratory-frame gas speed is ``v2 = u1 - u2``.
* ``velocity[2] = u5`` is the downstream gas speed in the reflected-shock frame.
  That gas is stationary at the wall, so ``u5`` also equals the reflected
  wave-speed magnitude in the laboratory frame (the wave travels back toward
  the unshocked gas). The upstream reflected-frame gas speed is
  ``u5_p_v2 = u5 + v2``. It is not ``u2``.

Thus the reflected mass balance is ``rho2*(v2 + u5) = rho5*u5`` and
``u5 = v2/(rho5/rho2 - 1)``. ``Mach`` is ``velocity/sonic_velocity`` in these
respective shock frames, not a vector of laboratory-frame Mach numbers.
The upstream reflected-shock Mach number would use ``(v2 + u5)/a2`` with
the appropriate upstream sound speed, not the reported downstream entry.

Sound speed is reported as ``sqrt(R*T*gamma_s/M)``. Initial and incident-frozen
states use frozen heat capacities; equilibrium states use the equilibrium
isentropic exponent. Reflected-frozen output retains the legacy incident-state
frozen heat-capacity basis and molecular weight, evaluated with the reflected
temperature in the sound-speed expression. For temperature-dependent heat
capacities this can differ from a locally evaluated frozen sound speed;
``Mach`` uses the reported convention. Retained last-valid states use the same
definition but remain unconverged and need not satisfy shock conservation.

``M21`` and ``M52`` are **molecular-weight ratios**, not Mach-number ratios.
The legacy ``M2/M1`` and ``M5/M2`` output labels and example variables retain
this meaning. Mach ratios, if needed, must be formed from ``Mach`` explicitly,
with due regard to the different reference frames.

.. toctree::
   :maxdepth: 1

   shock/example7
