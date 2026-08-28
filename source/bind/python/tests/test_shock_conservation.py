import numpy as np
import pytest


_CONSERVATION_RTOL = 5.0e-5


def _example7_case(cea):
    reac_names = ["H2", "O2", "Ar"]
    moles = np.array([0.05, 0.05, 0.9])
    p0 = cea.units.mmhg_to_bar(10.0)
    T0 = 300.0

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True)
    solver = cea.ShockSolver(prod, reactants=reac)
    soln = cea.ShockSolution(solver, reflected=True)
    weights = reac.moles_to_weights(moles)
    return solver, soln, weights, T0, p0


def _scaled_tolerance(lhs, rhs):
    return _CONSERVATION_RTOL * max(abs(lhs), abs(rhs), 1.0)


def _assert_conserved(lhs, rhs):
    residual = lhs - rhs
    assert abs(residual) <= _scaled_tolerance(lhs, rhs)


@pytest.mark.parametrize(
    ("mode", "u1", "kwargs"),
    [
        ("equilibrium", 1100.0, {"reflected": True}),
        ("equilibrium", 1400.0, {"reflected": True}),
        ("frozen", 1100.0, {"reflected": True, "incident_frozen": True, "reflected_frozen": True}),
        ("frozen", 1400.0, {"reflected": True, "incident_frozen": True, "reflected_frozen": True}),
        ("equilibrium/frozen", 1250.0, {"reflected": True, "reflected_frozen": True}),
        ("frozen/equilibrium", 1250.0, {"reflected": True, "incident_frozen": True}),
    ],
)
def test_shock_conservation_equations(cea_module, mode, u1, kwargs):
    cea = cea_module
    solver, soln, weights, T0, p0 = _example7_case(cea)

    solver.solve(soln, weights, T0, p0, u1=u1, **kwargs)

    assert soln.converged, f"{mode} shock solve failed to converge for u1={u1}"
    np.testing.assert_allclose(soln.Mach, soln.velocity / soln.sonic_velocity, rtol=1.0e-13)

    w1 = soln.velocity[0]
    w2 = soln.velocity[1]
    ur = soln.velocity[2]
    v2 = soln.v2

    rho1 = soln.density[0]
    rho2 = soln.density[1]
    rho5 = soln.density[2]

    p1 = soln.P[0] * 1.0e5
    p2 = soln.P[1] * 1.0e5
    p5 = soln.P[2] * 1.0e5

    h1 = soln.enthalpy[0]
    h2 = soln.enthalpy[1]
    h5 = soln.enthalpy[2]

    # Incident shock conservation in the incident shock frame.
    incident_mass_lhs = rho1 * w1
    incident_mass_rhs = rho2 * w2
    _assert_conserved(incident_mass_lhs, incident_mass_rhs)

    incident_momentum_lhs = p1 + rho1 * w1**2
    incident_momentum_rhs = p2 + rho2 * w2**2
    _assert_conserved(incident_momentum_lhs, incident_momentum_rhs)

    incident_energy_lhs = h1 + w1**2 / 2000.0
    incident_energy_rhs = h2 + w2**2 / 2000.0
    _assert_conserved(incident_energy_lhs, incident_energy_rhs)

    # Reflected shock conservation in the reflected-shock frame.
    reflected_upstream_speed = v2 + ur
    reflected_downstream_speed = ur

    reflected_mass_lhs = rho2 * reflected_upstream_speed
    reflected_mass_rhs = rho5 * reflected_downstream_speed
    _assert_conserved(reflected_mass_lhs, reflected_mass_rhs)

    reflected_momentum_lhs = p2 + rho2 * reflected_upstream_speed**2
    reflected_momentum_rhs = p5 + rho5 * reflected_downstream_speed**2
    _assert_conserved(reflected_momentum_lhs, reflected_momentum_rhs)

    reflected_energy_lhs = h2 + reflected_upstream_speed**2 / 2000.0
    reflected_energy_rhs = h5 + reflected_downstream_speed**2 / 2000.0
    _assert_conserved(reflected_energy_lhs, reflected_energy_rhs)

    assert soln.u5_p_v2 == pytest.approx(soln.v2 + soln.velocity[2], rel=1.0e-10, abs=1.0e-10)


@pytest.mark.parametrize("incident_frozen", [False, True])
@pytest.mark.parametrize("reflected,reflected_frozen", [(False, False), (True, False), (True, True)])
@pytest.mark.parametrize("use_mach", [False, True])
def test_shock_mach_reproducer(cea_module, incident_frozen, reflected, reflected_frozen, use_mach):
    cea = cea_module
    solver, _, weights, _, _ = _example7_case(cea)
    soln = cea.ShockSolution(solver, reflected=reflected)
    kwargs = dict(reflected=reflected, incident_frozen=incident_frozen, reflected_frozen=reflected_frozen)
    p0 = 5.0 * 1.01325 / 760.0
    solver.solve(soln, weights, 250.0, p0, u1=1250.0, **kwargs)
    expected_mach = soln.velocity / soln.sonic_velocity
    expected_T, expected_P = soln.T.copy(), soln.P.copy()

    # Reuse the same solution, using either input representation.
    speed = {"Mach1": soln.Mach[0]} if use_mach else {"u1": 1250.0}
    solver.solve(soln, weights, 250.0, p0, **speed, **kwargs)
    assert soln.converged
    np.testing.assert_allclose(soln.Mach, expected_mach, rtol=1.0e-10)
    np.testing.assert_allclose(soln.Mach, soln.velocity / soln.sonic_velocity, rtol=1.0e-13)
    np.testing.assert_allclose(soln.T, expected_T, rtol=1.0e-10)
    np.testing.assert_allclose(soln.P, expected_P, rtol=1.0e-10)
    assert soln.M21 == pytest.approx(soln.M[1] / soln.M[0], rel=1.0e-13)
    if reflected:
        assert soln.M52 == pytest.approx(soln.M[2] / soln.M[1], rel=1.0e-13)
        if reflected_frozen:
            # Finalization deliberately replaces the local iteration sound
            # speed with the legacy incident frozen-Cp basis. Mach must use it.
            cp2 = soln.cp_fr[1]
            gamma5 = cp2 / (cp2 - cea.R / (1000.0 * soln.M[1]))
            a5 = np.sqrt(cea.R * soln.T[2] * gamma5 / soln.M[1])
            assert soln.sonic_velocity[2] == pytest.approx(a5, rel=1.0e-12)
            assert soln.Mach[2] == pytest.approx(soln.velocity[2] / a5, rel=1.0e-12)
    if not incident_frozen and reflected and not reflected_frozen:
        np.testing.assert_allclose(soln.Mach, [4.1715730754, 0.6895985501, 0.6629497971], rtol=2.0e-8)


def _normal_shock(mach, gamma):
    # Calorically perfect normal shock, independent of CEA fields:
    # https://www.grc.nasa.gov/WWW/k-12/airplane/normal.html
    m2 = mach**2
    compression = (gamma + 1.0) * m2 / ((gamma - 1.0) * m2 + 2.0)
    pressure = (2.0 * gamma * m2 - (gamma - 1.0)) / (gamma + 1.0)
    downstream_mach = np.sqrt(((gamma - 1.0) * m2 + 2.0) / (2.0 * gamma * m2 - (gamma - 1.0)))
    return compression, pressure, downstream_mach


@pytest.mark.parametrize("incident_frozen", [False, True])
@pytest.mark.parametrize("reflected,reflected_frozen", [(False, False), (True, False), (True, True)])
@pytest.mark.parametrize("use_mach", [False, True])
def test_monatomic_normal_shock(cea_module, incident_frozen, reflected, reflected_frozen, use_mach):
    cea = cea_module
    argon = cea.Mixture(["Ar"])
    solver = cea.ShockSolver(argon, reactants=argon)
    soln = cea.ShockSolution(solver, reflected=reflected)
    weights = argon.moles_to_weights(np.array([1.0]))
    gamma, mach1, T1, p1 = 5.0 / 3.0, 3.0, 300.0, 1.0
    a1 = np.sqrt(gamma * cea.R * T1 / weights.sum())
    u1 = mach1 * a1
    r21, p21, mach2 = _normal_shock(mach1, gamma)
    T2 = T1 * p21 / r21
    a2 = a1 * np.sqrt(T2 / T1)
    v2 = u1 * (1.0 - 1.0 / r21)
    expected_mach = [mach1, mach2]
    expected_velocity = [u1, u1 / r21]
    expected_sound = [a1, a2]
    expected_T, expected_P = [T1, T2], [p1, p1 * p21]
    if reflected:
        # Station 5 is at rest at the wall. In the reflected frame,
        # w2r - w5r = v2, with w2r = Mr*a2 and w5r = w2r/r52.
        b = (gamma + 1.0) * v2 / (2.0 * a2)
        mr = (b + np.sqrt(b**2 + 4.0)) / 2.0
        r52, p52, mach5 = _normal_shock(mr, gamma)
        expected_mach.append(mach5)
        expected_velocity.append(mr * a2 / r52)
        expected_sound.append(a2 * np.sqrt(p52 / r52))
        expected_T.append(T2 * p52 / r52)
        expected_P.append(p1 * p21 * p52)

    speed = {"Mach1": mach1} if use_mach else {"u1": u1}
    solver.solve(soln, weights, T1, p1, **speed, reflected=reflected,
                 incident_frozen=incident_frozen, reflected_frozen=reflected_frozen)
    assert soln.converged
    for actual, expected in ((soln.Mach, expected_mach), (soln.velocity, expected_velocity),
                             (soln.sonic_velocity, expected_sound), (soln.T, expected_T), (soln.P, expected_P)):
        np.testing.assert_allclose(actual, expected, rtol=_CONSERVATION_RTOL)
    assert soln.M21 == pytest.approx(1.0, abs=1.0e-12)
    if reflected:
        assert soln.M52 == pytest.approx(1.0, abs=1.0e-12)
        assert soln.u5_p_v2 == pytest.approx(mr * a2, rel=_CONSERVATION_RTOL)
