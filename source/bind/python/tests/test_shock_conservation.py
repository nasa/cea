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
    ],
)
def test_shock_conservation_equations(cea_module, mode, u1, kwargs):
    cea = cea_module
    solver, soln, weights, T0, p0 = _example7_case(cea)

    solver.solve(soln, weights, T0, p0, u1=u1, **kwargs)

    assert soln.converged, f"{mode} shock solve failed to converge for u1={u1}"

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
