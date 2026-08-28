import re
import warnings

import numpy as np
import pytest


_LAST_VALID_U1_CANDIDATES = (
    1200.0,
    1300.0,
    1400.0,
    1500.0,
    1600.0,
    1700.0,
    1800.0,
    1900.0,
    2000.0,
    2100.0,
    2200.0,
    2300.0,
)


def _find_last_valid_incident_case(cea):
    reac = cea.Mixture(["N2 ", "CH4"], ions=True)
    prod = cea.Mixture(["N2 ", "CH4"], products_from_reactants=True, ions=True)
    solver = cea.ShockSolver(prod, reactants=reac, ions=True, transport=True, trace=1.0e-10)

    weights = reac.moles_to_weights(np.array([0.985, 0.015]))
    p0 = cea.units.mmhg_to_bar(0.224)
    T0 = 291.0
    outcomes = []

    for u1 in _LAST_VALID_U1_CANDIDATES:
        soln = cea.ShockSolution(solver, reflected=False)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            solver.solve(soln, weights, T0, p0, u1=u1, reflected=False)
        outcomes.append((u1, soln.last_error, soln.converged))
        if soln.last_error == cea.LAST_VALID_SOLUTION:
            return solver, soln, weights, T0, p0, u1, caught

    details = ", ".join(
        f"u1={u1:.0f}: last_error={last_error}, converged={converged}"
        for u1, last_error, converged in outcomes
    )
    pytest.fail(f"No LAST_VALID_SOLUTION case found across candidate velocities. {details}")


def test_shock_last_valid_solution_warns_and_preserves_state(cea_module):
    cea = cea_module
    solver, soln, weights, T0, p0, u1, caught = _find_last_valid_incident_case(cea)

    pattern = re.compile(
        r"CEA_LAST_VALID_SOLUTION.*last valid shock state retained.*converged is False"
    )
    assert any(
        isinstance(w.message, RuntimeWarning) and pattern.search(str(w.message))
        for w in caught
    )

    assert soln.last_error == cea.LAST_VALID_SOLUTION
    assert not soln.converged
    assert soln.P[1] > 0.0
    assert soln.T[1] > 0.0
    assert soln.velocity[1] > 0.0
    assert soln.sonic_velocity[1] > 0.0
    np.testing.assert_allclose(soln.Mach, soln.velocity / soln.sonic_velocity, rtol=1.0e-13)

    # A retained state is inspectable, not a converged normal shock. Do not
    # impose subsonic/conservation assertions on its unconverged downstream state.
    retained_mach = soln.Mach.copy()
    solver.solve(soln, weights, T0, p0, u1=1100.0, reflected=False, incident_frozen=True)
    assert soln.converged
    assert soln.last_error == cea.SUCCESS
    np.testing.assert_allclose(soln.Mach, soln.velocity / soln.sonic_velocity, rtol=1.0e-13)
    with pytest.warns(RuntimeWarning, match="CEA_LAST_VALID_SOLUTION"):
        solver.solve(soln, weights, T0, p0, u1=u1, reflected=False)
    assert not soln.converged
    assert soln.last_error == cea.LAST_VALID_SOLUTION
    np.testing.assert_array_equal(soln.Mach, retained_mach)


@pytest.mark.parametrize("frozen", [False, True])
@pytest.mark.parametrize("reflected", [False, True])
def test_shock_reuse_after_failure(cea_module, frozen, reflected):
    cea = cea_module
    argon = cea.Mixture(["Ar"])
    solver = cea.ShockSolver(argon, reactants=argon)
    soln = cea.ShockSolution(solver, reflected=reflected)
    weights = argon.moles_to_weights(np.array([1.0]))
    kwargs = dict(reflected=reflected, incident_frozen=frozen, reflected_frozen=frozen)

    solver.solve(soln, weights, 300.0, 1.0, Mach1=3.0, **kwargs)
    assert soln.converged
    expected = soln.Mach.copy()
    # Mach 12 fails at the reflected temperature cap; Mach 100 at the incident cap.
    failures = [(12.0, 2), (100.0, 1)] if reflected else [(100.0, 1)]
    for mach1, failed_index in failures:
        with pytest.warns(RuntimeWarning, match="CEA_NOT_CONVERGED"):
            solver.solve(soln, weights, 300.0, 1.0, Mach1=mach1, **kwargs)
        assert not soln.converged
        assert soln.last_error == cea.NOT_CONVERGED
        for values in (soln.Mach, soln.velocity, soln.sonic_velocity):
            np.testing.assert_array_equal(values[failed_index:], 0.0)
        np.testing.assert_allclose(
            soln.Mach[:failed_index],
            soln.velocity[:failed_index] / soln.sonic_velocity[:failed_index], rtol=1.0e-13,
        )
        solver.solve(soln, weights, 300.0, 1.0, Mach1=3.0, **kwargs)
        assert soln.converged
        assert soln.last_error == cea.SUCCESS
        np.testing.assert_array_equal(soln.Mach, expected)
