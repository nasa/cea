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
            return soln, caught

    details = ", ".join(
        f"u1={u1:.0f}: last_error={last_error}, converged={converged}"
        for u1, last_error, converged in outcomes
    )
    pytest.fail(f"No LAST_VALID_SOLUTION case found across candidate velocities. {details}")


def test_shock_last_valid_solution_warns_and_preserves_state(cea_module):
    cea = cea_module
    soln, caught = _find_last_valid_incident_case(cea)

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
