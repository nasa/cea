import numpy as np
import pytest


def test_shock_last_valid_solution_warns_and_preserves_state(cea_module):
    cea = cea_module

    reac = cea.Mixture(["N2 ", "CH4"], ions=True)
    prod = cea.Mixture(["N2 ", "CH4"], products_from_reactants=True, ions=True)
    solver = cea.ShockSolver(prod, reactants=reac, ions=True, transport=True, trace=1.0e-10)
    soln = cea.ShockSolution(solver, reflected=False)

    weights = reac.moles_to_weights(np.array([0.985, 0.015]))
    p0 = cea.units.mmhg_to_bar(0.224)
    T0 = 291.0

    with pytest.warns(
        RuntimeWarning,
        match=r"CEA_LAST_VALID_SOLUTION.*last valid shock state retained.*converged is False",
    ):
        solver.solve(soln, weights, T0, p0, u1=1900.0, reflected=False)

    assert soln.last_error == cea.LAST_VALID_SOLUTION
    assert not soln.converged
    assert soln.P[1] > 0.0
    assert soln.T[1] > 0.0
    assert soln.velocity[1] > 0.0
    assert soln.sonic_velocity[1] > 0.0
