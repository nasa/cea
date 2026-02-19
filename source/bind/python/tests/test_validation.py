import numpy as np
import pytest


def _basic_reactants(cea_module):
    reactants = ["H2", "O2"]
    reactants_mix = cea_module.Mixture(reactants)
    products_mix = cea_module.Mixture(reactants, products_from_reactants=True)
    weights = reactants_mix.moles_to_weights(np.array([2.0, 1.0], dtype=np.float64))
    return reactants_mix, products_mix, weights


def test_eqsolver_rejects_non_mixture_reactants(cea_module):
    products = cea_module.Mixture(["H2", "O2"], products_from_reactants=True)
    with pytest.raises(TypeError):
        cea_module.EqSolver(products, reactants=["H2", "O2"])


def test_rocket_solver_rejects_conflicting_hc_tc(cea_module):
    reactants_mix, products_mix, weights = _basic_reactants(cea_module)
    solver = cea_module.RocketSolver(products_mix, reactants=reactants_mix)
    soln = cea_module.RocketSolution(solver)
    with pytest.raises(ValueError):
        solver.solve(soln, weights, pc=10.0, pi_p=2.0, hc=100.0, tc=3000.0)


def test_shock_solver_requires_one_input(cea_module):
    reactants_mix, products_mix, weights = _basic_reactants(cea_module)
    solver = cea_module.ShockSolver(products_mix, reactants=reactants_mix)
    soln = cea_module.ShockSolution(solver, reflected=False)
    with pytest.raises(ValueError):
        solver.solve(soln, weights, T0=300.0, p0=1.0)
    with pytest.raises(ValueError):
        solver.solve(soln, weights, T0=300.0, p0=1.0, u1=1000.0, Mach1=2.0)
