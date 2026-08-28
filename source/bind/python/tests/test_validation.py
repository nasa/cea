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


def test_eqsolver_preserves_double_precision_states(cea_module):
    reactants, products, weights = _basic_reactants(cea_module)
    solver = cea_module.EqSolver(products, reactants=reactants)
    temperature = 2359.2343858155864
    pressure = 7.242735337192802
    soln = cea_module.EqSolution(solver)
    solver.solve(soln, cea_module.TP, temperature, pressure, weights)
    assert soln.converged
    assert soln.T == temperature
    assert soln.P == pressure


def test_rocket_solver_rejects_conflicting_hc_tc(cea_module):
    reactants_mix, products_mix, weights = _basic_reactants(cea_module)
    solver = cea_module.RocketSolver(products_mix, reactants=reactants_mix)
    soln = cea_module.RocketSolution(solver)
    with pytest.raises(ValueError):
        solver.solve(soln, weights, pc=10.0, pi_p=2.0, hc=100.0, tc=3000.0)


@pytest.mark.parametrize("iac", [True, False])
def test_rocket_solver_rejects_weights_for_omitted_reactants(cea_module, iac):
    _, products, weights = _basic_reactants(cea_module)
    solver = cea_module.RocketSolver(products)
    soln = cea_module.RocketSolution(solver)
    assert solver.num_reactants > len(weights)

    with pytest.raises(ValueError, match="reactants=reac") as exc:
        solver.solve(soln, weights, pc=10.0, pi_p=2.0, tc=3000.0, iac=iac, ac_at=2.0)
    assert f"expected {solver.num_reactants}, got {len(weights)}" in str(exc.value)


@pytest.mark.parametrize("iac", [True, False])
@pytest.mark.parametrize("weight_count", [0, 1, 3])
def test_rocket_solver_rejects_wrong_weight_count(cea_module, iac, weight_count):
    reactants, products, weights = _basic_reactants(cea_module)
    solver = cea_module.RocketSolver(products, reactants=reactants)
    soln = cea_module.RocketSolution(solver)
    invalid_weights = np.zeros(weight_count)
    count = min(weight_count, len(weights))
    invalid_weights[:count] = weights[:count]

    with pytest.raises(ValueError, match=f"expected 2, got {weight_count}"):
        solver.solve(soln, invalid_weights, pc=10.0, pi_p=2.0, tc=3000.0, iac=iac, ac_at=2.0)


@pytest.mark.parametrize("weights", [None, np.array(1.0), np.ones((1, 2)), np.ones((2, 1))])
def test_rocket_solver_rejects_non_vector_weights(cea_module, weights):
    reactants, products, _ = _basic_reactants(cea_module)
    solver = cea_module.RocketSolver(products, reactants=reactants)
    soln = cea_module.RocketSolution(solver)

    with pytest.raises(ValueError, match="weights must be a one-dimensional array"):
        solver.solve(soln, weights, pc=10.0, pi_p=2.0, tc=3000.0)


@pytest.mark.parametrize("iac", [True, False])
def test_rocket_solver_defaults_to_product_reactants(cea_module, iac):
    reactants, products, weights = _basic_reactants(cea_module)
    product_weights = np.zeros(products.num_species)
    for species, weight in zip(reactants.species_names, weights):
        product_weights[products.species_names.index(species)] = weight

    default_solver = cea_module.RocketSolver(products)
    explicit_solver = cea_module.RocketSolver(products, reactants=products)
    default_soln = cea_module.RocketSolution(default_solver)
    explicit_soln = cea_module.RocketSolution(explicit_solver)
    for solver, soln in ((default_solver, default_soln), (explicit_solver, explicit_soln)):
        solver.solve(soln, product_weights, pc=10.0, pi_p=2.0, tc=3000.0, iac=iac, ac_at=2.0)
        assert soln.converged

    for prop in ("T", "P", "Isp", "Mach", "nj"):
        actual = getattr(default_soln, prop)
        assert np.all(np.isfinite(actual))
        np.testing.assert_array_equal(actual, getattr(explicit_soln, prop))


def test_shock_solver_requires_one_input(cea_module):
    reactants_mix, products_mix, weights = _basic_reactants(cea_module)
    solver = cea_module.ShockSolver(products_mix, reactants=reactants_mix)
    soln = cea_module.ShockSolution(solver, reflected=False)
    with pytest.raises(ValueError):
        solver.solve(soln, weights, T0=300.0, p0=1.0)
    with pytest.raises(ValueError):
        solver.solve(soln, weights, T0=300.0, p0=1.0, u1=1000.0, Mach1=2.0)
