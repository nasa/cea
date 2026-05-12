import numpy as np
import pytest


def _wrapper_case_inputs():
    reactants = ["H2", "O2"]
    fuel_amounts = np.array([2.0, 0.0], dtype=np.float64)
    oxid_amounts = np.array([0.0, 1.0], dtype=np.float64)
    return reactants, fuel_amounts, oxid_amounts


def test_matlab_module_import_and_eq_solve(cea_module):
    import cea.matlab

    reactants, fuel_amounts, oxid_amounts = _wrapper_case_inputs()
    soln = cea.matlab.eq_solve(
        cea_module.TP,
        reactants,
        T=3000.0,
        P=cea_module.units.atm_to_bar(1.0),
        fuel_amounts=fuel_amounts,
        oxid_amounts=oxid_amounts,
        moles=True,
    )

    assert isinstance(soln, cea_module.EqSolution)
    assert soln.converged


def test_matlab_eq_solve_matches_direct_solver(cea_module):
    import cea.matlab

    reactants, fuel_amounts, oxid_amounts = _wrapper_case_inputs()
    wrapper_soln = cea.matlab.eq_solve(
        cea_module.TP,
        reactants,
        T=3000.0,
        P=cea_module.units.atm_to_bar(1.0),
        fuel_amounts=fuel_amounts,
        oxid_amounts=oxid_amounts,
        moles=True,
    )

    reactants_mix = cea_module.Mixture(reactants)
    products_mix = cea_module.Mixture(reactants, products_from_reactants=True)
    weights = reactants_mix.moles_to_weights(fuel_amounts)
    weights += reactants_mix.moles_to_weights(oxid_amounts)
    direct_solver = cea_module.EqSolver(products_mix, reactants=reactants_mix)
    direct_soln = cea_module.EqSolution(direct_solver)
    direct_solver.solve(
        direct_soln,
        cea_module.TP,
        3000.0,
        cea_module.units.atm_to_bar(1.0),
        weights,
    )

    assert wrapper_soln.converged == direct_soln.converged
    assert wrapper_soln.last_error == direct_soln.last_error
    np.testing.assert_allclose(wrapper_soln.T, direct_soln.T, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(wrapper_soln.P, direct_soln.P, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(wrapper_soln.MW, direct_soln.MW, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(
        wrapper_soln.nj,
        direct_soln.nj,
        rtol=0.0,
        atol=0.0,
    )


def test_root_eq_solve_shim_warns_and_forwards(cea_module):
    reactants, fuel_amounts, oxid_amounts = _wrapper_case_inputs()
    with pytest.warns(
        DeprecationWarning,
        match=r"cea\.eq_solve is deprecated; use cea\.matlab\.eq_solve instead\.",
    ):
        soln = cea_module.eq_solve(
            cea_module.TP,
            reactants,
            T=3000.0,
            P=cea_module.units.atm_to_bar(1.0),
            fuel_amounts=fuel_amounts,
            oxid_amounts=oxid_amounts,
            moles=True,
        )

    assert soln.converged
