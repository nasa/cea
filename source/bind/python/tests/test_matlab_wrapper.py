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

    assert not isinstance(soln, cea_module.EqSolution)
    assert not hasattr(soln, "solver")
    assert soln.converged


def test_libcea_eq_solve_returns_eqsolution(cea_module):
    from cea.lib import libcea

    reactants, fuel_amounts, oxid_amounts = _wrapper_case_inputs()
    soln = libcea.eq_solve(
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


def test_matlab_eq_solve_matches_compiled_eq_solve(cea_module):
    import cea.matlab
    from cea.lib import libcea

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

    compiled_soln = libcea.eq_solve(
        cea_module.TP,
        reactants,
        T=3000.0,
        P=cea_module.units.atm_to_bar(1.0),
        fuel_amounts=fuel_amounts,
        oxid_amounts=oxid_amounts,
        moles=True,
    )

    assert wrapper_soln.converged == compiled_soln.converged
    assert wrapper_soln.last_error == compiled_soln.last_error
    np.testing.assert_allclose(wrapper_soln.T, compiled_soln.T, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(wrapper_soln.P, compiled_soln.P, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(wrapper_soln.MW, compiled_soln.MW, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(
        wrapper_soln.enthalpy,
        compiled_soln.enthalpy,
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        wrapper_soln.cp_eq,
        compiled_soln.cp_eq,
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        wrapper_soln.nj,
        compiled_soln.nj,
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        wrapper_soln.ln_nj,
        compiled_soln.ln_nj,
        rtol=0.0,
        atol=0.0,
    )
    assert wrapper_soln.mass_fractions == compiled_soln.mass_fractions
    assert wrapper_soln.mole_fractions == compiled_soln.mole_fractions


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

    assert isinstance(soln, cea_module.EqSolution)
    assert soln.converged
