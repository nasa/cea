import numpy as np
import pytest

# Field names taken from the `property` blocks on each compiled solution
# class in CEA.pyx (not from cea.matlab's own field lists), so these checks
# don't just validate the wrapper's field lists against themselves.
_EQ_ALL_FIELDS = (
    "T", "P", "volume", "density", "M", "MW", "enthalpy", "energy", "entropy",
    "gibbs_energy", "gamma_s", "cp_fr", "cp_eq", "cp", "cv_fr", "cv_eq", "cv",
    "viscosity", "conductivity_fr", "conductivity_eq", "Pr_fr", "Pr_eq",
    "nj", "ln_nj", "n",
)

_ROCKET_ALL_FIELDS = (
    "T", "P", "volume", "density", "M", "MW", "enthalpy", "energy", "entropy",
    "gibbs_energy", "gamma_s", "cp_fr", "cp_eq", "cp", "cv_fr", "cv_eq", "cv",
    "Mach", "sonic_velocity", "ae_at", "c_star", "coefficient_of_thrust",
    "Isp", "Isp_vacuum", "viscosity", "conductivity_fr", "conductivity_eq",
    "Pr_fr", "Pr_eq", "nj", "ln_nj", "n",
)

# ShockSolution.nj/ln_nj/n are excluded: the compiled class's own properties
# raise AttributeError (ShockSolution._get_weights calls a _get_size method
# that only RocketSolution defines), so cea.matlab reconstructs them from
# mole_fractions instead (see _reconstruct_shock_amounts) and there is no
# working compiled-side value to compare against.
_SHOCK_ALL_FIELDS = (
    "T", "P", "velocity", "Mach", "sonic_velocity", "rho12", "rho52", "P21",
    "P52", "T21", "T52", "M21", "M52", "v2", "u5_p_v2", "volume", "density",
    "M", "MW", "enthalpy", "energy", "entropy", "gibbs_energy", "gamma_s",
    "cp_fr", "cp_eq", "cp", "cv_fr", "cv_eq", "cv", "viscosity",
    "conductivity_fr", "conductivity_eq", "Pr_fr", "Pr_eq",
)

_DETONATION_ALL_FIELDS = (
    "P1", "T1", "H1", "M1", "gamma1", "sonic_velocity1", "P", "T", "density",
    "enthalpy", "energy", "gibbs_energy", "entropy", "Mach", "velocity",
    "sonic_velocity", "gamma_s", "P_P1", "T_T1", "M_M1", "rho_rho1", "cp_fr",
    "cv_fr", "cp_eq", "cv_eq", "M", "MW", "viscosity", "conductivity_fr",
    "conductivity_eq", "Pr_fr", "Pr_eq", "nj", "ln_nj", "n",
)


def _assert_fraction_dicts_match(wrapper_fractions, compiled_fractions):
    assert wrapper_fractions.keys() == compiled_fractions.keys()
    for name, compiled_value in compiled_fractions.items():
        np.testing.assert_allclose(
            wrapper_fractions[name],
            compiled_value,
            rtol=0.0,
            atol=0.0,
            err_msg=f"fraction {name!r} did not match the compiled solution",
        )


def _assert_all_fields_match(wrapper_soln, compiled_soln, field_names):
    for name in field_names:
        np.testing.assert_allclose(
            getattr(wrapper_soln, name),
            getattr(compiled_soln, name),
            rtol=0.0,
            atol=0.0,
            err_msg=f"field {name!r} did not match the compiled solution",
        )


def _wrapper_case_inputs():
    reactants = ["H2", "O2"]
    fuel_amounts = np.array([2.0, 0.0], dtype=np.float64)
    oxid_amounts = np.array([0.0, 1.0], dtype=np.float64)
    return reactants, fuel_amounts, oxid_amounts


def _rocket_case_inputs():
    reactants = ["H2(L)", "O2(L)"]
    fuel_amounts = np.array([1.0, 0.0], dtype=np.float64)
    oxid_amounts = np.array([0.0, 1.0], dtype=np.float64)
    return reactants, fuel_amounts, oxid_amounts


def _shock_case_inputs():
    reactants = ["H2", "O2", "Ar"]
    fuel_amounts = np.array([0.05, 0.0, 0.0], dtype=np.float64)
    oxid_amounts = np.array([0.0, 0.05, 0.9], dtype=np.float64)
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
    _assert_all_fields_match(wrapper_soln, compiled_soln, _EQ_ALL_FIELDS)
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


def test_matlab_rocket_solve_matches_compiled_solver(cea_module):
    import cea.matlab

    reactants, fuel_amounts, oxid_amounts = _rocket_case_inputs()
    t_reac = np.array([20.27, 90.17], dtype=np.float64)
    wrapper_soln = cea.matlab.rocket_solve(
        reactants,
        pc=53.3172,
        pi_p=np.array([10.0, 100.0], dtype=np.float64),
        subar=np.array([1.58], dtype=np.float64),
        supar=np.array([25.0], dtype=np.float64),
        T_reac=t_reac,
        fuel_amounts=fuel_amounts,
        oxid_amounts=oxid_amounts,
        of_ratio=5.55157,
        iac=True,
    )

    reactants_mix = cea_module.Mixture(reactants)
    products_mix = cea_module.Mixture(reactants, products_from_reactants=True)
    compiled_weights = reactants_mix.of_ratio_to_weights(
        oxid_amounts,
        fuel_amounts,
        5.55157,
    )
    compiled_hc = reactants_mix.calc_property(
        cea_module.ENTHALPY,
        compiled_weights,
        t_reac,
    ) / cea_module.R
    compiled_solver = cea_module.RocketSolver(products_mix, reactants=reactants_mix)
    compiled_soln = cea_module.RocketSolution(compiled_solver)
    compiled_solver.solve(
        compiled_soln,
        compiled_weights,
        53.3172,
        pi_p=np.array([10.0, 100.0], dtype=np.float64),
        subar=np.array([1.58], dtype=np.float64),
        supar=np.array([25.0], dtype=np.float64),
        hc=compiled_hc,
        iac=True,
    )

    assert not isinstance(wrapper_soln, cea_module.RocketSolution)
    assert not hasattr(wrapper_soln, "solver")
    assert wrapper_soln.converged == compiled_soln.converged
    assert wrapper_soln.last_error == compiled_soln.last_error
    assert wrapper_soln.num_pts == compiled_soln.num_pts
    _assert_all_fields_match(wrapper_soln, compiled_soln, _ROCKET_ALL_FIELDS)
    _assert_fraction_dicts_match(wrapper_soln.mass_fractions, compiled_soln.mass_fractions)
    _assert_fraction_dicts_match(wrapper_soln.mole_fractions, compiled_soln.mole_fractions)


def test_matlab_shock_solve_matches_compiled_solver(cea_module):
    import cea.matlab

    reactants, fuel_amounts, oxid_amounts = _shock_case_inputs()
    wrapper_soln = cea.matlab.shock_solve(
        reactants,
        T0=300.0,
        p0=cea_module.units.mmhg_to_bar(10.0),
        u1=1100.0,
        fuel_amounts=fuel_amounts,
        oxid_amounts=oxid_amounts,
        moles=True,
        reflected=True,
    )

    reactants_mix = cea_module.Mixture(reactants)
    products_mix = cea_module.Mixture(reactants, products_from_reactants=True)
    compiled_weights = reactants_mix.moles_to_weights(fuel_amounts + oxid_amounts)
    compiled_solver = cea_module.ShockSolver(products_mix, reactants=reactants_mix)
    compiled_soln = cea_module.ShockSolution(compiled_solver, reflected=True)
    compiled_solver.solve(
        compiled_soln,
        compiled_weights,
        300.0,
        cea_module.units.mmhg_to_bar(10.0),
        u1=1100.0,
        reflected=True,
    )

    assert not isinstance(wrapper_soln, cea_module.ShockSolution)
    assert not hasattr(wrapper_soln, "solver")
    assert wrapper_soln.converged == compiled_soln.converged
    assert wrapper_soln.last_error == compiled_soln.last_error
    _assert_all_fields_match(wrapper_soln, compiled_soln, _SHOCK_ALL_FIELDS)
    _assert_fraction_dicts_match(wrapper_soln.mass_fractions, compiled_soln.mass_fractions)
    _assert_fraction_dicts_match(wrapper_soln.mole_fractions, compiled_soln.mole_fractions)


def test_matlab_detonation_solve_matches_compiled_solver(cea_module):
    import cea.matlab

    reactants, fuel_amounts, oxid_amounts = _wrapper_case_inputs()
    wrapper_soln = cea.matlab.detonation_solve(
        reactants,
        T1=298.15,
        p1=1.0,
        fuel_amounts=fuel_amounts,
        oxid_amounts=oxid_amounts,
        moles=False,
        r_eq=1.0,
    )

    reactants_mix = cea_module.Mixture(reactants)
    products_mix = cea_module.Mixture(reactants, products_from_reactants=True)
    of_ratio = reactants_mix.chem_eq_ratio_to_of_ratio(oxid_amounts, fuel_amounts, 1.0)
    compiled_weights = reactants_mix.of_ratio_to_weights(oxid_amounts, fuel_amounts, of_ratio)
    compiled_solver = cea_module.DetonationSolver(products_mix, reactants=reactants_mix)
    compiled_soln = cea_module.DetonationSolution(compiled_solver)
    compiled_solver.solve(compiled_soln, compiled_weights, 298.15, 1.0)

    assert not isinstance(wrapper_soln, cea_module.DetonationSolution)
    assert not hasattr(wrapper_soln, "solver")
    assert wrapper_soln.converged == compiled_soln.converged
    assert wrapper_soln.last_error == compiled_soln.last_error
    _assert_all_fields_match(wrapper_soln, compiled_soln, _DETONATION_ALL_FIELDS)
    assert wrapper_soln.mass_fractions == compiled_soln.mass_fractions
    assert wrapper_soln.mole_fractions == compiled_soln.mole_fractions
