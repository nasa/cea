import subprocess
import sys
import textwrap

import pytest


@pytest.mark.smoke
def test_rocket_abort_raises_runtimeerror_and_python_continues():
    code = textwrap.dedent(
        """
        import numpy as np
        import cea

        reac_names = ["H2(L)", "O2(L)"]
        reactant_temperatures = np.array([20.27, 90.17], dtype=np.float64)
        fuel_weights = np.array([1.0, 0.0], dtype=np.float64)
        oxidant_weights = np.array([0.0, 1.0], dtype=np.float64)

        reac = cea.Mixture(reac_names)
        prod = cea.Mixture(reac_names, products_from_reactants=True)
        solver = cea.RocketSolver(prod, reactants=reac)
        soln = cea.RocketSolution(solver)

        weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, 5.55157)
        hc = reac.calc_property(cea.ENTHALPY, weights, reactant_temperatures) / cea.R

        try:
            solver.solve(soln, weights, pc=53.3172, supar=[1.0], hc=hc, iac=True)
        except RuntimeError as exc:
            message = str(exc)
            assert "CEA_FORTRAN_ABORT" in message
            assert "Supersonic area ratio must be greater than 1.0" in message
            print("CAUGHT_RUNTIME_ERROR")
        else:
            raise AssertionError("Expected RocketSolver.solve to raise RuntimeError")

        print("CONTINUED_AFTER_ERROR")
        """
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        check=False,
        timeout=20,
    )
    combined = result.stdout + result.stderr

    assert result.returncode == 0, combined
    assert "CAUGHT_RUNTIME_ERROR" in result.stdout, combined
    assert "CONTINUED_AFTER_ERROR" in result.stdout, combined


@pytest.mark.smoke
def test_detonation_frozen_raises_runtimeerror_and_python_continues():
    code = textwrap.dedent(
        """
        import numpy as np
        import cea

        reac_names = ["H2", "O2"]
        reac = cea.Mixture(reac_names)
        prod = cea.Mixture(reac_names, products_from_reactants=True)
        solver = cea.DetonationSolver(prod, reactants=reac)
        soln = cea.DetonationSolution(solver)

        weights = np.array([0.0, 1.0], dtype=np.float64)

        try:
            solver.solve(soln, weights, 298.15, 1.0, frozen=True)
        except RuntimeError as exc:
            message = str(exc)
            assert "CEA_FORTRAN_ABORT" in message
            assert "frozen composition not supported yet" in message
            print("CAUGHT_RUNTIME_ERROR")
        else:
            raise AssertionError("Expected DetonationSolver.solve to raise RuntimeError")

        print("CONTINUED_AFTER_ERROR")
        """
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        check=False,
        timeout=20,
    )
    combined = result.stdout + result.stderr

    assert result.returncode == 0, combined
    assert "CAUGHT_RUNTIME_ERROR" in result.stdout, combined
    assert "CONTINUED_AFTER_ERROR" in result.stdout, combined
