import subprocess
import sys
import textwrap

import pytest


# ---------------------------------------------------------------------------
# EqSolver smoke tests
# ---------------------------------------------------------------------------

@pytest.mark.smoke
def test_smooth_truncation_eqsolver_smoke():
    code = textwrap.dedent(
        """
        import numpy as np
        import cea

        reac = cea.Mixture(["H2", "O2"])
        prod = cea.Mixture(["H2", "O2"], products_from_reactants=True)
        weights = reac.moles_to_weights(np.array([2.0, 1.0], dtype=np.float64))

        solver = cea.EqSolver(prod, reactants=reac, smooth_truncation=True)
        solution = cea.EqSolution(solver)

        h0 = reac.calc_property(cea.ENTHALPY, weights, 298.15) / cea.R
        p0 = 1.0  # bar

        solver.solve(solution, cea.HP, h0, p0, weights)
        assert solution.converged
        """
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr + result.stdout


@pytest.mark.smoke
def test_smooth_truncation_eqsolver_custom_width_smoke():
    code = textwrap.dedent(
        """
        import numpy as np
        import cea

        reac = cea.Mixture(["H2", "O2"])
        prod = cea.Mixture(["H2", "O2"], products_from_reactants=True)
        weights = reac.moles_to_weights(np.array([2.0, 1.0], dtype=np.float64))

        solver = cea.EqSolver(prod, reactants=reac, smooth_truncation=True, truncation_width=0.5)
        solution = cea.EqSolution(solver)

        h0 = reac.calc_property(cea.ENTHALPY, weights, 298.15) / cea.R
        p0 = 1.0  # bar

        solver.solve(solution, cea.HP, h0, p0, weights)
        assert solution.converged
        """
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr + result.stdout


# ---------------------------------------------------------------------------
# RocketSolver / ShockSolver / DetonationSolver — constructor acceptance
# ---------------------------------------------------------------------------

@pytest.mark.smoke
def test_smooth_truncation_rocket_solver_smoke():
    code = textwrap.dedent(
        """
        import numpy as np
        import cea

        reactants = ["H2", "O2"]
        reac = cea.Mixture(reactants)
        prod = cea.Mixture(reactants, products_from_reactants=True)
        weights = reac.moles_to_weights(np.array([2.0, 1.0], dtype=np.float64))

        solver = cea.RocketSolver(prod, reactants=reac, smooth_truncation=True)
        soln = cea.RocketSolution(solver)

        solver.solve(soln, weights, pc=10.0, pi_p=2.0, tc=3000.0)
        assert soln.num_pts > 0
        """
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr + result.stdout


@pytest.mark.smoke
def test_smooth_truncation_shock_solver_smoke():
    code = textwrap.dedent(
        """
        import numpy as np
        import cea

        reactants = ["H2", "O2", "Ar"]
        reac = cea.Mixture(reactants)
        prod = cea.Mixture(reactants, products_from_reactants=True)
        moles = np.array([0.05, 0.05, 0.9])
        weights = reac.moles_to_weights(moles)

        solver = cea.ShockSolver(prod, reactants=reac, smooth_truncation=True)
        soln = cea.ShockSolution(solver, reflected=False)

        p0 = cea.units.mmhg_to_bar(10.0)
        solver.solve(soln, weights, T0=300.0, p0=p0, u1=1200.0)
        """
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr + result.stdout


@pytest.mark.smoke
def test_smooth_truncation_detonation_solver_smoke():
    code = textwrap.dedent(
        """
        import numpy as np
        import cea

        reactants = ["H2", "O2"]
        reac = cea.Mixture(reactants)
        prod = cea.Mixture(reactants, products_from_reactants=True)
        weights = reac.moles_to_weights(np.array([2.0, 1.0], dtype=np.float64))

        solver = cea.DetonationSolver(prod, reactants=reac, smooth_truncation=True)
        soln = cea.DetonationSolution(solver)

        solver.solve(soln, weights, T1=298.15, p1=1.0)
        """
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr + result.stdout


# ---------------------------------------------------------------------------
# Validation: truncation_width=0.0 with smooth_truncation=True must abort
# ---------------------------------------------------------------------------

@pytest.mark.smoke
def test_smooth_truncation_zero_width_aborts():
    """truncation_width=0.0 is explicitly invalid and should abort the process."""
    code = textwrap.dedent(
        """
        import cea

        prod = cea.Mixture(["H2", "O2"], products_from_reactants=True)
        solver = cea.EqSolver(prod, smooth_truncation=True, truncation_width=0.0)
        """
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0, (
        "Expected non-zero exit when truncation_width=0.0 with smooth_truncation=True, "
        "but process exited successfully."
    )
