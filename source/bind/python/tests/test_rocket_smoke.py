import subprocess
import sys
import textwrap

import pytest


@pytest.mark.smoke
def test_rocket_solver_smoke():
    code = textwrap.dedent(
        """
        import numpy as np
        import cea

        reactants = ["H2", "O2"]
        reactants_mix = cea.Mixture(reactants)
        products_mix = cea.Mixture(reactants, products_from_reactants=True)
        weights = reactants_mix.moles_to_weights(np.array([2.0, 1.0], dtype=np.float64))
        solver = cea.RocketSolver(products_mix, reactants=reactants_mix)
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
def test_rocket_solver_smoke_without_exit_conditions():
    code = textwrap.dedent(
        """
        import numpy as np
        import cea

        reactants = ["H2", "O2"]
        reactants_mix = cea.Mixture(reactants)
        products_mix = cea.Mixture(reactants, products_from_reactants=True)
        weights = reactants_mix.moles_to_weights(np.array([2.0, 1.0], dtype=np.float64))
        solver = cea.RocketSolver(products_mix, reactants=reactants_mix)
        soln = cea.RocketSolution(solver)

        solver.solve(soln, weights, pc=10.0, tc=3000.0)
        """
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr + result.stdout
