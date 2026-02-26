import subprocess
import sys
import textwrap

import numpy as np
import pytest

import cea


@pytest.mark.smoke
def test_issue42_detonation_no_inner_eq_failure():
    code = textwrap.dedent(
        """
        import numpy as np
        import cea

        reac_names = ["O2", "C2H4"]
        reac = cea.Mixture(reac_names)
        prod = cea.Mixture(reac_names, products_from_reactants=True)

        wo = np.array((1.0, 0.0))
        wf = np.array((0.0, 1.0))
        of_ratio = reac.weight_eq_ratio_to_of_ratio(wo, wf, 0.4)
        weights = reac.of_ratio_to_weights(wo, wf, of_ratio)

        solver = cea.DetonationSolver(prod, reactants=reac)
        soln = cea.DetonationSolution(solver)
        solver.solve(soln, weights, T1=298.15, p1=4.0)
        assert soln.converged
        """
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        check=False,
    )
    combined = result.stdout + result.stderr
    assert result.returncode == 0, combined
    assert "DETON: Equilibrium solver did not converge" not in combined, combined


def test_issue42_reused_eqsolution_matches_fresh():
    reac_names = ["O2", "C2H4"]
    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True)

    p1 = 60.0
    h1 = 129.92208593485029
    p2 = 95.242559745142628
    t2 = 3585.9524742656135

    wo = np.array((1.0, 0.0))
    wf = np.array((0.0, 1.0))
    of_ratio = reac.weight_eq_ratio_to_of_ratio(wo, wf, 0.4)
    weights = reac.of_ratio_to_weights(wo, wf, of_ratio)

    solver = cea.EqSolver(prod, reactants=reac)

    reuse_soln = cea.EqSolution(solver)
    solver.solve(reuse_soln, cea.HP, h1, p1, weights)
    assert reuse_soln.converged
    solver.solve(reuse_soln, cea.TP, t2, p2, weights)
    assert reuse_soln.converged

    fresh_soln = cea.EqSolution(solver)
    solver.solve(fresh_soln, cea.TP, t2, p2, weights)
    assert fresh_soln.converged

    assert reuse_soln.T == pytest.approx(fresh_soln.T, rel=1.0e-6)
    assert reuse_soln.n == pytest.approx(fresh_soln.n, rel=1.0e-6)
    assert np.allclose(reuse_soln.nj, fresh_soln.nj, rtol=1.0e-6, atol=1.0e-12)
