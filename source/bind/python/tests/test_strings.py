import numpy as np
import pytest

@pytest.mark.smoke
def test_mixture_accepts_python_strings(cea_module) -> None:
    mix = cea_module.Mixture(["H2", "O2", "Ar"])
    weights = mix.moles_to_weights(np.array([2.0, 1.0, 0.5], dtype=np.float64))
    density = mix.calc_property(cea_module.DENSITY, weights, temperature=3000.0, pressure=1.0)
    assert np.isfinite(density)
