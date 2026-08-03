import numpy as np
import pytest


def test_species_names_round_trip(cea_module):
    mix = cea_module.Mixture(["H2", "O2", "Ar"])
    names = mix.species_names
    assert len(names) == 3
    for name in ("H2", "O2", "Ar"):
        assert name in names


def test_moles_weights_round_trip(cea_module):
    mix = cea_module.Mixture(["H2", "O2", "Ar"])
    moles = np.array([2.0, 1.0, 0.5], dtype=np.float64)
    weights = mix.moles_to_weights(moles)
    roundtrip = mix.weights_to_moles(weights)
    moles_norm = moles / np.sum(moles)
    roundtrip_norm = roundtrip / np.sum(roundtrip)
    assert np.allclose(roundtrip_norm, moles_norm, rtol=1e-12, atol=1e-12)


def test_calc_property_requires_pressure(cea_module):
    mix = cea_module.Mixture(["H2", "O2"])
    weights = np.array([0.4, 0.6], dtype=np.float64)
    with pytest.raises(ValueError):
        mix.calc_property(cea_module.ENTROPY, weights, temperature=3000.0)


def test_calc_property_rejects_unknown_type(cea_module):
    mix = cea_module.Mixture(["H2", "O2"])
    weights = np.array([0.4, 0.6], dtype=np.float64)
    with pytest.raises(ValueError):
        mix.calc_property(999999, weights, temperature=3000.0)


@pytest.mark.parametrize(
    ("species", "weights", "temperatures", "expected_density"),
    [
        (["O2(L)"], np.array([1.0]), np.array([90.0]), 0.2907792655383313),
        (
            ["CH4(L)", "O2(L)"],
            np.array([0.2, 0.8]),
            np.array([111.0, 90.0]),
            0.22505975471169806,
        ),
    ],
)
def test_calc_liquid_density_is_finite(
    cea_module, species, weights, temperatures, expected_density
):
    mix = cea_module.Mixture(species)

    density = mix.calc_property(
        cea_module.DENSITY, weights, temperature=temperatures, pressure=68.0
    )
    volume = mix.calc_property(
        cea_module.VOLUME, weights, temperature=temperatures, pressure=68.0
    )

    assert np.isfinite(density)
    assert density > 0.0
    assert density == pytest.approx(expected_density)
    assert density * volume == pytest.approx(1.0)


def test_multitemperature_density_matches_single_temperature(cea_module):
    mix = cea_module.Mixture(["O2(L)"])
    weights = np.array([1.0])

    density_single = mix.calc_property(
        cea_module.DENSITY, weights, temperature=90.0, pressure=68.0
    )
    density_multi = mix.calc_property(
        cea_module.DENSITY, weights, temperature=np.array([90.0]), pressure=68.0
    )

    assert density_multi == pytest.approx(density_single)
