import pytest


def test_reactant_reference_temperature_range_before_mixture(cea_module):
    reactant = cea_module.Reactant("(CH2)x(cr)")

    temperature_range = reactant.get_valid_temperature_range()

    assert temperature_range == pytest.approx((288.15, 308.15))
    assert all(isinstance(value, float) for value in temperature_range)


def test_reactant_fitted_temperature_range(cea_module):
    reactant = cea_module.Reactant("NH4CLO4(I)")

    assert reactant.get_valid_temperature_range() == pytest.approx((100.0, 513.15))


def test_reactant_temperature_range_rejects_unknown_name(cea_module):
    reactant = cea_module.Reactant("not-a-reactant")

    with pytest.raises(ValueError, match="not found in thermo database"):
        reactant.get_valid_temperature_range()
