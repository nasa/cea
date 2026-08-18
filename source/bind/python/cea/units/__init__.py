"""
Unit conversion helpers for the Python bindings.
"""

from __future__ import annotations

__all__ = [
    "celsius_to_kelvin",
    "kelvin_to_celsius",
    "fahrenheit_to_kelvin",
    "kelvin_to_fahrenheit",
    "rankine_to_kelvin",
    "kelvin_to_rankine",
    "atm_to_bar",
    "bar_to_atm",
    "psi_to_bar",
    "bar_to_psi",
    "mmhg_to_bar",
    "bar_to_mmhg",
    "g_per_cm3_to_kg_per_m3",
    "kg_per_m3_to_g_per_cm3",
    "cm3_per_g_to_m3_per_kg",
    "m3_per_kg_to_cm3_per_g",
    "m_per_s_to_ft_per_s",
    "ft_per_s_to_m_per_s",
    "joule_to_cal",
    "cal_to_joule",
    "kj_per_kg_to_cal_per_g",
    "cal_per_g_to_kj_per_kg",
]


# Temperature -----------------------------------------------------------------

def celsius_to_kelvin(value):
    """Convert °C to K."""

    return value + 273.15


def kelvin_to_celsius(value):
    """Convert K to °C."""

    return value - 273.15


def fahrenheit_to_kelvin(value):
    """Convert °F to K."""

    return (value - 32.0) / 1.8 + 273.15


def kelvin_to_fahrenheit(value):
    """Convert K to °F."""

    return (value - 273.15) * 1.8 + 32.0


def rankine_to_kelvin(value):
    """Convert °R to K."""

    return value / 1.8


def kelvin_to_rankine(value):
    """Convert K to °R."""

    return value * 1.8


# Pressure --------------------------------------------------------------------

_ATM_TO_BAR = 1.01325
# 1 atm = 14.695948775513 psi, derived from the exact SI base-unit definitions
# of the lbf, inch, and atm. Matches the Fortran core's psi_to_bar
# (source/units.f90) to ~9 significant figures.
_PSI_TO_BAR = _ATM_TO_BAR/14.695948775513
_MMHG_TO_BAR = _ATM_TO_BAR/760.0


def atm_to_bar(value):
    """Convert atmospheres to bar."""

    return value * _ATM_TO_BAR


def bar_to_atm(value):
    """Convert bar to atmospheres."""

    return value / _ATM_TO_BAR


def psi_to_bar(value):
    """Convert psi to bar."""

    return value * _PSI_TO_BAR


def bar_to_psi(value):
    """Convert bar to psi."""

    return value / _PSI_TO_BAR


def mmhg_to_bar(value):
    """Convert mmHg to bar."""

    return value * _MMHG_TO_BAR


def bar_to_mmhg(value):
    """Convert bar to mmHg."""

    return value / _MMHG_TO_BAR


# Density & specific volume ----------------------------------------------------

def g_per_cm3_to_kg_per_m3(value):
    """Convert g/cm³ to kg/m³."""

    return value * 1000.0


def kg_per_m3_to_g_per_cm3(value):
    """Convert kg/m³ to g/cm³."""

    return value / 1000.0


def cm3_per_g_to_m3_per_kg(value):
    """Convert cm³/g to m³/kg."""

    return value * 1.0e-3


def m3_per_kg_to_cm3_per_g(value):
    """Convert m³/kg to cm³/g."""

    return value * 1.0e3


# Velocity --------------------------------------------------------------------
_M_TO_FT = 3.28084


def m_per_s_to_ft_per_s(value):
    """Convert m/s to ft/s."""

    return value * _M_TO_FT


def ft_per_s_to_m_per_s(value):
    """Convert ft/s to m/s."""

    return value / _M_TO_FT


# Energy ----------------------------------------------------------------------
_JOULE_PER_CAL = 4.184


def joule_to_cal(value):
    """Convert joules to calories."""

    return value / _JOULE_PER_CAL


def cal_to_joule(value):
    """Convert calories to joules."""

    return value * _JOULE_PER_CAL


def kj_per_kg_to_cal_per_g(value):
    """Convert kJ/kg to cal/g."""

    return value / _JOULE_PER_CAL


def cal_per_g_to_kj_per_kg(value):
    """Convert cal/g to kJ/kg."""

    return value * _JOULE_PER_CAL
