"""
Accuracy regression tests for RP-1311 sample problems.

Each test mirrors the solver setup from the corresponding sample in
source/bind/python/cea/samples/rp1311/ and asserts key scalar outputs
against reference values obtained by running those samples.

Reference values were captured with rel=1e-5 tolerance in mind.
"""
import numpy as np
import pytest

# Omit list shared by Examples 3 and 4
_OMIT_34 = [
    "CCN", "CNC", "C2N2", "C2O",
    "C3H4,allene", "C3H4,propyne", "C3H4,cyclo-", "C3",
    "C3H5,allyl", "C3H6,propylene", "C3H6,cyclo-", "C3H3,propargyl",
    "C3H6O", "C3H7,n-propyl", "C3H7,i-propyl", "Jet-A(g)",
    "C3O2", "C4", "C4H2", "C3H8O,2propanol",
    "C4H4,1,3-cyclo-", "C4H6,butadiene", "C4H6,2-butyne", "C3H8O,1propanol",
    "C4H8,tr2-butene", "C4H8,isobutene", "C4H8,cyclo-", "C4H6,cyclo-",
    "(CH3COOH)2", "C4H9,n-butyl", "C4H9,i-butyl", "C4H8,1-butene",
    "C4H9,s-butyl", "C4H9,t-butyl", "C4H10,isobutane", "C4H8,cis2-buten",
    "C4H10,n-butane", "C4N2", "C5", "C3H8",
    "C5H6,1,3cyclo-", "C5H8,cyclo-", "C5H10,1-pentene", "C10H21,n-decyl",
    "C5H10,cyclo-", "C5H11,pentyl", "C5H11,t-pentyl", "C12H10,biphenyl",
    "C5H12,n-pentane", "C5H12,i-pentane", "CH3C(CH3)2CH3", "C12H9,o-bipheny",
    "C6H6", "C6H5OH,phenol", "C6H10,cyclo-", "C6H2",
    "C6H12,1-hexene", "C6H12,cyclo-", "C6H13,n-hexyl", "C6H5,phenyl",
    "C7H7,benzyl", "C7H8", "C7H8O,cresol-mx", "C6H5O,phenoxy",
    "C7H14,1-heptene", "C7H15,n-heptyl", "C7H16,n-heptane", "C10H8,azulene",
    "C8H8,styrene", "C8H10,ethylbenz", "C8H16,1-octene", "C10H8,napthlene",
    "C8H17,n-octyl", "C8H18,isooctane", "C8H18,n-octane", "C9H19,n-nonyl",
    "Jet-A(L)", "C6H6(L)", "H2O(s)", "H2O(L)",
]


@pytest.mark.rp1311
def test_example1_tp_equilibrium(cea_module):
    """RP-1311 Example 1: TP equilibrium, H2/Air at multiple phi, P, T."""
    cea = cea_module
    reac_names = ["H2", "Air"]
    prod_names = ["Ar", "C", "CO", "CO2", "H", "H2", "H2O", "HNO", "HO2", "HNO2",
                  "HNO3", "N", "NH", "NO", "N2", "N2O3", "O", "O2", "OH", "O3"]
    pressures = cea.units.atm_to_bar(np.array([1.0, 0.1, 0.01]))
    fuel_moles = np.array([1.0, 0.0])
    oxidant_moles = np.array([0.0, 1.0])

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(prod_names)
    solver = cea.EqSolver(prod, reactants=reac)
    soln = cea.EqSolution(solver)
    fw = reac.moles_to_weights(fuel_moles)
    ow = reac.moles_to_weights(oxidant_moles)

    # phi=1.0, P=1 atm, T=3000 K
    of1 = reac.chem_eq_ratio_to_of_ratio(ow, fw, 1.0)
    w1 = reac.of_ratio_to_weights(ow, fw, of1)
    solver.solve(soln, cea.TP, 3000.0, pressures[0], w1)
    assert soln.converged
    assert soln.density * 1e-3 == pytest.approx(9.178126199101762e-05, rel=1e-5)
    assert cea.units.joule_to_cal(soln.entropy) == pytest.approx(2.8794139352625625, rel=1e-5)
    assert soln.gamma_s == pytest.approx(1.1311926054937052, rel=1e-5)
    assert cea.units.joule_to_cal(soln.cp_eq) == pytest.approx(1.6816419961717988, rel=1e-5)

    # phi=1.0, P=1 atm, T=2000 K
    solver.solve(soln, cea.TP, 2000.0, pressures[0], w1)
    assert soln.converged
    assert soln.density * 1e-3 == pytest.approx(0.00014989459487915605, rel=1e-5)
    assert soln.gamma_s == pytest.approx(1.225763874534332, rel=1e-5)

    # phi=1.5, P=0.01 atm, T=2000 K  (last condition)
    of15 = reac.chem_eq_ratio_to_of_ratio(ow, fw, 1.5)
    w15 = reac.of_ratio_to_weights(ow, fw, of15)
    solver.solve(soln, cea.TP, 2000.0, pressures[2], w15)
    assert soln.converged
    assert soln.density * 1e-3 == pytest.approx(1.2929027883808208e-06, rel=1e-5)
    assert soln.gamma_s == pytest.approx(1.2051203751977477, rel=1e-5)


@pytest.mark.rp1311
def test_example2_tv_equilibrium_transport(cea_module):
    """RP-1311 Example 2: TV equilibrium + transport, H2/Air (phi=1.0)."""
    cea = cea_module
    reac_names = ["H2", "Air"]
    prod_names = ["Ar", "C", "CO", "CO2", "H", "H2", "H2O", "HNO", "HO2", "HNO2",
                  "HNO3", "N", "NH", "NO", "N2", "N2O3", "O", "O2", "OH", "O3"]
    densities = 1.0e3 * np.array([9.1864e-5, 8.0877e-6, 6.6054e-7])
    fuel_moles = np.array([1.0, 0.0])
    oxidant_moles = np.array([0.0, 1.0])

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(prod_names)
    solver = cea.EqSolver(prod, reactants=reac, transport=True)
    soln = cea.EqSolution(solver)
    fw = reac.moles_to_weights(fuel_moles)
    ow = reac.moles_to_weights(oxidant_moles)
    of_phi1 = reac.weight_eq_ratio_to_of_ratio(ow, fw, 1.0)
    w = reac.of_ratio_to_weights(ow, fw, of_phi1)

    # density[0], T=3000 K
    solver.solve(soln, cea.TV, 3000.0, 1.0 / densities[0], w)
    assert soln.converged
    assert cea.units.bar_to_atm(soln.P) == pytest.approx(1.000871440004418, rel=1e-5)
    assert cea.units.joule_to_cal(soln.entropy) == pytest.approx(2.8792839379777275, rel=1e-5)
    assert soln.viscosity == pytest.approx(0.9357605004867593, rel=1e-4)
    assert soln.Pr_eq == pytest.approx(0.3558398119457039, rel=1e-4)
    assert soln.gamma_s == pytest.approx(1.131199663437041, rel=1e-5)

    # density[2], T=3000 K (lowest density)
    solver.solve(soln, cea.TV, 3000.0, 1.0 / densities[2], w)
    assert soln.converged
    assert cea.units.bar_to_atm(soln.P) == pytest.approx(0.009988512678686615, rel=1e-5)
    assert cea.units.joule_to_cal(soln.entropy) == pytest.approx(4.009691725408506, rel=1e-5)


@pytest.mark.rp1311
def test_example3_hp_equilibrium_omit(cea_module):
    """RP-1311 Example 3: HP equilibrium + omit list, Air/C7H8(L)/C8H18(L)."""
    cea = cea_module
    reac_names = ["Air", "C7H8(L)", "C8H18(L),n-octa"]
    ow = np.array([1.0, 0.0, 0.0])
    fw = np.array([0.0, 0.4, 0.6])
    T_reac = np.array([700.0, 298.15, 298.15])
    of_ratio = 17.0

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True, omit=_OMIT_34)
    solver = cea.EqSolver(prod, reactants=reac, trace=1e-15)
    soln = cea.EqSolution(solver)
    w = reac.of_ratio_to_weights(ow, fw, of_ratio)
    h0 = reac.calc_property(cea.ENTHALPY, w, T_reac)

    # P=100 bar (highest pressure)
    solver.solve(soln, cea.HP, h0 / cea.R, 100.0, w)
    assert soln.converged
    assert soln.T == pytest.approx(2418.6601475091556, rel=1e-5)
    assert soln.entropy == pytest.approx(8.168090909551534, rel=1e-5)
    assert soln.gamma_s == pytest.approx(1.2256993099599578, rel=1e-5)
    assert soln.mole_fractions["CO2"] == pytest.approx(0.11536898106241796, rel=1e-4)

    # P=1 bar (lowest pressure)
    solver.solve(soln, cea.HP, h0 / cea.R, 1.0, w)
    assert soln.converged
    assert soln.T == pytest.approx(2338.839916883393, rel=1e-5)
    assert soln.gamma_s == pytest.approx(1.1799908487697623, rel=1e-5)


@pytest.mark.rp1311
def test_example4_uv_equilibrium_omit(cea_module):
    """RP-1311 Example 4: UV equilibrium + omit list, Air/C7H8(L)/C8H18(L)."""
    cea = cea_module
    reac_names = ["Air", "C7H8(L)", "C8H18(L),n-octa"]
    ow = np.array([1.0, 0.0, 0.0])
    fw = np.array([0.0, 0.4, 0.6])
    of_ratio = 17.0
    u_R = -45.1343
    density = 14.428  # kg/m^3

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True, omit=_OMIT_34)
    solver = cea.EqSolver(prod, reactants=reac, trace=1e-15)
    soln = cea.EqSolution(solver)
    w = reac.of_ratio_to_weights(ow, fw, of_ratio)

    solver.solve(soln, cea.UV, u_R, 1.0 / density, w)
    assert soln.converged
    assert soln.T == pytest.approx(2418.5382545643524, rel=1e-5)
    assert soln.P == pytest.approx(99.9736940627695, rel=1e-5)
    assert soln.entropy == pytest.approx(8.168086603243845, rel=1e-5)
    assert soln.gamma_s == pytest.approx(1.225707104935455, rel=1e-5)


@pytest.mark.rp1311
def test_example6_detonation_transport(cea_module):
    """RP-1311 Example 6: CJ detonation + transport, H2/O2 (phi=1.0)."""
    cea = cea_module
    reac_names = ["H2", "O2"]
    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True)
    solver = cea.DetonationSolver(prod, reactants=reac, transport=True)
    soln = cea.DetonationSolution(solver)
    ow = np.array([0.0, 1.0])
    fw = np.array([1.0, 0.0])
    of_ratio = reac.chem_eq_ratio_to_of_ratio(ow, fw, 1.0)
    w = reac.of_ratio_to_weights(ow, fw, of_ratio)

    # T1=298.15 K, p1=1 bar
    solver.solve(soln, w, 298.15, 1.0)
    assert soln.converged
    assert soln.T == pytest.approx(3674.2808344247833, rel=1e-5)
    assert soln.P_P1 == pytest.approx(18.768747358331048, rel=1e-5)
    assert soln.T_T1 == pytest.approx(12.323598556572396, rel=1e-5)
    assert soln.velocity == pytest.approx(2835.5301144596524, rel=1e-5)
    assert soln.Mach == pytest.approx(5.271794754549113, rel=1e-5)
    assert soln.viscosity == pytest.approx(1.1411610227886602, rel=1e-4)

    # T1=298.15 K, p1=20 bar
    solver.solve(soln, w, 298.15, 20.0)
    assert soln.converged
    assert soln.T == pytest.approx(4283.381543661223, rel=1e-5)
    assert soln.velocity == pytest.approx(2993.805224611246, rel=1e-5)

    # T1=500 K, p1=1 bar
    solver.solve(soln, w, 500.0, 1.0)
    assert soln.converged
    assert soln.T == pytest.approx(3600.031140211338, rel=1e-5)
    assert soln.velocity == pytest.approx(2778.311547817278, rel=1e-5)


@pytest.mark.rp1311
def test_example7_shock_incident_reflected(cea_module):
    """RP-1311 Example 7: Shock tube (incident + reflected), H2/O2/Ar."""
    cea = cea_module
    reac_names = ["H2", "O2", "Ar"]
    moles = np.array([0.05, 0.05, 0.9])
    p0 = cea.units.mmhg_to_bar(10.0)
    T0 = 300.0

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True)
    solver = cea.ShockSolver(prod, reactants=reac)
    soln = cea.ShockSolution(solver, reflected=True)
    w = reac.moles_to_weights(moles)

    # u1=1100 m/s (lowest velocity)
    solver.solve(soln, w, T0, p0, u1=1100.0, reflected=True)
    assert soln.converged
    assert soln.T[1] == pytest.approx(1528.5164120286374, rel=1e-5)
    assert soln.P[1] == pytest.approx(0.10932970678558974, rel=1e-5)
    assert soln.T[2] == pytest.approx(2108.3756721712966, rel=1e-5)
    assert soln.P21 == pytest.approx(8.20040255874062, rel=1e-5)
    assert soln.P52 == pytest.approx(2.3693284648883495, rel=1e-5)

    # u1=1400 m/s (highest velocity)
    solver.solve(soln, w, T0, p0, u1=1400.0, reflected=True)
    assert soln.converged
    assert soln.T[1] == pytest.approx(2249.060164401022, rel=1e-5)
    assert soln.T[2] == pytest.approx(3039.4461380251523, rel=1e-5)
    assert soln.P21 == pytest.approx(19.432335645509717, rel=1e-5)
    assert soln.P52 == pytest.approx(3.569192372318813, rel=1e-5)


@pytest.mark.rp1311
def test_example8_iac_rocket(cea_module):
    """RP-1311 Example 8: IAC rocket, H2(L)/O2(L), pi_p + subar + supar."""
    cea = cea_module
    reac_names = ["H2(L)", "O2(L)"]
    T_reac = np.array([20.27, 90.17])
    fw = np.array([1.0, 0.0])
    ow = np.array([0.0, 1.0])
    of_ratio = 5.55157
    pc = 53.3172
    pi_p = [10.0, 100.0, 1000.0]
    subar = [1.58]
    supar = [25.0, 50.0, 75.0]

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True)
    solver = cea.RocketSolver(prod, reactants=reac)
    soln = cea.RocketSolution(solver)
    w = reac.of_ratio_to_weights(ow, fw, of_ratio)
    hc = reac.calc_property(cea.ENTHALPY, w, T_reac) / cea.R
    solver.solve(soln, w, pc, pi_p, subar=subar, supar=supar, hc=hc, iac=True)

    assert soln.num_pts == 9
    assert soln.T[0] == pytest.approx(3383.844614137881, rel=1e-5)
    assert soln.c_star[0] == pytest.approx(2332.3358043354447, rel=1e-5)
    assert soln.Isp_vacuum[-1] == pytest.approx(4554.912803667864, rel=1e-5)
    assert soln.coefficient_of_thrust[-1] == pytest.approx(1.8861438785811016, rel=1e-5)
    assert soln.gamma_s[2] == pytest.approx(1.1724599771394038, rel=1e-5)


@pytest.mark.rp1311
def test_example9_fac_rocket(cea_module):
    """RP-1311 Example 9: FAC rocket, H2(L)/O2(L), ac_at=1.58."""
    cea = cea_module
    reac_names = ["H2(L)", "O2(L)"]
    T_reac = np.array([20.27, 90.17])
    fw = np.array([1.0, 0.0])
    ow = np.array([0.0, 1.0])
    of_ratio = 5.55157
    pc = 53.3172
    pi_p = [10.0, 100.0, 1000.0]
    supar = [25.0, 50.0, 75.0]

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True)
    solver = cea.RocketSolver(prod, reactants=reac)
    soln = cea.RocketSolution(solver)
    w = reac.of_ratio_to_weights(ow, fw, of_ratio)
    hc = reac.calc_property(cea.ENTHALPY, w, T_reac) / cea.R
    solver.solve(soln, w, pc, pi_p, supar=supar, ac_at=1.58, iac=False, hc=hc)

    assert soln.num_pts == 10
    assert soln.T[0] == pytest.approx(3383.844614137881, rel=1e-5)
    assert soln.c_star[0] == pytest.approx(2330.967589918211, rel=1e-5)
    assert soln.Isp_vacuum[-1] == pytest.approx(4554.315454546227, rel=1e-5)
    assert soln.coefficient_of_thrust[-1] == pytest.approx(1.8869095701543406, rel=1e-5)


@pytest.mark.rp1311
def test_example10_fac_rocket_mdot(cea_module):
    """RP-1311 Example 10: FAC rocket + mass flow rate, H2(L)/O2(L)."""
    cea = cea_module
    reac_names = ["H2(L)", "O2(L)"]
    T_reac = np.array([20.27, 90.17])
    fw = np.array([1.0, 0.0])
    ow = np.array([0.0, 1.0])
    of_ratio = 5.55157
    pc = 53.3172
    pi_p = [10.0, 100.0, 1000.0]
    supar = [25.0, 50.0, 75.0]

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True)
    solver = cea.RocketSolver(prod, reactants=reac)
    soln = cea.RocketSolution(solver)
    w = reac.of_ratio_to_weights(ow, fw, of_ratio)
    hc = reac.calc_property(cea.ENTHALPY, w, T_reac) / cea.R
    solver.solve(soln, w, pc, pi_p, supar=supar, mdot=1333.9, iac=False, hc=hc)

    assert soln.T[0] == pytest.approx(3383.844614137881, rel=1e-5)
    assert soln.Isp_vacuum[-1] == pytest.approx(4554.328423841097, rel=1e-5)
    assert soln.c_star[0] == pytest.approx(2330.9972184994995, rel=1e-5)


@pytest.mark.rp1311
def test_example11_iac_rocket_ions_transport(cea_module):
    """RP-1311 Example 11: IAC rocket + transport + ions, Li(cr)/F2(L)."""
    cea = cea_module
    reac_names = ["Li(cr)", "F2(L)"]
    T_reac = np.array([298.15, 85.02])
    fuel_moles = np.array([1.0, 0.0])
    oxidant_moles = np.array([0.0, 0.5556])
    pc = cea.units.psi_to_bar(1000.0)
    pi_p = [68.0457]
    subar = [10.0]
    supar = [10.0, 20.0, 100.0]

    reac = cea.Mixture(reac_names, ions=True)
    prod = cea.Mixture(reac_names, products_from_reactants=True, ions=True)
    solver = cea.RocketSolver(prod, reactants=reac, transport=True, ions=True)
    soln = cea.RocketSolution(solver)
    fw = reac.moles_to_weights(fuel_moles)
    ow = reac.moles_to_weights(oxidant_moles)
    w = fw + ow
    hc = reac.calc_property(cea.ENTHALPY, w, T_reac) / cea.R
    solver.solve(soln, w, pc, pi_p, subar=subar, supar=supar, iac=True, hc=hc)

    assert soln.T[0] == pytest.approx(5680.806086249066, rel=1e-5)
    assert soln.Isp_vacuum[-1] == pytest.approx(4504.508490838814, rel=1e-5)
    assert soln.c_star[0] == pytest.approx(2274.4011463000497, rel=1e-5)
    assert soln.viscosity[0] == pytest.approx(1.440860840873496, rel=1e-4)


@pytest.mark.rp1311
def test_example12_iac_rocket_frozen(cea_module):
    """RP-1311 Example 12: IAC rocket + frozen composition (n_frz=2), CH6N2(L)/N2O4(L)."""
    cea = cea_module
    reac_names = ["CH6N2(L)", "N2O4(L)"]
    fw = np.array([1.0, 0.0])
    ow = np.array([0.0, 1.0])
    of_ratio = 2.5
    pc = cea.units.psi_to_bar(1000)
    pi_p = 68.0457
    supar = [5.0, 10.0, 25.0, 50.0, 75.0, 100.0, 150.0, 200.0]

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(["CO", "CO2", "H", "HNO", "HNO2", "HO2",
                        "H2", "H2O", "H2O2", "N", "NO", "NO2",
                        "N2", "N2O", "O", "OH", "O2", "HCO", "NH",
                        "CH4", "NH2", "NH3", "H2O(L)", "C(gr)"])
    solver = cea.RocketSolver(prod, reactants=reac)
    soln = cea.RocketSolution(solver)
    w = reac.of_ratio_to_weights(ow, fw, of_ratio)
    hc = reac.calc_property(cea.ENTHALPY, w, 298.15) / cea.R
    solver.solve(soln, w, pc, pi_p, supar=supar, iac=True, hc=hc, n_frz=2)

    assert soln.T[0] == pytest.approx(3382.318113350732, rel=1e-5)
    assert soln.c_star[0] == pytest.approx(1707.9200870494435, rel=1e-5)
    assert soln.Isp_vacuum[-1] == pytest.approx(3287.030968566901, rel=1e-5)
    assert soln.coefficient_of_thrust[-1] == pytest.approx(1.884600866418363, rel=1e-5)


@pytest.mark.rp1311
def test_example13_iac_rocket_insert_trace(cea_module):
    """RP-1311 Example 13: IAC rocket + BeO(L) insert + trace, N2H4(L)/Be(a)/H2O2(L)."""
    cea = cea_module
    reac_names = ["N2H4(L)", "Be(a)", "H2O2(L)"]
    fw = np.array([0.8, 0.2, 0.0])
    ow = np.array([0.0, 0.0, 1.0])
    reac_T = np.array([298.15, 298.15, 298.15])
    pct_fuel = 67.0
    of_ratio = (100.0 - pct_fuel) / pct_fuel
    pc = cea.units.psi_to_bar(3000)
    pi_p = [3.0, 10.0, 30.0, 300.0]

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True)
    solver = cea.RocketSolver(prod, reactants=reac, trace=1e-10, insert=["BeO(L)"])
    soln = cea.RocketSolution(solver)
    w = reac.of_ratio_to_weights(ow, fw, of_ratio)
    hc = reac.calc_property(cea.ENTHALPY, w, reac_T) / cea.R
    solver.solve(soln, w, pc, pi_p, iac=True, hc=hc)

    assert soln.T[0] == pytest.approx(3002.53999543011, rel=1e-5)
    assert soln.Isp[-1] == pytest.approx(3575.6292392518762, rel=1e-5)
    assert soln.Isp_vacuum[-1] == pytest.approx(3770.3024620774268, rel=1e-5)
    assert soln.c_star[0] == pytest.approx(1946.6824479636061, rel=1e-5)
    assert soln.coefficient_of_thrust[-1] == pytest.approx(1.8367809516093834, rel=1e-5)


@pytest.mark.rp1311
def test_example14_tp_condensed(cea_module):
    """RP-1311 Example 14: TP equilibrium + condensed species, H2(L)/O2(L)."""
    cea = cea_module
    reac_names = ["H2(L)", "O2(L)"]
    p = cea.units.atm_to_bar(0.05)
    fuel_moles = np.array([100.0, 0.0])
    oxidant_moles = np.array([0.0, 60.0])

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True)
    solver = cea.EqSolver(prod, reactants=reac)
    soln = cea.EqSolution(solver)
    fw = reac.moles_to_weights(fuel_moles)
    ow = reac.moles_to_weights(oxidant_moles)
    w = fw + ow

    # T=1000 K (hot end, no condensation)
    solver.solve(soln, cea.TP, 1000.0, p, w)
    assert soln.converged
    assert soln.density == pytest.approx(0.011751775633966005, rel=1e-5)
    assert soln.gamma_s == pytest.approx(1.2566523234911122, rel=1e-5)
    assert soln.entropy == pytest.approx(13.535615036085526, rel=1e-5)

    # T=300 K (cold end, near condensation)
    solver.solve(soln, cea.TP, 300.0, p, w)
    assert soln.converged
    assert soln.density == pytest.approx(0.13035416916736045, rel=1e-5)
    assert soln.gamma_s == pytest.approx(1.0343906413887525, rel=1e-5)
