import numpy as np
import cea

"""
Example 5 from RP-1311
- HP equilibrium for a solid propellant blend
- Includes one custom reactant (CHOS-Binder) not present in thermo.lib
"""

pressures = np.array([34.473652, 17.236826, 8.618413, 3.447365, 0.344737])
weights = np.array([0.7206, 0.1858, 0.09, 0.002, 0.0016], dtype=np.float64)
T_reac = np.array([298.15, 298.15, 298.15, 298.15, 298.15], dtype=np.float64)

omit_names = [
    "COOH", "C2", "C2H", "CHCO,ketyl", "C2H2,vinylidene", "CH2CO,ketene", "C2H3,vinyl",
    "CH3CO,acetyl", "C2H4O,ethylen-o", "CH3CHO,ethanal", "CH3COOH", "(HCOOH)2",
    "C2H5", "C2H6", "CH3N2CH3", "CH3OCH3", "C2H5OH", "CCN", "CNC", "C2N2",
    "C2O", "C3", "C3H3,propargyl", "C3H4,allene", "C3H4,propyne", "C3H4,cyclo-",
    "C3H5,allyl", "C3H6,propylene", "C3H6,cyclo-", "C3H6O", "C3H7,n-propyl",
    "C3H7,i-propyl", "C3H8", "C3H8O,1propanol", "C3H8O,2propanol", "C3O2",
    "C4", "C4H2", "C4H4,1,3-cyclo-", "C4H6,butadiene", "C4H6,2-butyne", "C4H6,cyclo-",
    "C4H8,1-butene", "C4H8,cis2-buten", "C4H8,tr2-butene", "C4H8,isobutene", "C4H8,cyclo-",
    "(CH3COOH)2", "C4H9,n-butyl", "C4H9,i-butyl", "C4H9,s-butyl", "C4H9,t-butyl",
    "C4H10,isobutane", "C4H10,n-butane", "C4N2", "C5", "C5H6,1,3cyclo-", "C5H8,cyclo-",
    "C5H10,1-pentene", "C5H10,cyclo-", "C5H11,pentyl", "C5H11,t-pentyl", "C5H12,n-pentane",
    "C5H12,i-pentane", "CH3C(CH3)2CH3", "C6H2", "C6H5,phenyl", "C6H5O,phenoxy",
    "C6H6", "C6H5OH,phenol", "C6H10,cyclo-", "C6H12,1-hexene", "C6H12,cyclo-",
    "C6H13,n-hexyl", "C7H7,benzyl", "C7H8", "C7H8O,cresol-mx", "C7H14,1-heptene",
    "C7H15,n-heptyl", "C7H16,n-heptane", "C8H8,styrene", "C8H10,ethylbenz",
    "C8H16,1-octene", "C8H17,n-octyl", "C8H18,isooctane", "C8H18,n-octane",
    "C9H19,n-nonyl", "C10H8,naphthale", "C10H21,n-decyl", "C12H9,o-bipheny", "C12H10,biphenyl",
    "Jet-A(g)", "HNCO", "HNO", "HNO2", "HNO3", "HCCN", "HCHO,formaldehy", "HCOOH",
    "NH", "NH2", "NH2OH", "NCN", "N2H2", "NH2NO2", "N2H4", "H2O2",
    "(HCOOH)2", "C6H6(L)", "C7H8(L)", "C8H18(L),n-octa", "Jet-A(L)", "H2O(s)", "H2O(L)",
]

# Convert h,cal = -2999.082 cal/mol to SI J/kg for Python Reactant input.
chos_binder_mw_kg_per_mol = 14.6652984484e-3
chos_binder_h_si = cea.units.cal_to_joule(-2999.082) / chos_binder_mw_kg_per_mol

reactants = [
    "NH4CLO4(I)",
    cea.Reactant(
        name="CHOS-Binder",
        formula={"C": 1.0, "H": 1.86955, "O": 0.031256, "S": 0.008415},
        molecular_weight=14.6652984484,
        enthalpy=chos_binder_h_si,
        temperature=298.15,
    ),
    "AL(cr)",
    "MgO(cr)",
    "H2O(L)",
]

reac = cea.Mixture(reactants)
prod = cea.Mixture(reactants, products_from_reactants=True, omit=omit_names)
solver = cea.EqSolver(prod, reactants=reac)
solution = cea.EqSolution(solver)

h0 = reac.calc_property(cea.ENTHALPY, weights, T_reac)

print(f"{'P (bar)':>10s} {'T (K)':>10s} {'Gamma_s':>10s} {'Cp (cal/g-K)':>14s}")
for p in pressures:
    solver.solve(solution, cea.HP, h0 / cea.R, p, weights)
    cp_cal_per_gk = cea.units.joule_to_cal(solution.cp_eq)
    print(f"{p:10.6f} {solution.T:10.2f} {solution.gamma_s:10.4f} {cp_cal_per_gk:14.4f}")
