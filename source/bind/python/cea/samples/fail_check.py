import numpy as np

# importing cea prints the loaded thermo.lib and trans.lib
# setting cea.set_log_level(cea.LOG_NONE) stops these
# logs, but results in lots of other logs showing up
# maybe have the prints for this info be only for log level
# debug or higher? Basically, would like them not to be
# shown under normal import conditions. maybe have simple
# variables of the path to them for inspection in code?
import cea

reac_names = ["H2(L)", "O2(L)"]
T_reactant = np.array([20.27, 90.17])  # Reactant temperatures (K)
fuel_weights = np.array([1.0, 0.0])
oxidant_weights = np.array([0.0, 1.0])
of_ratio = 5.55157
pc = 53.3172  # Chamber pressure (bar)
supar = [1.0]  # Supersonic area ratio <-- this will error out in solver.solve()

reac = cea.Mixture(reac_names)
prod = cea.Mixture(reac_names, products_from_reactants=True)

solver = cea.RocketSolver(prod, reactants=reac)
solution = cea.RocketSolution(solver)

weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio)
hc = reac.calc_property(cea.ENTHALPY, weights, T_reactant)/cea.R

# the underlying compiled fortran library will call "stop" for a supersonic
# area ratio less than 1.0, but this also kills the python interpreter
try:
    solver.solve(solution, weights, pc, supar=supar, hc=hc, iac=True)
except:
    print("solve failed")

print("This never print!")