% Define the Python environment
pythonExe = '/Users/mleader/miniconda3/envs/cea/bin/python'; % Replace with your Python path
pyenv('Version', pythonExe);

clear; clc;

cea = py.importlib.import_module('cea');
ceam = py.importlib.import_module('cea.matlab');

reac_names = py.list({"H2", "O2"});
fuel_amounts = py.numpy.array([2.0, 0.0]);
oxid_amounts = py.numpy.array([0.0, 1.0]);

solution = ceam.detonation_solve(reac_names, 298.15, 1.0, ...
    fuel_amounts=fuel_amounts, oxid_amounts=oxid_amounts, ...
    r_eq=1.0);

fprintf('Converged: %d\n', logical(solution.converged));
fprintf('Detonation T [K]: %.3f\n', double(solution.T));
fprintf('Detonation P [bar]: %.6f\n', double(solution.P));
fprintf('Velocity [m/s]: %.3f\n', double(solution.velocity));
fprintf('P/P1: %.6f\n', double(solution.P_P1));
