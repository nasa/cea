% Define the Python environment
pythonExe = '/Users/mleader/miniconda3/envs/cea/bin/python'; % Replace with your Python path
pyenv('Version', pythonExe);

clear; clc;

cea = py.importlib.import_module('cea');
ceam = py.importlib.import_module('cea.matlab');

reac_names = py.list({"H2", "O2", "Ar"});
fuel_amounts = py.numpy.array([0.05, 0.0, 0.0]);
oxid_amounts = py.numpy.array([0.0, 0.05, 0.9]);
p0 = cea.units.mmhg_to_bar(10.0);

solution = ceam.shock_solve(reac_names, 300.0, p0, ...
    u1=1100.0, ...
    fuel_amounts=fuel_amounts, oxid_amounts=oxid_amounts, ...
    moles=true, reflected=true);

fprintf('Converged: %d\n', logical(solution.converged));
fprintf('Incident T [K]: %.3f\n', double(solution.T.item(int32(1))));
fprintf('Incident P [bar]: %.6f\n', double(solution.P.item(int32(1))));
fprintf('Reflected T [K]: %.3f\n', double(solution.T.item(int32(2))));
fprintf('P2/P1: %.6f\n', double(solution.P21));
fprintf('P5/P2: %.6f\n', double(solution.P52));
