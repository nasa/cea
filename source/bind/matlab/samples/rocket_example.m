% Define the Python environment
pythonExe = '/Users/mleader/miniconda3/envs/cea/bin/python'; % Replace with your Python path
pyenv('Version', pythonExe);

clear; clc;

cea = py.importlib.import_module('cea');
ceam = py.importlib.import_module('cea.matlab');

reac_names = py.list({"H2(L)", "O2(L)"});
t_reac = py.numpy.array([20.27, 90.17]);
fuel_amounts = py.numpy.array([1.0, 0.0]);
oxid_amounts = py.numpy.array([0.0, 1.0]);
pi_p = py.numpy.array([10.0, 100.0]);
subar = py.numpy.array([1.58]);
supar = py.numpy.array([25.0]);

solution = ceam.rocket_solve(reac_names, 53.3172, ...
    pi_p=pi_p, subar=subar, supar=supar, ...
    T_reac=t_reac, ...
    fuel_amounts=fuel_amounts, oxid_amounts=oxid_amounts, ...
    of_ratio=5.55157, iac=true);

fprintf('Converged: %d\n', logical(solution.converged));
fprintf('Stations: %d\n', double(solution.num_pts));
fprintf('Chamber T [K]: %.3f\n', double(solution.T.item(int32(0))));
fprintf('Exit P [bar]: %.6f\n', double(solution.P.item(int32(double(solution.num_pts) - 1))));
fprintf('Exit Isp_vac [m/s]: %.3f\n', ...
    double(solution.Isp_vacuum.item(int32(double(solution.num_pts) - 1))));
