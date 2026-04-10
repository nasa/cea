% Define the Python environment
pythonExe = '/Users/mleader/miniconda3/envs/cea/bin/python'; % Replace with your Python path
pyenv('Version', pythonExe);

clear; clc;

cea = py.importlib.import_module('cea');

atm_to_bar = 1.01325;
J_to_cal = 1./4.184;
R = 8314.51;

% Species
reac_names = py.list({"H2", "Air"});
prod_names = py.list({"Ar",   "C",   "CO",  "CO2", "H", ...
                      "H2",   "H2O", "HNO", "HO2", "HNO2", ...
                      "HNO3", "N",   "NH",  "NO",  "N2", ...
                      "N2O3", "O",   "O2",  "OH",  "O3"});

% Thermo states
pressures = atm_to_bar*[1.0, 0.1, 0.01];
temperatures = [3000.0, 2000.0];

% Mixture states
fuel_moles = [1.0, 0.0];
oxidant_moles = [0.0, 1.0];
chem_eq_ratios = [1.0, 1.5];

% Store output
n = length(chem_eq_ratios)*length(pressures)*length(temperatures);
of_ratio_out = zeros(n,1);
T_out = zeros(n,1);
P_out = zeros(n,1);
rho = zeros(n,1);
volume = zeros(n,1);
enthalpy = zeros(n,1);
energy = zeros(n,1);
gibbs = zeros(n,1);
entropy = zeros(n,1);
molecular_weight_M = zeros(n,1);
molecular_weight_MW = zeros(n,1);
gamma_s = zeros(n,1);
cp_eq = zeros(n,1);
cp_fr = zeros(n,1);
cv_eq = zeros(n,1);
cv_fr = zeros(n,1);
i = 1;

% Equilibrium solve
for r_eq = chem_eq_ratios
    for p = pressures
        for t = temperatures
            solution = cea.eq_solve(cea.TP, reac_names, T=t, P=p, ...
                                    fuel_amounts=py.numpy.array(fuel_moles), ...
                                    oxid_amounts=py.numpy.array(oxidant_moles), ...
                                    moles=true, ...
                                    r_eq=r_eq, only=prod_names);
            % Store the output
            T_out(i) = t;
            P_out(i) = p/atm_to_bar;
            if solution.converged
                rho(i) = solution.density*1.e-3;
                volume(i) = solution.volume*1.e3;
                enthalpy(i) = solution.enthalpy*J_to_cal;
                energy(i) = solution.energy*J_to_cal;
                gibbs(i) = solution.gibbs_energy*J_to_cal;
                entropy(i) = solution.entropy*J_to_cal;
                molecular_weight_M(i) = solution.M;
                molecular_weight_MW(i) = solution.MW;
                gamma_s(i) = solution.gamma_s;
                cp_eq(i) = solution.cp_eq*J_to_cal;
                cp_fr(i) = solution.cp_fr*J_to_cal;
                cv_eq(i) = solution.cv_eq*J_to_cal;
                cv_fr(i) = solution.cv_fr*J_to_cal;

            i = i + 1;
            end
        end
    end
end

fprintf('Done: %d/%d\n', i-1, n);