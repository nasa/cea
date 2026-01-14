import os
import warnings
import pandas as pd

# This program will run the tests in `test_names`, executing the `.inp` files
# with the CEA main interface, and compare the output files generated in `test_dir`
# against reference output files located in `reference_dir`.

# from .test_output.get_output import run_tests
from parse_output import parse_output

# Parameter inputs
print_all = False  # if True, print all comparisons; if False, print only incorrect values
rtol = 1e-2  # Relative acceptance tolerance for all values
round_vals = True  # Attempt to round the test values to the precision of the reference value
test_names = [
    "example1", "example2", "example3", "example4", "example5", "example14",   # Equilibrium problems
    "example8", "example9", "example10", "example11", "example12", "example13", # Rocket problems
    "example7",  # Shock problem
    "example6"   # Deton problem
    ]
reference_dir = "./reference_output/"
test_dir = "./test_output/"

# Initialize values
test_passed = True  # Flag for each individual test case
passed_count = 0  # Count the number of tests that fully passed
num_tests = len(test_names)

# Convert new variable names to the old format for comparison
thermo_val_to_test = {"P":"P", "T":"T", "RHO":"Density", "H":"H", "U":"U", "G":"G", "S":"S",
                      "M":"M", "(dLV/dLP)t":"(dln(V)/dln(P))t", "(dLV/dLT)p":"(dln(V)/dln(T))p",
                      "Cp":"Cp", "GAMMAs":"Gamma_s", "SON VEL":"Son. Vel.", "Pinj/P":"Pinj/P", "Pinf/P":"Pinf/P"}
trans_val_to_test = {"VISC":"Visc", "Cp_fr":"Cp_fr", "CONDUCTIVITY_fr":"Conductivity_fr", "PRANDTL NUMBER_fr":"Prandtl Number_fr",
                     "Cp_eq":"Cp_eq", "CONDUCTIVITY_eq":"Conductivity_eq", "PRANDTL NUMBER_eq":"Prandtl Number_eq"}
rocket_val_to_test = {"Ae/At":"Ae/At", "CSTAR":"C*", "CF":"Cf", "Ivac":"Ivac", "Isp":"Isp"}
shock_val_to_test = thermo_val_to_test | {"MACH NUMBER1":"Mach1", "U1":"u1",
                      "P2/P1":"P2/P1", "T2/T1":"T2/T1", "M2/M1":"M2/M1", "RHO2/RHO1":"rho2/rho1", "V2":"v2",
                      "U5":"u5",
                      "P5/P2":"P5/P2", "T5/T2":"T5/T2", "M5/M2":"M5/M2", "RHO5/RHO2":"rho5/rho2", "U5+V2":"u5+u2"}
deton_val_to_test = thermo_val_to_test | {"P1":"P1", "T1":"T1", "H1":"H1", "M1":"M1", "GAMMA1":"Gamma1", "SON VEL1":"Son. Vel.1",
                                          "P/P1":"P/P1", "T/T1":"T/T1", "M/M1":"M/M1", "RHO/RHO1":"rho/rho1",
                                          "DET MACH NUMBER":"Det. Mach Number", "DET VEL":"Det. Vel."}

def ref_round(ref_val, test_val):
    # Figure out the rounding precision of the reference value
    max_digits = 5
    round_digits = max_digits
    rounded = False
    for i in range(max_digits, -1, -1):
        if ref_val == round(ref_val, i):
            rounded = True
            round_digits = i
        else:
            break

    if rounded:
        return round(test_val, round_digits)

    return test_val

def run_tests(test_names):

    run_dir = "~/git/cea/build-dev/source"

    for test in test_names:
        # Execute the code on the input file
        print(f"Running {test}")
        subprocess.run(run_dir+"/cea"+f" {test}", shell=False, check=True)
        subprocess.run(f"mv {test}.out {test_dir}", shell=False, check=True)
        print()

    return

# Run the tests
run_tests(test_names)

# Initialize arrays for csv output
test_names_csv = []
variable = []
value = []
reference = []
value_type = []
abs_error = []
rel_error = []

for test in test_names:
    print(f"Starting test case: {test}")
    print("----------------------------")
    print()
    test_passed = True

    # Get the validation output
    thermo_ref, amounts_ref, transport_ref, rocket_ref, shock_ref, deton_ref = parse_output(reference_dir+f"{test}.out")

    # Get the test output
    thermo, amounts, transport, rocket, shock, deton = parse_output(test_dir+f"{test}.out")

    # Compare thermo output
    # ---------------------
    for var in thermo_ref:
        if var not in thermo_val_to_test:
            continue
        # Make sure the reference and test arrays are the same length
        ref_len = len(thermo_ref[var]["vals"])
        test_len = len(thermo[thermo_val_to_test[var]]["vals"])
        if ref_len != test_len:
            test_passed = False
            warnings.warn(f"Property {var} has reference length of {ref_len}; test array has length of {test_len}. SKIPPING.")
            print("Reference: ", thermo_ref[var]["vals"])
            print("Test:      ", thermo[thermo_val_to_test[var]]["vals"])
            continue

        for i in range(ref_len):
            ref_val = thermo_ref[var]["vals"][i]
            test_val = thermo[thermo_val_to_test[var]]["vals"][i]

            # Round the test value to the number of digits in the reference value
            if round_vals:
                test_val = ref_round(ref_val, test_val)

            # Compute the absolute and relative error
            abs_err = abs(test_val - ref_val)
            rel_err = abs_err#/abs(ref_val)
            if abs(ref_val) > 1e-20:
                rel_err /= abs(ref_val)

            # Store the values for csv output
            test_names_csv.append(test)
            variable.append(var)
            value.append(test_val)
            reference.append(ref_val)
            value_type.append("thermo")
            abs_error.append(abs_err)
            rel_error.append(rel_err)

            if (rel_err > rtol) or print_all:
                if test_passed:
                    # Print the error headers
                    print("                    Reference    | Test         | Rel. Error ")
                    print("--------------------------------------------------------------")
                test_passed = False
                if (abs(ref_val) > 1e-12) and ((abs(ref_val) < 1e-3) or (abs(ref_val) > 1e6)):
                    if (abs(test_val) > 1e-12) and ((abs(test_val) < 1e-3) or (abs(test_val) > 1e6)):
                        print(f"{var:18s}: {ref_val:12.4e} | {test_val:12.4e} | {100*rel_err:11.3f}%")
                    else:
                        print(f"{var:18s}: {ref_val:12.4e} | {test_val:12.4f} | {100*rel_err:11.3f}%")
                else:
                    if (abs(test_val) > 1e-12) and ((abs(test_val) < 1e-3) or (abs(test_val) > 1e6)):
                        print(f"{var:18s}: {ref_val:12.4f} | {test_val:12.4e} | {100*rel_err:11.3f}%")
                    else:
                        print(f"{var:18s}: {ref_val:12.4f} | {test_val:12.4f} | {100*rel_err:11.3f}%")

    # Compare transport output
    # ------------------------
    for var in transport_ref:
        # Make sure the reference and test arrays are the same length
        ref_len = len(transport_ref[var]["vals"])
        test_len = len(transport[trans_val_to_test[var]]["vals"])
        if ref_len != test_len:
            test_passed = False
            warnings.warn(f"Property {var} has reference length of {ref_len}; test array has length of {test_len}. SKIPPING.")
            print("Reference: ", transport_ref[var]["vals"])
            print("Test:      ", transport[trans_val_to_test[var]]["vals"])
            continue

        for i in range(ref_len):
            ref_val = transport_ref[var]["vals"][i]
            test_val = transport[trans_val_to_test[var]]["vals"][i]

            # Round the test value to the number of digits in the reference value
            if round_vals:
                test_val = ref_round(ref_val, test_val)

            # Compute the absolute and relative error
            abs_err = abs(test_val - ref_val)
            rel_err = abs_err#/abs(ref_val)
            if abs(ref_val) > 1e-20:
                rel_err /= abs(ref_val)

            # Store the values for csv output
            test_names_csv.append(test)
            variable.append(var)
            value.append(test_val)
            reference.append(ref_val)
            value_type.append("transport")
            abs_error.append(abs_err)
            rel_error.append(rel_err)

            if (rel_err > rtol) or print_all:
                if test_passed:
                    # Print the error headers
                    print("                    Reference    | Test         | Rel. Error ")
                    print("--------------------------------------------------------------")
                test_passed = False
                if (abs(ref_val) > 1e-12) and ((abs(ref_val) < 1e-3) or (abs(ref_val) > 1e6)):
                    if (abs(test_val) > 1e-12) and ((abs(test_val) < 1e-3) or (abs(test_val) > 1e6)):
                        print(f"{var:18s}: {ref_val:12.4e} | {test_val:12.4e} | {100*rel_err:11.3f}%")
                    else:
                        print(f"{var:18s}: {ref_val:12.4e} | {test_val:12.4f} | {100*rel_err:11.3f}%")
                else:
                    if (abs(test_val) > 1e-12) and ((abs(test_val) < 1e-3) or (abs(test_val) > 1e6)):
                        print(f"{var:18s}: {ref_val:12.4f} | {test_val:12.4e} | {100*rel_err:11.3f}%")
                    else:
                        print(f"{var:18s}: {ref_val:12.4f} | {test_val:12.4f} | {100*rel_err:11.3f}%")

    # Compare rocket output
    # ---------------------
    for var in rocket_ref:
        if var not in rocket_val_to_test:
            continue
        # Make sure the reference and test arrays are the same length
        ref_len = len(rocket_ref[var]["vals"])
        test_len = len(rocket[rocket_val_to_test[var]]["vals"])
        if ref_len != test_len:
            test_passed = False
            warnings.warn(f"Property {var} has reference length of {ref_len}; test array has length of {test_len}. SKIPPING.")
            print("Reference: ", rocket_ref[var]["vals"])
            print("Test:      ", rocket[rocket_val_to_test[var]]["vals"])
            continue

        for i in range(ref_len):
            ref_val = rocket_ref[var]["vals"][i]
            test_val = rocket[rocket_val_to_test[var]]["vals"][i]

            # Round the test value to the number of digits in the reference value
            if round_vals:
                test_val = ref_round(ref_val, test_val)

            # Compute the absolute and relative error
            abs_err = abs(test_val - ref_val)
            rel_err = abs_err#/abs(ref_val)
            if abs(ref_val) > 1e-20:
                rel_err /= abs(ref_val)

            # Store the values for csv output
            test_names_csv.append(test)
            variable.append(var)
            value.append(test_val)
            reference.append(ref_val)
            value_type.append("rocket")
            abs_error.append(abs_err)
            rel_error.append(rel_err)

            if (rel_err > rtol) or print_all:
                if test_passed:
                    # Print the error headers
                    print("                    Reference    | Test         | Rel. Error ")
                    print("--------------------------------------------------------------")
                test_passed = False
                if (abs(ref_val) > 1e-12) and ((abs(ref_val) < 1e-3) or (abs(ref_val) > 1e6)):
                    if (abs(test_val) > 1e-12) and ((abs(test_val) < 1e-3) or (abs(test_val) > 1e6)):
                        print(f"{var:18s}: {ref_val:12.4e} | {test_val:12.4e} | {100*rel_err:11.3f}%")
                    else:
                        print(f"{var:18s}: {ref_val:12.4e} | {test_val:12.4f} | {100*rel_err:11.3f}%")
                else:
                    if (abs(test_val) > 1e-12) and ((abs(test_val) < 1e-3) or (abs(test_val) > 1e6)):
                        print(f"{var:18s}: {ref_val:12.4f} | {test_val:12.4e} | {100*rel_err:11.3f}%")
                    else:
                        print(f"{var:18s}: {ref_val:12.4f} | {test_val:12.4f} | {100*rel_err:11.3f}%")

    # Compare species output
    # ----------------------
    for name in amounts_ref:
        # Check that this species exists in the test output
        if not name in amounts:
            test_passed = False
            warnings.warn(f"Species {name} not found in test output.")
            continue

        # Make sure the reference and test arrays are the same length
        ref_len = len(amounts_ref[name])
        test_len = len(amounts[name])
        if ref_len != test_len:
            test_passed = False
            warnings.warn(f"Species {name} has reference length of {ref_len}; test array has length of {test_len}. SKIPPING.")
            print("Reference: ", amounts_ref[name])
            print("Test:      ", amounts[name])
            continue

        for i in range(ref_len):
            ref_val = amounts_ref[name][i]
            test_val = amounts[name][i]

            # Round the test value to the number of digits in the reference value
            if round_vals:
                test_val = ref_round(ref_val, test_val)

            abs_err = abs(test_val - ref_val)
            rel_err = abs_err#/abs(ref_val)
            if abs(ref_val) > 1e-20:
                rel_err /= abs(ref_val)

            # Store the values for csv output
            test_names_csv.append(test)
            variable.append(name)
            value.append(test_val)
            reference.append(ref_val)
            value_type.append("species")
            abs_error.append(abs_err)
            rel_error.append(rel_err)

            if (rel_err > rtol) or print_all:
                if test_passed:
                    # Print the error headers
                    print("                    Reference    | Test         | Rel. Error ")
                    print("--------------------------------------------------------------")
                test_passed = False
                if (abs(ref_val) > 1e-12) and ((abs(ref_val) < 1e-3) or (abs(ref_val) > 1e6)):
                    if (abs(test_val) > 1e-12) and ((abs(test_val) < 1e-3) or (abs(test_val) > 1e6)):
                        print(f"{name:18s}: {ref_val:12.4e} | {test_val:12.4e} | {100*rel_err:11.3f}%")
                    else:
                        print(f"{name:18s}: {ref_val:12.4e} | {test_val:12.4f} | {100*rel_err:11.3f}%")
                else:
                    if (abs(test_val) > 1e-12) and ((abs(test_val) < 1e-3) or (abs(test_val) > 1e6)):
                        print(f"{name:18s}: {ref_val:12.4f} | {test_val:12.4e} | {100*rel_err:11.3f}%")
                    else:
                        print(f"{name:18s}: {ref_val:12.4f} | {test_val:12.4f} | {100*rel_err:11.3f}%")

    if test_passed:
        passed_count += 1
    print()

print(f"------- {passed_count}/{num_tests} tests passed. -------")

# Save the results to a CSV file
df = pd.DataFrame({
    "test_name": test_names_csv,
    "value_type": value_type,
    "variable": variable,
    "value": value,
    "reference": reference,
    "abs_error": abs_error,
    "rel_error": rel_error
})
df.to_csv("test_results.csv", index=False)