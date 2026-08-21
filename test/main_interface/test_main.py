import os
import shutil
import subprocess
import warnings
import csv
from pathlib import Path
import numpy as np

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
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
reference_dir = SCRIPT_DIR / "reference_output"
test_dir = SCRIPT_DIR / "test_output"


def resolve_cea_exe(repo_root: Path) -> str:
    cea_exe = os.environ.get("CEA_EXE")
    if cea_exe is not None:
        return cea_exe

    candidate_exes = []

    run_dir = os.environ.get("CEA_RUN_DIR")
    if run_dir is not None:
        candidate_exes.append(Path(run_dir).expanduser() / "cea")

    # Prefer binaries built from this checkout to avoid stale external builds.
    candidate_exes.extend(sorted(repo_root.glob("build*/source/cea")))

    # Keep the historical default location as a final fallback.
    candidate_exes.append(Path("~/git/cea/build-dev/source/cea").expanduser())

    existing_exes = [path for path in candidate_exes if path.exists()]
    if existing_exes:
        return str(max(existing_exes, key=lambda path: path.stat().st_mtime))

    cea_exe = shutil.which("cea")
    if cea_exe is not None:
        return cea_exe

    raise FileNotFoundError(
        "Could not locate `cea` executable. Set CEA_EXE/CEA_RUN_DIR or build `build*/source/cea` in this repository."
    )

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
shock_val_to_test = {
    "MACH NUMBER1":"Mach1", "U1":"u1", "P":"P", "T":"T", "H":"H", "U":"U", "G":"G", "S":"S",
    "M":"M", "(dLV/dLP)t":"(dln(V)/dln(P))t", "(dLV/dLT)p":"(dln(V)/dln(T))p",
    "Cp":"Cp", "GAMMAs":"Gamma_s", "SON VEL":"Son. Vel.", "U2":"u2",
    "P2/P1":"P2/P1", "T2/T1":"T2/T1", "M2/M1":"M2/M1", "RHO2/RHO1":"rho2/rho1", "V2":"v2",
    "U5":"u5", "P5/P2":"P5/P2", "T5/T2":"T5/T2", "M5/M2":"M5/M2",
    "RHO5/RHO2":"rho5/rho2", "U5+V2":"u5+v2",
}
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


def resolve_dataset_key(dataset, key):
    if key in dataset:
        return key

    fallback_keys = {
        "Pinf/P": "Pinj/P",
        "Pinj/P": "Pinf/P",
    }
    fallback_key = fallback_keys.get(key)
    if fallback_key in dataset:
        return fallback_key

    return None


def find_block_starts(vals):
    if len(vals) < 3:
        return [0]

    starts = [0]
    first_val = vals[0]
    second_val = vals[1]
    for i in range(1, len(vals) - 1):
        if np.isclose(vals[i], first_val, rtol=1e-6, atol=1e-12) and np.isclose(vals[i + 1], second_val, rtol=1e-6, atol=1e-12):
            starts.append(i)

    return starts


def last_block_length(dataset, control_keys):
    for control_key in control_keys:
        dataset_key = resolve_dataset_key(dataset, control_key)
        if dataset_key is None:
            continue
        starts = find_block_starts(dataset[dataset_key]["vals"])
        if len(starts) > 1:
            return len(dataset[dataset_key]["vals"]) - starts[-1]

    return None


def block_sizes(dataset, control_keys):
    for control_key in control_keys:
        dataset_key = resolve_dataset_key(dataset, control_key)
        if dataset_key is None:
            continue
        starts = find_block_starts(dataset[dataset_key]["vals"])
        if len(starts) == 1:
            continue
        sizes = []
        for i, start in enumerate(starts):
            end = starts[i + 1] if i + 1 < len(starts) else len(dataset[dataset_key]["vals"])
            sizes.append(end - start)
        return sizes

    return []


def alignment_score(ref_vals, test_vals):
    score = 0.0
    for ref_val, test_val in zip(ref_vals, test_vals):
        if round_vals:
            test_val = ref_round(ref_val, test_val)
        abs_err = abs(test_val - ref_val)
        rel_err = abs_err
        if abs(ref_val) > 1e-20:
            rel_err /= abs(ref_val)
        score += rel_err
    return score


def select_aligned_series(ref_vals, test_vals, candidate_lengths=None):
    if candidate_lengths is None:
        candidate_lengths = []

    valid_lengths = []
    for candidate_len in candidate_lengths:
        if candidate_len is None:
            continue
        if candidate_len <= len(ref_vals) and candidate_len <= len(test_vals):
            valid_lengths.append(candidate_len)

    if len(ref_vals) != len(test_vals):
        valid_lengths.append(min(len(ref_vals), len(test_vals)))

    valid_lengths = sorted(set(valid_lengths), reverse=True)
    if not valid_lengths:
        return ref_vals, test_vals

    best = None
    for compare_len in valid_lengths:
        for ref_start in range(len(ref_vals) - compare_len + 1):
            ref_window = ref_vals[ref_start:ref_start + compare_len]
            for test_start in range(len(test_vals) - compare_len + 1):
                test_window = test_vals[test_start:test_start + compare_len]
                score = alignment_score(ref_window, test_window)
                candidate = (score, -compare_len, ref_start + test_start, ref_window, test_window)
                if best is None or candidate < best:
                    best = candidate

    return best[3], best[4]

def run_tests(test_names):
    cea_exe = resolve_cea_exe(REPO_ROOT)

    test_dir.mkdir(parents=True, exist_ok=True)

    for test in test_names:
        # Execute the code on the input file
        print(f"Running {test}")
        input_base = SCRIPT_DIR / test
        subprocess.run([cea_exe, str(input_base)], cwd=REPO_ROOT, check=True)
        out_file = input_base.with_suffix(".out")
        if not out_file.exists():
            raise FileNotFoundError(f"Expected output file not found: {out_file}")
        shutil.move(str(out_file), str(test_dir / out_file.name))
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
    thermo_ref, amounts_ref, transport_ref, rocket_ref, shock_ref, deton_ref = parse_output(str(reference_dir / f"{test}.out"))

    # Get the test output
    thermo, amounts, transport, rocket, shock, deton = parse_output(str(test_dir / f"{test}.out"))

    thermo_block_len_ref = last_block_length(thermo_ref, ["Pinf/P", "Pinj/P"])
    thermo_block_len_test = last_block_length(thermo, ["Pinf/P", "Pinj/P"])
    rocket_block_len_ref = last_block_length(rocket_ref, ["Ae/At"])
    rocket_block_len_test = last_block_length(rocket, ["Ae/At"])
    thermo_block_sizes_ref = block_sizes(thermo_ref, ["Pinf/P", "Pinj/P"])
    thermo_block_sizes_test = block_sizes(thermo, ["Pinf/P", "Pinj/P"])
    rocket_block_sizes_ref = block_sizes(rocket_ref, ["Ae/At"])
    rocket_block_sizes_test = block_sizes(rocket, ["Ae/At"])

    # Compare thermo output
    # ---------------------
    for var in thermo_ref:
        if var not in thermo_val_to_test:
            continue
        test_key = resolve_dataset_key(thermo, thermo_val_to_test[var])
        if test_key is None:
            test_passed = False
            warnings.warn(f"Property {var} not found in test output. SKIPPING.")
            continue

        ref_vals, test_vals = select_aligned_series(
            thermo_ref[var]["vals"],
            thermo[test_key]["vals"],
            candidate_lengths=thermo_block_sizes_ref + thermo_block_sizes_test + [thermo_block_len_ref, thermo_block_len_test],
        )

        # Make sure the reference and test arrays are the same length
        ref_len = len(ref_vals)
        test_len = len(test_vals)
        if ref_len != test_len:
            test_passed = False
            warnings.warn(f"Property {var} has reference length of {ref_len}; test array has length of {test_len}. SKIPPING.")
            print("Reference: ", ref_vals)
            print("Test:      ", test_vals)
            continue

        for i in range(ref_len):
            ref_val = ref_vals[i]
            test_val = test_vals[i]

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
        test_key = resolve_dataset_key(rocket, rocket_val_to_test[var])
        if test_key is None:
            test_passed = False
            warnings.warn(f"Property {var} not found in test output. SKIPPING.")
            continue

        ref_vals, test_vals = select_aligned_series(
            rocket_ref[var]["vals"],
            rocket[test_key]["vals"],
            candidate_lengths=rocket_block_sizes_ref + rocket_block_sizes_test + [rocket_block_len_ref, rocket_block_len_test],
        )
        # Make sure the reference and test arrays are the same length
        ref_len = len(ref_vals)
        test_len = len(test_vals)
        if ref_len != test_len:
            test_passed = False
            warnings.warn(f"Property {var} has reference length of {ref_len}; test array has length of {test_len}. SKIPPING.")
            print("Reference: ", ref_vals)
            print("Test:      ", test_vals)
            continue

        for i in range(ref_len):
            ref_val = ref_vals[i]
            test_val = test_vals[i]

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

    # Compare shock output
    # --------------------
    for ref_key in shock_ref:
        mode, state, ref_var = ref_key.split(":", 2)
        if ref_var == "RHO":
            test_var = "rho1" if state == "initial" else "rho"
        else:
            test_var = shock_val_to_test.get(ref_var)
        if test_var is None:
            test_passed = False
            warnings.warn(f"Shock property {ref_key} has no comparison mapping. SKIPPING.")
            continue

        test_key = f"{mode}:{state}:{test_var}"
        if test_key not in shock:
            test_passed = False
            warnings.warn(f"Shock property {ref_key} not found as {test_key} in test output. SKIPPING.")
            continue

        ref_vals = shock_ref[ref_key]["vals"]
        test_vals = shock[test_key]["vals"]
        ref_len = len(ref_vals)
        test_len = len(test_vals)
        if ref_len != test_len:
            test_passed = False
            warnings.warn(
                f"Shock property {ref_key} has reference length of {ref_len}; "
                f"test array has length of {test_len}. SKIPPING."
            )
            continue

        for i in range(ref_len):
            ref_val = ref_vals[i]
            test_val = test_vals[i]
            if round_vals:
                test_val = ref_round(ref_val, test_val)

            abs_err = abs(test_val - ref_val)
            rel_err = abs_err
            if abs(ref_val) > 1e-20:
                rel_err /= abs(ref_val)

            test_names_csv.append(test)
            variable.append(ref_key)
            value.append(test_val)
            reference.append(ref_val)
            value_type.append("shock")
            abs_error.append(abs_err)
            rel_error.append(rel_err)

            if (rel_err > rtol) or print_all:
                if test_passed:
                    print("                    Reference    | Test         | Rel. Error ")
                    print("--------------------------------------------------------------")
                test_passed = False
                print(f"{ref_key:36s}: {ref_val:12.4e} | {test_val:12.4e} | {100*rel_err:11.3f}%")

    # Compare species output
    # ----------------------
    for name in amounts_ref:
        # Check that this species exists in the test output
        if not name in amounts:
            test_passed = False
            warnings.warn(f"Species {name} not found in test output.")
            continue

        ref_vals, test_vals = select_aligned_series(
            amounts_ref[name],
            amounts[name],
            candidate_lengths=thermo_block_sizes_ref + thermo_block_sizes_test + [thermo_block_len_ref, thermo_block_len_test],
        )

        # Make sure the reference and test arrays are the same length
        ref_len = len(ref_vals)
        test_len = len(test_vals)
        if ref_len != test_len:
            test_passed = False
            warnings.warn(f"Species {name} has reference length of {ref_len}; test array has length of {test_len}. SKIPPING.")
            print("Reference: ", ref_vals)
            print("Test:      ", test_vals)
            continue

        for i in range(ref_len):
            ref_val = ref_vals[i]
            test_val = test_vals[i]

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

# Save the results to a CSV file without requiring pandas.
with open("test_results.csv", "w", newline="") as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(["test_name", "value_type", "variable", "value", "reference", "abs_error", "rel_error"])
    for row in zip(test_names_csv, value_type, variable, value, reference, abs_error, rel_error):
        writer.writerow(row)
