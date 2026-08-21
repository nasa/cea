import numpy as np
from enum import Enum
import os

class Parse(Enum):
    NONE = 0
    THERMO = 1
    SPECIES = 2
    TRANSPORT = 3
    ROCKET = 4
    SHOCK = 5
    DETON = 6

def parse_numeric(value_str):

    value_str = value_str.strip()

    # Handle standard E notation first
    if 'E' in value_str.upper():
        return float(value_str)

    # Split the string into mantissa and exponent
    parts = value_str.replace(' ', '').split('-')

    if len(parts) == 1:  # No negative exponent
        return float(parts[0])

    if len(parts) == 2:
        if (parts[0] == ""): # Negative float
            return -float(parts[1])
        else: # Simple case like "0.8085-02"
            return float(parts[0]) * (10 ** -float(parts[1]))

    if len(parts) == 3 and not parts[0]:  # Negative mantissa case like "-0.8085-02"
        return -float(parts[1]) * (10 ** -float(parts[2]))

    raise ValueError(f"Unable to parse scientific notation: {value_str}")

def normalize_array(array):
    normed_array = []

    idx = 0
    for i, val in enumerate(array):
        if len(val) > 2:
            normed_array.append(val)
            idx += 1
        elif i > 0:
            if "E" in val:
                normed_array[idx-1] += val
            else:
                normed_array[idx-1] += "E"+val

    return normed_array

# TODO: Read the output ".out" files and convert them to a dict
def parse_output(fname):

    # Example output:
    # thermo = {"P":{"units":"bar", vals:np.array([1.0, 0.1, 0.01])},
    #           "T":{"units":"K",   vals:np.array([3000.0, 3200.0, 3400.0])}, ...}
    # amounts = {"H":np.array([0.3232, 0.1333, 0.00422]),
    #            "O":np.array([0.0562, 0.2367, 0.51123]), ...}

    print(f"-- Parsing {fname} --")

    thermo = {}
    amounts = {}
    transport = {}
    rocket = {}
    shock = {}
    deton = {}

    # Initalize flags
    parse_flag = Parse.NONE
    trace_next = False  # Trace species get parsed on the next line
    frozen = False  # Transport property type
    count_blank = 0  # Number of blank lines in a row
    shock_mode = "unknown"
    shock_state = "unknown"

    # Parse the file
    f = open(fname, "r")
    for line in f:
        vals = line.split()

        # Skip blank lines and comments
        if len(vals) == 0:
            count_blank += 1
            if count_blank >= 2:
                parse_flag = Parse.NONE
                trace_next = False
            continue
        else:
            count_blank = 0

        if vals[0] == "!":
            continue

        if vals[0].startswith("WARNING") or vals[0] == "ANSWERS":
            continue

        upper_line = line.upper()
        if "EQUILIBRIUM COMPOSITION FOR INCIDENT SHOCKED" in upper_line:
            shock_mode = "equilibrium"
            continue
        if "FROZEN COMPOSITION FOR INCIDENT SHOCKED" in upper_line:
            shock_mode = "frozen"
            continue

        # Check for data headings
        if vals[0] == "THERMODYNAMIC":
            if len(vals) > 1:
                if vals[1] == "PROPERTIES":
                    parse_flag = Parse.THERMO
                    continue

        if (vals[0] == "CHAMBER") or (vals[0] == "INJECTOR"):
            parse_flag = Parse.THERMO
            continue

        if vals[0] == "INITIAL":
            parse_flag = Parse.SHOCK
            shock_state = "initial"
            continue

        if vals[0] == "SHOCKED":
            parse_flag = Parse.SHOCK
            shock_state = "reflected" if "REFLECTED" in upper_line else "incident"
            if "EQUILIBRIUM" in upper_line:
                shock_mode = "equilibrium"
            elif "FROZEN" in upper_line:
                shock_mode = "frozen"
            continue

        if (vals[0] == "UNBURNED") or (vals[0] == "BURNED"):
            parse_flag = Parse.DETON
            continue

        if (vals[0] == "DETONATION"):
            if len(vals) > 1:
                if (vals[1] == "PARAMETERS"):
                    parse_flag = Parse.DETON
                    continue

        if vals[0] == "MOLE":
            if len(vals) > 1:
                if vals[1] == "FRACTIONS":
                    parse_flag = Parse.SPECIES
                    trace_next = False
                    continue

        if vals[0] == "MASS":
            if len(vals) > 1:
                if vals[1] == "FRACTIONS":
                    parse_flag = Parse.SPECIES
                    trace_next = False
                    continue

        if vals[0] == "TRANSPORT":
            parse_flag = Parse.TRANSPORT
            continue

        if vals[0] == "PERFORMANCE":
            parse_flag = Parse.ROCKET
            continue

        # Parse mixture thermo properties
        if parse_flag == Parse.THERMO:
            # Parse the name and the units
            name_and_units = line[:17].strip()
            property_name = name_and_units.split(",")[0]
            unit_name = ""
            if len(name_and_units.split(",")) > 1:
                unit_name = name_and_units.split(",")[1].strip()

            # Remove the property name and units, and split the line by values
            vals = normalize_array(line[17:].split())

            # Convert the vals to numeric
            for i in range(len(vals)):
                vals[i] = parse_numeric(vals[i])

            # Check if these values are already in the results dict; we are just adding more points
            if (property_name in thermo):  # Adding new entries to existing values
                thermo[property_name]["vals"] = np.append(thermo[property_name]["vals"], np.array(vals))
            else:  # Adding a new property type
                thermo[property_name] = {"units":unit_name, "vals":np.array(vals)}

        # Parse rocket performance parameters
        if parse_flag == Parse.ROCKET:
            # Parse the name and the units
            name_and_units = line[:17].strip()
            property_name = name_and_units.split(",")[0]
            unit_name = ""
            if len(name_and_units.split(",")) > 1:
                unit_name = name_and_units.split(",")[1].strip()

            # Remove the property name and units, and split the line by values
            vals = normalize_array(line[17:].split())

            # Convert the vals to numeric
            for i in range(len(vals)):
                vals[i] = parse_numeric(vals[i])

            # Check if these values are already in the results dict; we are just adding more points
            if (property_name in rocket):  # Adding new entries to existing values
                rocket[property_name]["vals"] = np.append(rocket[property_name]["vals"], np.array(vals))
            else:  # Adding a new property type
                rocket[property_name] = {"units":unit_name, "vals":np.array(vals)}

        if parse_flag == Parse.SHOCK:
            # Parse the name and the units
            name_and_units = line[:17].strip()
            property_name = name_and_units.split(",")[0]
            unit_name = ""
            if len(name_and_units.split(",")) > 1:
                unit_name = name_and_units.split(",")[1].strip()

            # Remove the property name and units, and split the line by values
            vals = normalize_array(line[17:].split())

            # Convert the vals to numeric
            for i in range(len(vals)):
                vals[i] = parse_numeric(vals[i])

            property_key = f"{shock_mode}:{shock_state}:{property_name}"

            # Check if these values are already in the results dict; we are just adding more points
            if (property_key in shock):  # Adding new entries to existing values
                shock[property_key]["vals"] = np.append(shock[property_key]["vals"], np.array(vals))
            else:  # Adding a new property type
                shock[property_key] = {"units":unit_name, "vals":np.array(vals)}

        if parse_flag == Parse.DETON:
            # Parse the name and the units
            name_and_units = line[:17].strip()
            property_name = name_and_units.split(",")[0]
            unit_name = ""
            if len(name_and_units.split(",")) > 1:
                unit_name = name_and_units.split(",")[1].strip()

            # Remove the property name and units, and split the line by values
            vals = normalize_array(line[17:].split())

            # Convert the vals to numeric
            for i in range(len(vals)):
                vals[i] = parse_numeric(vals[i])

            # Check if these values are already in the results dict; we are just adding more points
            if (property_name in deton):  # Adding new entries to existing values
                deton[property_name]["vals"] = np.append(deton[property_name]["vals"], np.array(vals))
            else:  # Adding a new property type
                deton[property_name] = {"units":unit_name, "vals":np.array(vals)}

        # Parse species concentrations
        if parse_flag == Parse.SPECIES:

            # Add trace species as zeros
            if line.strip() == "* THERMODYNAMIC PROPERTIES FITTED TO 20000.K":
                continue
            elif "PRODUCTS WHICH WERE CONSIDERED" in line:
                continue
            elif "NOTE." in line:
                parse_flag = Parse.NONE
                continue
            elif "WERE LESS THAN" in line:
                trace_next = True
                continue
            elif trace_next:
                name_vals = line.split()
                if (name_vals[0] == "O/F") or (name_vals[0] == "NOTE"):
                    trace_next = False
                    parse_flag = Parse.NONE
                    continue
                elif (name_vals[0] == "SHOCKED"):
                    trace_next = False
                    parse_flag = Parse.SHOCK
                    upper_line = line.upper()
                    shock_state = "reflected" if "REFLECTED" in upper_line else "incident"
                    if "EQUILIBRIUM" in upper_line:
                        shock_mode = "equilibrium"
                    elif "FROZEN" in upper_line:
                        shock_mode = "frozen"
                    continue
                else:
                    for species_name in name_vals:
                        trim_name = species_name.strip().strip("*")
                        if (trim_name in amounts):  # Adding new entries to existing values
                            amounts[trim_name] = np.append(amounts[trim_name], np.zeros(npts))
                        else:  # Adding a new species
                            amounts[trim_name] = np.zeros(npts)
                    continue

            line_vals = line.split()
            num_vals = []
            species_name = ""
            for val in line_vals:
                if (val[0].isalpha() or val[0]=="*"):
                    # Process the previous species data
                    if len(num_vals) > 0:
                        vals = normalize_array(num_vals)

                        # Convert the vals to numeric
                        npts = len(vals)
                        for i in range(npts):
                            vals[i] = parse_numeric(vals[i])

                        # Check if these values are already in the results dict; we are just adding more points
                        if (species_name in amounts):  # Adding new entries to existing values
                            amounts[species_name] = np.append(amounts[species_name], np.array(vals))
                        else:  # Adding a new species
                            amounts[species_name] = np.array(vals)

                    # Save the new species name
                    species_name = val.strip().strip("*")
                    num_vals = []
                    continue
                else:
                    num_vals.append(val)

            # Process the last species data
            vals = normalize_array(num_vals)

            # Convert the vals to numeric
            npts = len(vals)
            for i in range(npts):
                vals[i] = parse_numeric(vals[i])

            # Check if these values are already in the results dict; we are just adding more points
            if (species_name in amounts):  # Adding new entries to existing values
                amounts[species_name] = np.append(amounts[species_name], np.array(vals))
            else:  # Adding a new species
                amounts[species_name] = np.array(vals)

        # Parse transport properties
        if parse_flag == Parse.TRANSPORT:

            if "CONDUCTIVITY IN UNITS OF" in line:
                continue
            elif "WITH EQUILIBRIUM REACTIONS" in line:
                frozen = False
                continue
            elif "WITH FROZEN REACTIONS" in line:
                frozen = True
                continue

            # Parse the name and the units
            name_and_units = line[:17].strip()
            property_name = name_and_units.split(",")[0]
            unit_name = ""
            if len(name_and_units.split(",")) > 1:
                unit_name = name_and_units.split(",")[1].strip()

            # Check for valid transport property names
            if property_name in ["Visc", "VISC", "Cp", "Conductivity", "CONDUCTIVITY", "Prandtl Number", "PRANDTL NUMBER"]:

                # Append the suffix to the property name
                if property_name not in ["Visc", "VISC"]:
                    if frozen:
                        property_name += "_fr"
                    else:
                        property_name += "_eq"

                # Remove the property name and units, and split the line by values
                vals = normalize_array(line[17:].split())

                # Convert the vals to numeric
                for i in range(len(vals)):
                    vals[i] = parse_numeric(vals[i])

                # Check if these values are already in the results dict; we are just adding more points
                if (property_name in transport):  # Adding new entries to existing values
                    transport[property_name]["vals"] = np.append(transport[property_name]["vals"], np.array(vals))
                else:  # Adding a new property type
                    transport[property_name] = {"units":unit_name, "vals":np.array(vals)}

                # Frozen transport is last:
                if frozen and property_name in ["Prandtl Number", "PRANDTL NUMBER"]:
                    parse_flag = Parse.NONE

    f.close()

    return thermo, amounts, transport, rocket, shock, deton


if __name__ == "__main__":

    # Save all of the reference results as csv files
    out_dir = "./reference_output"
    for filename in os.listdir(out_dir):
        file_path = os.path.join(out_dir, filename)

        # Check if it's a file (not a directory)
        if os.path.isfile(file_path):
            # Do something with the file
            print("Starting new file: ",filename)
            thermo, amounts, transport, rocket, shock, deton = parse_output(file_path)
            for key in thermo:
                print(key)
                print(thermo[key])
                print()

            for key in amounts:
                print(key, ":", amounts[key])
