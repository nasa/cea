Legacy Interface
****************

The CEA legacy interface is designed for backward compatibility with older versions of the software (ie CEA2). The CEA2 interface uses the text-based file format for input and output described in RP-1311 Part 2 [1]_.

Usage
-----

To use the legacy interface, call the compiled cea module from the command line with the following syntax:
```
./cea <input_file>
```
where `<input_file>` is the path to the input file in the legacy format, without the `.inp` extension. The output will be written to a file with the same name as the input file with a `.out` extension.
By default, the compiled cea module will be placed in the `build` directory after compiling the program.

.. [1] McBride, B.J., Gordon, S., "Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications II. Users Manual and Program Description: Users Manual and Program Description - 2",
    NASA RP-1311, 1996. [NTRS](https://ntrs.nasa.gov/citations/19960044559)
