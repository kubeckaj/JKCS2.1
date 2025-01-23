import sys
import os


def display_help():
    print("This script converts the coordinates of an xyz file from Bohr to Angstroms.")
    print("Conversion factor: 0.529177")
    print("Input file: xyz")
    sys.exit(0)


def convert_xyz(input_file):
    conversion_factor = 0.529177
    output_file = input_file.replace(".xyz", "_converted.xyz")

    try:
        with open(input_file, "r") as infile, open(output_file, "w") as outfile:
            converted_lines = 0
            while True:
                # Read number of atoms
                num_atoms_line = infile.readline()
                if not num_atoms_line:
                    break  # End of file
                try:
                    num_atoms = int(num_atoms_line.strip())
                except ValueError:
                    print(f"Error: Unexpected format for number of atoms in '{input_file}'.")
                    return

                # Read the comment line
                comment_line = infile.readline()
                if not comment_line:
                    print(f"Error: Missing comment line in '{input_file}'.")
                    return

                # Write the first two lines to the output
                outfile.write(num_atoms_line)
                outfile.write(comment_line)

                # Convert and write the atom coordinates
                for _ in range(num_atoms):
                    line = infile.readline()
                    if not line:
                        print(f"Error: Missing atom coordinate lines in '{input_file}'.")
                        return
                    parts = line.split()
                    if len(parts) != 4:
                        print(f"Error: Invalid line format in '{input_file}': {line.strip()}")
                        return
                    atom, x, y, z = parts
                    x_angstrom = float(x) * conversion_factor
                    y_angstrom = float(y) * conversion_factor
                    z_angstrom = float(z) * conversion_factor
                    outfile.write(
                        f"{atom:<2} {x_angstrom:>12.8f} {y_angstrom:>12.8f} {z_angstrom:>12.8f}\n"
                    )

                    converted_lines += 1

            if converted_lines == 0:
                print(f"Error: No coordinates were converted in '{input_file}'. Please check the input file.")
                return

        print(f"Conversion complete for '{input_file}'. Output written to '{output_file}'.")

    except Exception as e:
        print(f"Error processing '{input_file}': {e}")


def main():
    # Check for input arguments
    if len(sys.argv) == 1:
        print("Error: Please provide one or more xyz files as 'JKBohr2A.py file1.xyz file2.xyz ...'")
        sys.exit(1)

    if sys.argv[1] in ["-help", "--help"]:
        display_help()

    # Process each input file
    for input_file in sys.argv[1:]:
        if not input_file.endswith(".xyz"):
            print(f"Error: '{input_file}' is not an .xyz file. Skipping.")
            continue
        if not os.path.isfile(input_file):
            print(f"Error: File '{input_file}' does not exist. Skipping.")
            continue

        convert_xyz(input_file)


if __name__ == "__main__":
    main()
