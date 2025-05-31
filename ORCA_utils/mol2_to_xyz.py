import argparse
import os # For path manipulation in the example
import datetime # For adding date to XYZ header

# A utility to convert a mol2 file to an xyz file for use in ORCA.
def mol2_to_xyz(mol2_file, xyz_file):
    """
    Convert MOL2 file to XYZ format for ORCA calculations
    """
    with open(mol2_file, 'r') as f:
        lines = f.readlines()

    # Attempt to parse molecule name for the XYZ comment
    molecule_name = "molecule" # Default name
    if lines[0].startswith("@<TRIPOS>MOLECULE"):
        if len(lines) > 1:
            name_line = lines[1].strip()
            if name_line: # Ensure the line is not empty
                molecule_name = name_line.split()[0] 
    
    # Find the @<TRIPOS>ATOM section
    atom_start = None
    atom_end = None
    
    for i, line in enumerate(lines):
        if '@<TRIPOS>ATOM' in line:
            atom_start = i + 1
        elif atom_start is not None and line.startswith('@<TRIPOS>'):
            atom_end = i
            break
    
    if atom_start is None:
        raise ValueError(f"Could not find @<TRIPOS>ATOM section in {mol2_file}")
    
    if atom_end is None: # If no other @<TRIPOS> section found after ATOM
        # Search for the end of atom lines by checking for non-empty lines
        # that don't conform to atom line structure (e.g. @<TRIPOS>BOND)
        # or simply end of file.
        for i in range(atom_start, len(lines)):
            if lines[i].startswith('@<TRIPOS>') or not lines[i].strip():
                atom_end = i
                break
        if atom_end is None: # Reached end of file
            atom_end = len(lines)

    # Parse atoms
    atoms = []
    for i in range(atom_start, atom_end):
        line_content = lines[i].strip()
        if line_content:  # Skip empty lines
            parts = line_content.split()
            if len(parts) >= 6: # Basic check for enough columns
                # atom_id = parts[0] # Not used in XYZ
                # atom_name = parts[1] # Not used in XYZ
                try:
                    x = float(parts[2])
                    y = float(parts[3])
                    z = float(parts[4])
                    atom_type = parts[5]
                    # Get element from atom type (e.g., C.3 -> C, h -> H)
                    element = atom_type.split('.')[0].capitalize()
                    atoms.append((element, x, y, z))
                except ValueError:
                    print(f"Warning: Skipping malformed atom line: {line_content}")
            # else: # Optional: warn about lines that don't have enough parts
                # print(f"Warning: Skipping line with insufficient parts in ATOM section: {line_content}")
    
    if not atoms:
        raise ValueError(f"No atoms found or parsed from {mol2_file}")

    # Write XYZ file
    current_date = datetime.date.today().isoformat()
    # Get absolute path for the input file for the header
    abs_mol2_file_path = os.path.abspath(mol2_file)
    comment_line = (
        f"{molecule_name} | Source: {abs_mol2_file_path} | "
        f"Converted: {current_date} | Script Author: Markus Williams"
    )
    with open(xyz_file, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"{comment_line}\n")
        for element, x, y, z in atoms:
            f.write(f"{element:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n")
    
    # This print is useful for command-line operation, handled in main()
    # print(f"Converted {len(atoms)} atoms from {mol2_file} to {xyz_file}") 
    return len(atoms)

def main():
    parser = argparse.ArgumentParser(
        description="Convert a MOL2 file to an XYZ file for ORCA.",
        formatter_class=argparse.RawTextHelpFormatter # Allows for better formatting in help text
    )
    parser.add_argument(
        "-i", "--input",
        dest="mol2_file",
        help="Path to the input MOL2 file."
    )
    parser.add_argument(
        "-o", "--output",
        dest="xyz_file",
        help="Path to the output XYZ file."
    )
    parser.add_argument(
        "--run_example",
        action="store_true",
        help="Run with an internal example MOL2 content.\n"
             "If this flag is used, -i/--input and -o/--output are ignored."
    )

    args = parser.parse_args()

    if args.run_example:
        # Define sample mol2 data for the example
        sample_mol2_content = """@<TRIPOS>MOLECULE
example_mol
3 2 0 0 0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 C1 0.0 0.0 0.0 C.3 1 LIG 0.0
2 H1 1.0 0.0 0.0 H 1 LIG 0.0
3 O1 0.0 1.0 0.0 O.2 1 LIG 0.0
@<TRIPOS>BOND
1 1 2 1
2 1 3 2
"""
        example_mol2_file = "example_temp.mol2"
        example_xyz_file = "example_temp.xyz"
        with open(example_mol2_file, 'w') as f:
            f.write(sample_mol2_content)
        
        print(f"Running example with internal data, input: {example_mol2_file}, output: {example_xyz_file}")
        try:
            num_atoms = mol2_to_xyz(example_mol2_file, example_xyz_file)
            print(f"Successfully converted {num_atoms} atoms from {example_mol2_file} to {example_xyz_file}")
            with open(example_xyz_file, 'r') as f_xyz:
                print("\nGenerated XYZ file content:")
                print(f_xyz.read())
        except Exception as e:
            print(f"Error during example processing: {e}")
        finally:
            # Clean up temporary example files
            if os.path.exists(example_mol2_file):
                os.remove(example_mol2_file)
            if os.path.exists(example_xyz_file):
                os.remove(example_xyz_file)
    else:
        if not args.mol2_file or not args.xyz_file:
            parser.error(
                "the following arguments are required: -i/--input, -o/--output\n"
                "(unless --run_example is specified)"
            )
        try:
            num_atoms = mol2_to_xyz(args.mol2_file, args.xyz_file)
            print(f"Successfully converted {num_atoms} atoms from {args.mol2_file} to {args.xyz_file}")
        except FileNotFoundError:
            print(f"Error: Input MOL2 file not found at {args.mol2_file}")
        except ValueError as ve:
            print(f"Error processing MOL2 file: {ve}")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()

# The previous test code block is now integrated into the main() function
# when --run_example is used, or directly uses command-line arguments.
# Removed the old hardcoded test block from here.
