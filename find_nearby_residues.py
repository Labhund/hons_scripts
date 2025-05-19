# find_nearby_residues.py
# PyMOL script to find residues within specified distances of selections
# and write the results to separate text files.

from pymol import cmd, stored

def find_nearby_residues(structure_name, selection_name, distances, output_filename):
    """
    Finds residues within specified distances of a selection in a structure
    and writes the information (Chain, ResID, ResName) to a text file.

    Args:
        structure_name (str): The name of the structure object in PyMOL (e.g., "1r1k").
        selection_name (str): The name of the initial selection (e.g., "P1A").
                              This selection must already exist or be selectable (e.g., by resn).
        distances (list): A list of distance cutoffs in Angstroms.
        output_filename (str): The name for the output text file.
    """
    
    # Define the full source selection string
    full_selection_query = f"({selection_name} and {structure_name})"
    
    # Check if the selection exists and has atoms
    if selection_name not in cmd.get_names("selections"):
        print(f"Selection '{selection_name}' not found. Trying 'resn {selection_name}'...")
        cmd.select(selection_name, f"resn {selection_name} and {structure_name}")
        
    if cmd.count_atoms(full_selection_query) == 0:
        print(f"Error: Selection '{selection_name}' in structure '{structure_name}' contains no atoms. Skipping.")
        return

    print(f"Processing: Structure={structure_name}, Selection={selection_name}")
    
    with open(output_filename, 'w') as f:
        f.write(f"# Residues near selection '{selection_name}' in structure '{structure_name}'\n")
        f.write(f"# Source selection query: {full_selection_query}\n")
        f.write("# Format: Chain, ResidueNumber, ResidueName\n")
        f.write("=" * 60 + "\n\n")

        for dist in distances:
            # Define a unique name for the temporary selection at this distance
            nearby_sel_name = f"nearby_{structure_name}_{selection_name}_{str(dist).replace('.', '_')}A"
            
            # Select residues where any atom is within 'dist' Angstroms of the source selection
            # 'byres' ensures whole residues are selected
            # 'all' refers to any atoms in the structure
            cmd.select(nearby_sel_name, f"byres (all within {dist} of {full_selection_query})")
            
            # Optional: Exclude the original selection itself if it's made of residues you don't want listed
            # cmd.select(nearby_sel_name, f"{nearby_sel_name} and not ({full_selection_query})")

            # Use iterate_state to get info for the selected residues
            # Using a set ensures uniqueness based on (chain, resi, resn)
            stored.residues_info = set() 
            cmd.iterate_state(1, nearby_sel_name, "stored.residues_info.add((chain, resi, resn))")

            f.write(f"## Residues within {dist} Å:\n")
            if stored.residues_info:
                # Sort residues for readability (by chain, then residue number)
                sorted_residues = sorted(list(stored.residues_info), key=lambda x: (x[0], int(x[1])))
                for chain, resi, resn in sorted_residues:
                    f.write(f"{chain}, {resi}, {resn}\n")
            else:
                f.write("None\n")
            f.write("\n") # Add a blank line for separation
            
            # Clean up the temporary selection for this distance
            cmd.delete(nearby_sel_name)
            
    print(f"Output for {structure_name}/{selection_name} written to {output_filename}")

# --- Configuration ---

# List of target structures, their initial selections, and desired output filenames
targets = [
    {"struct": "1r1k", "sel": "P1A", "outfile": "1r1k_P1A_nearby_residues.txt"},
    {"struct": "2r40", "sel": "20E", "outfile": "2r40_20E_nearby_residues.txt"},
    {"struct": "1r20", "sel": "HWG", "outfile": "1r20_HWG_nearby_residues.txt"}
]

# List of distance cutoffs in Angstroms (Corrected your list slightly)
distances_A = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5] 

# --- Script Execution ---

print("Starting residue proximity analysis...")

# Optional: Load structures if they aren't already loaded
# Make sure the object names match the 'struct' keys in the targets list
# cmd.fetch("1r1k", async_=0)
# cmd.fetch("2r40", async_=0)
# cmd.fetch("1r20", async_=0)

# Optional: Define your initial selections if they aren't already defined
# Example: Select ligand P1A in object 1r1k
# cmd.select("P1A", "resn P1A and 1r1k") 
# cmd.select("20E", "resn 20E and 2r40") 
# cmd.select("HWG", "resn HWG and 1r20") 

# Iterate through each target configuration and run the analysis
for target in targets:
    find_nearby_residues(target["struct"], 
                         target["sel"], 
                         distances_A, 
                         target["outfile"])

print("Script finished.")
cmd.refresh() # Refresh PyMOL display (optional)


