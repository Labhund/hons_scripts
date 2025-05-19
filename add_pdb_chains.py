def assign_chains_to_pdb(input_pdb_path, output_pdb_path, chain_assignments_list, default_chain=' '):
    """
    Assigns chain IDs to residues in a PDB file based on residue number ranges.

    Args:
        input_pdb_path (str): Path to the input PDB file.
        output_pdb_path (str): Path to save the modified PDB file.
        chain_assignments_list (list): A list of tuples. Each tuple should be
                                     (start_res_num, end_res_num, chain_id_to_assign).
                                     Ranges are inclusive.
                                     Example: [(1, 270, 'A'), (271, 500, 'B')]
        default_chain (str): Chain ID to use if a residue doesn't fall into any defined range.
    """
    try:
        with open(input_pdb_path, 'r') as infile, open(output_pdb_path, 'w') as outfile:
            for line in infile:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        res_num_str = line[22:26].strip()
                        res_num = int(res_num_str)
                        
                        assigned_chain_this_residue = None
                        for start_res, end_res, chain_id in chain_assignments_list:
                            if start_res <= res_num <= end_res:
                                assigned_chain_this_residue = chain_id
                                break
                        
                        if assigned_chain_this_residue:
                            # PDB chain ID is at column 22 (index 21)
                            modified_line = line[:21] + assigned_chain_this_residue + line[22:]
                            outfile.write(modified_line)
                        else:
                            # If no specific assignment, use default or keep original
                            modified_line = line[:21] + default_chain + line[22:]
                            outfile.write(modified_line)
                    except ValueError:
                        # If residue number is not an integer, write line as is or with default chain
                        # This might happen for some HETATM or non-standard lines
                        modified_line = line[:21] + default_chain + line[22:]
                        outfile.write(modified_line)
                        print(f"Warning: Could not parse residue number for line: {line.strip()}")
                else:
                    outfile.write(line)
        print(f"Successfully processed PDB. Output saved to: {output_pdb_path}")
    except FileNotFoundError:
        print(f"Error: Input PDB file not found at {input_pdb_path}")
    except Exception as e:
        print(f"An error occurred: {e}")

# --- How to use the script ---
# 1. Define the path to your input PDB (the AF2/H++ processed one)
#    and where you want to save the output.
# input_pdb = "../EcR_USP_fromAf2_notails_H++_pH7.4_0.15mM.pdb" # Your actual path
# output_pdb = "../EcR_USP_fromAf2_notails_H++_pH7.4_0.15mM_chains.pdb" # Your desired output path

# 2. VERY IMPORTANT: Define your chain assignments.
#    You need to know the residue number ranges for each protein/chain
#    in your AF2/H++ PDB file.
#    For example, if residues 1-300 are EcR (assign chain 'D')
#    and residues 301-550 are USP (assign chain 'E'):
#
#    my_chain_assignments = [
#        (1, 300, 'D'),    # All residues from 1 to 300 will be chain D
#        (301, 550, 'E')  # All residues from 301 to 550 will be chain E
#    ]
#    Adjust these numbers and chain IDs based on YOUR specific structure.
#    The residue numbers here refer to the numbering in your input AF2/H++ PDB.

# 3. Call the function:
# assign_chains_to_pdb(input_pdb, output_pdb, my_chain_assignments)

# 4. If you have residues that don't fall into these ranges (e.g. ligands, ions)
#    and you want them to have a specific chain ID or be blank, adjust the
#    `default_chain` argument (e.g., default_chain='L' for a ligand).
#    Currently, it defaults to a space.
def assign_chains_to_pdb(input_pdb_path, output_pdb_path, chain_assignments_list, default_chain=' '):
    """
    Assigns chain IDs to residues in a PDB file based on residue number ranges.

    Args:
        input_pdb_path (str): Path to the input PDB file.
        output_pdb_path (str): Path to save the modified PDB file.
        chain_assignments_list (list): A list of tuples. Each tuple should be
                                     (start_res_num, end_res_num, chain_id_to_assign).
                                     Ranges are inclusive.
                                     Example: [(1, 270, 'A'), (271, 500, 'B')]
        default_chain (str): Chain ID to use if a residue doesn't fall into any defined range.
    """
    try:
        with open(input_pdb_path, 'r') as infile, open(output_pdb_path, 'w') as outfile:
            for line in infile:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        res_num_str = line[22:26].strip()
                        res_num = int(res_num_str)
                        
                        assigned_chain_this_residue = None
                        for start_res, end_res, chain_id in chain_assignments_list:
                            if start_res <= res_num <= end_res:
                                assigned_chain_this_residue = chain_id
                                break
                        
                        if assigned_chain_this_residue:
                            # PDB chain ID is at column 22 (index 21)
                            modified_line = line[:21] + assigned_chain_this_residue + line[22:]
                            outfile.write(modified_line)
                        else:
                            # If no specific assignment, use default or keep original
                            modified_line = line[:21] + default_chain + line[22:]
                            outfile.write(modified_line)
                    except ValueError:
                        # If residue number is not an integer, write line as is or with default chain
                        # This might happen for some HETATM or non-standard lines
                        modified_line = line[:21] + default_chain + line[22:]
                        outfile.write(modified_line)
                        print(f"Warning: Could not parse residue number for line: {line.strip()}")
                else:
                    outfile.write(line)
        print(f"Successfully processed PDB. Output saved to: {output_pdb_path}")
    except FileNotFoundError:
        print(f"Error: Input PDB file not found at {input_pdb_path}")
    except Exception as e:
        print(f"An error occurred: {e}")

# --- How to use the script ---
# 1. Define the path to your input PDB (the AF2/H++ processed one)
#    and where you want to save the output.
# input_pdb = "../EcR_USP_fromAf2_notails_H++_pH7.4_0.15mM.pdb" # Your actual path
# output_pdb = "../EcR_USP_fromAf2_notails_H++_pH7.4_0.15mM_chains.pdb" # Your desired output path

# 2. VERY IMPORTANT: Define your chain assignments.
#    You need to know the residue number ranges for each protein/chain
#    in your AF2/H++ PDB file.
#    For example, if residues 1-300 are EcR (assign chain 'D')
#    and residues 301-550 are USP (assign chain 'E'):
#
#    my_chain_assignments = [
#        (1, 300, 'D'),    # All residues from 1 to 300 will be chain D
#        (301, 550, 'E')  # All residues from 301 to 550 will be chain E
#    ]
#    Adjust these numbers and chain IDs based on YOUR specific structure.
#    The residue numbers here refer to the numbering in your input AF2/H++ PDB.

# 3. Call the function:
# assign_chains_to_pdb(input_pdb, output_pdb, my_chain_assignments)

# 4. If you have residues that don't fall into these ranges (e.g. ligands, ions)
#    and you want them to have a specific chain ID or be blank, adjust the
#    `default_chain` argument (e.g., default_chain='L' for a ligand).
#    Currently, it defaults to a space.
# --- How to use the script ---
# 1. Define the path to your input PDB (the AF2/H++ processed one)
#    and where you want to save the output.
input_pdb = "EcR_USP_fromAf2_notails_H++_pH7.4_0.15mM.pdb" # Make sure this path is correct
output_pdb = "EcR_USP_fromAf2_notails_H++_pH7.4_0.15mM_chains.pdb" # Desired output path

# 2. Define your chain assignments based on your observation in PyMOL:
my_chain_assignments = [
    (4, 264, 'A'),    # Residues 4 through 264 will be assigned to chain A
    (315, 551, 'D')  # Residues 315 through 551 will be assigned to chain D
]

# 3. Call the function:
assign_chains_to_pdb(input_pdb, output_pdb, my_chain_assignments)

print(f"PDB processing complete. Check '{output_pdb}' for the result.")

