import argparse
import sys
from collections import defaultdict

# Dictionary mapping three-letter amino acid codes to one-letter codes
three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def extract_pdb_sequence_data(file_path):
    """
    Extracts protein sequence data (one-letter code, residue number) per chain
    from a PDB file.

    Args:
        file_path (str): The path to the PDB file.

    Returns:
        dict: A dictionary where keys are chain IDs and values are lists of
              tuples (one-letter code, residue number).
              Returns an empty dict if the file cannot be processed or is empty.
    """
    # Use defaultdict to easily append to lists for each chain
    chain_data = defaultdict(list)
    last_residue_id = None

    try:
        with open(file_path, 'r') as f:
            print(f"Processing file: {file_path}", file=sys.stderr)
            for line in f:
                if line.startswith("ATOM"):
                    try:
                        chain_id = line[21]
                        res_num_str = line[22:26].strip()
                        res_name_3 = line[17:20].strip()

                        if not res_num_str:
                            continue

                        res_num = int(res_num_str)
                        current_residue_id = (chain_id, res_num)

                        if current_residue_id != last_residue_id:
                            res_name_1 = three_to_one.get(res_name_3, 'X')
                            # Append (one-letter code, residue number) to the list for this chain
                            chain_data[chain_id].append((res_name_1, res_num))
                            last_residue_id = current_residue_id
                    except (IndexError, ValueError) as e:
                        print(f"Skipping line due to parsing error: {line.strip()} - Error: {e}", file=sys.stderr)
                        continue
        return dict(chain_data) # Convert back to a regular dict for clarity

    except FileNotFoundError:
        print(f"Error: File not found at {file_path}", file=sys.stderr)
        return {}
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        return {}

def main():
    parser = argparse.ArgumentParser(
        description="Extract protein sequence (one-letter codes), start/end residue numbers, and residue information per chain from a PDB file.",
        epilog="Example usage: python extract_pdb_seq.py -f ./1r1k_apo_full.pdb"
    )
    parser.add_argument(
        "-f", "--file",
        required=True,
        help="Path to the input PDB file.",
        metavar="FILE_PATH"
    )
    args = parser.parse_args()

    # Extract sequence data grouped by chain
    extracted_data = extract_pdb_sequence_data(args.file)

    if extracted_data:
        print("\nSequence String (per chain):")
        # Sort chains by ID for consistent output
        for chain_id in sorted(extracted_data.keys()):
            residues = extracted_data[chain_id]
            if not residues: # Skip if a chain somehow ended up empty
                continue

            # Extract residue numbers and one-letter codes
            residue_numbers = [res_num for _, res_num in residues]
            sequence_string = "".join([res_code for res_code, _ in residues])

            # Find start and end residue numbers
            start_res = min(residue_numbers)
            end_res = max(residue_numbers)

            # Print the results for the current chain
            print(f"Chain {chain_id} (Res {start_res}-{end_res}): {sequence_string}")

    else:
        print("No sequence data extracted.", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()

