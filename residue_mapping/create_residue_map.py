import argparse
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio import pairwise2
# from Bio.pairwise2 import format_alignment

# Global variable to hold the conversion dictionary and its source
conversion_dict = None
conversion_dict_source = "None"
# To store our defined function, initialized to None
_three_to_one_actual_func = None

try:
    from Bio.SeqUtils import three_to_one as biopython_sequtils_three_to_one
    # Test if Bio.SeqUtils.three_to_one is functional
    try:
        biopython_sequtils_three_to_one("ALA") # Test with a standard code
        def _temp_func_sequtils(res_name):
            # Bio.SeqUtils.three_to_one usually expects uppercase or handles case.
            cleaned_res_name = res_name.strip().upper()
            return biopython_sequtils_three_to_one(cleaned_res_name)
        _three_to_one_actual_func = _temp_func_sequtils
        conversion_dict_source = "Bio.SeqUtils.three_to_one"
        print("Successfully configured to use functional Bio.SeqUtils.three_to_one.")
    except Exception as e:
        print(f"Bio.SeqUtils.three_to_one imported but non-functional ({e}). Falling back.")
        _three_to_one_actual_func = None # Explicitly mark as not usable

except ImportError:
    print("Could not import three_to_one from Bio.SeqUtils. Trying Bio.Data.IUPACData...")
    _three_to_one_actual_func = None # Ensure it's None if import fails

if _three_to_one_actual_func is None: # Fallback if Bio.SeqUtils.three_to_one was not set
    try:
        from Bio.Data.IUPACData import protein_letters_3to1
        conversion_dict = protein_letters_3to1 # This is a direct reference to the dict
        conversion_dict_source = "Bio.Data.IUPACData.protein_letters_3to1"
        print(f"Successfully imported protein_letters_3to1 from Bio.Data.IUPACData. Type: {type(conversion_dict)}")
        if isinstance(conversion_dict, dict) and conversion_dict:
            sample_key_example = next(iter(conversion_dict.keys()))
            print(f"Sample keys from protein_letters_3to1 (e.g., '{sample_key_example}'): {list(conversion_dict.keys())[:10]}")
        else:
            print("Warning: protein_letters_3to1 is not a dictionary or is empty!")
            raise ImportError("protein_letters_3to1 unusable")

        # Define the converter for Bio.Data.IUPACData path
        def _temp_func_iupacdata(res_name):
            global conversion_dict # Uses the IUPACData dict
            # ***** THIS IS THE KEY CORRECTION *****
            # Keys are TitleCase like 'Ala', 'Val' in Biopython 1.85
            cleaned_res_name = res_name.strip().capitalize()
            if not isinstance(conversion_dict, dict): # Should not happen if import was successful
                raise TypeError(f"{conversion_dict_source} is not a dictionary at point of use.")
            return conversion_dict[cleaned_res_name]
        _three_to_one_actual_func = _temp_func_iupacdata
        print(f"Configured to use {conversion_dict_source} (keys are TitleCase).")

    except ImportError:
        print("Could not import or use protein_letters_3to1. Using a hardcoded fallback.")
        _conversion_dict_hardcoded = { # Hardcoded uses UPPERCASE keys
            "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
            "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
            "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
            "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
            "SEC": "U", "PYL": "O", "ASX": "B", "GLX": "Z", "XAA": "X", "UNK": "X",
        }
        conversion_dict = _conversion_dict_hardcoded # Set global conversion_dict for this path
        conversion_dict_source = "Hardcoded Fallback"
        print(f"Using hardcoded fallback dictionary. Type: {type(conversion_dict)}")
        def _temp_func_hardcoded(res_name):
            global conversion_dict # Uses the hardcoded one
            cleaned_res_name = res_name.strip().upper() # Match hardcoded keys
            return conversion_dict[cleaned_res_name]
        _three_to_one_actual_func = _temp_func_hardcoded
        print("Configured to use hardcoded fallback (keys are UPPERCASE).")

# The single, definitive three_to_one function for the script
def three_to_one(res_name):
    if _three_to_one_actual_func is None:
        # This should not be reached if the logic above is correct
        raise RuntimeError("Amino acid conversion function was not properly initialized.")
    return _three_to_one_actual_func(res_name)

def get_sequence_and_res_info(pdb_file_path):
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("PDB_STRUCTURE", pdb_file_path)
    except Exception as e:
        print(f"Error parsing PDB file {pdb_file_path}: {e}")
        return None, None

    residue_info_list = []
    sequence_chars = []
    processed_count = 0
    error_count = 0

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ' and is_aa(residue.get_resname(), standard=True):
                    resname_3letter = residue.get_resname()
                    resnum = residue.get_id()[1]
                    processed_count +=1
                    
                    try:
                        reschar_1letter = three_to_one(resname_3letter)
                        residue_info_list.append((resname_3letter.strip(), resnum, reschar_1letter))
                        sequence_chars.append(reschar_1letter)
                    except KeyError:
                        error_count += 1
                        # Adjust how cleaned_name is generated for the warning based on the source
                        original_cleaned_name_for_warning = ""
                        if conversion_dict_source == "Bio.Data.IUPACData.protein_letters_3to1":
                            original_cleaned_name_for_warning = resname_3letter.strip().capitalize()
                        else: # For Bio.SeqUtils or hardcoded, .upper() was the target
                            original_cleaned_name_for_warning = resname_3letter.strip().upper()
                        
                        is_in_dict_debug = False
                        # Ensure conversion_dict is not None and is a dict before trying to access it
                        if conversion_dict is not None and isinstance(conversion_dict, dict):
                             is_in_dict_debug = original_cleaned_name_for_warning in conversion_dict

                        print(f"Warning: Skipped residue '{original_cleaned_name_for_warning}' (Original: '{resname_3letter}') resnum {resnum} in {pdb_file_path} "
                              f"during 1-letter conversion (KeyError using {conversion_dict_source}). "
                              f"Debug: '{original_cleaned_name_for_warning}' in dict at error point? {is_in_dict_debug}")
                        # if error_count < 5 and conversion_dict is not None and isinstance(conversion_dict, dict) : # More debug for first few errors
                        #     print(f"   First 10 keys of current conversion_dict: {list(conversion_dict.keys())[:10]}")


    print(f"Total residues processed by is_aa for {pdb_file_path}: {processed_count}, KeyErrors: {error_count}")
    if not sequence_chars:
        print(f"Warning: No standard amino acid sequence extracted from {pdb_file_path}.")
        return None, None
        
    return "".join(sequence_chars), residue_info_list

def main():
    parser = argparse.ArgumentParser(
        description="Aligns sequences from two PDB files and generates a residue mapping file."
    )
    parser.add_argument("-r", "--ref", required=True, help="Path to the reference PDB file.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input PDB file.")
    parser.add_argument("-o", "--output", required=True, help="Name for the output residue mapping file.")

    args = parser.parse_args()

    print(f"Processing reference PDB: {args.ref}")
    ref_seq, ref_info_list = get_sequence_and_res_info(args.ref)
    if not ref_seq or not ref_info_list:
        print(f"Could not process reference PDB file: {args.ref}. Exiting.")
        return

    print(f"Processing input PDB: {args.input}")
    input_seq, input_info_list = get_sequence_and_res_info(args.input)
    if not input_seq or not input_info_list:
        print(f"Could not process input PDB file: {args.input}. Exiting.")
        return

    print(f"Reference sequence length: {len(ref_seq)} residues")
    print(f"Input sequence length: {len(input_seq)} residues")

    if not ref_seq or not input_seq:
        print("One or both sequences are empty. Cannot perform alignment. Exiting.")
        return

    print("Performing global sequence alignment...")
    alignments = pairwise2.align.globalms(ref_seq, input_seq, 2, -1, -0.5, -0.1, one_alignment_only=True)

    if not alignments:
        print("Error: Could not generate an alignment. Exiting.")
        return

    alignment = alignments[0]
    # print(pairwise2.format_alignment(*alignment)) # For debugging alignment

    output_lines = ["#RefResName\tRefResNum\tInputResNum\n"]
    ref_list_idx = 0
    input_list_idx = 0
    mapped_pairs_count = 0
    
    aligned_ref_seq, aligned_input_seq, score, begin, end = alignment

    for i in range(len(aligned_ref_seq)):
        ref_char = aligned_ref_seq[i]
        input_char = aligned_input_seq[i]

        if ref_char != '-' and input_char != '-': # Aligned pair
            if ref_list_idx < len(ref_info_list) and input_list_idx < len(input_info_list):
                ref_resname, ref_resnum, _ = ref_info_list[ref_list_idx]
                _, input_resnum, _ = input_info_list[input_list_idx]
                output_lines.append(f"{ref_resname}\t{ref_resnum}\t{input_resnum}\n")
                mapped_pairs_count +=1
            else: # Should not happen if sequence extraction and alignment are correct
                print("Warning: Alignment produced indices out of bounds for residue info lists. Stopping.")
                break 
            ref_list_idx += 1
            input_list_idx += 1
        elif ref_char != '-': # Gap in input_seq
            if ref_list_idx < len(ref_info_list): ref_list_idx += 1
            else: print("Warning: Ran out of reference residues during gap. Stopping."); break
        elif input_char != '-': # Gap in ref_seq
            if input_list_idx < len(input_info_list): input_list_idx += 1
            else: print("Warning: Ran out of input residues during gap. Stopping."); break
            
    try:
        with open(args.output, 'w') as f:
            f.writelines(output_lines)
        print(f"\nSuccessfully wrote {mapped_pairs_count} mapped residue pairs to: {args.output}")
        if mapped_pairs_count == 0 and (ref_seq and input_seq):
            print("Warning: No residue pairs were mapped. Check PDB files and alignment parameters.")
    except IOError as e:
        print(f"Error writing output file {args.output}: {e}")

if __name__ == "__main__":
    main()

