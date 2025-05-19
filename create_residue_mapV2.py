#!/usr/bin/env python3

import argparse
from Bio.PDB import PDBParser, Superimposer, Selection
from Bio.PDB.Polypeptide import PPBuilder
import numpy as np
from pathlib import Path
import sys

def get_ca_atoms(chain):
    """Extracts C-alpha atoms from a chain."""
    return [atom for atom in chain.get_atoms() if atom.name == 'CA']

def get_residues_from_atoms(ca_atoms):
    """Gets unique parent residues from a list of C-alpha atoms."""
    residues = []
    seen_residues = set()
    for atom in ca_atoms:
        residue = atom.get_parent()
        if residue.id not in seen_residues:
            residues.append(residue)
            seen_residues.add(residue.id)
    return residues

def align_chains_and_get_map(ref_chain, input_chain, gap_open_penalty, gap_extend_penalty):
    """
    Aligns two chains based on C-alpha atoms using sequence alignment of their residues
    and returns the residue mapping.
    """
    
    ref_ca_atoms = get_ca_atoms(ref_chain)
    input_ca_atoms = get_ca_atoms(input_chain)

    if not ref_ca_atoms or not input_ca_atoms:
        # print(f"Warning: No C-alpha atoms found in ref chain {ref_chain.id} or input chain {input_chain.id}. Skipping alignment.", file=sys.stderr)
        return [], 0.0, 0

    # Build sequences of 3-letter codes for alignment
    ppb = PPBuilder()
    try:
        ref_seq_pp_list = ppb.build_peptides(ref_chain)
        input_seq_pp_list = ppb.build_peptides(input_chain)
    except Exception as e:
        # print(f"Warning: Could not build peptide for ref chain {ref_chain.id} or input chain {input_chain.id} using PPBuilder. Error: {e}. Skipping.", file=sys.stderr)
        return [], 0.0, 0


    if not ref_seq_pp_list or not input_seq_pp_list:
        # print(f"Warning: Could not build peptide sequence for ref chain {ref_chain.id} or input chain {input_chain.id}. Skipping.", file=sys.stderr)
        return [], 0.0, 0

    # Take the first polypeptide chain if multiple are returned (e.g. due to chain breaks or multiple polypeptides in one PDB 'chain' entity)
    ref_sequence_str = str(ref_seq_pp_list[0].get_sequence())
    input_sequence_str = str(input_seq_pp_list[0].get_sequence())
    
    ref_residues_list = Selection.unfold_entities(ref_seq_pp_list[0], 'R')
    input_residues_list = Selection.unfold_entities(input_seq_pp_list[0], 'R')

    from Bio import pairwise2
    from Bio.pairwise2 import format_alignment

    alignments = pairwise2.align.globalms(ref_sequence_str, input_sequence_str, 2, -1, gap_open_penalty, gap_extend_penalty, one_alignment_only=True)

    if not alignments:
        # print(f"Warning: Could not align sequences for ref chain {ref_chain.id} and input chain {input_chain.id}.", file=sys.stderr)
        return [], 0.0, 0

    alignment = alignments[0]
    aligned_ref_seq, aligned_input_seq, score, begin, end = alignment
    
    aligned_residue_pairs = []
    ref_res_idx = 0
    input_res_idx = 0
    num_aligned_residues = 0

    for i in range(len(aligned_ref_seq)):
        ref_char = aligned_ref_seq[i]
        input_char = aligned_input_seq[i]

        current_ref_res = None
        current_input_res = None

        if ref_char != '-':
            if ref_res_idx < len(ref_residues_list):
                current_ref_res = ref_residues_list[ref_res_idx]
            ref_res_idx += 1
        
        if input_char != '-':
            if input_res_idx < len(input_residues_list):
                current_input_res = input_residues_list[input_res_idx]
            input_res_idx += 1

        if current_ref_res and current_input_res and ref_char != '-' and input_char != '-':
            aligned_residue_pairs.append(
                (current_ref_res, current_input_res)
            )
            num_aligned_residues +=1
            
    return aligned_residue_pairs, score, num_aligned_residues


def main():
    parser = argparse.ArgumentParser(
        description="Create a residue map file (.crmp) between a reference PDB and an input PDB. "
                    "The map explicitly states reference and input chain IDs for each mapped residue pair. "
                    "Handles multi-chain PDBs by attempting to map each reference chain to each input chain "
                    "and selecting the best alignment (highest number of aligned residues, then best score).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--ref', required=True, help="Reference PDB file path (e.g., experimental structure like 1R1K).")
    parser.add_argument('--input', required=True, help="Input PDB file path (e.g., simulation starting structure, AlphaFold model).")
    parser.add_argument('--output', required=True, help="Output residue map file path (.crmp).")
    parser.add_argument('--ref_chains', nargs='+', help="Specific chain IDs from the reference PDB to process (e.g., A D). If None, all chains are processed.")
    parser.add_argument('--input_chains', nargs='+', help="Specific chain IDs from the input PDB to consider for mapping (e.g., A B). If None, all chains are considered.")
    parser.add_argument('--align_gap_open', type=float, default=-10, help="Gap open penalty for sequence alignment.")
    parser.add_argument('--align_gap_extend', type=float, default=-0.5, help="Gap extend penalty for sequence alignment.")
    parser.add_argument('--default_input_chain_id', type=str, default="X", 
                        help="Default chain ID to use for input PDB chains if their ID is blank/empty. Must be a single character.")

    args = parser.parse_args()

    if not Path(args.ref).exists():
        print(f"Error: Reference PDB file not found: {args.ref}", file=sys.stderr)
        sys.exit(1)
    if not Path(args.input).exists():
        print(f"Error: Input PDB file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    if len(args.default_input_chain_id) != 1:
        print(f"Error: --default_input_chain_id ('{args.default_input_chain_id}') must be a single character.", file=sys.stderr)
        sys.exit(1)

    pdb_parser = PDBParser(QUIET=True)
    try:
        ref_structure = pdb_parser.get_structure("reference", args.ref)
        input_structure = pdb_parser.get_structure("input", args.input)
    except Exception as e:
        print(f"Error parsing PDB files: {e}", file=sys.stderr)
        sys.exit(1)


    output_map_file = Path(args.output)
    output_map_file.parent.mkdir(parents=True, exist_ok=True)

    processed_ref_chains = set()

    with open(output_map_file, 'w') as f:
        f.write("# RefChainID\tRefResName\tRefResNum\tInputChainID\tInputResName\tInputResNum\n") # Header
        
        ref_chains_to_process = []
        if args.ref_chains:
            for chain_id in args.ref_chains:
                if chain_id in ref_structure[0]:
                    ref_chains_to_process.append(ref_structure[0][chain_id])
                else:
                    print(f"Warning: Reference chain {chain_id} not found in {args.ref}. Skipping.", file=sys.stderr)
        else:
            ref_chains_to_process = list(ref_structure[0]) 

        input_chains_to_consider = []
        if args.input_chains:
            for chain_id in args.input_chains:
                if chain_id in input_structure[0]:
                    input_chains_to_consider.append(input_structure[0][chain_id])
                else:
                    print(f"Warning: Input chain {chain_id} not found in {args.input}. Skipping.", file=sys.stderr)
        else:
            input_chains_to_consider = list(input_structure[0])


        for ref_chain in ref_chains_to_process:
            ref_chain_id_for_log = ref_chain.id.strip()
            if not ref_chain_id_for_log: 
                print(f"Warning: Reference chain with blank ID found in {args.ref}. Skipping.", file=sys.stderr)
                continue
            
            print(f"\nProcessing reference chain: {ref_chain_id_for_log}")
            best_input_chain_obj = None
            best_alignment_map = []
            best_num_aligned = -1  
            best_score = -float('inf') 

            for input_chain_obj in input_chains_to_consider:
                # Determine the input chain ID for logging, using default if blank
                current_input_chain_id_for_log = input_chain_obj.id.strip()
                if not current_input_chain_id_for_log:
                    current_input_chain_id_for_log = f"'{args.default_input_chain_id}' (defaulted from blank)"
                else:
                    current_input_chain_id_for_log = f"'{current_input_chain_id_for_log}'"

                print(f"  Aligning ref chain '{ref_chain_id_for_log}' with input chain {current_input_chain_id_for_log}...")
                
                aligned_pairs, score, num_aligned = align_chains_and_get_map(
                    ref_chain, input_chain_obj, args.align_gap_open, args.align_gap_extend
                )
                print(f"    -> Found {num_aligned} aligned residues with score {score:.2f}.")

                if num_aligned > best_num_aligned:
                    best_num_aligned = num_aligned
                    best_score = score
                    best_alignment_map = aligned_pairs
                    best_input_chain_obj = input_chain_obj
                elif num_aligned == best_num_aligned and score > best_score: 
                    best_score = score
                    best_alignment_map = aligned_pairs
                    best_input_chain_obj = input_chain_obj

            if best_input_chain_obj and best_alignment_map:
                ref_chain_id_to_write = ref_chain_id_for_log # Already stripped

                # --- KEY FIX IS HERE ---
                # Use default if input chain ID is blank, otherwise use its actual ID (stripped)
                input_chain_id_to_write = best_input_chain_obj.id.strip()
                if not input_chain_id_to_write:
                    input_chain_id_to_write = args.default_input_chain_id
                # --- END OF KEY FIX ---

                final_input_chain_id_for_log = input_chain_id_to_write
                if best_input_chain_obj.id.strip() == "":
                     final_input_chain_id_for_log = f"'{input_chain_id_to_write}' (defaulted from blank)"
                else:
                     final_input_chain_id_for_log = f"'{input_chain_id_to_write}'"


                print(f"  Best match for ref chain '{ref_chain_id_to_write}': input chain {final_input_chain_id_for_log} ({best_num_aligned} residues, score {best_score:.2f})")
                
                for ref_res, input_res in best_alignment_map:
                    ref_res_num = ref_res.id[1]
                    input_res_num = input_res.id[1]
                    ref_res_name = ref_res.get_resname().strip()
                    input_res_name = input_res.get_resname().strip()
                    
                    f.write(f"{ref_chain_id_to_write}\t{ref_res_name}\t{ref_res_num}\t"
                            f"{input_chain_id_to_write}\t{input_res_name}\t{input_res_num}\n")
                processed_ref_chains.add(ref_chain.id)
            else:
                print(f"  Warning: No suitable alignment found for reference chain '{ref_chain_id_for_log}' with any input chains.", file=sys.stderr)
    
    print(f"\nResidue map written to: {output_map_file}")
    if not processed_ref_chains:
        print("Warning: No reference chains were successfully processed or mapped.", file=sys.stderr)

if __name__ == '__main__':
    main()

