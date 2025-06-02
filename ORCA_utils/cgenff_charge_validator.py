import re
import argparse
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple, Any
import csv
import sys # For exiting early if critical info for offset is missing

@dataclass
class CGenFFAtom:
    name: str
    atom_type: str
    charge: float
    penalty: float
    element: str
    original_index: int # Its 0-based index in the STR file

@dataclass
class OrcaCharge:
    index: int # 0-based index from ORCA output
    element: str
    charge: float

def get_element_from_cgenff_atom_name(atom_name: str) -> str:
    """Extracts element symbol from CGenFF atom name (e.g., C1->C, H01->H)."""
    match = re.match(r"([A-Za-z]+)", atom_name)
    if match:
        if len(match.group(1)) == 2 and match.group(1).upper() in ["CL", "BR", "NA", "MG", "ZN", "FE", "CU", "MN"]:
            return match.group(1).upper()
        element_symbol = match.group(1)[0].upper()
        if element_symbol in ['H', 'C', 'N', 'O', 'S', 'P', 'F', 'I', 'K', 'B']:
            return element_symbol
    return atom_name[0].upper()

def parse_cgenff_str(str_filepath: str) -> List[CGenFFAtom]:
    """Parses ATOM lines from a CGenFF .str file."""
    atoms: List[CGenFFAtom] = []
    atom_pattern = re.compile(r"^ATOM\s+(\S+)\s+(\S+)\s+([-\d.]+)\s+!\s+([-\d.]+).*")
    try:
        with open(str_filepath, 'r') as f:
            for i, line in enumerate(f):
                match = atom_pattern.match(line)
                if match:
                    name, atom_type, charge_str, penalty_str = match.groups()
                    element = get_element_from_cgenff_atom_name(name)
                    atoms.append(CGenFFAtom(
                        name=name,
                        atom_type=atom_type,
                        charge=float(charge_str),
                        penalty=float(penalty_str),
                        element=element,
                        original_index=len(atoms)
                    ))
    except FileNotFoundError:
        print(f"Error: CGenFF file not found at {str_filepath}")
        return []
    return atoms

def parse_orca_out(orca_filepath: str) -> Dict[str, List[OrcaCharge]]:
    """Parses Loewdin and Mulliken charges from an ORCA .out file."""
    charges: Dict[str, List[OrcaCharge]] = {"loewdin": [], "mulliken": []}
    
    loewdin_pattern = re.compile(r"^\s*(\d+)\s+([A-Za-z]+)\s*:\s*([-\d.]+)")
    mulliken_pattern = re.compile(r"^\s*(\d+)\s+([A-Za-z]+)\s+\S+\s+\S+\s+([-\d.]+).*")

    try:
        with open(orca_filepath, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: ORCA output file not found at {orca_filepath}")
        return charges

    current_section = None 
    mayer_header_skipped = False

    for line in lines:
        if "LOEWDIN ATOMIC CHARGES" in line:
            current_section = "loewdin"
            charges["loewdin"] = [] 
            continue
        elif "MAYER POPULATION ANALYSIS" in line:
            current_section = "mulliken"
            charges["mulliken"] = [] 
            mayer_header_skipped = False 
            continue
        
        if current_section == "loewdin":
            if "LOEWDIN REDUCED ORBITAL CHARGES" in line:
                current_section = None
                continue
        elif current_section == "mulliken":
            if mayer_header_skipped: 
                if "Sum of atomic valences" in line or \
                   "Mayer bond orders" in line:
                    current_section = None
                    continue
        
        if current_section and not line.strip(): 
            if (current_section == "loewdin" and charges["loewdin"]) or \
               (current_section == "mulliken" and charges["mulliken"] and mayer_header_skipped): 
                current_section = None
                continue

        if current_section == "loewdin":
            if "----" in line.strip() and not charges["loewdin"]:
                continue
            match = loewdin_pattern.match(line)
            if match:
                idx, element, charge_str = match.groups()
                charges["loewdin"].append(OrcaCharge(
                    index=int(idx),
                    element=element.strip(), 
                    charge=float(charge_str)
                ))
        
        elif current_section == "mulliken":
            if not mayer_header_skipped:
                if "ATOM       NA         ZA         QA" in line:
                    mayer_header_skipped = True
                continue 
            
            match = mulliken_pattern.match(line)
            if match:
                idx, element, qa_charge_str = match.groups()
                charges["mulliken"].append(OrcaCharge(
                    index=int(idx), 
                    element=element.strip(), 
                    charge=float(qa_charge_str)
                ))
    return charges

def get_orca_atom_index_offset_from_mol2(mol2_filepath: str, target_residue_name: str) -> Optional[int]:
    """
    Parses a MOL2 file to find the minimum atom ID for a target residue name.
    Returns the 0-based offset (min_atom_id - 1) or None if residue not found.
    """
    min_atom_id_for_residue: Optional[int] = None
    in_atom_section = False
    found_residue_atoms = False

    try:
        with open(mol2_filepath, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.strip() == "@<TRIPOS>ATOM":
                    in_atom_section = True
                    continue
                if line.strip() == "@<TRIPOS>BOND" or \
                   (in_atom_section and line.startswith("@<TRIPOS>")): # End of ATOM section or another TRIPOS section starts
                    in_atom_section = False
                    break 
                
                if in_atom_section:
                    parts = line.split()
                    # MOL2: atom_id atom_name x y z atom_type [subst_id subst_name [charge]]
                    if len(parts) >= 8: # Need at least up to subst_name
                        try:
                            atom_id = int(parts[0])
                            res_name_in_mol2 = parts[7] # subst_name
                            
                            if res_name_in_mol2.strip().upper() == target_residue_name.strip().upper():
                                found_residue_atoms = True
                                if min_atom_id_for_residue is None or atom_id < min_atom_id_for_residue:
                                    min_atom_id_for_residue = atom_id
                        except ValueError:
                            print(f"Warning: Could not parse atom ID on line {line_num} in '{mol2_filepath}': {line.strip()}")
                            continue
    except FileNotFoundError:
        print(f"Error: MOL2 file '{mol2_filepath}' not found.")
        return None # Indicate critical failure
    
    if not found_residue_atoms:
        print(f"Error: Target residue '{target_residue_name}' not found in '{mol2_filepath}'. Cannot determine ORCA atom index offset.")
        return None # Indicate critical failure
        
    if min_atom_id_for_residue is not None:
        offset = min_atom_id_for_residue - 1
        print(f"Found target residue '{target_residue_name}' in '{mol2_filepath}'. Min atom ID: {min_atom_id_for_residue}. Calculated ORCA atom index offset: {offset}")
        return offset
    else:
        # This case should ideally be caught by found_residue_atoms, but as a fallback:
        print(f"Error: Target residue '{target_residue_name}' was matched but no valid atom IDs found in '{mol2_filepath}'.")
        return None

def generate_comparison_data(
    cgenff_atoms: List[CGenFFAtom],
    orca_data: Dict[str, List[OrcaCharge]],
    penalty_threshold: float,
    orca_atom_index_offset: int = 0 # Default to 0 if not provided
) -> Tuple[List[str], List[Dict[str, Any]]]:
    """Generates headers and data rows for the comparison report."""
    
    headers = [
        "Atom Name", "CGenFF Type", "CGenFF Q", "Penalty", "ORCA Calc Idx",
        "Loewdin Q", "Mulliken Q", "Mean QM Q",
        "Abs_Diff", "Pct_Diff", "Sym_Std_Diff"
    ]
    data_rows: List[Dict[str, Any]] = []

    loewdin_map = {oc.index: oc for oc in orca_data.get('loewdin', [])}
    mulliken_map = {oc.index: oc for oc in orca_data.get('mulliken', [])}

    for atom in cgenff_atoms:
        if atom.penalty > penalty_threshold:
            loewdin_q_val: Optional[float] = None
            mulliken_q_val: Optional[float] = None
            mean_qm_q: Optional[float] = None
            abs_diff: Optional[float] = None
            pct_diff: Optional[float] = None
            sym_std_diff: Optional[float] = None

            effective_orca_index = atom.original_index + orca_atom_index_offset

            if effective_orca_index in loewdin_map:
                orca_loewdin_atom = loewdin_map[effective_orca_index]
                if orca_loewdin_atom.element.upper() == atom.element.upper():
                    loewdin_q_val = orca_loewdin_atom.charge
                else:
                    print(f"Warning: Element mismatch for Loewdin. CGenFF: {atom.name}({atom.element}, STR idx {atom.original_index}), Target ORCA: ({orca_loewdin_atom.element}, ORCA idx {effective_orca_index})")
            else:
                print(f"Warning: No Loewdin charge found for CGenFF atom {atom.name} (STR idx {atom.original_index}) at calculated ORCA index {effective_orca_index}.")


            if effective_orca_index in mulliken_map:
                orca_mulliken_atom = mulliken_map[effective_orca_index]
                if orca_mulliken_atom.element.upper() == atom.element.upper():
                    mulliken_q_val = orca_mulliken_atom.charge
                else:
                    print(f"Warning: Element mismatch for Mulliken. CGenFF: {atom.name}({atom.element}, STR idx {atom.original_index}), Target ORCA: ({orca_mulliken_atom.element}, ORCA idx {effective_orca_index})")
            else:
                print(f"Warning: No Mulliken charge found for CGenFF atom {atom.name} (STR idx {atom.original_index}) at calculated ORCA index {effective_orca_index}.")

            if loewdin_q_val is not None and mulliken_q_val is not None:
                mean_qm_q = (loewdin_q_val + mulliken_q_val) / 2.0
                abs_diff = atom.charge - mean_qm_q
                
                if abs(mean_qm_q) < 1e-6: # Using a small epsilon for near-zero check
                    pct_diff = None 
                else:
                    pct_diff = ((atom.charge - mean_qm_q) / mean_qm_q) * 100.0
                
                denominator_sym_std = (abs(atom.charge) + abs(mean_qm_q))
                if denominator_sym_std < 1e-9: # Avoid division by zero or very small numbers
                    sym_std_diff = None
                else:
                    sym_std_diff = ((atom.charge - mean_qm_q) / denominator_sym_std) * 100.0
            
            row_data: Dict[str, Any] = {
                "Atom Name": atom.name,
                "CGenFF Type": atom.atom_type,
                "CGenFF Q": atom.charge,
                "Penalty": atom.penalty,
                "ORCA Calc Idx": effective_orca_index, # Store the calculated ORCA index
                "Loewdin Q": loewdin_q_val,
                "Mulliken Q": mulliken_q_val,
                "Mean QM Q": mean_qm_q,
                "Abs_Diff": abs_diff,
                "Pct_Diff": pct_diff,
                "Sym_Std_Diff": sym_std_diff,
            }
            data_rows.append(row_data)
            
    return headers, data_rows

def format_report_for_console(headers: List[str], data_rows: List[Dict[str, Any]]) -> str:
    """Formats the report data for console printing."""
    report_lines: List[str] = []
    
    # Adjust column widths, adding ORCA Calc Idx
    # Headers: Atom Name, CGenFF Type, CGenFF Q, Penalty, ORCA Calc Idx, Loewdin Q, Mulliken Q, Mean QM Q, Abs_Diff, Pct_Diff, Sym_Std_Diff
    # Index:      0            1           2         3            4             5           6           7          8         9            10
    header_format = (
        f"{headers[0]:<10} | {headers[1]:<12} | {headers[2]:>10} | {headers[3]:>10} | "
        f"{headers[4]:>13} | {headers[5]:>10} | {headers[6]:>10} | {headers[7]:>10} | "
        f"{headers[8]:>10} | {headers[9]:>10} | {headers[10]:>12}"
    )
    report_lines.append(header_format)
    line_length = len(header_format)
    report_lines.append("-" * line_length)

    if not data_rows:
        report_lines.append("No atoms found with penalty above the threshold (or matching other criteria).")
    else:
        for row in data_rows:
            loewdin_q_str = f"{row['Loewdin Q']:.6f}" if row['Loewdin Q'] is not None else "N/A"
            mulliken_q_str = f"{row['Mulliken Q']:.6f}" if row['Mulliken Q'] is not None else "N/A"
            mean_qm_q_str = f"{row['Mean QM Q']:.6f}" if row['Mean QM Q'] is not None else "N/A"
            abs_diff_str = f"{row['Abs_Diff']:.6f}" if row['Abs_Diff'] is not None else "N/A"
            pct_diff_str = f"{row['Pct_Diff']:.2f}%" if row['Pct_Diff'] is not None else "N/A"
            sym_std_diff_str = f"{row['Sym_Std_Diff']:.2f}%" if row['Sym_Std_Diff'] is not None else "N/A"
            
            report_lines.append(
                f"{row['Atom Name']:<10} | {row['CGenFF Type']:<12} | {row['CGenFF Q']:>10.6f} | {row['Penalty']:>10.3f} | "
                f"{str(row['ORCA Calc Idx']):>13} | {loewdin_q_str:>10} | {mulliken_q_str:>10} | {mean_qm_q_str:>10} | "
                f"{abs_diff_str:>10} | {pct_diff_str:>10} | {sym_std_diff_str:>12}"
            )
    return "\n".join(report_lines)

def write_report_to_csv(filepath: str, headers: List[str], data_rows: List[Dict[str, Any]]):
    """Writes the report data to a CSV file."""
    try:
        with open(filepath, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=headers)
            writer.writeheader()
            for row in data_rows:
                csv_row = {}
                for header in headers:
                    value = row[header]
                    if value is None:
                        csv_row[header] = "N/A"
                    elif isinstance(value, float):
                        if header.endswith(" Q") or header == "Abs_Diff":
                            csv_row[header] = f"{value:.6f}"
                        elif header == "Penalty":
                            csv_row[header] = f"{value:.3f}"
                        elif header == "Pct_Diff" or header == "Sym_Std_Diff":
                             csv_row[header] = f"{value:.2f}%"
                        else:
                            csv_row[header] = value 
                    else: # Handles ORCA Calc Idx (int) and other string types
                        csv_row[header] = value
                writer.writerow(csv_row)
        print(f"Report successfully written to {filepath}")
    except IOError:
        print(f"Error: Could not write CSV file to {filepath}")

def main():
    parser = argparse.ArgumentParser(
        description="Compare CGenFF charges with ORCA Loewdin and Mulliken charges for atoms with high penalty. "
                    "Can apply an atom index offset for ORCA output from cluster calculations using a MOL2 file."
    )
    parser.add_argument("str_file", help="Path to the CGenFF .str file (for the molecule of interest).")
    parser.add_argument("orca_out_file", help="Path to the ORCA .out file (can be for a larger cluster).")
    parser.add_argument(
        "--penalty_threshold",
        type=float,
        default=10.0,
        help="CH_PENALTY threshold for reporting (default: 10.0)."
    )
    parser.add_argument(
        "--output_csv",
        help="Optional path to save the report as a CSV file."
    )
    parser.add_argument(
        "--cluster_mol2_file",
        help="Optional path to the MOL2 file of the full cluster system used for the ORCA calculation. "
             "Required if --target_residue_name is used."
    )
    parser.add_argument(
        "--target_residue_name",
        help="Residue name of the molecule of interest within the --cluster_mol2_file "
             "to determine the ORCA atom index offset."
    )
    args = parser.parse_args()

    if args.target_residue_name and not args.cluster_mol2_file:
        parser.error("--target_residue_name requires --cluster_mol2_file to be specified.")
        sys.exit(1)

    cgenff_atoms = parse_cgenff_str(args.str_file)
    if not cgenff_atoms:
        print("No CGenFF atoms parsed from .str file. Exiting.")
        return

    orca_charges = parse_orca_out(args.orca_out_file)
    if not orca_charges.get('loewdin') and not orca_charges.get('mulliken'):
        print("No Loewdin or Mulliken charges parsed from ORCA .out file. Exiting.")
        return
        
    orca_atom_index_offset = 0
    if args.cluster_mol2_file and args.target_residue_name:
        print(f"Attempting to determine ORCA atom index offset using MOL2 file: '{args.cluster_mol2_file}' for residue: '{args.target_residue_name}'")
        offset = get_orca_atom_index_offset_from_mol2(args.cluster_mol2_file, args.target_residue_name)
        if offset is None:
            # Error messages are printed within get_orca_atom_index_offset_from_mol2
            print("Exiting due to failure in determining ORCA atom index offset.")
            sys.exit(1)
        orca_atom_index_offset = offset
    elif args.target_residue_name and not args.cluster_mol2_file:
        # This case is caught by parser.error, but as a safeguard:
        print("Error: --target_residue_name was provided without --cluster_mol2_file. Exiting.")
        sys.exit(1)


    print(f"CGenFF atoms parsed: {len(cgenff_atoms)}")
    print(f"ORCA Loewdin charges parsed: {len(orca_charges.get('loewdin',[]))}")
    print(f"ORCA Mulliken charges parsed: {len(orca_charges.get('mulliken',[]))}")
    if args.cluster_mol2_file and args.target_residue_name:
        print(f"Using ORCA atom index offset: {orca_atom_index_offset}")

    headers, data_rows = generate_comparison_data(
        cgenff_atoms, 
        orca_charges, 
        args.penalty_threshold,
        orca_atom_index_offset # Pass the calculated offset
    )
    
    console_report = format_report_for_console(headers, data_rows)
    print("\n--- CHARGE COMPARISON REPORT ---")
    print(console_report)

    if args.output_csv:
        write_report_to_csv(args.output_csv, headers, data_rows)

if __name__ == "__main__":
    main()
