import re
import argparse
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple, Any
import csv # Added for CSV output

@dataclass
class CGenFFAtom:
    name: str
    atom_type: str
    charge: float
    penalty: float
    element: str
    original_index: int # Its index in the full list of atoms from STR file

@dataclass
class OrcaCharge:
    index: int # 0-based index from ORCA output
    element: str
    charge: float

def get_element_from_cgenff_atom_name(atom_name: str) -> str:
    """Extracts element symbol from CGenFF atom name (e.g., C1->C, H01->H)."""
    match = re.match(r"([A-Za-z]+)", atom_name)
    if match:
        # Prioritize two-letter symbols if the full matched string is a known two-letter element
        if len(match.group(1)) == 2 and match.group(1).upper() in ["CL", "BR", "NA", "MG", "ZN", "FE", "CU", "MN"]: # Add more as needed
            return match.group(1).upper()
        # Fallback to first letter for single-letter elements or others
        element_symbol = match.group(1)[0].upper()
        if element_symbol in ['H', 'C', 'N', 'O', 'S', 'P', 'F', 'I', 'K', 'B']: # Added I, K, B
            return element_symbol
    return atom_name[0].upper() # Fallback to first letter if regex fails or no specific match

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
                        original_index=len(atoms) # 0-based index as parsed
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
        # Detect section starts
        if "LOEWDIN ATOMIC CHARGES" in line:
            current_section = "loewdin"
            charges["loewdin"] = [] 
            continue
        elif "MAYER POPULATION ANALYSIS" in line:
            current_section = "mulliken"
            charges["mulliken"] = [] 
            mayer_header_skipped = False # Reset for the start of this specific section
            continue
        
        # Detect section ends by specific keywords or change of major section
        if current_section == "loewdin":
            if "LOEWDIN REDUCED ORBITAL CHARGES" in line:
                current_section = None
                continue
        elif current_section == "mulliken":
            # Termination conditions for Mulliken section should only apply AFTER its header is parsed.
            if mayer_header_skipped: # <<< KEY CHANGE: Only check if header was already skipped
                if "Sum of atomic valences" in line or \
                   "Mayer bond orders" in line:
                    current_section = None
                    # mayer_header_skipped = False; # Not strictly needed to reset, as section is ending
                    continue
        
        # If a blank line is encountered after some data has been collected for the current section,
        # it can also signify the end of the data block.
        if current_section and not line.strip(): # Blank line check
            if (current_section == "loewdin" and charges["loewdin"]) or \
               (current_section == "mulliken" and charges["mulliken"] and mayer_header_skipped): # For Mulliken, also ensure header was passed
                current_section = None
                continue

        # Process lines based on current section
        if current_section == "loewdin":
            # Skip the "----" separator line often found after the section header, if no data yet
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
            # Skip known header lines for Mulliken section before attempting to parse data
            if not mayer_header_skipped:
                if "ATOM       NA         ZA         QA" in line:
                    mayer_header_skipped = True
                # This continue is vital: it skips the current line (e.g. '****' border, description, or 'ATOM NA...')
                # from falling into the data parsing logic below in this same iteration.
                # It also prevents such lines from being evaluated by the termination logic above if mayer_header_skipped is false.
                continue 
            
            # After header is skipped (mayer_header_skipped is True), and if not a termination line, try to match data.
            match = mulliken_pattern.match(line)
            if match:
                idx, element, qa_charge_str = match.groups()
                charges["mulliken"].append(OrcaCharge(
                    index=int(idx), 
                    element=element.strip(), 
                    charge=float(qa_charge_str)
                ))
                
    return charges

def generate_comparison_data(
    cgenff_atoms: List[CGenFFAtom],
    orca_data: Dict[str, List[OrcaCharge]],
    penalty_threshold: float
) -> Tuple[List[str], List[Dict[str, Any]]]:
    """Generates headers and data rows for the comparison report."""
    
    headers = [
        "Atom Name", "CGenFF Type", "CGenFF Q", "Penalty", 
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

            orca_atom_idx_to_match = atom.original_index 

            if orca_atom_idx_to_match in loewdin_map:
                orca_loewdin_atom = loewdin_map[orca_atom_idx_to_match]
                if orca_loewdin_atom.element.upper() == atom.element.upper():
                    loewdin_q_val = orca_loewdin_atom.charge
                else:
                    print(f"Warning: Element mismatch for Loewdin. CGenFF: {atom.name}({atom.element}, idx {atom.original_index}), ORCA: ({orca_loewdin_atom.element}, idx {orca_atom_idx_to_match})")
            
            if orca_atom_idx_to_match in mulliken_map:
                orca_mulliken_atom = mulliken_map[orca_atom_idx_to_match]
                if orca_mulliken_atom.element.upper() == atom.element.upper():
                    mulliken_q_val = orca_mulliken_atom.charge
                else:
                    print(f"Warning: Element mismatch for Mulliken. CGenFF: {atom.name}({atom.element}, idx {atom.original_index}), ORCA: ({orca_mulliken_atom.element}, idx {orca_atom_idx_to_match})")

            # Calculate new metrics if QM charges are available
            if loewdin_q_val is not None and mulliken_q_val is not None:
                mean_qm_q = (loewdin_q_val + mulliken_q_val) / 2.0
                
                abs_diff = atom.charge - mean_qm_q
                
                if abs(mean_qm_q) < 0.01:
                    pct_diff = None # Will be formatted as N/A
                else:
                    pct_diff = ((atom.charge - mean_qm_q) / mean_qm_q) * 100.0
                
                denominator_sym_std = (abs(atom.charge) + abs(mean_qm_q)) + 1e-9
                sym_std_diff = ((atom.charge - mean_qm_q) / denominator_sym_std) * 100.0
            
            row_data: Dict[str, Any] = {
                "Atom Name": atom.name,
                "CGenFF Type": atom.atom_type,
                "CGenFF Q": atom.charge,
                "Penalty": atom.penalty,
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
    
    # Adjust column widths as needed
    header_format = (
        f"{headers[0]:<10} | {headers[1]:<12} | {headers[2]:>10} | {headers[3]:>10} | "
        f"{headers[4]:>10} | {headers[5]:>10} | {headers[6]:>10} | " # Mean QM Q
        f"{headers[7]:>10} | {headers[8]:>10} | {headers[9]:>12}"  # Abs_Diff, Pct_Diff, Sym_Std_Diff
    )
    report_lines.append(header_format)
    line_length = len(header_format)
    report_lines.append("-" * line_length)

    if not data_rows:
        report_lines.append("No atoms found with penalty above the threshold.")
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
                f"{loewdin_q_str:>10} | {mulliken_q_str:>10} | {mean_qm_q_str:>10} | "
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
                # Prepare row for DictWriter, handling None for N/A and formatting floats
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
                             csv_row[header] = f"{value:.2f}%" # Adding % for CSV as well
                        else:
                            csv_row[header] = value 
                    else:
                        csv_row[header] = value
                writer.writerow(csv_row)
        print(f"Report successfully written to {filepath}")
    except IOError:
        print(f"Error: Could not write CSV file to {filepath}")

def main():
    parser = argparse.ArgumentParser(
        description="Compare CGenFF charges with ORCA Loewdin and Mulliken charges for atoms with high penalty."
    )
    parser.add_argument("str_file", help="Path to the CGenFF .str file.")
    parser.add_argument("orca_out_file", help="Path to the ORCA .out file.")
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
    args = parser.parse_args()

    cgenff_atoms = parse_cgenff_str(args.str_file)
    if not cgenff_atoms:
        return

    orca_charges = parse_orca_out(args.orca_out_file)
    
    print(f"CGenFF atoms parsed: {len(cgenff_atoms)}")
    print(f"ORCA Loewdin charges parsed: {len(orca_charges.get('loewdin',[]))}")
    print(f"ORCA Mulliken charges parsed: {len(orca_charges.get('mulliken',[]))}")

    headers, data_rows = generate_comparison_data(cgenff_atoms, orca_charges, args.penalty_threshold)
    
    console_report = format_report_for_console(headers, data_rows)
    print(console_report)

    if args.output_csv:
        write_report_to_csv(args.output_csv, headers, data_rows)

if __name__ == "__main__":
    main()