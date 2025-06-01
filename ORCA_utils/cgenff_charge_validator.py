import re
import argparse
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple

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
        element_symbol = match.group(1)[0].upper()
        # Basic check for common elements, can be expanded
        if element_symbol in ['H', 'C', 'N', 'O', 'S', 'P', 'F']:
            return element_symbol
        # For two-letter symbols like Cl, Br, if atom names are e.g. CL1, BR1
        if len(match.group(1)) > 1 and match.group(1)[:2].upper() in ["CL", "BR"]:
            return match.group(1)[:2].upper()
    return atom_name[0].upper() # Fallback to first letter

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

def generate_comparison_report(
    cgenff_atoms: List[CGenFFAtom],
    orca_data: Dict[str, List[OrcaCharge]],
    penalty_threshold: float
) -> str:
    """Generates a formatted string table comparing charges."""
    report_lines: List[str] = []
    
    headers = ["Atom Name", "CGenFF Type", "CGenFF Q", "Penalty", "Loewdin Q", "Mulliken Q"]
    report_lines.append(f"{headers[0]:<10} | {headers[1]:<12} | {headers[2]:>12} | {headers[3]:>10} | {headers[4]:>12} | {headers[5]:>12}")
    report_lines.append("-" * (10 + 3 + 12 + 3 + 12 + 3 + 10 + 3 + 12 + 3 + 12))

    # Create dictionaries for quick lookup of ORCA charges by index
    loewdin_map = {oc.index: oc for oc in orca_data.get('loewdin', [])}
    mulliken_map = {oc.index: oc for oc in orca_data.get('mulliken', [])}

    atoms_reported_count = 0
    for atom in cgenff_atoms:
        if atom.penalty > penalty_threshold:
            atoms_reported_count += 1
            loewdin_q_str = "N/A"
            mulliken_q_str = "N/A"

            # CGenFF original_index is 0-based, ORCA index from output is also 0-based.
            orca_atom_idx_to_match = atom.original_index 

            # Match Loewdin charge
            if orca_atom_idx_to_match in loewdin_map:
                orca_loewdin_atom = loewdin_map[orca_atom_idx_to_match]
                # Element check for safety, case-insensitive
                if orca_loewdin_atom.element.upper() == atom.element.upper():
                    loewdin_q_str = f"{orca_loewdin_atom.charge:.6f}"
                else:
                    print(f"Warning: Element mismatch for Loewdin. CGenFF: {atom.name}({atom.element}, idx {atom.original_index}), ORCA: ({orca_loewdin_atom.element}, idx {orca_atom_idx_to_match})")
            
            # Match Mulliken charge
            if orca_atom_idx_to_match in mulliken_map:
                orca_mulliken_atom = mulliken_map[orca_atom_idx_to_match]
                # Element check for safety, case-insensitive
                if orca_mulliken_atom.element.upper() == atom.element.upper():
                    mulliken_q_str = f"{orca_mulliken_atom.charge:.6f}"
                else:
                    print(f"Warning: Element mismatch for Mulliken. CGenFF: {atom.name}({atom.element}, idx {atom.original_index}), ORCA: ({orca_mulliken_atom.element}, idx {orca_atom_idx_to_match})")

            report_lines.append(
                f"{atom.name:<10} | {atom.atom_type:<12} | {atom.charge:>12.6f} | {atom.penalty:>10.3f} | "
                f"{loewdin_q_str:>12} | {mulliken_q_str:>12}"
            )
            
    if atoms_reported_count == 0:
        report_lines.append("No atoms found with penalty above the threshold.")
        
    return "\n".join(report_lines)

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
    args = parser.parse_args()

    cgenff_atoms = parse_cgenff_str(args.str_file)
    if not cgenff_atoms:
        return

    orca_charges = parse_orca_out(args.orca_out_file)
    
    print(f"CGenFF atoms parsed: {len(cgenff_atoms)}")
    print(f"ORCA Loewdin charges parsed: {len(orca_charges.get('loewdin',[]))}")
    print(f"ORCA Mulliken charges parsed: {len(orca_charges.get('mulliken',[]))}")


    report = generate_comparison_report(cgenff_atoms, orca_charges, args.penalty_threshold)
    print(report)

if __name__ == "__main__":
    main()