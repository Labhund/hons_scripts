# ORCA_utils

This directory contains utilities for working with molecular files and charge validation in computational chemistry workflows, especially for ORCA and CGenFF.

## Scripts

### 1. `mol2_to_xyz.py`

**Purpose:**  
Convert a MOL2 file to an XYZ file for use in ORCA calculations.

**Usage:**
```bash
python mol2_to_xyz.py [MOL2_FILE] [-i MOL2_FILE] [-o XYZ_FILE | --deffnm] [--run_example]
```

**Arguments:**
- `MOL2_FILE` (positional): Path to the input MOL2 file.
- `-i`, `--input`: Path to the input MOL2 file (overrides positional).
- `-o`, `--output`: Path to the output XYZ file.
- `--deffnm`: Automatically set output filename to match input (e.g., `input.mol2` → `input.xyz`).
- `--run_example`: Run an internal example conversion.

**Example:**
```bash
python mol2_to_xyz.py ligand.mol2 -o ligand.xyz
python mol2_to_xyz.py -i ligand.mol2 --deffnm
python mol2_to_xyz.py --run_example
```

---

### 2. `cgenff_charge_validator.py`

**Purpose:**  
Compare CGenFF charges (from a `.str` file) with ORCA Loewdin and Mulliken charges (from an ORCA `.out` file) for atoms with high penalty scores. Supports atom index offsetting for cluster calculations using a MOL2 file.

**Usage:**
```bash
python cgenff_charge_validator.py STR_FILE ORCA_OUT_FILE [--penalty_threshold THRESH] [--output_csv CSV_FILE] [--cluster_mol2_file MOL2_FILE --target_residue_name RESNAME]
```

**Arguments:**
- `STR_FILE`: Path to the CGenFF `.str` file.
- `ORCA_OUT_FILE`: Path to the ORCA `.out` file.
- `--penalty_threshold`: Only report atoms with penalty above this value (default: 10.0).
- `--output_csv`: Write the report to a CSV file.
- `--cluster_mol2_file`: MOL2 file for the full cluster (needed for offset calculation).
- `--target_residue_name`: Residue name in the MOL2 file to determine atom index offset.

**Example:**
```bash
python cgenff_charge_validator.py ligand.str orca.out --penalty_threshold 15 --output_csv report.csv
python cgenff_charge_validator.py ligand.str orca_cluster.out --cluster_mol2_file cluster.mol2 --target_residue_name LIG
```

---

## Requirements

- Python 3.7+
- No external dependencies (uses only standard library).

---

## Notes

- Both scripts print warnings and errors for malformed input or mismatches.
- `mol2_to_xyz.py` can run a built-in example for demonstration.
- `cgenff_charge_validator.py` can handle atom index offsets for cluster calculations if provided with the appropriate MOL2 and residue name.

---

**Author:** Markus