# PDB Sequence Extractor (`extract_pdb_seq.py`)

## Project Overview

`extract_pdb_seq.py` is a Python script designed to parse Protein Data Bank (PDB) files and extract protein sequence information. It identifies unique amino acid residues for each chain present in the PDB file, converts their three-letter codes to one-letter codes, and outputs the sequence along with the start and end residue numbers for each chain.

## Core Functionalities

*   **PDB Parsing:** Reads ATOM records from a given PDB file.
*   **Chain Identification:** Groups residues by their chain ID.
*   **Residue Information Extraction:** Extracts residue name (three-letter code) and residue number.
*   **Sequence Conversion:** Converts standard three-letter amino acid codes to their corresponding one-letter codes (e.g., ALA to A). Unknown or non-standard residues are represented by 'X'.
*   **Uniqueness:** Ensures that each residue (defined by chain ID and residue number) is processed only once, even if multiple atoms define it.
*   **Output:** For each chain, the script prints:
    *   The chain ID.
    *   The start and end residue numbers observed in the PDB file for that chain.
    *   The complete protein sequence for that chain in one-letter code format.

## Dependencies

*   **Python 3.x**
*   Standard Python libraries:
    *   `argparse` (for command-line argument parsing)
    *   `sys` (for system-specific parameters and functions, like `stderr`)
    *   `collections.defaultdict` (for convenient data structuring)

No external packages need to be installed beyond a standard Python distribution.

## Script Usage

The script is run from the command line and requires the path to a PDB file as input.

```bash
python extract_pdb_seq.py -f <path_to_pdb_file>
```

### Arguments:

*   `-f FILE_PATH`, `--file FILE_PATH`: (Required) Specifies the path to the input PDB file.

### Example:

```bash
python extract_pdb_seq.py -f ./my_protein.pdb
```

Or using a PDB file from your repository context:
```bash
python extract_pdb_seq.py -f ./1r1k_apo_full.pdb
```

## Output Format

The script first prints status messages or errors to `stderr`. The primary output (sequence information) is printed to `stdout`.

**Example Output (to `stdout`):**

```
Sequence String (per chain):
Chain A (Res 10-250): MGSDKIHVNETAGGTSIVYRGTYIGH...
Chain B (Res 5-190): MKIHSLEFGKSLVTYRGFHILK...
```

*(The sequences above are illustrative examples)*

## Error Handling

*   **File Not Found:** If the specified PDB file does not exist, an error message is printed to `stderr`, and the script exits.
*   **Parsing Errors:** If a line in the PDB file cannot be parsed as expected (e.g., due to formatting issues or missing data in ATOM records), a message is printed to `stderr` indicating the line and error, and that line is skipped.
*   **No Data:** If no sequence data can be extracted (e.g., empty file, no valid ATOM records), a message is printed to `stderr`, and the script exits with a non-zero status.
