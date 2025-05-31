# Honours Year Computational Chemistry Scripts (`hons_scripts`)

## About This Repository

This repository, `hons_scripts`, is a collection of Python scripts and utilities developed during my Honours year in Biochemistry and Computational Chemistry at The University of Sydney. These tools were created to assist with various tasks in the analysis of Molecular Dynamics (MD) simulations, protein structure manipulation, and data visualization.

The primary focus of these scripts is to streamline common workflows encountered in computational chemistry research, particularly related to protein dynamics, structural analysis, and Markov State Model (MSM) construction.

## Repository Structure and Contents

The repository is organized into several subdirectories, each containing specific tools or related files:

*   **`/msmanalysis/`**:
    *   Contains `msmanalysis.py`, a comprehensive Python script for performing Markov State Model (MSM) analysis on MD trajectory data. This script guides users through featurization, dimensionality reduction (TICA), clustering, MSM estimation, validation (Chapman-Kolmogorov), and metastable state identification (PCCA++), including representative structure extraction.
    *   See the `msmanalysis/README.md` for detailed usage instructions and dependencies.

*   **`/extract_pdb_seq/`**:
    *   Contains `extract_pdb_seq.py`, a Python script for parsing Protein Data Bank (PDB) files to extract protein sequences for each chain, converting three-letter amino acid codes to one-letter codes.
    *   See the `extract_pdb_seq/README.md` for detailed usage instructions.

*   **`/xvg_plotter/`**:
    *   Contains scripts and potentially notebooks for plotting data from `.xvg` files, which are commonly output by GROMACS.
    *   *(You might want to add a specific script name here if there's a primary one, or briefly describe its contents based on its current README).*
    *   See the `xvg_plotter/README.md` for more details.

*   **`/residue_mapping/`**:
    *   Contains scripts or utilities related to mapping residue numbers or schemes, potentially between different PDB structures or sequence numbering conventions.
    *   *(Similarly, add a brief description based on its current README).*
    *   See the `residue_mapping/README.md` for more details.

*   **`/miscellaneous/`**:
    *   A collection of smaller utility scripts or one-off analysis tools that don't fit into the more structured categories above. Currently includes:
        *   `analyse_md.sh`: *(Briefly describe its purpose if you recall or can quickly check)*
        *   `compare_xvg_plots.py`: *(Briefly describe its purpose)*
    *   These scripts may have minimal documentation; inspect them directly for usage.

## General Usage

Each primary script or toolset is located within its respective subdirectory, which also contains a dedicated `README.md` file. These individual READMEs provide detailed information on:
*   The purpose of the script(s).
*   Specific dependencies and installation instructions.
*   Command-line arguments and usage examples.
*   Expected input and output formats.

Please navigate to the relevant subdirectory for the tool you are interested in and consult its `README.md` for comprehensive guidance.

## Dependencies

Dependencies are generally listed in the `README.md` file for each specific script or subdirectory, as they can vary. Common libraries used across several Python scripts include:
*   Python 3.x
*   NumPy
*   MDTraj
*   PyEMMA
*   Matplotlib

Ensure you have these installed in your Python environment if you intend to use the relevant scripts. A `requirements.txt` file may be added in the future for easier environment setup.

## Contributing

While this repository primarily documents scripts from my Honours project, suggestions or identified issues are welcome via the "Issues" tab.

---

*This README was last updated: May 23, 2025.*
