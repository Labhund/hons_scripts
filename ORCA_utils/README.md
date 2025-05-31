#### ORCA Utils ####

This directory hosts a collection of Python scripts designed to streamline workflows and assist with quantum chemistry calculations using the [ORCA Quantum Chemistry Program](https://orcaforum.kofo.mpg.de/). While currently a work in progress, the aim is to provide a suite of helpful tools for pre-processing, post-processing, and analysis of ORCA jobs.

---

### Current Utilities

#### 1. `mol2_to_xyz.py`

*   **Purpose**: Converts molecular structure files from MOL2 format to XYZ format. XYZ is a simple and widely supported format, often required as input for ORCA and other quantum chemistry software.
*   **Features**:
    *   Parses atom coordinates and element types from MOL2 files.
    *   Extracts the molecule name from the MOL2 file for use in the XYZ comment line.
    *   Generates an XYZ file with the atom count on the first line, a detailed comment line on the second (including source file, conversion date, and script author), followed by atom coordinates.
    *   Command-line interface for easy conversion.
    *   Can be imported as a module into other Python scripts.
*   **Usage**:
    *   **Command-line**:
        ```bash
        python mol2_to_xyz.py -i <input_mol2_file> -o <output_xyz_file>
        ```
        Example:
        ```bash
        python mol2_to_xyz.py -i my_molecule.mol2 -o my_molecule.xyz
        ```
    *   **Run internal example**:
        ```bash
        python mol2_to_xyz.py --run_example
        ```
    *   **As a module**:
        ```python
        from mol2_to_xyz import mol2_to_xyz

        try:
            num_atoms = mol2_to_xyz("input.mol2", "output.xyz")
            print(f"Successfully converted {num_atoms} atoms.")
        except Exception as e:
            print(f"An error occurred: {e}")
        ```

---

### Future Plans

This collection is expected to grow. Potential future utilities include:

*   **Input File Generation**: Scripts to help generate ORCA input files from templates or other formats, potentially automating the setup of common calculation types (e.g., geometry optimizations, frequency calculations, TD-DFT).
*   **Output File Parsing**: Tools to extract key information from ORCA output files (`.out`, `.gbw`, `.prop`, etc.), such as:
    *   Optimized geometries (and conversion to various formats).
    *   Energies (SCF, ZPE, thermal corrections).
    *   Vibrational frequencies and normal modes.
    *   UV-Vis or ECD spectra data from TD-DFT calculations.
    *   Molecular orbital energies and occupancies.
*   **Batch Processing**: Utilities to run a series of ORCA calculations or processing steps on multiple molecules.
*   **Visualization Aids**: Scripts to generate data in formats suitable for common molecular visualization software or plotting libraries (e.g., VMD, ChimeraX, Matplotlib).
*   **Job Management**: Simple tools for submitting or monitoring ORCA jobs, especially in cluster environments (though this might be limited by system-specific queuing systems).
*   **Property Calculation Helpers**: Scripts to assist in calculating specific molecular properties derived from ORCA outputs.

Contributions and suggestions are welcome as the project develops.