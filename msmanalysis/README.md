# MSM Analysis Pipeline (`msmanalysis.py`)

## Project Overview

`msmanalysis.py` is a comprehensive Python script for performing Markov State Model (MSM) analysis on molecular dynamics (MD) trajectory data. It leverages the PyEMMA and MDTraj libraries to guide the user through a typical MSM construction workflow, from raw trajectory data to the identification and characterization of metastable states.

The pipeline includes:
1.  **Featurization**: Defining and extracting relevant features from the MD trajectory (e.g., backbone torsions).
2.  **Dimensionality Reduction**: Applying Time-lagged Independent Component Analysis (TICA) to reduce the dimensionality of the feature space.
3.  **Clustering**: Discretizing the TICA-projected data into microstates using KMeans clustering.
4.  **Implied Timescales (ITS) Estimation**: Calculating and plotting implied timescales to help choose an appropriate lag time for MSM construction.
5.  **MSM Estimation**: Building the Markov State Model from the discretized trajectories.
6.  **Metastable State Identification**: Using PCCA++ (Perron-Cluster Cluster Analysis) to group microstates into a user-defined number of metastable (macro) states.
7.  **MSM Validation**: Performing the Chapman-Kolmogorov test to assess the quality and Markovianity of the MSM.
8.  **Representative Structure Extraction**: Saving representative PDB structures for each identified metastable state.
9.  **Output Management**: Saving all key data, plots, and the final MSM object to a structured output directory, including a detailed log file.

## Core Workflow Steps

The script executes the following major steps sequentially:

1.  **Argument Parsing & Configuration**:
    *   Sets up analysis parameters from command-line arguments and hardcoded defaults.
    *   Dynamically creates an output directory (e.g., `msm_analysis_output_k<clusters>_lag<msm_lag>_stride<stride>`), backing up any existing directory with the same name.
    *   Initializes a logger to write all `print` statements to both the console and a log file within the output directory.

2.  **Featurizer Initialization**:
    *   Initializes a PyEMMA featurizer using a provided topology file (e.g., PDB).
    *   By default, adds backbone torsions (phi, psi angles) with cosine/sine transformation.

3.  **Feature Extraction**:
    *   Loads the MD trajectory (e.g., XTC, DCD) using MDTraj via PyEMMA's `coordinates.load`.
    *   Applies the defined features during trajectory loading.
    *   Supports striding to reduce data size and processing time.
    *   Utilizes multiple CPU cores for parallel processing if available (`n_jobs`).

4.  **Dimensionality Reduction (TICA)**:
    *   Performs TICA on the extracted features.
    *   Key parameters: `tica_lag` (lag time for TICA), `tica_dim` (number of TICA dimensions to retain).
    *   Outputs:
        *   `tica_projection.png`: A 2D hexbin plot of the trajectory projected onto the first two TICs.
        *   `tica_timescales.png`: A plot of the TICA timescales.
    *   The raw feature trajectory data is deleted after this step to save memory.

5.  **Clustering (KMeans)**:
    *   Clusters the TICA-projected data into a specified number of microstates (`n_clusters`) using KMeans.
    *   Outputs:
        *   `dtrajs.npy`: The discretized trajectories (sequences of cluster labels).
        *   `cluster_populations.png`: A bar chart showing the fractional population of each microstate.
    *   The TICA output data is deleted after this step to save memory.

6.  **Implied Timescales (ITS) Test**:
    *   Calculates implied timescales over a range of lag times (`its_lags`).
    *   Outputs:
        *   `implied_timescales.png`: A plot of the ITS. This plot is crucial for selecting an appropriate `msm_lag` where the timescales appear to converge.
    *   **User Action Required**: Examine this plot to choose a suitable `msm_lag`.

7.  **Markov State Model (MSM) Estimation**:
    *   Estimates an MSM from the discretized trajectories using the chosen `msm_lag`.
    *   Calculates and prints MSM properties: number of connected states, active state fraction, VAMP-2 score, and MSM timescales.

8.  **Metastable State Identification (PCCA++)**:
    *   If the MSM is successfully built and has sufficient connected states, performs PCCA++ to coarse-grain the MSM into `n_metastable` states.
    *   Prints properties of each metastable state (e.g., probability, constituent microstates).
    *   Outputs:
        *   `metastable_states_free_energy.png`: A plot of the free energy surface (projected onto the first two TICs) with the centers of the identified metastable states highlighted.

9.  **MSM Validation (Chapman-Kolmogorov Test)**:
    *   If PCCA++ was successful, performs the Chapman-Kolmogorov (CK) test to validate the MSM.
    *   Outputs:
        *   `chapman_kolmogorov_test_manual.png`: A plot showing the CK test results. For a valid MSM, the model's predictions should match the direct estimations of state populations over various lag times.
    *   **User Action Required**: Examine this plot to assess MSM validity.
    *   *Note: The script uses a manual plotting routine for CK results as there were previous issues with PyEMMA's direct plotting function with error bars.*

10. **Representative Structure Extraction**:
    *   If PCCA++ was successful, extracts a specified number of representative PDB structures (`n_samples_per_state`) for each metastable state.
    *   Maps frames from the analysis (strided) trajectory back to the original, unstrided trajectory.
    *   Loads the full, unstrided trajectory using MDTraj.
    *   Saves multi-model PDB files (e.g., `metastable_state_0_samples.pdb`) into a `representative_structures` subdirectory.
    *   Prints the 0-based frame indices (from the original trajectory) of the extracted structures.

11. **Save MSM Object**:
    *   Saves the estimated PyEMMA MSM object to a file (`msm_model.pyemma`) for later analysis or reloading.

12. **Final Summary**:
    *   Prints total execution time and a summary of key parameters and outcomes.

## Key Configuration Parameters

Many parameters can be controlled via command-line arguments. Key internal (default) and CLI-configurable parameters include:

*   **Input Files (Hardcoded Defaults - modify script if needed):**
    *   `topology_file = "reference.pdb"`
    *   `trajectory_file = "protein_only_stride10_truncated.xtc"`
*   **Output Directory:**
    *   `base_output_dir`: Base name, dynamically generated as `msm_analysis_output_k<n_clusters>_lag<msm_lag>_stride<stride>`.
*   **TICA Parameters (Hardcoded Defaults):**
    *   `tica_lag = 50` (frames)
    *   `tica_dim = 10`
*   **Clustering Parameters:**
    *   `n_clusters = 100` (Default, CLI: `-k` or `--clusters`)
*   **ITS Parameters (Hardcoded Default):**
    *   `its_lags = range(1, 101, 5)` (frames)
*   **MSM Parameters:**
    *   `msm_lag = 50` (frames, Default, CLI: `-l` or `--lag`). **Crucial to set based on ITS plot.**
*   **PCCA++ Parameters:**
    *   `n_metastable = 4` (Default, CLI: `-n` or `--meta` or `--n_metastable`)
*   **Trajectory Loading:**
    *   `stride = 1` (Default, CLI: `-s` or `--stride`). Stride for *analysis loading*, not necessarily the stride of the input trajectory file itself.
*   **Parallelization (Hardcoded Default):**
    *   `n_jobs = -1` (Use all available CPU cores)
*   **Structure Extraction (Hardcoded Default):**
    *   `n_samples_per_state = 5`

## Dependencies

*   **Python 3.x**
*   **MDTraj**: For trajectory manipulation and PDB file handling.
*   **PyEMMA**: For MSM construction and analysis (includes TICA, clustering, ITS, MSM, PCCA++, CK test).
*   **NumPy**: For numerical operations.
*   **Matplotlib**: For plotting.
*   **tqdm**: For progress bars (though not explicitly used for user-facing progress bars in the provided snippet, PyEMMA might use it internally).
*   **multiprocess**: Used for `freeze_support()` which is essential for robust multiprocessing, especially on Windows or when using 'spawn' as the start method.
*   Standard Python libraries: `os`, `sys`, `shutil`, `traceback`, `time`, `argparse`.

**Installation:**
You can install the primary dependencies using pip:
```bash
pip install mdtraj pyemma numpy matplotlib tqdm multiprocess
```

## Script Usage

The script is run from the command line.

```bash
python msmanalysis.py [OPTIONS]
```

### Key Command-Line Arguments:

*   `-k CLUSTERS`, `--clusters CLUSTERS`: Number of clusters (microstates) for KMeans. Default: `100`.
*   `-l LAG`, `--lag LAG`: MSM lag time in frames (should be chosen based on ITS plot). Default: `50`.
*   `-s STRIDE`, `--stride STRIDE`: Stride to use when loading trajectory features for analysis. Default: `1`.
*   `-n N_METASTABLE`, `--meta N_METASTABLE`: Number of metastable states to identify with PCCA++. Default: `4`.

*(Note: The script also has commented-out arguments for `tica_lag` and `tica_dim` which could be easily enabled if needed.)*

### Example:

```bash
# Make sure reference.pdb and protein_only_stride10_truncated.xtc are in the same directory
# or adjust the hardcoded paths in the script.

# Run with 150 clusters, an MSM lag of 60 frames (chosen after viewing ITS),
# an analysis stride of 2, and aiming for 5 metastable states.
python msmanalysis.py -k 150 -l 60 -s 2 -n 5
```

## Output

The script generates a primary output directory named according to the pattern `msm_analysis_output_k<k_val>_lag<lag_val>_stride<stride_val>`. If such a directory already exists, it is backed up by renaming it with a `_bak` suffix (or `_bakN` if multiple backups exist).

Inside this directory, you will find:

*   `msmanalysis.log`: A comprehensive log file containing all console output, timestamps, and configuration details.
*   **Plots (PNG files):**
    *   `tica_projection.png`: 2D projection of data onto the first two TICs.
    *   `tica_timescales.png`: TICA timescales.
    *   `cluster_populations.png`: Population of each microstate.
    *   `implied_timescales.png`: Implied timescales plot (critical for choosing `msm_lag`).
    *   `metastable_states_free_energy.png`: Free energy landscape with metastable states.
    *   `chapman_kolmogorov_test_manual.png`: CK test results.
*   **Data Files:**
    *   `dtrajs.npy`: Discretized trajectories (NumPy array).
    *   `msm_model.pyemma`: The saved PyEMMA MSM object.
*   **Representative Structures:**
    *   A subdirectory `representative_structures/` containing PDB files (e.g., `metastable_state_0_samples.pdb`) for the sampled conformations from each metastable state.

## Important Notes & User Actions

*   **Multiprocessing:** The script includes `freeze_support()` from the `multiprocess` library, which is crucial for stable parallel execution, particularly on Windows systems or when Python uses the 'spawn' start method for new processes.
*   **ITS Plot for MSM Lag:** The most critical manual step is to **examine the `implied_timescales.png` plot**. You need to identify a lag time where the implied timescales (especially the slower ones) appear to converge or plateau. This lag time should then be used for the `-l` (or `--lag`) argument when running the MSM estimation step. The script's default `msm_lag` is a guess and likely needs adjustment.
*   **CK Test for Validation:** After running the full pipeline, **examine the `chapman_kolmogorov_test_manual.png` plot**. For a valid MSM, the predicted probabilities (lines) should closely match the estimated probabilities (dots) for each metastable state across different lag times.
*   **Memory Management:** The script attempts to free memory by deleting large intermediate data structures (like raw feature trajectories and TICA output) once they are no longer needed. However, MSM analysis, especially with long trajectories or many features/clusters, can be memory-intensive.
*   **File Paths:** The paths to the `topology_file` and `trajectory_file` are hardcoded. You will need to either place your files with these exact names in the script's directory or modify these paths within the script.

## Error Handling and Logging

*   The script includes extensive print statements that are logged to both the console and the `msmanalysis.log` file.
*   It attempts to catch common errors (e.g., `FileNotFoundError`, issues during PyEMMA calculations) and print informative messages, often including tracebacks.
*   In case of critical errors in early stages (e.g., file loading, TICA, clustering, ITS), the script will typically exit to prevent further issues.
*   The `prepare_output_directory` function includes logic to back up existing output directories to prevent accidental data loss.
