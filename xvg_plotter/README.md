# XVG Plotter (`pxvg`)

## Project Overview

This Python script, `xvg_plotter.py` (aliased as `pxvg`), is a versatile tool for plotting data from GROMACS `.xvg` files. XVG files are structured text files containing simulation analysis data, often with metadata headers (lines starting with `@` or `#`). This plotter leverages Python's `matplotlib` library to generate various types of plots, including line plots, 2D histograms, and secondary structure evolution maps.

A key feature is its ability to overlay crystallographic B-factors (from PDB files) onto line plots (e.g., RMSF data), potentially using residue mapping files to correlate simulation residue numbers with PDB residue numbers.

## Core Functionalities

1.  **Line Plots (Default Mode)**:
    *   Reads and plots one or more XVG files as standard line graphs.
    *   Automatically attempts to parse titles, axis labels, and series legends from XVG metadata.
    *   Supports plotting multiple Y-columns from a single XVG file (e.g., from `gmx sscount`).
    *   Handles multiple input XVG files, plotting them on the same axes with distinct labels.
    *   Customizable X and Y column selection (though often auto-detected for GROMACS files).

2.  **B-Factor Overlay (for Line Plots)**:
    *   Downloads PDB files from the RCSB PDB.
    *   Extracts B-factors for specified chains using an accompanying `extract_bfactors.awk` script.
    *   Can map B-factor residue numbers to simulation residue numbers using `.rmp` (residue map) files (compatible with the output of `create_residue_map.py` from the `residue_mapping` toolkit).
    *   Plots B-factors on a secondary Y-axis against the primary line plot data (e.g., RMSF vs. B-factors).

3.  **2D Histograms (`--plot_2d_histogram`)**:
    *   Combines X-Y data from one or more input XVG files to generate a single 2D histogram (heatmap).
    *   Useful for visualizing density of states, like Ramachandran plots.
    *   Supports customizable bin numbers, colormaps, and logarithmic color scaling.

4.  **Secondary Structure Evolution Maps (`--plot_ss_map`)**:
    *   Visualizes the evolution of secondary structure over time (or frames).
    *   Input is typically a file where rows are frames/time and columns are residues, with characters representing SS types (e.g., output from `gmx dssp` with `do_dssp -om`).
    *   Uses customizable character-to-integer mappings and color/label definitions for different secondary structures.
    *   Can use time values from a separate XVG file (e.g., `ss_count.xvg`) for the X-axis.

## Dependencies

*   **Python 3.x**
*   **NumPy**: For numerical operations.
*   **Matplotlib**: For plotting.
*   **AWK**: `gawk` or `mawk` must be installed and in the system PATH for B-factor extraction. The script `extract_bfactors.awk` (expected to be in `--awk_script_dir`) is used.

Install Python dependencies using pip:
```bash
pip install numpy matplotlib
```

## Script Usage (`xvg_plotter.py` or `pxvg`)

The script is command-line driven with numerous options. Here's a basic overview:

```bash
python xvg_plotter.py --inputs <file1.xvg> [<file2.xvg>...] --output <plot.png> [OPTIONS]
```

### Key Arguments:

**General:**
*   `--inputs <file(s)>`: (Required) One or more input XVG or data files.
*   `--output <filename>`: Output plot filename (e.g., `plot.png`). Default: `plot.png`.
*   `--title <text>`: Overall plot title.
*   `--xlabel <text>`, `--ylabel <text>`: Axis labels (often auto-detected from XVG).
*   `--figsize <W,H>`: Figure size in inches (e.g., `18,10`).
*   `--dpi <int>`: DPI for raster output (e.g., `300`).

**Line Plot Specific:**
*   `--x_col <int>`: Index of X column (0-based). Default: `0`.
*   `--y_col <int>`: Index of Y column (0-based) for simple XVG. Auto-detected for GROMACS files with `@ sN` legends. Default: `None` (uses column after X).
*   `--legend_fontsize <size>`: Font size for legend. Default: `x-small`.

**B-Factor Overlay (used with Line Plots):**
*   `--bfac_pdb_id <PDB_ID>`: PDB ID for B-factor extraction (e.g., `1R1K`).
*   `--bfac_chains <ChainID(s)>`: Chain ID(s) for B-factors (e.g., `A` or `A B`).
*   `--residue_map_files <mapfile(s).rmp>`: Path to residue map file(s). One file for all chains or one per chain.
*   `--awk_script_dir <path>`: Directory containing `extract_bfactors.awk`. Default: `.`.
*   `--pdb_download_dir <path>`: Directory to download/store PDB files. Default: `.`.

**2D Histogram:**
*   `--plot_2d_histogram`: Activates 2D histogram mode.
*   `--hist_bins <int>`: Number of bins. Default: `100`.
*   `--hist_cmap <colormap>`: Matplotlib colormap. Default: `viridis`.
*   `--hist_log_scale`: Use logarithmic color scale for counts.

**Secondary Structure Map:**
*   `--plot_ss_map`: Activates SS map mode.
*   `--ss_map_char_to_int "<char>:<int>,..."`: Mapping for SS characters to integers (e.g., `"~:0,S:1,H:2"`). Defaults to a GROMACS DSSP-like set.
*   `--ss_map_definitions "<id>:<label>:<color>,..."`: Definitions for SS colorbar (e.g., `"0:Loop:gray,1:Bend:yellow"`). Defaults to a pre-defined set.
*   `--time_values_xvg_file <file.xvg>`: XVG file to get time values for SS map X-axis (e.g., `ss_count.xvg`).
*   `--residue_offset <int>`: Starting residue number for Y-axis. Default: `1`.

**For detailed help on all arguments:**
```bash
python xvg_plotter.py --help
```

### `extract_bfactors.awk`

This AWK script is essential for the B-factor plotting functionality. It parses a PDB file and extracts residue numbers and B-factor values for C-alpha atoms of a specified chain.

**Purpose:**
*   Input: PDB file path, Target Chain ID.
*   Output: Lines of `ResidueNumber BFactorValue` for the target chain.

It should be located in the directory specified by `--awk_script_dir` (or the current directory by default).

## Examples

**1. Simple Line Plot:**
```bash
python xvg_plotter.py --inputs rmsd.xvg --output rmsd_plot.png --title "RMSD Evolution"
```

**2. Plot RMSF with B-Factor Overlay (using a residue map):**
```bash
python xvg_plotter.py --inputs rmsf_byres.xvg \
                      --output rmsf_bfac.png \
                      --title "RMSF vs B-Factors for 1R1K Chain A" \
                      --xlabel "Residue Number (Simulation)" \
                      --ylabel "RMSF (nm)" \
                      --bfac_pdb_id 1R1K \
                      --bfac_chains A \
                      --residue_map_files 1r1k_sim_to_pdb.rmp \
                      --awk_script_dir /path/to/your/scripts/
```

**3. Ramachandran Plot (2D Histogram):**
```bash
python xvg_plotter.py --inputs ramachandran.xvg \
                      --output rama_plot.png \
                      --plot_2d_histogram \
                      --xlabel "Phi Angle" --ylabel "Psi Angle" \
                      --title "Ramachandran Plot" \
                      --hist_log_scale
```

**4. Secondary Structure Evolution Map:**
```bash
python xvg_plotter.py --inputs ss_by_residue.dat \
                      --output ss_evolution.png \
                      --plot_ss_map \
                      --time_values_xvg_file ss_count.xvg \
                      --title "Secondary Structure Evolution" \
                      --ss_map_definitions "0:Coil:grey,4:Strand:red,6:Helix:blue"
                      # (Customize ss_map_char_to_int if your .dat uses different characters)
```

## Current Status & Notes

*   The script is feature-rich and handles many common GROMACS plotting tasks.
*   Error handling and warnings are included to guide users.
*   Default settings are provided for many options to simplify common use cases.
*   The B-factor mapping relies on correctly formatted PDB files and accurate residue map files if sequence/numbering differs significantly between simulation and PDB.
*   The `extract_bfactors.awk` script is a critical component for B-factor functionality.

## Repository Structure (Inferred)

Located within the `hons_scripts/xvg_plotter/` directory:
*   `README.md` (This file)
*   `xvg_plotter.py`: The main plotting script.
*   `extract_bfactors.awk`: AWK script for PDB B-factor extraction.
*   `old_ver/`: Contains older versions of the script.

---
*This README was generated based on the provided script and GitHub context.*
