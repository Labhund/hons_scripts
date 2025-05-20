#!/usr/bin/env python3

import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors # For LogNorm and colormaps
import matplotlib.rcParams['agg.path.chunksize'] = 50000
from pathlib import Path
import sys
import itertools # For color cycling
import urllib.request

# --- Utility and Data Handling Functions --- (Keep these as they are from your V8/previous version)
def download_pdb_file(pdb_id, pdb_dir="."):
    pdb_dir_path = Path(pdb_dir)
    pdb_dir_path.mkdir(parents=True, exist_ok=True)
    pdb_file = pdb_dir_path / f"{pdb_id.lower()}.pdb"
    
    if pdb_file.exists():
        # print(f"PDB file {pdb_file} already exists.") # Can be noisy
        return pdb_file

    print(f"Attempting to download {pdb_id.upper()}.pdb to {pdb_file}...")
    try:
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        urllib.request.urlretrieve(url, pdb_file)
        print(f"Successfully downloaded {pdb_file}")
        return pdb_file
    except Exception as e:
        print(f"Error: Could not download PDB file {pdb_id}.pdb: {e}", file=sys.stderr)
        print(f"Please ensure {pdb_id.lower()}.pdb is manually placed in {pdb_dir_path} or check internet connection.", file=sys.stderr)
        return None

def get_bfactors(pdb_id, chain_id, awk_script_path_str, pdb_download_dir="."):
    pdb_file_path = download_pdb_file(pdb_id, pdb_download_dir)
    if not pdb_file_path:
        return None, None

    awk_script_path = Path(awk_script_path_str)
    if not awk_script_path.exists():
        print(f"Error: AWK script not found at {awk_script_path}", file=sys.stderr)
        return None, None

    cmd = [
        "awk",
        "-v", f"TARGET_CHAIN={chain_id}",
        "-f", str(awk_script_path),
        str(pdb_file_path)
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=120)
        output = result.stdout.strip()

        if not output:
            # print(f"Warning: No B-factor data extracted for PDB {pdb_id}, chain {chain_id}.") # Can be noisy
            return np.array([]), np.array([])

        data = np.array([line.split() for line in output.splitlines() if line.strip()], dtype=float)
        
        if data.size == 0 or data.ndim != 2 or data.shape[1] != 2:
            # print(f"Warning: B-factor data for {pdb_id} chain {chain_id} is empty or not in expected 2-column format.") # Can be noisy
            return np.array([]), np.array([])
            
        return data[:, 0].astype(int), data[:, 1]
    except subprocess.CalledProcessError as e:
        print(f"Error executing AWK script for {pdb_id} chain {chain_id}: {e}", file=sys.stderr)
        print(f"AWK stderr: {e.stderr}", file=sys.stderr)
        return None, None
    except subprocess.TimeoutExpired:
        print(f"Timeout executing AWK script for {pdb_id} chain {chain_id}.", file=sys.stderr)
        return None, None
    except Exception as e:
        print(f"An unexpected error occurred in get_bfactors for {pdb_id} chain {chain_id}: {e}", file=sys.stderr)
        return None, None

def parse_residue_map_rmp(map_file_path_str):
    mapping = {}
    map_file_path = Path(map_file_path_str)
    if not map_file_path.exists():
        print(f"Error: Residue map file not found: {map_file_path}", file=sys.stderr)
        return None
    
    try:
        with open(map_file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                parts = line.strip().split()
                if not parts or line.startswith("#"):
                    continue
                if len(parts) == 3: 
                    try:
                        pdb_res_num = int(parts[1]) 
                        sim_res_num = int(parts[2]) 
                        mapping[pdb_res_num] = sim_res_num
                    except ValueError:
                        print(f"Warning: Skipping malformed numeric data in map file {map_file_path} at line {line_num}: {line.strip()}")
                else:
                    print(f"Warning: Skipping malformed line (expected 3 columns, got {len(parts)}) in map file {map_file_path} at line {line_num}: {line.strip()}")
    except Exception as e:
        print(f"Error reading or parsing residue map file {map_file_path}: {e}", file=sys.stderr)
        return None
    
    if not mapping:
        print(f"Warning: No valid mappings found in {map_file_path}.")
    # else: # Can be noisy
        # print(f"Successfully parsed {len(mapping)} mappings from {map_file_path}.")
    return mapping

def read_ss_map_data(filepath, char_to_int_map):
    """Reads a character matrix from filepath and converts it using char_to_int_map."""
    matrix_data = []
    unknown_chars = set()
    try:
        with open(filepath, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("@"): # Skip empty/comment/metadata lines
                    continue

                int_row = []
                for char_idx, char in enumerate(line):
                    val = char_to_int_map.get(char)
                    if val is None:
                        if char not in unknown_chars:
                            # print(f"Warning: Unknown character '{char}' in {filepath} at line {line_num}, char {char_idx+1}. Mapping to -1.", file=sys.stderr)
                            unknown_chars.add(char)
                        int_row.append(-1) # Default for unknown characters, e.g. -1
                    else:
                        int_row.append(val)
                if int_row: # Only add if row is not empty after processing
                    matrix_data.append(int_row)

        if unknown_chars:
            print(f"Warning: Unknown characters encountered in {filepath}: {sorted(list(unknown_chars))}. These were mapped to -1 (or your default for unknown).", file=sys.stderr)

        if not matrix_data:
            print(f"Warning: No data read from {filepath}. Returning empty array.", file=sys.stderr)
            return np.array([[]]) # Return an empty 2D array

        # Ensure all rows have the same length by padding shorter rows or truncating longer ones
        # This is important for creating a valid NumPy array for pcolormesh
        if matrix_data:
            max_len = max(len(row) for row in matrix_data)
            for row_idx, row in enumerate(matrix_data):
                if len(row) < max_len:
                    matrix_data[row_idx].extend([-1] * (max_len - len(row))) # Pad with -1
                elif len(row) > max_len:
                    matrix_data[row_idx] = row[:max_len] # Truncate

        return np.array(matrix_data, dtype=int)

    except Exception as e:
        print(f"Error reading or processing SS map file {filepath}: {e}", file=sys.stderr)
        return np.array([[]]) # Return empty 2D array on error


def parse_ss_definitions(ss_def_str):
    """Parses SS definitions string (e.g., "0:Coil:gray,1:Sheet:red") into dict and lists."""
    ss_map_info = {} # To store {'id': {'label': label, 'color': color}}
    colors = []
    labels = []
    ids = []
    if not ss_def_str: return {}, [], [], []
    try:
        definitions = ss_def_str.split(',')
        parsed_defs = []
        for def_item in definitions:
            parts = def_item.split(':')
            if len(parts) == 3:
                id_val, label, color = int(parts[0]), parts[1], parts[2]
                parsed_defs.append({'id': id_val, 'label': label, 'color': color})
            else:
                print(f"Warning: Skipping malformed SS definition: {def_item}", file=sys.stderr)

        # Sort by ID for consistent colormap and legend
        parsed_defs.sort(key=lambda x: x['id'])

        for item in parsed_defs:
            ss_map_info[item['id']] = {'label': item['label'], 'color': item['color']}
            ids.append(item['id'])
            colors.append(item['color'])
            labels.append(item['label'])

        return ss_map_info, colors, labels, ids
    except Exception as e:
        print(f"Error parsing SS definitions '{ss_def_str}': {e}", file=sys.stderr)
        return {}, [], [], []

def generate_ss_map_plot(ax, data_matrix, time_coords, residue_offset, ss_colors, ss_labels, ss_ids, title, xlabel, ylabel):
    """Generates a 2D heatmap for secondary structure data."""
    if data_matrix.size == 0 or data_matrix.shape[1] == 0: # Check if matrix has columns
        print("Warning: SS data matrix is empty or has no columns. Skipping plot.", file=sys.stderr)
        ax.text(0.5, 0.5, "No SS data to display", ha='center', va='center', transform=ax.transAxes)
        return

    if not ss_colors or not ss_labels or not ss_ids:
        print("Warning: SS definitions (colors, labels, ids) are missing. Cannot create SS map.", file=sys.stderr)
        ax.text(0.5, 0.5, "SS definitions missing", ha='center', va='center', transform=ax.transAxes)
        return

    cmap = mcolors.ListedColormap(ss_colors)

    # Create boundaries for the colormap to ensure discrete colors
    # Boundaries should be halfway between the integer values.
    # e.g., if ss_ids are [0, 1, 2, 5], boundaries are [-0.5, 0.5, 1.5, 3.5, 5.5]
    # We should also handle the -1 case for unknown characters if present in data_matrix

    present_ids_in_data = sorted(list(np.unique(data_matrix)))
    all_defined_ids = sorted(list(set(ss_ids))) # Ensure ss_ids are unique and sorted

    # Create a combined list of all relevant IDs for boundaries: defined SS IDs and any IDs actually in the data
    boundary_basis_ids = sorted(list(set(all_defined_ids + present_ids_in_data)))

    if not boundary_basis_ids: # If no IDs at all (empty data, no defs)
        print("Warning: No IDs available for SS map boundaries.", file=sys.stderr)
        ax.text(0.5, 0.5, "No IDs for SS map", ha='center', va='center', transform=ax.transAxes)
        return

    boundaries = [bid - 0.5 for bid in boundary_basis_ids]
    boundaries.append(boundary_basis_ids[-1] + 0.5)

    norm = mcolors.BoundaryNorm(boundaries, cmap.N)

    n_frames, n_residues = data_matrix.shape

    if time_coords is not None and len(time_coords) == n_frames:
        x_coords = time_coords
    else:
        if time_coords is not None and n_frames > 0 : # Mismatch
            print(f"Warning: Time coordinates length ({len(time_coords)}) doesn't match data frames ({n_frames}). Using frame numbers for X-axis.", file=sys.stderr)
        x_coords = np.arange(n_frames)
        if "Time" in xlabel: xlabel = "Frame Number" # Adjust label

    # For pcolormesh, X and Y define the corners of the cells.
    # So, if x_coords has N elements, x_edges should have N+1.
    if len(x_coords) > 1:
        x_edges = np.concatenate([ [x_coords[0] - (x_coords[1]-x_coords[0])/2.0 if len(x_coords)>1 else x_coords[0]-0.5],
                                   (x_coords[:-1] + x_coords[1:])/2.0,
                                   [x_coords[-1] + (x_coords[-1]-x_coords[-2])/2.0 if len(x_coords)>1 else x_coords[-1]+0.5] ])
    elif len(x_coords) == 1: # Single frame
        x_edges = np.array([x_coords[0]-0.5, x_coords[0]+0.5])
    else: # No frames
        print("Warning: No x_coords (time/frames) for pcolormesh. SS Map cannot be drawn.", file=sys.stderr)
        return

    y_edges = np.arange(residue_offset - 0.5, residue_offset + n_residues + 0.5)

    # data_matrix is (frames, residues). pcolormesh expects Z(num_Y_edges-1, num_X_edges-1)
    # So we need to transpose data_matrix: data_matrix.T will be (residues, frames)
    mesh = ax.pcolormesh(x_edges, y_edges, data_matrix.T, cmap=cmap, norm=norm, shading='auto')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10) # Smaller title for subplot

    ax.set_xlim(x_edges[0], x_edges[-1])
    ax.set_ylim(y_edges[0], y_edges[-1]) # pcolormesh sets this, but good to be explicit

    # Colorbar ticks should correspond to the centers of the colors (i.e., the ss_ids themselves)
    if ss_labels and ss_ids:
        # Filter ss_ids and ss_labels to only those present in boundary_basis_ids for relevant ticks
        relevant_ticks = [tick_id for tick_id in all_defined_ids if tick_id in boundary_basis_ids]
        relevant_labels_for_ticks = [ss_labels[all_defined_ids.index(tick_id)] for tick_id in relevant_ticks]

        if relevant_ticks:
            cbar = plt.colorbar(mesh, ax=ax, ticks=relevant_ticks, spacing='proportional')
            cbar.ax.set_yticklabels(relevant_labels_for_ticks)
            cbar.set_label("Secondary Structure")
        else:
            cbar = plt.colorbar(mesh, ax=ax) # Default colorbar if no relevant ticks
            cbar.set_label("Secondary Structure ID")

    # ax.invert_yaxis() # Typically residue numbers increase upwards, pcolormesh might handle this. Check output.
    # If Y axis is inverted by default by pcolormesh and you want it standard, remove this.
    # If it's not inverted and you want residue 1 at the bottom, this is fine.
    # Let's assume standard: residue 1 at bottom. If pcolormesh does it wrong, then invert.
    if residue_offset < residue_offset + n_residues -1 : # only if more than one residue
        ax.set_ylim(min(y_edges), max(y_edges)) # Ensure correct limits
    # If residue numbers are plotted "backwards" (e.g. N at top, 1 at bottom), then:
    # ax.invert_yaxis()

def read_xvg(filename, x_col_index=0, y_col_index=None):
    """
    Reads an XVG file, attempting to parse metadata for title, labels, and legends.
    Can handle single or multiple Y-series based on @ sN legend lines.
    """
    x_data_list, y_data_lists = [], []
    title_from_file, xlabel_from_file, ylabel_from_file = Path(filename).stem, "X-axis", "Y-axis"
    legend_labels = []
    y_column_indices_from_legends = []
    temp_legends = {} # Stores @sN legends before data block

    lines = []
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        print(f"Error reading file {filename}: {e}", file=sys.stderr)
        return None, None, None, title_from_file, xlabel_from_file, ylabel_from_file

    # First pass: get metadata and legend labels
    for line_idx, line in enumerate(lines):
        line = line.strip()
        if line.startswith('@'):
            parts = line.split('"')
            if "title" in line and len(parts) > 1: title_from_file = parts[1]
            elif "xaxis  label" in line and len(parts) > 1: xlabel_from_file = parts[1]
            elif "yaxis  label" in line and len(parts) > 1: ylabel_from_file = parts[1]
            elif line.startswith("@ s") and "legend" in line and len(parts) > 1:
                try:
                    s_index = int(line.split()[1][1:]) # Get N from @ sN
                    # Handle special characters like '⏨' for 3₁₀-Helices
                    legend_text = parts[1]
                    if "3⏨-Helices" in legend_text:
                        legend_text = legend_text.replace("3⏨-Helices", "$3_{10}$-Helices")
                    temp_legends[s_index] = legend_text
                except (ValueError, IndexError):
                    print(f"Warning: Could not parse legend line: {line} in {filename}", file=sys.stderr)
            continue

    # Determine y_column_indices and legend_labels based on sorted temp_legends
    if temp_legends:
        # GROMACS XVG: col 0 is X, col 1 is s0, col 2 is s1, ...
        # We assume x_col_index is 0 for files with @ sN legends.
        sorted_s_indices = sorted(temp_legends.keys())
        legend_labels = [temp_legends[s_idx] for s_idx in sorted_s_indices]
        # y_column_indices are 1, 2, 3,... corresponding to s0, s1, s2...
        # these are relative to the start of data columns, so add x_col_index
        y_column_indices_from_legends = [s_idx + 1 + x_col_index for s_idx in sorted_s_indices]


    # Prepare for data reading
    current_x_data, current_y_data_cols_list = [], []
    if y_column_indices_from_legends: # Multiple Y series from legends
        for _ in y_column_indices_from_legends: current_y_data_cols_list.append([])
    elif y_col_index is not None: # User specified a single Y column explicitly
         current_y_data_cols_list.append([])
         if not legend_labels: legend_labels.append(f"Column {y_col_index}") # Use if no @sN legends
    else: # Default to first column after X if no legends and no y_col_index specified
        current_y_data_cols_list.append([])
        # Default legend label will be added later if still empty
        if not legend_labels: legend_labels.append("Y-data")


    # Second pass: read data
    for line_idx, line in enumerate(lines):
        line = line.strip()
        if line.startswith('@') or line.startswith('#') or not line:
            continue # Skip metadata/comments

        cols = line.split()
        try:
            x_val = float(cols[x_col_index])
            current_x_data.append(x_val)

            if y_column_indices_from_legends:
                for i, y_col_idx in enumerate(y_column_indices_from_legends):
                    current_y_data_cols_list[i].append(float(cols[y_col_idx]))
            elif y_col_index is not None and len(cols) > y_col_index: # Single specified Y column
                 current_y_data_cols_list[0].append(float(cols[y_col_index]))
            elif y_col_index is None and len(cols) > x_col_index + 1: # Default Y column (first after X)
                current_y_data_cols_list[0].append(float(cols[x_col_index + 1]))
            else: # Not enough columns for Y data or missing data
                # Pad with NaNs to maintain array structure if a value is missing
                if y_column_indices_from_legends:
                    for i in range(len(y_column_indices_from_legends)): current_y_data_cols_list[i].append(np.nan)
                elif current_y_data_cols_list: # For single Y column cases
                    current_y_data_cols_list[0].append(np.nan)

        except (ValueError, IndexError) as e:
            # print(f"Warning: Skipping malformed data line {line_idx+1} in {filename}: '{line}'. Error: {e}", file=sys.stderr)
            # Pad with NaNs if a row is malformed to maintain array shapes
            if y_column_indices_from_legends:
                for i in range(len(y_column_indices_from_legends)): current_y_data_cols_list[i].append(np.nan)
            elif current_y_data_cols_list:
                current_y_data_cols_list[0].append(np.nan)
            continue
    
    if current_x_data:
        x_data_list = np.array(current_x_data, dtype=float)
        # Filter out empty lists from current_y_data_cols_list before converting to np.array
        y_data_lists = [np.array(y_col, dtype=float) for y_col in current_y_data_cols_list if y_col] 

    # If no legends were found (e.g. simple 2-column XVG) and we defaulted, ensure a label exists
    if not legend_labels and y_data_lists: # This condition might need refinement
        legend_labels.append(f"{Path(filename).stem} Y-data")
    elif not y_data_lists and y_col_index is not None: # Case where y_col specified but no data found
        pass # legend_labels might have "Column X"
    elif not y_data_lists: # No data at all
        legend_labels = []


    return x_data_list, y_data_lists, legend_labels, title_from_file, xlabel_from_file, ylabel_from_file

# --- Main Plotting Logic ---
# --- Main Plotting Logic ---
def main():
    parser = argparse.ArgumentParser(
        description="Plot XVG files as line plots (optionally with B-factors), combined 2D histograms, or secondary structure maps.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # General arguments
    parser.add_argument('--inputs', required=True, type=str, nargs='+',
                        help="Input file(s). For SS map, only the first is used for map data. "
                             "For 2D histogram, data from all files is combined. "
                             "For line plots, all files are plotted.")
    parser.add_argument('--output', type=str, default="plot.png", help="Output plot file name.")
    parser.add_argument('--title', type=str, help="Overall plot title. If not set, a default is generated based on input/mode.")
    parser.add_argument('--xlabel', type=str, help="X-axis label. Uses label from XVG/time file if not provided.")
    parser.add_argument('--ylabel', type=str, help="Y-axis label. Uses label from XVG if not provided.")
    parser.add_argument('--figsize', type=str, default="18,10", help="Figure size (width,height in inches).")
    parser.add_argument('--dpi', type=int, default=300, help="DPI for raster output.")

    # Arguments for Line Plots
    parser.add_argument('--x_col', type=int, default=0,
                        help="Index of X column in XVG (0-based). For files with @ sN legends (e.g. ss_count.xvg), "
                             "this is typically 0 and data columns are auto-detected.")
    parser.add_argument('--y_col', type=int, default=None,
                        help="Index of Y column in basic XVG (0-based). If not set for basic XVG, uses first column after X. "
                             "Ignored for XVG files with @ sN legends where Y columns are auto-detected.")
    parser.add_argument('--legend_fontsize', type=str, default="x-small", help="Legend font size (for line plots).")

    # Arguments for B-factor overlay (for Line Plots)
    parser.add_argument('--bfac_pdb_id', type=str, help="PDB ID for B-factor extraction (e.g., 1R1K). Used with line plots.")
    parser.add_argument('--bfac_chains', type=str, nargs='+', help="Chain IDs for B-factors. Used with line plots.")
    parser.add_argument('--residue_map_files', type=str, nargs='+',
                        help="Residue map file(s) (.rmp) for B-factors. Used with line plots.")
    parser.add_argument('--awk_script_dir', type=str, default=".", help="Directory of extract_bfactors.awk. Used with line plots.")
    parser.add_argument('--pdb_download_dir', type=str, default=".", help="Directory for PDB files. Used with line plots.")

    # Arguments for 2D Histogram
    parser.add_argument('--plot_2d_histogram', action='store_true',
                        help="Generate a single 2D histogram combining data from all input files. Ignores B-factor and SS map arguments.")
    parser.add_argument('--hist_bins', type=int, default=100, help="Number of bins for the 2D histogram.")
    parser.add_argument('--hist_cmap', type=str, default='viridis', help="Colormap for the 2D histogram.")
    parser.add_argument('--hist_log_scale', action='store_true', help="Use a logarithmic scale for the 2D histogram colorbar counts.")

    # Arguments for Secondary Structure Map
    parser.add_argument('--plot_ss_map', action='store_true',
                        help="Generate a 2D map for secondary structure (e.g., from ss_residue.dat). "
                             "Ignores --plot_2d_histogram and B-factor args.")
    parser.add_argument('--ss_map_char_to_int', type=str, default="",
                        help="Comma-separated char:int map for SS map, e.g., '~:0,S:1,H:2'. "
                             "If empty, uses a GROMACS DSSP default.")
    parser.add_argument('--ss_map_definitions', type=str, default="",
                        help="Comma-separated id:label:color map for SS map colorbar, e.g., '0:Loop:gray,1:Bend:yellow,2:Alpha Helix:red'. "
                             "If empty, uses a default set.")
    parser.add_argument('--time_values_xvg_file', type=str,
                        help="XVG file (e.g., ss_count.xvg) to get time values for SS map X-axis. "
                             "If not provided, frame numbers are used.")
    parser.add_argument('--residue_offset', type=int, default=1, help="Starting residue number for SS map Y-axis.")

    args = parser.parse_args()

    try:
        figsize_vals = tuple(map(float, args.figsize.split(',')))
        if len(figsize_vals) != 2: raise ValueError("Figsize must be two comma-separated floats.")
    except ValueError as e:
        print(f"Invalid figsize format '{args.figsize}'. Using default (18,10). Error: {e}", file=sys.stderr)
        figsize_vals = (18.0, 10.0)

    fig, ax1 = plt.subplots(figsize=figsize_vals)
    final_plot_title = args.title # Will be refined per mode if None

    # --- Plotting Mode Logic ---

    if args.plot_ss_map:
        print(f"\n--- Generating Secondary Structure Map ---")
        if not args.inputs:
            print("Error: --plot_ss_map requires an input data file for the map.", file=sys.stderr)
            plt.close(fig)
            sys.exit(1)

        ss_map_filepath = args.inputs[0]
        if len(args.inputs) > 1:
            print(f"Warning: --plot_ss_map uses only the first input file '{ss_map_filepath}' for map data. Others ignored.", file=sys.stderr)

        char_to_int_map_actual = {}
        if args.ss_map_char_to_int:
            try:
                for item in args.ss_map_char_to_int.split(','):
                    char, val_str = item.split(':', 1)
                    char_to_int_map_actual[char.strip()] = int(val_str.strip())
            except Exception as e:
                print(f"Error parsing --ss_map_char_to_int '{args.ss_map_char_to_int}': {e}. Using default.", file=sys.stderr)
                char_to_int_map_actual = {}

        if not char_to_int_map_actual:
            char_to_int_map_actual = {
                '~': 0, ' ': 0, '.':0, 'C':0, # Coil/Loop (added C for Coil)
                'S': 1,                       # Bend
                'T': 2,                       # Turn
                'B': 3,                       # Bridge (Beta-bridge)
                'E': 4,                       # Strand (Beta-strand)
                'G': 5,                       # 3-10 Helix
                'H': 6,                       # Alpha Helix
                'I': 7,                       # Pi Helix
                'P': 8,                       # PP Helix (Polyproline II)
                '=': 9                        # Break (custom)
            }
            print(f"Using default char_to_int_map: {char_to_int_map_actual}", file=sys.stderr)

        ss_data_matrix = read_ss_map_data(ss_map_filepath, char_to_int_map_actual)

        if ss_data_matrix.size == 0 or (ss_data_matrix.ndim == 2 and ss_data_matrix.shape[1] == 0):
             print(f"Error: No data loaded or processed from SS map file {ss_map_filepath}. Cannot generate plot.", file=sys.stderr)
             ax1.text(0.5, 0.5, f"Failed to load data from\n{Path(ss_map_filepath).name}", ha='center', va='center', color='red', transform=ax1.transAxes)
             final_plot_title = final_plot_title or f"Error Loading {Path(ss_map_filepath).name}"
        else:
            time_coords = None
            current_xlabel_ss = args.xlabel
            if args.time_values_xvg_file:
                x_time, _, _, _, xvg_xlabel, _ = read_xvg(args.time_values_xvg_file, x_col_index=0)
                if x_time is not None and x_time.size > 0:
                    if ss_data_matrix.shape[0] == x_time.size:
                        time_coords = x_time
                        if not current_xlabel_ss: current_xlabel_ss = xvg_xlabel
                    else:
                        print(f"Warning: Time coordinates length ({x_time.size}) from {args.time_values_xvg_file} "
                              f"doesn't match data frames ({ss_data_matrix.shape[0]}) from {ss_map_filepath}. Using frame numbers.", file=sys.stderr)
                else:
                    print(f"Warning: Could not read valid time values from {args.time_values_xvg_file}. Using frame numbers.", file=sys.stderr)

            if not current_xlabel_ss : current_xlabel_ss = "Time (ns)" if time_coords is not None else "Frame Number"

            ss_def_map_actual, ss_colors_actual, ss_labels_actual, ss_ids_actual = parse_ss_definitions(args.ss_map_definitions)
            if not ss_colors_actual:
                default_defs = [
                    (0, 'Loop/Coil', 'silver'), (1, 'Bend', 'yellow'), (2, 'Turn', 'cyan'),
                    (3, 'β-Bridge', 'orange'), (4, 'β-Strand', 'red'), (5, '$3_{10}$-Helix', 'magenta'),
                    (6, 'α-Helix', 'blue'), (7, 'π-Helix', 'pink'), (8, 'PP-Helix', 'green'),
                    (9, 'Break', 'black'), (-1, 'Unknown', 'darkgrey')
                ]
                default_def_str = ",".join([f"{id_}:{label}:{color}" for id_, label, color in default_defs])
                ss_def_map_actual, ss_colors_actual, ss_labels_actual, ss_ids_actual = parse_ss_definitions(default_def_str)
                print(f"Using default ss_map_definitions.", file=sys.stderr)

            final_plot_title = final_plot_title or f"Secondary Structure: {Path(ss_map_filepath).name}"
            current_ylabel_ss = args.ylabel or f"Residue Number (Offset: {args.residue_offset})"

            generate_ss_map_plot(ax1, ss_data_matrix, time_coords, args.residue_offset,
                                 ss_colors_actual, ss_labels_actual, ss_ids_actual,
                                 final_plot_title, current_xlabel_ss, current_ylabel_ss)
        ax1.set_title(final_plot_title, fontsize='medium')

    elif args.plot_2d_histogram:
        if not args.inputs:
            print("Error: --plot_2d_histogram requires at least one input file.", file=sys.stderr)
            plt.close(fig)
            sys.exit(1)

        print(f"\n--- Generating Combined 2D Histogram from {len(args.inputs)} file(s) ---")
        all_x_data, all_y_data = [], []
        current_xlabel_hist, current_ylabel_hist = args.xlabel, args.ylabel
        hist_title_set = False

        for i, xvg_filepath_str in enumerate(args.inputs):
            xvg_filepath = Path(xvg_filepath_str)
            print(f"  Reading data for histogram from: {xvg_filepath}")
            x_data, y_data_list, _, file_title, file_xlabel, file_ylabel = read_xvg(xvg_filepath_str, args.x_col, args.y_col)

            if x_data is not None and y_data_list and y_data_list[0] is not None and \
               x_data.size > 0 and y_data_list[0].size > 0:
                all_x_data.append(x_data)
                all_y_data.append(y_data_list[0])
                if i == 0:
                    if not current_xlabel_hist: current_xlabel_hist = file_xlabel
                    if not current_ylabel_hist: current_ylabel_hist = file_ylabel
                    if not final_plot_title:
                        final_plot_title = f"2D Histogram of {file_title}" if len(args.inputs) == 1 else "Combined 2D Histogram"
                        hist_title_set = True
            else:
                print(f"  Warning: No valid X/Y data read from {xvg_filepath_str} for histogram. Skipping.", file=sys.stderr)

        if not all_x_data or not all_y_data:
            print(f"Error: No valid data to plot combined 2D histogram. Exiting.", file=sys.stderr)
            ax1.text(0.5, 0.5, "No data for 2D histogram", ha='center', va='center', transform=ax1.transAxes)
            if not hist_title_set: final_plot_title = final_plot_title or "2D Histogram Error"
        else:
            combined_x = np.concatenate(all_x_data)
            combined_y = np.concatenate(all_y_data)

            if combined_x.size == 0 or combined_y.size == 0:
                print(f"Error: Combined data for 2D histogram is empty. Exiting.", file=sys.stderr)
                ax1.text(0.5, 0.5, "Empty combined data", ha='center', va='center', transform=ax1.transAxes)
                if not hist_title_set: final_plot_title = final_plot_title or "2D Histogram Error"
            else:
                norm_choice = matplotlib.colors.LogNorm() if args.hist_log_scale else None
                counts, xedges, yedges, im = ax1.hist2d(
                    combined_x, combined_y, bins=args.hist_bins, cmap=args.hist_cmap, norm=norm_choice, cmin=1 if args.hist_log_scale else None
                )
                cbar = fig.colorbar(im, ax=ax1)
                cbar.set_label('Counts' if not args.hist_log_scale else 'Log Counts')

        ax1.set_xlabel(current_xlabel_hist or "X-axis")
        ax1.set_ylabel(current_ylabel_hist or "Y-axis")
        if args.title and len(args.inputs) > 1 and not hist_title_set :
             final_plot_title = f"{args.title} (Combined)"
        ax1.set_title(final_plot_title or "2D Histogram", fontsize='medium')

        if args.bfac_pdb_id or args.plot_ss_map: # Check this condition
            print("Warning: B-factor and SS Map arguments are ignored when --plot_2d_histogram is active.", file=sys.stderr)

    else: # --- Line Plotting Logic (Default Mode) ---
        print(f"\n--- Processing {len(args.inputs)} XVG input file(s) for line plot ---")
        lines_for_legend, labels_for_legend = [], []
        color_cycler_lines = itertools.cycle(plt.cm.tab10.colors)

        first_file_processed_line = False
        current_xlabel_line, current_ylabel_line = args.xlabel, args.ylabel
        line_plot_title_set = False

        for i, xvg_filepath_str in enumerate(args.inputs):
            xvg_filepath = Path(xvg_filepath_str)
            x_data, y_data_list, legend_labels_from_file, file_title, file_xlabel, file_ylabel = \
                read_xvg(xvg_filepath_str, args.x_col, args.y_col)

            if x_data is None or not y_data_list or x_data.size == 0 or not any(y is not None and y.size > 0 for y in y_data_list):
                print(f"Warning: No valid data read from {xvg_filepath_str}. Skipping.", file=sys.stderr)
                continue

            if not first_file_processed_line:
                if not current_xlabel_line: current_xlabel_line = file_xlabel
                if not current_ylabel_line: current_ylabel_line = file_ylabel
                if not final_plot_title:
                    final_plot_title = file_title
                    line_plot_title_set = True
                first_file_processed_line = True

            for y_idx, y_data_series in enumerate(y_data_list):
                if y_data_series is None or y_data_series.size == 0:
                    continue

                plot_label_parts = []
                # Try to get a meaningful prefix from the directory structure if multiple inputs
                if len(args.inputs) > 1 and xvg_filepath.parent.name not in [".", Path.cwd().name, str(Path.cwd())]:
                     plot_label_parts.append(xvg_filepath.parent.name)

                stem_name = xvg_filepath.stem.replace("_per_res", "").replace("_calpha"," C\u03B1").replace("ramachandran", "Rama.")
                plot_label_parts.append(stem_name)

                # Use legend from file if available and specific
                if legend_labels_from_file and y_idx < len(legend_labels_from_file) and \
                   legend_labels_from_file[y_idx] not in ["Y-data", f"Column {args.y_col or (args.x_col + 1)}"]:
                    plot_label_parts.append(legend_labels_from_file[y_idx])
                elif len(y_data_list) > 1 : # If multiple y-series from one file but no specific legend for this series
                    plot_label_parts.append(f"Series {y_idx+1}")

                plot_label = " - ".join(plot_label_parts)


                x_to_plot, y_to_plot = x_data, y_data_series
                if x_to_plot.size > 0:
                    sort_indices = np.argsort(x_to_plot)
                    x_to_plot = x_to_plot[sort_indices]
                    y_to_plot = y_to_plot[sort_indices]

                # Optional: NaN insertion for gaps (if you want to re-enable)
                # if x_to_plot.size > 1:
                #     processed_x, processed_y = [x_to_plot[0]], [y_to_plot[0]]
                #     # Define residue_gap_threshold if re-enabling
                #     for k_idx in range(1, len(x_to_plot)):
                #         # if (x_to_plot[k_idx] - x_to_plot[k_idx-1]) > residue_gap_threshold:
                #         #     processed_x.append(np.nan); processed_y.append(np.nan)
                #         processed_x.append(x_to_plot[k_idx]); processed_y.append(y_to_plot[k_idx])
                #     x_to_plot, y_to_plot = np.array(processed_x), np.array(processed_y)

                line, = ax1.plot(x_to_plot, y_to_plot, color=next(color_cycler_lines), label=plot_label, alpha=0.8, linewidth=1.5)
                lines_for_legend.append(line)
                labels_for_legend.append(plot_label)

        if not lines_for_legend:
            print("Error: No valid XVG data could be plotted for line plot.", file=sys.stderr)
            ax1.text(0.5, 0.5, "No line data to plot", ha='center', va='center', transform=ax1.transAxes)
            if not line_plot_title_set: final_plot_title = final_plot_title or "Line Plot Error"
        else:
            ax1.set_xlabel(current_xlabel_line or "X-axis")
            ax1.set_ylabel(current_ylabel_line or "Y-axis", color='black')
            ax1.tick_params(axis='y', labelcolor='black')
            ax1.grid(True, alpha=0.3, axis='x', linestyle=':')
            ax1.grid(True, alpha=0.3, axis='y', linestyle=':', color='grey')

            ax2 = None
            any_bfactors_plotted = False
            bfactor_info_for_title_suffix = []

            if args.bfac_pdb_id and args.bfac_chains:
                shared_map_data = None
                use_shared_map = False
                use_individual_maps = False

                if args.residue_map_files:
                    if len(args.residue_map_files) == 1:
                        shared_map_data = parse_residue_map_rmp(args.residue_map_files[0])
                        if shared_map_data is None:
                            print(f"Warning: Could not parse shared map {args.residue_map_files[0]}. B-factors use original PDB numbers.", file=sys.stderr)
                        else: use_shared_map = True
                    elif len(args.residue_map_files) == len(args.bfac_chains):
                        use_individual_maps = True
                    else:
                        print("Error: Number of --residue_map_files must be 1 (shared) or match --bfac_chains count. B-factors use original PDB numbers.", file=sys.stderr)
                # else: # Already printed by your original script if no map files
                    # print("No --residue_map_files provided. B-factors (if plotted) will use original PDB residue numbers.")

                awk_script_full_path = Path(args.awk_script_dir) / "extract_bfactors.awk"
                bfactor_styles = [{'color': 'dimgray', 'linestyle': ':', 'marker': '.'}, {'color': 'darkslateblue', 'linestyle': '--', 'marker': '.'},
                                  {'color': 'darkgreen', 'linestyle': '-.', 'marker': '.'}, {'color': 'maroon', 'linestyle': ':', 'marker': '.'}]
                common_style_params = {'markersize': 3.5, 'alpha': 0.7, 'linewidth': 1.2}
                style_cycler_bfac = itertools.cycle(bfactor_styles)

                for chain_idx, chain_id_val in enumerate(args.bfac_chains):
                    original_bfac_res, original_bfac_vals = get_bfactors(
                        args.bfac_pdb_id, chain_id_val, str(awk_script_full_path), args.pdb_download_dir
                    )

                    if original_bfac_res is None or original_bfac_vals is None or original_bfac_res.size == 0:
                        print(f"Could not retrieve/process B-factors for {args.bfac_pdb_id} Chain {chain_id_val}.")
                        continue

                    bfactor_x_for_plotting = original_bfac_res
                    bfactor_y_for_plotting = original_bfac_vals
                    chain_mapped_status = "orig.PDB"
                    current_map_to_use = None
                    map_file_name_for_log = "N/A"

                    if use_shared_map and shared_map_data is not None:
                        current_map_to_use = shared_map_data
                        map_file_name_for_log = Path(args.residue_map_files[0]).name
                    elif use_individual_maps and chain_idx < len(args.residue_map_files):
                        map_file_path_str = args.residue_map_files[chain_idx]
                        map_file_name_for_log = Path(map_file_path_str).name
                        current_map_to_use = parse_residue_map_rmp(map_file_path_str)
                        if current_map_to_use is None:
                            print(f"Warning: Could not parse map {map_file_path_str} for chain {chain_id_val}. Plotting with original residue numbers.")

                    if current_map_to_use:
                        mapped_res_for_xvg_axis, mapped_vals, unmapped_count = [], [], 0
                        for res_idx_map, res_num_from_bfac_pdb in enumerate(original_bfac_res):
                            if res_num_from_bfac_pdb in current_map_to_use:
                                mapped_res_for_xvg_axis.append(current_map_to_use[res_num_from_bfac_pdb])
                                mapped_vals.append(original_bfac_vals[res_idx_map])
                            else: unmapped_count +=1

                        if mapped_res_for_xvg_axis:
                            bfactor_x_for_plotting = np.array(mapped_res_for_xvg_axis)
                            bfactor_y_for_plotting = np.array(mapped_vals)
                            chain_mapped_status = f"map:{map_file_name_for_log}"
                            if unmapped_count > 0:
                                print(f"Note: {unmapped_count} B-factor residues from Ch.{chain_id_val} not in map {map_file_name_for_log}.")
                        else:
                            print(f"Warning: No B-factors mapped for Ch.{chain_id_val} with map {map_file_name_for_log}. Using original PDB numbers.")

                    bfactor_legend_label = f"Bfac {args.bfac_pdb_id} Ch.{chain_id_val} ({chain_mapped_status})"
                    bfactor_info_for_title_suffix.append(f"Ch.{chain_id_val}({chain_mapped_status})")

                    if bfactor_x_for_plotting.size > 0:
                        if not ax2:
                            ax2 = ax1.twinx()
                            ax2_ylabel_base = f'B-factor (\u00C5\u00B2)' # Angstrom squared
                            ax2.set_ylabel(ax2_ylabel_base)
                            ax2.grid(True, alpha=0.2, axis='y', linestyle=':', color='lightgrey')

                        current_style = next(style_cycler_bfac)
                        line_bfac, = ax2.plot(bfactor_x_for_plotting, bfactor_y_for_plotting,
                                              label=bfactor_legend_label,
                                              **current_style, **common_style_params)
                        lines_for_legend.append(line_bfac)
                        labels_for_legend.append(bfactor_legend_label)
                        any_bfactors_plotted = True

                if ax2 and any_bfactors_plotted:
                     ax2.tick_params(axis='y', labelcolor='black')

                if any_bfactors_plotted and bfactor_info_for_title_suffix:
                    title_suffix_str = f" (Bfac: {args.bfac_pdb_id} " + ", ".join(bfactor_info_for_title_suffix) + ")"
                    if final_plot_title: # Append if title already exists
                        final_plot_title += title_suffix_str
                    else: # Set if title was empty
                        final_plot_title = title_suffix_str.strip() # Remove leading space if it's the only title part

            if lines_for_legend:
                sorted_legend_elements = sorted(zip(labels_for_legend, lines_for_legend), key=lambda x: "Bfac" in x[0])
                labels_for_legend_sorted = [elem[0] for elem in sorted_legend_elements]
                lines_for_legend_sorted = [elem[1] for elem in sorted_legend_elements]
                ax1.legend(lines_for_legend_sorted, labels_for_legend_sorted, loc='best', fontsize=args.legend_fontsize)

        ax1.set_title(final_plot_title or "Line Plot", fontsize='medium')


    # --- Final Plot Adjustments and Saving (Common for all modes) ---
    plt.tight_layout()
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        plt.savefig(output_path, dpi=args.dpi, bbox_inches='tight')
        print(f"\nPlot saved to {output_path}")
    except Exception as e:
        print(f"Error saving plot to {output_path}: {e}", file=sys.stderr)
    finally:
        plt.close(fig) # Close the figure to free memory

if __name__ == '__main__':
    main()

