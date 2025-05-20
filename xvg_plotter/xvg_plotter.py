#!/usr/bin/env python3

import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors # For LogNorm
matplotlib.rcParams['agg.path.chunksize'] = 50000
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

def read_xvg(filename, x_col=0, y_col=1):
    x_data, y_data = [], []
    # Initialize with generic labels; specific labels from XVG will override if present
    title_from_file, xlabel_from_file, ylabel_from_file = Path(filename).stem, "X-axis", "Y-axis"
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('@'):
                    parts = line.split('"')
                    if "title" in line and len(parts) > 1: title_from_file = parts[1]
                    if "xaxis  label" in line and len(parts) > 1: xlabel_from_file = parts[1]
                    if "yaxis  label" in line and len(parts) > 1: ylabel_from_file = parts[1]
                    continue
                if line.startswith('#') or not line:
                    continue
                
                cols = line.split()
                try:
                    x_data.append(float(cols[x_col]))
                    y_data.append(float(cols[y_col]))
                except (IndexError, ValueError):
                    print(f"Warning: Skipping malformed data line in {filename}: {line}")
                    continue
        
        if not x_data or not y_data:
             print(f"Error: No data columns could be read from {filename}.", file=sys.stderr)
             return None, None, title_from_file, xlabel_from_file, ylabel_from_file

        return np.array(x_data), np.array(y_data), title_from_file, xlabel_from_file, ylabel_from_file
    except FileNotFoundError:
        print(f"Error: XVG file not found: {filename}", file=sys.stderr)
        return None, None, title_from_file, xlabel_from_file, ylabel_from_file
    except Exception as e:
        print(f"Error reading XVG file {filename}: {e}", file=sys.stderr)
        return None, None, title_from_file, xlabel_from_file, ylabel_from_file

# --- Main Plotting Logic ---
def main():
    parser = argparse.ArgumentParser(
        description="Plot XVG files as line plots (optionally with B-factors) or as a combined 2D histogram.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--inputs', required=True, type=str, nargs='+', help="One or more input XVG file paths. For 2D histogram, data from all files will be combined.")
    parser.add_argument('--output', type=str, default="plot.png", help="Output plot file name.")
    
    parser.add_argument('--title', type=str, help="Overall plot title. If not set, a default is generated.") # MODIFIED: Default handled later
    parser.add_argument('--xlabel', type=str, help="X-axis label. Uses label from the first XVG if not provided.")
    parser.add_argument('--ylabel', type=str, help="Y-axis label. Uses label from the first XVG if not provided.")
    
    parser.add_argument('--bfac_pdb_id', type=str, help="PDB ID for B-factor extraction (e.g., 1R1K). Ignored for 2D histograms.")
    parser.add_argument('--bfac_chains', type=str, nargs='+', help="Chain IDs for B-factors. Ignored for 2D histograms.")
    parser.add_argument('--residue_map_files', type=str, nargs='+', 
                        help="Residue map file(s) (.rmp). Ignored for 2D histograms.")
    parser.add_argument('--awk_script_dir', type=str, default=".", help="Directory of extract_bfactors.awk. Ignored for 2D histograms.")
    parser.add_argument('--pdb_download_dir', type=str, default=".", help="Directory for PDB files. Ignored for 2D histograms.")

    parser.add_argument('--figsize', type=str, default="18,10", help="Figure size (width,height).")
    parser.add_argument('--dpi', type=int, default=300, help="DPI for raster output.")
    parser.add_argument('--legend_fontsize', type=str, default="x-small", help="Legend font size (for line plots).")

    parser.add_argument('--plot_2d_histogram', action='store_true', help="Generate a single 2D histogram combining data from all input files. B-factor arguments are ignored.")
    parser.add_argument('--hist_bins', type=int, default=100, help="Number of bins for the 2D histogram.")
    parser.add_argument('--hist_cmap', type=str, default='viridis', help="Colormap for the 2D histogram.")
    parser.add_argument('--hist_log_scale', action='store_true', help="Use a logarithmic scale for the 2D histogram colorbar counts.")

    args = parser.parse_args()

    try:
        figsize_vals = tuple(map(int, args.figsize.split(',')))
        if len(figsize_vals) != 2: raise ValueError("Figsize must be two comma-separated integers.")
    except ValueError as e:
        print(f"Invalid figsize format '{args.figsize}'. Using default (18,10). Error: {e}", file=sys.stderr)
        figsize_vals = (18, 10)

    fig, ax1 = plt.subplots(figsize=figsize_vals)
    
    # Initialize plot labels from the first file, can be overridden by args
    first_xvg_xlabel_file, first_xvg_ylabel_file, first_xvg_title_file = "X-axis", "Y-axis", "XVG Data Plot"
    if args.inputs:
        _, _, temp_title, temp_xlabel, temp_ylabel = read_xvg(args.inputs[0]) # Read first file for default labels
        if temp_xlabel: first_xvg_xlabel_file = temp_xlabel
        if temp_ylabel: first_xvg_ylabel_file = temp_ylabel
        if temp_title: first_xvg_title_file = temp_title


    final_plot_xlabel = args.xlabel or first_xvg_xlabel_file
    final_plot_ylabel = args.ylabel or first_xvg_ylabel_file
    final_plot_title = args.title # Will be set below if None

    if args.plot_2d_histogram:
        if not args.inputs:
            print("Error: --plot_2d_histogram requires at least one input file.", file=sys.stderr)
            sys.exit(1)
        
        print(f"\n--- Generating Combined 2D Histogram from {len(args.inputs)} file(s) ---")
        
        all_x_data = []
        all_y_data = []
        
        for i, xvg_filepath_str in enumerate(args.inputs):
            xvg_filepath = Path(xvg_filepath_str)
            print(f"  Reading data for histogram from: {xvg_filepath}")
            # Use default x_col=0, y_col=1 for hist2d
            x_data, y_data, xvg_title_file, xvg_xlabel_file, xvg_ylabel_file = read_xvg(xvg_filepath_str) 

            if x_data is not None and y_data is not None and x_data.size > 0 and y_data.size > 0:
                all_x_data.append(x_data)
                all_y_data.append(y_data)
                if i == 0: # Use labels from the first file if not overridden by args
                    final_plot_xlabel = args.xlabel or xvg_xlabel_file
                    final_plot_ylabel = args.ylabel or xvg_ylabel_file
                    if not args.title: # If user didn't provide a title
                        if len(args.inputs) == 1:
                            final_plot_title = f"2D Histogram of {xvg_filepath.stem}"
                        else:
                            final_plot_title = f"Combined 2D Histogram" # Generic for multiple files
            else:
                print(f"  Warning: No valid data read from {xvg_filepath_str}. Skipping for histogram.", file=sys.stderr)

        if not all_x_data or not all_y_data:
            print(f"Error: No valid data to plot combined 2D histogram. Exiting.", file=sys.stderr)
            sys.exit(1)

        combined_x = np.concatenate(all_x_data)
        combined_y = np.concatenate(all_y_data)

        if combined_x.size == 0 or combined_y.size == 0:
            print(f"Error: Combined data for 2D histogram is empty. Exiting.", file=sys.stderr)
            sys.exit(1)
            
        norm_choice = None
        if args.hist_log_scale:
            norm_choice = matplotlib.colors.LogNorm()

        counts, xedges, yedges, im = ax1.hist2d(
            combined_x, combined_y, 
            bins=args.hist_bins, 
            cmap=args.hist_cmap, 
            norm=norm_choice
        )
        cbar = fig.colorbar(im, ax=ax1)
        cbar.set_label('Counts')
        
        ax1.set_xlabel(final_plot_xlabel)
        ax1.set_ylabel(final_plot_ylabel)
        
        # If args.title was provided, it's already in final_plot_title.
        # If not, it was set based on single/multiple files.
        if args.title and len(args.inputs) > 1: # User provided title, and we combined files
             final_plot_title = f"{args.title} (Combined)"


        if args.bfac_pdb_id:
            print("Warning: B-factor plotting arguments are ignored when --plot_2d_histogram is active.", file=sys.stderr)

    else: # --- Original Line Plotting Logic ---
        final_plot_title = args.title or first_xvg_title_file # Default title from first file if not user-set
        color_cycler_rmsf = itertools.cycle(plt.cm.tab10.colors) 
        lines_for_legend, labels_for_legend = [], []
        
        print(f"\n--- Processing {len(args.inputs)} XVG input file(s) for line plot ---")
        for i, xvg_filepath_str in enumerate(args.inputs):
            xvg_filepath = Path(xvg_filepath_str)
            # print(f"Reading: {xvg_filepath}") # Can be noisy
            xvg_x, xvg_y, _, current_xvg_xlabel, current_xvg_ylabel = read_xvg(xvg_filepath_str)

            if xvg_x is None or xvg_y is None:
                print(f"Skipping {xvg_filepath_str} due to read error.", file=sys.stderr)
                continue
            
            # Apply sorting patch
            if xvg_x.size > 0:
                sort_indices = np.argsort(xvg_x)
                xvg_x = xvg_x[sort_indices]
                xvg_y = xvg_y[sort_indices]
            
            # Apply NaN insertion patch for gaps
            xvg_x_to_plot = xvg_x
            xvg_y_to_plot = xvg_y
            if xvg_x.size > 1:
                processed_x = [xvg_x[0]]
                processed_y = [xvg_y[0]]
                residue_gap_threshold = 1.5 
                for j in range(1, len(xvg_x)):
                    if (xvg_x[j] - xvg_x[j-1]) > residue_gap_threshold:
                        processed_x.append(np.nan)
                        processed_y.append(np.nan)
                    processed_x.append(xvg_x[j])
                    processed_y.append(xvg_y[j])
                xvg_x_to_plot = np.array(processed_x)
                xvg_y_to_plot = np.array(processed_y)

            # Set overall plot labels from the first file if not overridden by args
            # (already handled by final_plot_xlabel/ylabel initialization)

            label_parts = []
            # Try to get a meaningful prefix from the directory structure if multiple inputs
            if len(args.inputs) > 1 and xvg_filepath.parent.name not in [".", "analysis_output_V3"]:
                 label_parts.append(xvg_filepath.parent.name)
            
            stem_name = xvg_filepath.stem.replace("_per_res", "").replace("_calpha"," Cα").replace("ramachandran", "Rama.") # Shorten common names
            label_parts.append(stem_name)
            plot_label = " ".join(label_parts) if label_parts else xvg_filepath.name


            line, = ax1.plot(xvg_x_to_plot, xvg_y_to_plot, color=next(color_cycler_rmsf), label=plot_label, alpha=0.8, linewidth=1.5)        
            lines_for_legend.append(line)
            labels_for_legend.append(plot_label)
        
        if not lines_for_legend:
            print("Error: No valid XVG data could be plotted for line plot. Exiting.", file=sys.stderr)
            sys.exit(1)

        ax1.set_xlabel(final_plot_xlabel)
        ax1.set_ylabel(final_plot_ylabel, color='black')
        ax1.tick_params(axis='y', labelcolor='black')
        ax1.grid(True, alpha=0.3, axis='x', linestyle=':') 
        ax1.grid(True, alpha=0.3, axis='y', linestyle=':', color='grey') 

        ax2 = None 
        any_bfactors_plotted = False
        bfactor_info_for_title_suffix = [] # Renamed to avoid clash
        
        shared_map_data = None
        use_shared_map = False
        use_individual_maps = False

        if args.bfac_pdb_id and args.bfac_chains: # Only proceed if core B-factor args are present
            if args.residue_map_files: # Mapping is optional
                if len(args.residue_map_files) == 1:
                    # print(f"Using shared residue map file: {args.residue_map_files[0]} for all specified chains.") # Noisy
                    shared_map_data = parse_residue_map_rmp(args.residue_map_files[0])
                    if shared_map_data is None:
                        print(f"Warning: Could not parse shared map {args.residue_map_files[0]}. B-factors use original PDB numbers.", file=sys.stderr)
                    else:
                        use_shared_map = True
                elif len(args.residue_map_files) == len(args.bfac_chains):
                    # print("Using individual residue map file for each specified chain.") # Noisy
                    use_individual_maps = True
                else:
                    print("Error: Number of --residue_map_files must be 1 (shared) or match --bfac_chains count.", file=sys.stderr)
                    print("B-factors will use original PDB residue numbers if possible.", file=sys.stderr)
            else:
                print("No --residue_map_files provided. B-factors (if plotted) will use original PDB residue numbers.")


            # print(f"\n--- Processing B-factors for PDB {args.bfac_pdb_id}, Chains: {', '.join(args.bfac_chains)} ---") # Noisy
            awk_script_full_path = Path(args.awk_script_dir) / "extract_bfactors.awk"
            
            bfactor_styles = [
                {'color': 'dimgray', 'linestyle': ':', 'marker': '.'}, 
                {'color': 'darkslateblue', 'linestyle': '--', 'marker': '.'},
                {'color': 'darkgreen', 'linestyle': '-.', 'marker': '.'},
                {'color': 'maroon', 'linestyle': ':', 'marker': '.'},
            ] 
            common_style_params = {'markersize': 3.5, 'alpha': 0.7, 'linewidth': 1.2}
            style_cycler_bfac = itertools.cycle(bfactor_styles)

            for i, chain_id in enumerate(args.bfac_chains):
                # print(f"-- Processing Chain {chain_id} --") # Noisy
                original_bfac_res, original_bfac_vals = get_bfactors(
                    args.bfac_pdb_id, chain_id, str(awk_script_full_path), args.pdb_download_dir
                )

                if original_bfac_res is None or original_bfac_vals is None or original_bfac_res.size == 0:
                    print(f"Could not retrieve/process B-factors for {args.bfac_pdb_id} Chain {chain_id}.")
                    continue

                bfactor_x_for_plotting = original_bfac_res
                bfactor_y_for_plotting = original_bfac_vals
                chain_mapped_status = "orig.PDB" # Default status
                current_map_to_use = None
                map_file_name_for_log = "N/A"

                if use_shared_map and shared_map_data is not None:
                    current_map_to_use = shared_map_data
                    map_file_name_for_log = Path(args.residue_map_files[0]).name
                elif use_individual_maps:
                    map_file_path_str = args.residue_map_files[i]
                    map_file_name_for_log = Path(map_file_path_str).name
                    # print(f"Using individual map file: {map_file_path_str} for chain {chain_id}") # Noisy
                    current_map_to_use = parse_residue_map_rmp(map_file_path_str)
                    if current_map_to_use is None:
                        print(f"Warning: Could not parse map {map_file_path_str} for chain {chain_id}. Plotting with original residue numbers.")

                if current_map_to_use: # Attempt mapping only if a map was successfully loaded
                    mapped_res_for_xvg_axis, mapped_vals, unmapped_count = [], [], 0
                    for res_idx, res_num_from_bfac_pdb in enumerate(original_bfac_res):
                        if res_num_from_bfac_pdb in current_map_to_use:
                            mapped_res_for_xvg_axis.append(current_map_to_use[res_num_from_bfac_pdb])
                            mapped_vals.append(original_bfac_vals[res_idx])
                        else:
                            unmapped_count +=1
                    
                    if mapped_res_for_xvg_axis:
                        bfactor_x_for_plotting = np.array(mapped_res_for_xvg_axis)
                        bfactor_y_for_plotting = np.array(mapped_vals)
                        chain_mapped_status = f"map:{map_file_name_for_log}"
                        # print(f"Successfully mapped {len(mapped_vals)} B-factors for chain {chain_id} using map {map_file_name_for_log}.") # Noisy
                        if unmapped_count > 0:
                            print(f"Note: {unmapped_count} B-factor residues from Ch.{chain_id} not in map {map_file_name_for_log}.")
                    else:
                        print(f"Warning: No B-factors mapped for Ch.{chain_id} with map {map_file_name_for_log}. Using original PDB numbers.")
                
                bfactor_legend_label = f"Bfac {args.bfac_pdb_id} Ch.{chain_id} ({chain_mapped_status})"
                bfactor_info_for_title_suffix.append(f"Ch.{chain_id}({chain_mapped_status})")


                if bfactor_x_for_plotting.size > 0:
                    if not ax2: 
                        ax2 = ax1.twinx()
                        ax2_ylabel_base = f'B-factor (\u00C5\u00B2)' # Å²
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
                 ax2.tick_params(axis='y', labelcolor='black') # Keep B-factor y-axis ticks black

            if any_bfactors_plotted and bfactor_info_for_title_suffix:
                title_suffix_str = f" (Bfac: {args.bfac_pdb_id} " + ", ".join(bfactor_info_for_title_suffix) + ")"
                final_plot_title += title_suffix_str # Append B-factor info to the existing title
        
        if lines_for_legend:
            # Sort legend entries to group XVG lines first, then B-factors
            sorted_legend_elements = sorted(zip(labels_for_legend, lines_for_legend), key=lambda x: "Bfac" in x[0])
            labels_for_legend_sorted = [elem[0] for elem in sorted_legend_elements]
            lines_for_legend_sorted = [elem[1] for elem in sorted_legend_elements]
            ax1.legend(lines_for_legend_sorted, labels_for_legend_sorted, loc='best', fontsize=args.legend_fontsize)

    # --- Final Plot Adjustments and Saving (Common for both modes) ---
    ax1.set_title(final_plot_title, fontsize='medium')
    
    plt.tight_layout()
    
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        plt.savefig(output_path, dpi=args.dpi, bbox_inches='tight')
        print(f"\nPlot saved to {output_path}")
    except Exception as e:
        print(f"Error saving plot to {output_path}: {e}", file=sys.stderr)

if __name__ == '__main__':
    main()

