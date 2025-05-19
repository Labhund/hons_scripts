#!/usr/bin/env python3

import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import itertools # For color cycling

# --- Utility and Data Handling Functions --- (Mostly same as V6)

def download_pdb_file(pdb_id, pdb_dir="."):
    pdb_dir_path = Path(pdb_dir)
    pdb_dir_path.mkdir(parents=True, exist_ok=True)
    pdb_file = pdb_dir_path / f"{pdb_id.lower()}.pdb" 
    
    if pdb_file.exists():
        print(f"PDB file {pdb_file} already exists.")
        return pdb_file

    print(f"Attempting to download {pdb_id.upper()}.pdb to {pdb_file}...")
    try:
        import urllib.request
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        urllib.request.urlretrieve(url, pdb_file)
        print(f"Successfully downloaded {pdb_file}")
        return pdb_file
    except Exception as e:
        print(f"Error: Could not download PDB file {pdb_id}.pdb: {e}")
        print(f"Please ensure {pdb_id.lower()}.pdb is manually placed in {pdb_dir_path} or check internet connection.")
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
        print(f"Executing: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=120)
        output = result.stdout.strip()

        if not output:
            print(f"Warning: No B-factor data extracted for PDB {pdb_id}, chain {chain_id}.")
            return np.array([]), np.array([])

        data = np.array([line.split() for line in output.splitlines() if line.strip()], dtype=float)
        
        if data.size == 0 or data.ndim != 2 or data.shape[1] != 2:
            print(f"Warning: B-factor data for {pdb_id} chain {chain_id} is empty or not in expected 2-column format.")
            return np.array([]), np.array([])
            
        print(f"Successfully extracted {data.shape[0]} B-factors for {pdb_id} chain {chain_id}.")
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
                        pdb1_res_num = int(parts[1])
                        pdb2_res_num = int(parts[2])
                        mapping[pdb1_res_num] = pdb2_res_num
                    except ValueError:
                        print(f"Warning: Skipping malformed numeric data in map file {map_file_path} at line {line_num}: {line.strip()}")
                else:
                    print(f"Warning: Skipping malformed line (expected 3 columns, got {len(parts)}) in map file {map_file_path} at line {line_num}: {line.strip()}")
    except Exception as e:
        print(f"Error reading or parsing residue map file {map_file_path}: {e}", file=sys.stderr)
        return None
    
    if not mapping:
        print(f"Warning: No valid mappings found in {map_file_path}.")
    else:
        print(f"Successfully parsed {len(mapping)} mappings from {map_file_path}.")
    return mapping

def read_xvg(filename, x_col=0, y_col=1):
    x_data, y_data = [], []
    title, xlabel, ylabel = Path(filename).stem, "X-axis", "Y-axis"
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('@'):
                    parts = line.split('"')
                    if "title" in line and len(parts) > 1: title = parts[1]
                    if "xaxis  label" in line and len(parts) > 1: xlabel = parts[1]
                    if "yaxis  label" in line and len(parts) > 1: ylabel = parts[1]
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
             return None, None, title, xlabel, ylabel

        return np.array(x_data), np.array(y_data), title, xlabel, ylabel
    except FileNotFoundError:
        print(f"Error: XVG file not found: {filename}", file=sys.stderr)
        return None, None, title, xlabel, ylabel # Return None for data if file not found
    except Exception as e:
        print(f"Error reading XVG file {filename}: {e}", file=sys.stderr)
        return None, None, title, xlabel, ylabel


# --- Main Plotting Logic ---
def main():
    parser = argparse.ArgumentParser(
        description="Plot multiple XVG files on the primary Y-axis, optionally overlaying B-factors (mapped) on a secondary Y-axis.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--inputs', required=True, type=str, nargs='+', help="One or more input XVG file paths.")
    parser.add_argument('--output', type=str, default="plot.png", help="Output plot file name (e.g., plot.png, plot.svg).")
    
    parser.add_argument('--title', type=str, default="XVG Data Plot", help="Overall plot title.")
    parser.add_argument('--xlabel', type=str, help="X-axis label. If not provided, uses label from the first XVG file.")
    parser.add_argument('--ylabel', type=str, help="Y-axis label (for primary data - RMSF). If not provided, uses label from the first XVG file.")
    
    parser.add_argument('--bfac_pdb_id', type=str, help="PDB ID for B-factor extraction (e.g., 1R1K).")
    parser.add_argument('--bfac_chain', type=str, help="Chain ID in the PDB for B-factor extraction (e.g., A).")
    parser.add_argument('--residue_map_file', type=str, help="Path to the .rmp residue map file (e.g., 1r1k_A_to_sim_A_map.rmp).")
    
    parser.add_argument('--awk_script_dir', type=str, default=".", help="Directory containing extract_bfactors.awk.")
    parser.add_argument('--pdb_download_dir', type=str, default=".", help="Directory to store/find PDB files.")

    parser.add_argument('--figsize', type=str, default="15,8", help="Figure size (width,height) in inches, e.g., '12,7'.")
    parser.add_argument('--dpi', type=int, default=300, help="DPI for raster image output.")
    parser.add_argument('--legend_fontsize', type=str, default="small", help="Font size for the legend (e.g., 'xx-small', 'small', 'medium').")


    args = parser.parse_args()

    # --- Setup Plot ---
    try:
        figsize_vals = tuple(map(int, args.figsize.split(',')))
        if len(figsize_vals) != 2: raise ValueError("Figsize must be two comma-separated integers.")
    except ValueError as e:
        print(f"Invalid figsize format '{args.figsize}'. Using default (15,8). Error: {e}")
        figsize_vals = (15, 8)

    fig, ax1 = plt.subplots(figsize=figsize_vals)
    
    # Define a color cycle for the primary plots
    # Using 'tab10' which has 10 distinct colors. For more, consider 'tab20' or other colormaps.
    color_cycler = itertools.cycle(plt.cm.tab10.colors) 

    lines_for_legend = []
    labels_for_legend = []
    
    first_xvg_xlabel = None
    first_xvg_ylabel = None

    # --- Read and Plot XVG Data ---
    print(f"\n--- Processing {len(args.inputs)} XVG input file(s) ---")
    for i, xvg_filepath_str in enumerate(args.inputs):
        xvg_filepath = Path(xvg_filepath_str)
        print(f"Reading: {xvg_filepath}")
        xvg_x, xvg_y, xvg_title, xvg_xlabel_file, xvg_ylabel_file = read_xvg(xvg_filepath_str)

        if xvg_x is None or xvg_y is None:
            print(f"Skipping {xvg_filepath_str} due to read error.", file=sys.stderr)
            continue
        
        if i == 0: # Store labels from the first file as potential defaults
            first_xvg_xlabel = xvg_xlabel_file
            first_xvg_ylabel = xvg_ylabel_file

        # Attempt to create a sensible label from the path components
        # e.g., g1_rep1/analysis_output/rmsf_calpha_per_res.xvg -> g1_rep1 rmsf_calpha
        label_parts = []
        if xvg_filepath.parent.name != "analysis_output":
             label_parts.append(xvg_filepath.parent.name) # e.g. g1_rep1
        label_parts.append(xvg_filepath.stem.replace("_per_res", "").replace("_calpha"," Cα")) # e.g. rmsf Cα
        plot_label = " ".join(label_parts) if label_parts else xvg_filepath.name


        current_color = next(color_cycler)
        line, = ax1.plot(xvg_x, xvg_y, color=current_color, label=plot_label, alpha=0.8, linewidth=1.5)
        lines_for_legend.append(line)
        labels_for_legend.append(plot_label)
    
    if not lines_for_legend:
        print("No valid XVG data was plotted. Exiting.", file=sys.stderr)
        sys.exit(1)

    # Set primary axis labels
    ax1.set_xlabel(args.xlabel if args.xlabel else first_xvg_xlabel if first_xvg_xlabel else "Residue Number")
    ax1.set_ylabel(args.ylabel if args.ylabel else first_xvg_ylabel if first_xvg_ylabel else "RMSF (nm)", color='black') # Main Y-axis color
    ax1.tick_params(axis='y', labelcolor='black')
    ax1.grid(True, alpha=0.3, axis='x', linestyle=':') 
    ax1.grid(True, alpha=0.3, axis='y', linestyle=':', color='grey') 


    # --- B-factor Processing ---
    bfactor_x_for_plotting = None
    bfactor_y_for_plotting = None
    bfactor_legend_label = ""

    if args.bfac_pdb_id and args.bfac_chain:
        print(f"\n--- Processing B-factors for PDB {args.bfac_pdb_id}, Chain {args.bfac_chain} ---")
        awk_script_full_path = Path(args.awk_script_dir) / "extract_bfactors.awk"
        
        original_bfac_res, original_bfac_vals = get_bfactors(
            args.bfac_pdb_id, 
            args.bfac_chain,
            str(awk_script_full_path),
            args.pdb_download_dir
        )

        if original_bfac_res is not None and original_bfac_vals is not None and original_bfac_res.size > 0:
            bfactor_legend_label = f"B-factors ({args.bfac_pdb_id} Ch.{args.bfac_chain}"
            
            if args.residue_map_file:
                print(f"Using residue map file: {args.residue_map_file}")
                res_map_bfac_pdb_to_xvg = parse_residue_map_rmp(args.residue_map_file)

                if res_map_bfac_pdb_to_xvg:
                    mapped_res_for_xvg_axis = []
                    mapped_vals = []
                    unmapped_count = 0
                    for i, res_num_from_bfac_pdb in enumerate(original_bfac_res):
                        if res_num_from_bfac_pdb in res_map_bfac_pdb_to_xvg:
                            xvg_equiv_res_num = res_map_bfac_pdb_to_xvg[res_num_from_bfac_pdb]
                            mapped_res_for_xvg_axis.append(xvg_equiv_res_num)
                            mapped_vals.append(original_bfac_vals[i])
                        else:
                            unmapped_count +=1
                    
                    if mapped_res_for_xvg_axis:
                        bfactor_x_for_plotting = np.array(mapped_res_for_xvg_axis)
                        bfactor_y_for_plotting = np.array(mapped_vals)
                        print(f"Successfully mapped {len(mapped_vals)} B-factors using {args.residue_map_file}.")
                        if unmapped_count > 0:
                            print(f"Note: {unmapped_count} B-factor residues from {args.bfac_pdb_id} Ch.{args.bfac_chain} were not found in the map file.")
                        bfactor_legend_label += " mapped"
                    else:
                        print(f"Warning: No B-factors could be mapped using {args.residue_map_file}.")
                else:
                    print(f"Warning: Could not parse or found no data in residue map file. B-factors will not be mapped.")
            else: 
                print("No residue map file provided. Using original B-factor residue numbers for B-factor plot.")
                bfactor_x_for_plotting = original_bfac_res
                bfactor_y_for_plotting = original_bfac_vals
            
            if bfactor_legend_label and not bfactor_legend_label.endswith(")"): bfactor_legend_label += ")"
        else:
            print(f"Could not retrieve or process B-factors for {args.bfac_pdb_id} Chain {args.bfac_chain}.")
        print("--- End B-factor Processing ---\n")

    # Plot B-factors on secondary axis if available
    if bfactor_x_for_plotting is not None and bfactor_y_for_plotting is not None and bfactor_x_for_plotting.size > 0:
        ax2 = ax1.twinx()
        color_ax2 = 'dimgray' # A distinct color for B-factors
        # Using Unicode for Angstrom squared: Å²
        ax2_ylabel = f'B-factor (\u00C5\u00B2)'
        if bfactor_legend_label and "mapped" in bfactor_legend_label:
             ax2_ylabel += " mapped"

        ax2.set_ylabel(ax2_ylabel, color=color_ax2)
        line_bfac, = ax2.plot(bfactor_x_for_plotting, bfactor_y_for_plotting, color=color_ax2, linestyle=':', marker='o', markersize=2.5, label=bfactor_legend_label, alpha=0.7, linewidth=1.2)
        ax2.tick_params(axis='y', labelcolor=color_ax2)
        ax2.grid(True, alpha=0.2, axis='y', linestyle=':', color=color_ax2)
        
        lines_for_legend.append(line_bfac)
        labels_for_legend.append(bfactor_legend_label or "B-factors")
        
    # Finalize plot
    plot_title = args.title
    if args.bfac_pdb_id and "mapped" in bfactor_legend_label:
        plot_title += f" (B-factors mapped from {args.bfac_pdb_id} Ch.{args.bfac_chain})"
    elif args.bfac_pdb_id:
         plot_title += f" (B-factors from {args.bfac_pdb_id} Ch.{args.bfac_chain})"


    ax1.set_title(plot_title, fontsize='large', pad=15)

    if lines_for_legend:
        ax1.legend(lines_for_legend, labels_for_legend, loc='best', fontsize=args.legend_fontsize)

    fig.tight_layout()

    try:
        plt.savefig(args.output, dpi=args.dpi, bbox_inches='tight')
        print(f"Plot saved to {args.output}")
    except Exception as e:
        print(f"Error saving plot to {args.output}: {e}", file=sys.stderr)

    # plt.show() # Uncomment to display plot interactively

if __name__ == "__main__":
    main()

