#!/usr/bin/env python3

import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import itertools # For color cycling
import urllib.request

# --- Utility and Data Handling Functions ---

def download_pdb_file(pdb_id, pdb_dir="."):
    pdb_dir_path = Path(pdb_dir)
    pdb_dir_path.mkdir(parents=True, exist_ok=True)
    pdb_file = pdb_dir_path / f"{pdb_id.lower()}.pdb"
    
    if pdb_file.exists():
        print(f"INFO: PDB file {pdb_file} already exists.")
        return pdb_file

    print(f"INFO: Attempting to download {pdb_id.upper()}.pdb to {pdb_file}...")
    try:
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        urllib.request.urlretrieve(url, pdb_file)
        print(f"INFO: Successfully downloaded {pdb_file}")
        return pdb_file
    except Exception as e:
        print(f"ERROR: Could not download PDB file {pdb_id}.pdb: {e}", file=sys.stderr)
        print(f"Please ensure {pdb_id.lower()}.pdb is manually placed in {pdb_dir_path} or check internet connection.", file=sys.stderr)
        return None

def get_bfactors(pdb_id, chain_id, awk_script_path_str, pdb_download_dir="."):
    """
    Extracts B-factors for a specific chain from a PDB file.
    Returns original PDB residue numbers for that chain and their B-factors.
    """
    pdb_file_path = download_pdb_file(pdb_id, pdb_download_dir)
    if not pdb_file_path:
        return None, None

    awk_script_path = Path(awk_script_path_str)
    if not awk_script_path.exists():
        print(f"ERROR: AWK script not found at {awk_script_path}", file=sys.stderr)
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
            print(f"WARNING: No B-factor data extracted for PDB {pdb_id}, chain {chain_id}.")
            return np.array([]), np.array([])

        data = np.array([line.split() for line in output.splitlines() if line.strip()], dtype=float)
        
        if data.size == 0 or data.ndim != 2 or data.shape[1] != 2:
            print(f"WARNING: B-factor data for {pdb_id} chain {chain_id} is empty or not in expected 2-column format.")
            return np.array([]), np.array([])
            
        print(f"INFO: Successfully extracted {data.shape[0]} B-factors for {pdb_id} chain {chain_id} (original PDB numbering).")
        return data[:, 0].astype(int), data[:, 1] # Col 0: ResNum, Col 1: B-factor
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Error executing AWK script for {pdb_id} chain {chain_id}: {e}", file=sys.stderr)
        print(f"AWK stderr: {e.stderr}", file=sys.stderr)
        return None, None
    except subprocess.TimeoutExpired:
        print(f"ERROR: Timeout executing AWK script for {pdb_id} chain {chain_id}.", file=sys.stderr)
        return None, None
    except Exception as e:
        print(f"ERROR: An unexpected error occurred in get_bfactors for {pdb_id} chain {chain_id}: {e}", file=sys.stderr)
        return None, None

def parse_residue_map_crmp(map_file_path_str):
    """
    Parses a .crmp (chain-explicit residue map) file.
    The map file format is: RefChainID RefResName RefResNum InputChainID InputResName InputResNum
    Returns a dictionary: {(RefChainID, RefPDBResNum): (InputChainID, SimResNum)}
    """
    mapping = {}
    map_file_path = Path(map_file_path_str)
    if not map_file_path.exists():
        print(f"ERROR: Residue map file not found: {map_file_path}", file=sys.stderr)
        return None
    
    try:
        with open(map_file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                
                parts = line.split()
                if len(parts) == 6: 
                    try:
                        ref_chain_id = parts[0]
                        # ref_res_name = parts[1] # Not strictly needed for mapping key/value
                        ref_res_num = int(parts[2])
                        input_chain_id = parts[3]
                        # input_res_name = parts[4] # Not strictly needed
                        sim_res_num = int(parts[5])
                        
                        mapping[(ref_chain_id, ref_res_num)] = (input_chain_id, sim_res_num)
                    except ValueError:
                        print(f"WARNING: Skipping malformed numeric data in map file {map_file_path} at line {line_num}: '{line}'", file=sys.stderr)
                else:
                    print(f"WARNING: Skipping malformed line (expected 6 columns, got {len(parts)}) in map file {map_file_path} at line {line_num}: '{line}'", file=sys.stderr)
    except Exception as e:
        print(f"ERROR: Error reading or parsing residue map file {map_file_path}: {e}", file=sys.stderr)
        return None
    
    if not mapping:
        print(f"WARNING: No valid mappings found in {map_file_path}.")
    else:
        print(f"INFO: Successfully parsed {len(mapping)} chain-explicit mappings from {map_file_path}.")
    return mapping

def read_xvg(filename, x_col=0, y_col=1):
    x_data, y_data = [], []
    title, xlabel, ylabel = Path(filename).stem, "X-axis", "Y-axis" # Defaults
    
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
                    print(f"WARNING: Skipping malformed data line in {filename}: '{line}'", file=sys.stderr)
                    continue
        
        if not x_data or not y_data:
             print(f"ERROR: No data columns could be read from {filename}.", file=sys.stderr)
             return None, None, title, xlabel, ylabel

        return np.array(x_data), np.array(y_data), title, xlabel, ylabel
    except FileNotFoundError:
        print(f"ERROR: XVG file not found: {filename}", file=sys.stderr)
        return None, None, title, xlabel, ylabel
    except Exception as e:
        print(f"ERROR: Error reading XVG file {filename}: {e}", file=sys.stderr)
        return None, None, title, xlabel, ylabel

# --- Main Plotting Logic ---
def main():
    parser = argparse.ArgumentParser(
        description="Plot multiple XVG files, optionally overlaying B-factors from PDB chains on a secondary Y-axis. "
                    "Uses chain-explicit .crmp residue map files for aligning B-factor residue numbers to simulation numbering.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--inputs', required=True, type=str, nargs='+', help="One or more input XVG file paths.")
    parser.add_argument('--output', type=str, default="plot.png", help="Output plot file name.")
    
    parser.add_argument('--title', type=str, default="XVG Data Plot", help="Overall plot title.")
    parser.add_argument('--xlabel', type=str, help="X-axis label. Uses label from the first XVG if not provided.")
    parser.add_argument('--ylabel', type=str, help="Y-axis label (primary data - e.g., RMSF). Uses label from the first XVG if not provided.")
    
    parser.add_argument('--bfac_pdb_id', type=str, help="PDB ID for B-factor extraction (e.g., 1R1K). This is the 'reference' PDB.")
    parser.add_argument('--bfac_chains', type=str, nargs='+', 
                        help="One or more chain IDs from the --bfac_pdb_id (reference PDB) for B-factor extraction (e.g., A D).")
    parser.add_argument('--residue_map_files', type=str, nargs='+', 
                        help="Residue map file(s) in .crmp format. Provide one map file to be used for all specified --bfac_chains, "
                             "or one map file per chain (must match order and number of --bfac_chains). "
                             "If omitted, B-factors use original PDB residue numbers from --bfac_pdb_id.")
    
    parser.add_argument('--awk_script_dir', type=str, default=".", help="Directory of extract_bfactors.awk.")
    parser.add_argument('--pdb_download_dir', type=str, default=".", help="Directory for PDB files.")

    parser.add_argument('--figsize', type=str, default="18,10", help="Figure size (width,height), e.g., '12,7'.")
    parser.add_argument('--dpi', type=int, default=300, help="DPI for raster output.")
    parser.add_argument('--legend_fontsize', type=str, default="x-small", help="Legend font size.")

    args = parser.parse_args()

    try:
        figsize_vals = tuple(map(int, args.figsize.split(',')))
        if len(figsize_vals) != 2: raise ValueError("Figsize must be two comma-separated integers.")
    except ValueError as e:
        print(f"Invalid figsize format '{args.figsize}'. Using default (18,10). Error: {e}", file=sys.stderr)
        figsize_vals = (18, 10)

    fig, ax1 = plt.subplots(figsize=figsize_vals)
    color_cycler_rmsf = itertools.cycle(plt.cm.tab10.colors) 
    lines_for_legend, labels_for_legend = [], []
    first_xvg_xlabel, first_xvg_ylabel = None, None

    print(f"\n--- Processing {len(args.inputs)} XVG input file(s) ---")
    for i, xvg_filepath_str in enumerate(args.inputs):
        xvg_filepath = Path(xvg_filepath_str)
        xvg_x, xvg_y, _, xvg_xlabel_file, xvg_ylabel_file = read_xvg(xvg_filepath_str)

        if xvg_x is None or xvg_y is None:
            print(f"Skipping {xvg_filepath_str} due to read error.", file=sys.stderr)
            continue
        
        if i == 0: 
            first_xvg_xlabel = xvg_xlabel_file
            first_xvg_ylabel = xvg_ylabel_file

        label_parts = [xvg_filepath.parent.name] if xvg_filepath.parent.name != "." and xvg_filepath.parent.name != Path.cwd().name else []
        stem_name = xvg_filepath.stem.replace("_per_res", "").replace("_calpha"," Cα") # Using actual alpha
        label_parts.append(stem_name)
        plot_label = " ".join(label_parts) if label_parts else xvg_filepath.name

        line, = ax1.plot(xvg_x, xvg_y, color=next(color_cycler_rmsf), label=plot_label, alpha=0.8, linewidth=1.5)
        lines_for_legend.append(line)
        labels_for_legend.append(plot_label)
    
    if not lines_for_legend:
        print("ERROR: No valid XVG data could be plotted. Exiting.", file=sys.stderr)
        sys.exit(1)

    ax1.set_xlabel(args.xlabel or first_xvg_xlabel or "Residue Number (Simulation)")
    ax1.set_ylabel(args.ylabel or first_xvg_ylabel or "RMSF (nm)", color='black')
    ax1.tick_params(axis='y', labelcolor='black')
    ax1.grid(True, alpha=0.3, axis='x', linestyle=':') 
    ax1.grid(True, alpha=0.3, axis='y', linestyle=':', color='grey') 

    ax2 = None 
    any_bfactors_plotted = False
    plot_title = args.title
    bfactor_info_for_title = []
    
    # --- B-factor and Mapping Setup ---
    shared_map_data = None
    use_shared_map = False
    use_individual_maps = False

    if args.bfac_pdb_id and args.bfac_chains and args.residue_map_files:
        if len(args.residue_map_files) == 1:
            print(f"INFO: Using shared residue map file: {args.residue_map_files[0]} for all specified reference chains.")
            shared_map_data = parse_residue_map_crmp(args.residue_map_files[0]) # Use new parser
            if shared_map_data is None:
                print(f"WARNING: Could not parse the shared map file {args.residue_map_files[0]}. B-factors will use original PDB residue numbers if plotted.", file=sys.stderr)
            else:
                use_shared_map = True
        elif len(args.residue_map_files) == len(args.bfac_chains):
            print("INFO: Using individual residue map file for each specified reference chain.")
            use_individual_maps = True
        else:
            print("ERROR: Number of --residue_map_files must be 1 (for a shared map) or match the number of --bfac_chains.", file=sys.stderr)
            print("B-factors will be plotted with original PDB residue numbers if possible.", file=sys.stderr)

    if args.bfac_pdb_id and args.bfac_chains:
        print(f"\n--- Processing B-factors for Reference PDB {args.bfac_pdb_id}, Chains: {', '.join(args.bfac_chains)} ---")
        awk_script_full_path = Path(args.awk_script_dir) / "extract_bfactors.awk"
        
        bfactor_styles = [ # Define styles for B-factor plots
            {'color': 'dimgray', 'linestyle': ':', 'marker': '.'}, 
            {'color': 'darkslateblue', 'linestyle': '--', 'marker': '.'},
            {'color': 'darkgreen', 'linestyle': '-.', 'marker': '.'},
            {'color': 'maroon', 'linestyle': ':', 'marker': '.'},
        ] 
        common_style_params = {'markersize': 3.5, 'alpha': 0.7, 'linewidth': 1.2}
        style_cycler_bfac = itertools.cycle(bfactor_styles)

        for i, ref_chain_id in enumerate(args.bfac_chains): # ref_chain_id is from the reference PDB
            print(f"-- Processing Reference Chain {ref_chain_id} --")
            original_ref_pdb_res_nums, original_bfac_vals = get_bfactors(
                args.bfac_pdb_id, ref_chain_id, str(awk_script_full_path), args.pdb_download_dir
            )

            if original_ref_pdb_res_nums is None or original_bfac_vals is None or original_ref_pdb_res_nums.size == 0:
                print(f"WARNING: Could not retrieve or process B-factors for {args.bfac_pdb_id} Chain {ref_chain_id}.")
                continue

            bfactor_x_for_plotting = original_ref_pdb_res_nums # Default to original if no map or mapping fails
            bfactor_y_for_plotting = original_bfac_vals
            chain_mapped_status = "orig. PDB res." 
            current_map_to_use = None
            map_file_name_for_log = "N/A (no map provided or applicable)"

            if use_shared_map and shared_map_data is not None:
                current_map_to_use = shared_map_data
                map_file_name_for_log = args.residue_map_files[0]
            elif use_individual_maps:
                map_file_path_str = args.residue_map_files[i]
                map_file_name_for_log = map_file_path_str
                print(f"INFO: Using individual map file: {map_file_path_str} for reference chain {ref_chain_id}")
                current_map_to_use = parse_residue_map_crmp(map_file_path_str) # Use new parser
                if current_map_to_use is None:
                    print(f"WARNING: Could not parse map file {map_file_path_str} for reference chain {ref_chain_id}. Plotting with original residue numbers.")
            
            if current_map_to_use:
                mapped_sim_res_nums, mapped_bfac_vals, unmapped_count = [], [], 0
                for res_idx, ref_res_num_from_pdb in enumerate(original_ref_pdb_res_nums):
                    # Key for crmp map is (RefChainID, RefPDBResNum)
                    map_key = (ref_chain_id, ref_res_num_from_pdb) 
                    mapped_info = current_map_to_use.get(map_key)
                    
                    if mapped_info:
                        _input_chain_id, sim_res_num = mapped_info # unpack, we need sim_res_num for x-axis
                        mapped_sim_res_nums.append(sim_res_num)
                        mapped_bfac_vals.append(original_bfac_vals[res_idx])
                    else:
                        unmapped_count +=1
                
                if mapped_sim_res_nums:
                    bfactor_x_for_plotting = np.array(mapped_sim_res_nums)
                    bfactor_y_for_plotting = np.array(mapped_bfac_vals)
                    chain_mapped_status = "mapped to sim."
                    print(f"INFO: Successfully mapped {len(mapped_bfac_vals)} B-factors for ref chain {ref_chain_id} using map from {Path(map_file_name_for_log).name}.")
                    if unmapped_count > 0:
                        print(f"NOTE: {unmapped_count} B-factor residues from ref chain {ref_chain_id} were not found in map from {Path(map_file_name_for_log).name}.")
                else:
                    print(f"WARNING: No B-factors could be mapped for ref chain {ref_chain_id} using map from {Path(map_file_name_for_log).name}. Plotting with original PDB residue numbers.")
            
            bfactor_legend_label = f"B-factors ({args.bfac_pdb_id} Ch.{ref_chain_id}"
            if args.residue_map_files : 
                 bfactor_legend_label += f" {chain_mapped_status}"
            bfactor_legend_label += ")"
            
            bfactor_info_for_title.append(f"Ch.{ref_chain_id} ({chain_mapped_status if args.residue_map_files else 'orig. PDB res.'})")

            if bfactor_x_for_plotting.size > 0:
                if not ax2: 
                    ax2 = ax1.twinx()
                    ax2_ylabel_base = f'B-factor (\u00C5\u00B2)' # Å²
                    ax2.set_ylabel(ax2_ylabel_base) 
                    ax2.grid(True, alpha=0.2, axis='y', linestyle=':', color='lightgrey')

                current_style = next(style_cycler_bfac)
                # ax2.tick_params(axis='y', labelcolor=current_style['color']) # Color y-axis ticks per B-factor line can be messy

                line_bfac, = ax2.plot(bfactor_x_for_plotting, bfactor_y_for_plotting, 
                                      label=bfactor_legend_label, 
                                      **current_style, **common_style_params)
                lines_for_legend.append(line_bfac)
                labels_for_legend.append(bfactor_legend_label)
                any_bfactors_plotted = True
        
        if any_bfactors_plotted:
            if bfactor_info_for_title:
                title_suffix = f" (B-factors: {args.bfac_pdb_id} " + ", ".join(bfactor_info_for_title) + ")"
                # Avoid appending if default title is used and it already contains "(B-factors..."
                if not (args.title == parser.get_default("title") and "(B-factors" in args.title):
                     plot_title += title_suffix
                elif args.title != parser.get_default("title"): # User provided a custom title
                     plot_title += title_suffix
            if ax2: # Ensure ax2 was created before setting its tick color (if desired for the last plotted B-factor line)
                ax2.tick_params(axis='y', labelcolor='dimgray') # General color for B-factor axis


    ax1.set_title(plot_title, fontsize='medium')
    
    if lines_for_legend:
        ax1.legend(lines_for_legend, labels_for_legend, loc='best', fontsize=args.legend_fontsize)
    
    plt.tight_layout()
    
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        plt.savefig(output_path, dpi=args.dpi, bbox_inches='tight')
        print(f"\nPlot saved to {output_path}")
    except Exception as e:
        print(f"ERROR: Error saving plot to {output_path}: {e}", file=sys.stderr)

if __name__ == '__main__':
    main()

