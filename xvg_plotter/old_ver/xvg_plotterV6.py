#!/usr/bin/env python3

import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# --- Utility and Data Handling Functions ---

def download_pdb_file(pdb_id, pdb_dir="."):
    """
    Downloads a PDB file from RCSB PDB if it doesn't already exist.
    Basic implementation, consider using 'requests' for robustness in a real scenario.
    """
    pdb_dir_path = Path(pdb_dir)
    pdb_dir_path.mkdir(parents=True, exist_ok=True)
    # PDB files are often named in lowercase, but IDs are uppercase.
    # RCSB download links use uppercase PDB IDs.
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
    """
    Extracts B-factors for a given PDB ID and chain using an AWK script.
    Returns (residue_numbers, bfactor_values) or (None, None) on error.
    """
    pdb_file_path = download_pdb_file(pdb_id, pdb_download_dir)
    if not pdb_file_path:
        return None, None

    awk_script_path = Path(awk_script_path_str)
    if not awk_script_path.exists():
        print(f"Error: AWK script not found at {awk_script_path}")
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
        return data[:, 0].astype(int), data[:, 1]  # Residue numbers as int, B-factors as float
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
    """
    Parses a .rmp file (ResName PDB1_ResNum PDB2_ResNum) into a dictionary
    mapping PDB1_ResNum (from B-factor PDB) -> PDB2_ResNum (from XVG/simulation PDB).
    """
    mapping = {}
    map_file_path = Path(map_file_path_str)
    if not map_file_path.exists():
        print(f"Error: Residue map file not found: {map_file_path}", file=sys.stderr)
        return None
    
    try:
        with open(map_file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                parts = line.strip().split()
                if not parts or line.startswith("#"):  # Skip empty lines and comments
                    continue
                if len(parts) == 3: # ResName PDB1_ResNum PDB2_ResNum
                    try:
                        pdb1_res_num = int(parts[1]) # Key: residue number from B-factor PDB
                        pdb2_res_num = int(parts[2]) # Value: residue number for XVG x-axis
                        mapping[pdb1_res_num] = pdb2_res_num
                    except ValueError:
                        print(f"Warning: Skipping malformed numeric data in map file {map_file_path} at line {line_num}: {line.strip()}")
                else:
                    print(f"Warning: Skipping malformed line (expected 3 columns, got {len(parts)}) in map file {map_file_path} at line {line_num}: {line.strip()}")
    except Exception as e:
        print(f"Error reading or parsing residue map file {map_file_path}: {e}", file=sys.stderr)
        return None
    
    if not mapping:
        print(f"Warning: No valid mappings found in {map_file_path}. The file might be empty or incorrectly formatted.")
    else:
        print(f"Successfully parsed {len(mapping)} mappings from {map_file_path}.")
    return mapping


def read_xvg(filename, x_col=0, y_col=1):
    """
    Basic XVG reader. Assumes comments start with '#' or '@'.
    Returns x_data, y_data, title, xlabel, ylabel.
    """
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
                    print(f"Warning: Skipping malformed data line in {filename}: {line}")
                    continue
        
        if not x_data or not y_data:
             print(f"Error: No data columns could be read from {filename}.", file=sys.stderr)
             return None, None, title, xlabel, ylabel

        return np.array(x_data), np.array(y_data), title, xlabel, ylabel
    except FileNotFoundError:
        print(f"Error: XVG file not found: {filename}", file=sys.stderr)
        sys.exit(1) # Exit if primary input file is not found
    except Exception as e:
        print(f"Error reading XVG file {filename}: {e}", file=sys.stderr)
        sys.exit(1)


# --- Main Plotting Logic ---
def main():
    parser = argparse.ArgumentParser(
        description="Plot XVG files, optionally overlaying B-factors with residue mapping.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input', required=True, type=str, help="Input XVG file path.")
    parser.add_argument('--output', type=str, default="plot.png", help="Output plot file name (e.g., plot.png, plot.svg).")
    
    parser.add_argument('--title', type=str, help="Overall plot title. Overrides XVG title if provided.")
    parser.add_argument('--xlabel', type=str, help="X-axis label. Overrides XVG label if provided.")
    parser.add_argument('--ylabel', type=str, help="Y-axis label (for primary data). Overrides XVG label if provided.")
    
    parser.add_argument('--bfac_pdb_id', type=str, help="PDB ID for B-factor extraction (e.g., 1R1K).")
    parser.add_argument('--bfac_chain', type=str, help="Chain ID in the PDB for B-factor extraction (e.g., A).")
    parser.add_argument('--residue_map_file', type=str, help="Path to the .rmp residue map file (e.g., 1r1k_A_to_sim_A_map.rmp).")
    
    parser.add_argument('--awk_script_dir', type=str, default=".", help="Directory containing extract_bfactors.awk.")
    parser.add_argument('--pdb_download_dir', type=str, default=".", help="Directory to store/find PDB files.")

    parser.add_argument('--figsize', type=str, default="12,7", help="Figure size (width,height) in inches, e.g., '10,6'.")
    parser.add_argument('--dpi', type=int, default=300, help="DPI for raster image output.")

    args = parser.parse_args()

    # --- Read XVG Data ---
    xvg_x, xvg_y, xvg_title, xvg_xlabel, xvg_ylabel = read_xvg(args.input)
    if xvg_x is None:
        print(f"Could not read data from {args.input}. Exiting.", file=sys.stderr)
        sys.exit(1)

    # Override XVG metadata if command line args are given
    plot_title = args.title if args.title else xvg_title
    plot_xlabel = args.xlabel if args.xlabel else xvg_xlabel
    plot_ylabel = args.ylabel if args.ylabel else xvg_ylabel

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
                # The map should map B-factor PDB residue numbers (keys) to XVG residue numbers (values)
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
                        print(f"Warning: No B-factors could be mapped using {args.residue_map_file}. Check map content and PDB residue numbers.")
                else:
                    print(f"Warning: Could not parse or found no data in residue map file. B-factors will not be mapped.")
            else: # No map file provided, use original B-factor residue numbers (might not align with XVG x-axis)
                print("No residue map file provided. Using original B-factor residue numbers. This might not align correctly with the XVG x-axis if residue numberings differ.")
                bfactor_x_for_plotting = original_bfac_res
                bfactor_y_for_plotting = original_bfac_vals
            
            if bfactor_legend_label: bfactor_legend_label += ")"
        else:
            print(f"Could not retrieve or process B-factors for {args.bfac_pdb_id} Chain {args.bfac_chain}.")
        print("--- End B-factor Processing ---\n")

    # --- Plotting ---
    try:
        figsize_vals = tuple(map(int, args.figsize.split(',')))
        if len(figsize_vals) != 2: raise ValueError("Figsize must be two comma-separated integers.")
    except ValueError as e:
        print(f"Invalid figsize format '{args.figsize}'. Using default (12,7). Error: {e}")
        figsize_vals = (12, 7)

    fig, ax1 = plt.subplots(figsize=figsize_vals)

    # Plot primary XVG data
    color_ax1 = 'tab:red'
    ax1.set_xlabel(plot_xlabel)
    ax1.set_ylabel(plot_ylabel, color=color_ax1)
    line1, = ax1.plot(xvg_x, xvg_y, color=color_ax1, label=plot_ylabel or Path(args.input).stem)
    ax1.tick_params(axis='y', labelcolor=color_ax1)
    ax1.grid(True, alpha=0.3, axis='x', linestyle=':') 
    ax1.grid(True, alpha=0.3, axis='y', linestyle=':', color=color_ax1) 

    lines_for_legend = [line1]
    labels_for_legend = [line1.get_label()]

    # Plot B-factors on secondary axis if available
    if bfactor_x_for_plotting is not None and bfactor_y_for_plotting is not None and bfactor_x_for_plotting.size > 0:
        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
        color_ax2 = 'tab:blue'
        # Using Unicode for Angstrom squared: Å²
        ax2.set_ylabel(f'B-factor (\u00C5\u00B2){" - " + bfactor_legend_label if bfactor_legend_label else ""}', color=color_ax2)
        line2, = ax2.plot(bfactor_x_for_plotting, bfactor_y_for_plotting, color=color_ax2, linestyle='--', marker='.', markersize=3, label=bfactor_legend_label)
        ax2.tick_params(axis='y', labelcolor=color_ax2)
        ax2.grid(True, alpha=0.3, axis='y', linestyle=':', color=color_ax2)
        
        lines_for_legend.append(line2)
        labels_for_legend.append(line2.get_label())
        
        # Adjust title if B-factors are plotted
        if args.bfac_pdb_id:
             plot_title += f" with B-factors ({args.bfac_pdb_id} Ch.{args.bfac_chain}{' mapped' if args.residue_map_file and 'mapped' in bfactor_legend_label else ''})"


    ax1.set_title(plot_title, fontsize='large', pad=15) # Add some padding to title

    # Unified legend
    if lines_for_legend:
        # Position legend to avoid overlap as much as possible
        ax1.legend(lines_for_legend, labels_for_legend, loc='best', fontsize='small')

    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    try:
        plt.savefig(args.output, dpi=args.dpi, bbox_inches='tight')
        print(f"Plot saved to {args.output}")
    except Exception as e:
        print(f"Error saving plot to {args.output}: {e}", file=sys.stderr)

    # plt.show() # Uncomment to display plot interactively

if __name__ == "__main__":
    main()

