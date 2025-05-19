#!/usr/bin/env python3
"""
XVG File Plotter - A utility to visualize GROMACS .xvg files

Usage:
    python xvg_plotter.py --input path/to/file1.xvg [path/to/file2.xvg ...] [--noline] [--gaussian] [--no-show] [--output filename.png]
"""

import sys
import re
import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path
from scipy import stats

def parse_xvg_file(filepath):
    """
    Parse an XVG file and extract data and metadata.
    Handles lines with trailing non-numeric text (e.g., Ramachandran plots).

    Args:
        filepath: Path to the .xvg file

    Returns:
        tuple: (data_array, y_legends, title, xlabel, ylabel, filename_stem)
               y_legends is a list of legends for the Y-axes.
    """
    filepath_obj = Path(filepath)
    filename_stem = filepath_obj.stem
    with open(filepath_obj, 'r') as f:
        lines = f.readlines()

    raw_data = []
    y_legends_dict = {}
    title = "GROMACS Data"
    xlabel = "X-axis" # Default X-axis label
    ylabel = "Y-axis" # Default Y-axis label

    for line in lines:
        line = line.strip()
        if not line:
            continue

        if line.startswith('@'):
            if line.startswith('@    title'):
                title = line.split('"')[1] if '"' in line else line.split()[-1]
            elif line.startswith('@    xaxis  label'):
                xlabel = line.split('"')[1] if '"' in line else line.split()[-1]
            elif line.startswith('@    yaxis  label'):
                ylabel = line.split('"')[1] if '"' in line else line.split()[-1]
            else:
                legend_match = re.match(r'@\s+s(\d+)\s+legend\s+"([^"]+)"', line)
                if legend_match:
                    s_index = int(legend_match.group(1))
                    legend_text = legend_match.group(2)
                    y_legends_dict[s_index] = legend_text
        elif not line.startswith('#'):
            parts = line.split()
            numeric_values_for_line = []
            for part in parts:
                try:
                    numeric_values_for_line.append(float(part))
                except ValueError:
                    break
            if numeric_values_for_line:
                raw_data.append(numeric_values_for_line)

    if not raw_data:
        print(f"Warning: No data lines found or parsed successfully in {filepath}.")
        return np.array([]), [], title, xlabel, ylabel, filename_stem

    num_cols_first_row = len(raw_data[0])
    consistent_raw_data = [row for row in raw_data if len(row) == num_cols_first_row]
    if len(consistent_raw_data) != len(raw_data):
        print(f"Warning: Some data lines in {filepath} had an inconsistent number of numeric columns. "
              f"Using {len(consistent_raw_data)} rows out of {len(raw_data)} that matched the first row's column count ({num_cols_first_row}).")
    raw_data = consistent_raw_data

    if not raw_data:
        print(f"Warning: No consistent data lines found after filtering for column consistency in {filepath}.")
        return np.array([]), [], title, xlabel, ylabel, filename_stem

    data_array = np.array(raw_data)

    final_y_legends = []
    if data_array.ndim == 2 and data_array.shape[1] > 1:
        for y_col_idx_in_data in range(1, data_array.shape[1]): # For XVG, data usually starts with X, then Y1, Y2...
            legend_key_for_dict = y_col_idx_in_data - 1 # s0 usually corresponds to the first Y column
            fallback_legend_name = f"Y-Col {y_col_idx_in_data} (s{legend_key_for_dict})"
            final_y_legends.append(y_legends_dict.get(legend_key_for_dict, fallback_legend_name))
    elif data_array.ndim == 1 and data_array.size > 0 : # Should not happen for typical XVG with X, Y data
        print(f"Warning: Data in {filepath} is 1D after parsing. Shape: {data_array.shape}. Cannot plot as X,Y.")
        return np.array([]), [], title, xlabel, ylabel, filename_stem
    elif data_array.size == 0 :
        print(f"Warning: No valid data points resulted after parsing {filepath}. Final data_array shape: {data_array.shape}")
        return np.array([]), [], title, xlabel, ylabel, filename_stem
    
    if not final_y_legends and data_array.ndim == 2 and data_array.shape[1] > 1: # if no s_legends, create one for the first Y column
        final_y_legends.append(f"{filename_stem} Y-data")


    return data_array, final_y_legends, title, xlabel, ylabel, filename_stem

def analyze_gaussian(data_column):
    if len(data_column) < 2:
        return ({'mu': np.nan, 'sigma': np.nan},
                {'mu': np.nan, 'sigma': np.nan},
                {'mu': np.nan, 'sigma': np.nan})
    mid_point = len(data_column) // 2
    first_half = data_column[:mid_point] if mid_point > 0 else np.array([])
    second_half = data_column[mid_point:] if mid_point > 0 else np.array([])
    mu1, sigma1 = stats.norm.fit(first_half) if len(first_half) > 1 else (np.nan, np.nan)
    mu2, sigma2 = stats.norm.fit(second_half) if len(second_half) > 1 else (np.nan, np.nan)
    mu_total, sigma_total = stats.norm.fit(data_column) if len(data_column) > 1 else (np.nan, np.nan)
    return ({'mu': mu1, 'sigma': sigma1}, {'mu': mu2, 'sigma': sigma2}, {'mu': mu_total, 'sigma': sigma_total})

def plot_with_gaussian(data, y_legends, title, xlabel, ylabel, input_filename_stem):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=False, height_ratios=[2, 1])
    for y_col_idx in range(1, data.shape[1]):
        legend_idx = y_col_idx - 1
        legend_label = y_legends[legend_idx] if legend_idx < len(y_legends) else f"Y-Col {y_col_idx}"
        ax1.plot(data[:, 0], data[:, y_col_idx], label=legend_label)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_title(f"{title} ({input_filename_stem}) - Time Series")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    for y_col_idx in range(1, data.shape[1]):
        column_data = data[:, y_col_idx]
        legend_idx = y_col_idx - 1
        legend_label = y_legends[legend_idx] if legend_idx < len(y_legends) else f"Y-Col {y_col_idx}"
        if len(column_data) < 2: continue
        stats_first, stats_second, stats_total = analyze_gaussian(column_data)
        ax2.hist(column_data, bins=50, density=True, alpha=0.3, label=f'Dist: {legend_label}')
        min_val, max_val = np.min(column_data), np.max(column_data)
        x_smooth = np.linspace(min_val, max_val, 200) if min_val != max_val else np.array([min_val])
        if x_smooth.size > 0:
            if not np.isnan(stats_first['mu']) and not np.isnan(stats_first['sigma']) and stats_first['sigma'] > 1e-6:
                 ax2.plot(x_smooth, stats.norm.pdf(x_smooth, stats_first['mu'], stats_first['sigma']), '--',
                         label=f'{legend_label} First Half (μ={stats_first["mu"]:.2f}, σ={stats_first["sigma"]:.2f})')
            if not np.isnan(stats_second['mu']) and not np.isnan(stats_second['sigma']) and stats_second['sigma'] > 1e-6:
                ax2.plot(x_smooth, stats.norm.pdf(x_smooth, stats_second['mu'], stats_second['sigma']), ':',
                         label=f'{legend_label} Second Half (μ={stats_second["mu"]:.2f}, σ={stats_second["sigma"]:.2f})')
    ax2.set_xlabel(ylabel) # X-axis of histogram is the value from Y-axis of top plot
    ax2.set_ylabel('Probability Density')
    ax2.set_title('Distribution Analysis')
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
    ax2.grid(True, alpha=0.3)
    plt.tight_layout(rect=[0, 0, 0.9, 1])
    return fig

def plot_ramachandran(data, title, xlabel, ylabel, input_filename_path, output_path=None, no_show=False):
    if data.shape[1] < 2: # Need at least Phi and Psi columns
        print(f"Error: Ramachandran plot requires at least 2 data columns (e.g., Phi, Psi), got {data.shape[1]}. Skipping {input_filename_path.name}")
        return
    # Assuming first column is phi, second is psi if more than 2 cols (e.g. X, phi, psi)
    # Or if just 2 cols, then they are phi, psi
    phi_idx, psi_idx = (0, 1) if data.shape[1] == 2 else (data.shape[1]-2, data.shape[1]-1) # Heuristic: take last two if more than 2
    if "phi" in xlabel.lower() and "psi" in ylabel.lower() and data.shape[1] > 2: # If labels are explicit and there's an X col
         phi_idx, psi_idx = 0, 1 # If xvg has x, phi, psi, but labels are for phi, psi
    
    # If data has 3 columns (time, phi, psi), use last two. If 2 (phi, psi) use those.
    col_offset = 1 if data.shape[1] > 2 else 0 # If there's a time/index column
    phi_angles = data[:, col_offset]
    psi_angles = data[:, col_offset+1]

    plt.figure(figsize=(8, 7))
    final_xlabel = xlabel if "phi" in xlabel.lower() else "phi (degrees)"
    final_ylabel = ylabel if "psi" in ylabel.lower() else "psi (degrees)"
    
    counts, xedges, yedges, im = plt.hist2d(phi_angles, psi_angles, bins=120, cmap='viridis', cmin=1, range=[[-180, 180], [-180, 180]])
    plt.colorbar(im, label='Density', shrink=0.8)
    plt.xlabel(final_xlabel)
    plt.ylabel(final_ylabel)
    plot_title = f"Ramachandran Plot: {input_filename_path.name}"
    if title and title != "GROMACS Data" and "Ramachandran" not in title : # Use file's title if more specific
        plot_title = f"{title} ({input_filename_path.name})"
    plt.title(plot_title)
    plt.xlim(-180, 180); plt.ylim(-180, 180)
    plt.axhline(0, color='grey', lw=0.5, linestyle='--'); plt.axvline(0, color='grey', lw=0.5, linestyle='--')
    plt.xticks(np.arange(-180, 181, 60)); plt.yticks(np.arange(-180, 181, 60))
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout(rect=[0, 0, 0.95, 1]) # Adjust rect to make space for colorbar legend
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Ramachandran plot saved to {output_path}")
    if not no_show:
        plt.show()
    plt.close()

def main():
    parser = argparse.ArgumentParser(
        description='Plot XVG files with optional Gaussian analysis, Ramachandran handling, and multi-file scatter overlays.',
        formatter_class=argparse.RawTextHelpFormatter # To allow newlines in help text
    )
    parser.add_argument('--input', '-i', type=str, nargs='+', help='Path(s) to the XVG file(s). \nMultiple files can be specified for scatter plots to overlay them.', required=False)
    parser.add_argument('--gaussian', '-g', action='store_true', help='Enable Gaussian distribution analysis (for single time series data)')
    parser.add_argument('--noline', action='store_true', help='Generate a scatter plot. If multiple inputs, overlays them.')
    parser.add_argument('--output', '-o', type=str, help='Output file path (optional). \nFor multi-file scatter, this is the direct output name.')
    parser.add_argument('--no-show', action='store_true', help='Do not display the plot (save only)')
    parser.add_argument('filepaths', nargs='*', type=str, help=argparse.SUPPRESS) # For positional arguments

    args = parser.parse_args()

    input_files = []
    if args.input:
        input_files.extend(args.input)
    if args.filepaths: # Positional arguments
        input_files.extend(args.filepaths)
    
    if not input_files:
        parser.error("No input files provided. Use --input or positional arguments.")
        return

    # --- Multi-file scatter plot (--noline and >1 input) ---
    if args.noline and len(input_files) > 1:
        print(f"Generating combined scatter plot for: {', '.join(input_files)}")
        plt.figure(figsize=(10, 7)) # Adjusted size for potentially more legend items
        
        first_file_xlabel = "Principal Component 1" # Default
        first_file_ylabel = "Principal Component 2" # Default
        first_file_title = "Combined Projection Plot"
        
        for i, file_path_str in enumerate(input_files):
            data, y_legends, title, xlabel, ylabel, filename_stem = parse_xvg_file(file_path_str)
            if data.size == 0 or data.ndim < 2 or data.shape[1] < 2:
                print(f"Skipping {file_path_str} for combined plot: No valid 2D data (X, Y). Shape: {data.shape}")
                continue
            
            if i == 0: # Use metadata from the first file for combined plot
                first_file_xlabel = xlabel if "PC" in xlabel or "Principal Component" in xlabel else "PC1" # Make it more generic for PCA
                first_file_ylabel = ylabel if "PC" in ylabel or "Principal Component" in ylabel else "PC2"
                first_file_title = f"Overlay: {Path(input_files[0]).stem} & {Path(input_files[1]).stem}" if len(input_files) == 2 else "Combined PCA Projection"


            # For 2D projections (e.g. PC1 vs PC2), X is data[:,0], Y is data[:,1]
            # The y_legends from parse_xvg_file might be for multiple Y columns if the file was e.g. time vs multiple energies.
            # For PCA 2D projections, we typically have X=PC1, Y=PC2.
            # So we plot data[:,0] vs data[:,1]
            plot_label = filename_stem # Use filename stem for legend in combined plots
            plt.scatter(data[:, 0], data[:, 1], label=plot_label, s=10, alpha=0.7)

        plt.xlabel(first_file_xlabel)
        plt.ylabel(first_file_ylabel)
        plt.title(first_file_title)
        plt.legend(fontsize='small', loc='best')
        plt.grid(True, alpha=0.5)
        plt.tight_layout()

        if args.output:
            output_path_combined = Path(args.output)
        else:
            # Create a default name for combined plot
            base_name = Path(input_files[0]).stem
            output_path_combined = Path(input_files[0]).with_name(f"combined_scatter_{base_name}.png")
        
        plt.savefig(output_path_combined, dpi=300, bbox_inches='tight')
        print(f"Combined scatter plot saved to {output_path_combined}")
        if not args.no_show:
            plt.show()
        plt.close()
        return # Handled combined plot, so exit

    # --- Process single files or other plot types ---
    # (This part remains largely the same, but iterates if multiple inputs were given without --noline for combined scatter)
    for file_path_str in input_files:
        input_file_path = Path(file_path_str)
        data, y_legends, title, xlabel, ylabel, input_filename_stem = parse_xvg_file(file_path_str)

        if data.size == 0 or data.ndim < 2 : # For most plots, we need at least 2D array (e.g. X vs Y)
            print(f"No valid 2D data to plot for {file_path_str}. Exiting or skipping.")
            continue # Skip this file if multiple inputs, or exit if single.

        is_ramachandran_file = "ramachandran" in input_filename_stem.lower() or \
                               ("phi" in xlabel.lower() and "psi" in ylabel.lower()) or \
                               ("phi" in title.lower() and "psi" in title.lower())


        default_suffix = ".png"
        if args.gaussian and not is_ramachandran_file:
            default_suffix = f"_{input_filename_stem}_gaussian.png"
        elif is_ramachandran_file:
             default_suffix = f"_{input_filename_stem}_ramachandran.png"
        elif args.noline: # Single file scatter
            default_suffix = f"_{input_filename_stem}_scatter.png"
        else:
             default_suffix = f"_{input_filename_stem}.png"

        if args.output and len(input_files) == 1: # Only use -o directly if one input file (and not handled by combined scatter)
            output_path = Path(args.output)
        else: # Default naming for single files or series of plots
            output_path = input_file_path.with_name(input_filename_stem + default_suffix)
            if args.output and len(input_files) > 1: # Using -o as a prefix/directory
                 output_dir = Path(args.output)
                 if output_dir.suffix: # if user gave a file name, use its parent as dir
                     output_dir = output_dir.parent
                 output_dir.mkdir(parents=True, exist_ok=True)
                 output_path = output_dir / (input_filename_stem + default_suffix)


        if is_ramachandran_file and data.shape[1] >= 2: # Needs at least 2 columns (phi, psi)
            print(f"Detected Ramachandran data for {file_path_str}. Generating 2D histogram.")
            plot_ramachandran(data, title, xlabel, ylabel, input_file_path, output_path=output_path, no_show=args.no_show)
        elif data.shape[1] < 2: # Not enough columns for X, Y plot
            print(f"Warning: Data in {file_path_str} has shape {data.shape}, not suitable for standard X,Y plotting. Skipping.")
            continue
        elif args.gaussian:
            print(f"Generating Gaussian analysis plot for {file_path_str}")
            fig = plot_with_gaussian(data, y_legends, title, xlabel, ylabel, input_filename_stem)
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {output_path}")
            if not args.no_show:
                plt.show()
            plt.close(fig)
        else: # Default plotting (line or scatter for single file)
            plot_type_msg = "scatter plot (due to --noline)" if args.noline else "standard time series plot"
            print(f"Generating {plot_type_msg} for {file_path_str}")
            plt.figure(figsize=(10, 6))
            
            # Standard plot: data[:,0] is X, subsequent columns are Y series
            for y_col_idx in range(1, data.shape[1]):
                legend_idx = y_col_idx - 1
                legend_label = y_legends[legend_idx] if legend_idx < len(y_legends) else f"Y-Col {y_col_idx}"
                
                if args.noline: # This is for single file scatter
                    plt.scatter(data[:, 0], data[:, y_col_idx], label=legend_label, s=10)
                else:
                    plt.plot(data[:, 0], data[:, y_col_idx], label=legend_label)

            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            final_title = f"{title} ({input_filename_stem})" if title != "GROMACS Data" else f"{input_filename_stem}"
            plt.title(final_title)
            
            # Improved legend display logic
            if data.shape[1] > 2 or (y_legends and any(not leg.startswith("Y-Col") for leg in y_legends)):
                 plt.legend()

            plt.grid(True, alpha=0.5)
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {output_path}")
            if not args.no_show:
                plt.show()
            plt.close()

if __name__ == "__main__":
    main()

