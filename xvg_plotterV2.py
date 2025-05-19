#!/usr/bin/env python3
"""
XVG File Plotter - A utility to visualize GROMACS .xvg files

Usage:
    python xvg_plotter.py --input path/to/file.xvg [--gaussian] [--no-show] [--output filename.png]
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
        tuple: (data_array, y_legends, title, xlabel, ylabel)
               y_legends is a list of legends for the Y-axes.
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    raw_data = []
    y_legends_dict = {}  # Stores legends like {0: "Potential", 1: "Kinetic En."}
    title = "GROMACS Data"
    xlabel = "Time" 
    ylabel = "Value" 

    for line in lines:
        line = line.strip()
        if not line:
            continue

        if line.startswith('@    title'):
            title = line.split('"')[1] if '"' in line else line.split()[-1]
        elif line.startswith('@    xaxis  label'):
            xlabel = line.split('"')[1] if '"' in line else line.split()[-1]
        elif line.startswith('@    yaxis  label'):
            ylabel = line.split('"')[1] if '"' in line else line.split()[-1]
        elif line.startswith('@    s') and 'legend' in line: # Match GROMACS legend lines
            legend_match = re.match(r'@\s+s(\d+)\s+legend\s+"([^"]+)"', line)
            if legend_match:
                s_index = int(legend_match.group(1)) # GROMACS s-index (s0, s1, ...)
                legend_text = legend_match.group(2)
                y_legends_dict[s_index] = legend_text
        elif not line.startswith('@') and not line.startswith('#'):
            parts = line.split()
            numeric_values_for_line = []
            for part in parts:
                try:
                    numeric_values_for_line.append(float(part))
                except ValueError:
                    # Stop reading this line if a non-float part is encountered
                    # This handles cases like "phi psi RES-ID" in Ramachandran plots
                    break 
            
            if numeric_values_for_line: # Only add if we extracted some numbers
                raw_data.append(numeric_values_for_line)
            # else:
                # Optional: print a debug message if a non-comment line yielded no numeric data
                # print(f"Debug: Line did not yield leading numeric data: '{line}' in {filepath}")
    
    # print(f"Debug: y_legends_dict for {filepath}: {y_legends_dict}") # Uncomment for testing legends

    if not raw_data: 
        print(f"Warning: No data lines found or parsed successfully in {filepath}.")
        return np.array([]), [], title, xlabel, ylabel

    # Check for consistent number of columns based on the first data row
    if raw_data:
        num_cols_first_row = len(raw_data[0])
        consistent_raw_data = [row for row in raw_data if len(row) == num_cols_first_row]
        if len(consistent_raw_data) != len(raw_data):
            print(f"Warning: Some data lines in {filepath} had an inconsistent number of numeric columns. "
                  f"Using {len(consistent_raw_data)} rows out of {len(raw_data)} that matched the first row's column count ({num_cols_first_row}).")
        raw_data = consistent_raw_data

    if not raw_data: # If filtering made it empty
        print(f"Warning: No consistent data lines found after filtering for column consistency in {filepath}.")
        return np.array([]), [], title, xlabel, ylabel

    data_array = np.array(raw_data)

    final_y_legends = []
    if data_array.ndim == 2 and data_array.shape[1] > 1: 
        # data_array[:,0] is X. data_array[:,1] is the 1st Y-series, data_array[:,2] is the 2nd Y-series ...
        # y_col_idx_in_data ranges from 1 up to (num_columns - 1)
        for y_col_idx_in_data in range(1, data_array.shape[1]): 
            # The legend for the 1st Y-series (data_array[:,1]) should be from s0.
            # The legend for the 2nd Y-series (data_array[:,2]) should be from s1.
            # So, the key for y_legends_dict is (y_col_idx_in_data - 1)
            legend_key_for_dict = y_col_idx_in_data - 1 
            
            # Fallback name if legend not found in dictionary
            fallback_legend_name = f"Y-Col {y_col_idx_in_data} (s{legend_key_for_dict})"
            
            final_y_legends.append(y_legends_dict.get(legend_key_for_dict, fallback_legend_name))
            
    elif data_array.ndim == 1 and data_array.size > 0: 
        print(f"Warning: Data in {filepath} is 1D after parsing. Cannot plot as time series. Shape: {data_array.shape}")
        return np.array([]), [], title, xlabel, ylabel
    elif data_array.size == 0 : 
        print(f"Warning: No valid data points resulted after parsing {filepath}. Final data_array shape: {data_array.shape}")
        return np.array([]), [], title, xlabel, ylabel
    
    # print(f"Debug: final_y_legends for {filepath}: {final_y_legends}") # Uncomment for testing legends
        
    return data_array, final_y_legends, title, xlabel, ylabel

def analyze_gaussian(data_column):
    """
    Perform Gaussian analysis on a single data column.
    """
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
    
    return (
        {'mu': mu1, 'sigma': sigma1},
        {'mu': mu2, 'sigma': sigma2},
        {'mu': mu_total, 'sigma': sigma_total}
    )

def plot_with_gaussian(data, y_legends, title, xlabel, ylabel, input_filename_stem):
    """
    Create a plot with Gaussian distribution analysis.
    y_legends is a list, where y_legends[k] is for data[:, k+1]
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=False, height_ratios=[2, 1])

    # Plot original time series in top panel
    for y_col_idx in range(1, data.shape[1]): # y_col_idx is 1 for data[:,1], 2 for data[:,2], etc.
        legend_idx = y_col_idx - 1 # y_legends is 0-indexed for Y columns
        legend_label = y_legends[legend_idx] if legend_idx < len(y_legends) and y_legends[legend_idx] else f"Series Y{y_col_idx}"
        ax1.plot(data[:, 0], data[:, y_col_idx], label=legend_label)

    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_title(f"{title} ({input_filename_stem}) - Time Series")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot Gaussian distribution in bottom panel
    for y_col_idx in range(1, data.shape[1]):
        column_data = data[:, y_col_idx]
        legend_idx = y_col_idx - 1
        legend_label = y_legends[legend_idx] if legend_idx < len(y_legends) and y_legends[legend_idx] else f"Series Y{y_col_idx}"

        if len(column_data) < 2: 
            continue

        stats_first, stats_second, stats_total = analyze_gaussian(column_data)

        hist, bins, _ = ax2.hist(column_data, bins=50, density=True, alpha=0.3,
                                 label=f'Dist: {legend_label}')
        
        min_val, max_val = np.min(column_data), np.max(column_data)
        if min_val == max_val: 
            x_smooth = np.array([min_val])
        else:
            x_smooth = np.linspace(min_val, max_val, 200)
        
        if x_smooth.size > 0:
            if not np.isnan(stats_first['mu']) and not np.isnan(stats_first['sigma']) and stats_first['sigma'] > 1e-6:
                 ax2.plot(x_smooth, stats.norm.pdf(x_smooth, stats_first['mu'], stats_first['sigma']), '--',
                         label=f'{legend_label} First Half (μ={stats_first["mu"]:.2f}, σ={stats_first["sigma"]:.2f})')
            if not np.isnan(stats_second['mu']) and not np.isnan(stats_second['sigma']) and stats_second['sigma'] > 1e-6:
                ax2.plot(x_smooth, stats.norm.pdf(x_smooth, stats_second['mu'], stats_second['sigma']), ':',
                         label=f'{legend_label} Second Half (μ={stats_second["mu"]:.2f}, σ={stats_second["sigma"]:.2f})')

    ax2.set_xlabel(ylabel) 
    ax2.set_ylabel('Probability Density')
    ax2.set_title('Distribution Analysis')
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout(rect=[0, 0, 0.9, 1]) 
    return fig

def plot_ramachandran(data, title, xlabel, ylabel, input_filename_path, output_path=None, no_show=False):
    """
    Create a Ramachandran plot (2D histogram).
    Assumes data has at least two columns: phi and psi.
    """
    if data.shape[1] < 2:
        print(f"Error: Ramachandran plot requires at least 2 columns of data, got {data.shape[1]}. Skipping {input_filename_path.name}")
        return

    phi_angles = data[:, 0]
    psi_angles = data[:, 1]

    plt.figure(figsize=(8, 7))
    
    final_xlabel = xlabel if "phi" in xlabel.lower() or "φ" in xlabel else "φ (degrees)"
    final_ylabel = ylabel if "psi" in ylabel.lower() or "ψ" in ylabel else "ψ (degrees)"
    
    counts, xedges, yedges, im = plt.hist2d(phi_angles, psi_angles, bins=120, cmap='viridis', cmin=1, range=[[-180, 180], [-180, 180]])
    
    plt.colorbar(im, label='Density', shrink=0.8)

    plt.xlabel(final_xlabel)
    plt.ylabel(final_ylabel)
    
    plot_title = f"Ramachandran Plot: {input_filename_path.name}"
    if title and title != "GROMACS Data" and title != "Ramachandran Plot": # Avoid redundant titles
        plot_title = f"{title} ({input_filename_path.name})"
    plt.title(plot_title)

    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.axhline(0, color='grey', lw=0.5, linestyle='--')
    plt.axvline(0, color='grey', lw=0.5, linestyle='--')
    plt.xticks(np.arange(-180, 181, 60))
    plt.yticks(np.arange(-180, 181, 60))
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.gca().set_aspect('equal', adjustable='box')

    plt.tight_layout(rect=[0, 0, 0.95, 1]) # Adjust rect to make space for colorbar

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Ramachandran plot saved to {output_path}")
    if not no_show:
        plt.show()
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Plot XVG files with optional Gaussian analysis and Ramachandran handling.')
    parser.add_argument('--input', '-i', type=str, help='Path to the XVG file', required=False)
    parser.add_argument('--gaussian', '-g', action='store_true', help='Enable Gaussian distribution analysis (for time series data)')
    parser.add_argument('--output', '-o', type=str, help='Output file path (optional)')
    parser.add_argument('--no-show', action='store_true', help='Do not display the plot (save only)')
    
    parser.add_argument('filepath', nargs='?', type=str, help=argparse.SUPPRESS) 
    
    args = parser.parse_args()
    
    if args.filepath and not args.input:
        args.input = args.filepath
    elif not args.input and not args.filepath:
        parser.error("the following arguments are required: --input or positional filepath")
        return

    input_file_path = Path(args.input)
    input_filename_stem = input_file_path.stem

    data, y_legends, title, xlabel, ylabel = parse_xvg_file(args.input)

    if data.size == 0 or data.ndim < 2:
        print(f"No valid 2D data to plot for {args.input}. Exiting.")
        return
    
    if args.output:
        output_path = Path(args.output)
    else:
        suffix = "_gaussian.png" if args.gaussian and "ramachandran" not in input_filename_stem.lower() else ".png"
        output_path = input_file_path.with_name(input_filename_stem + suffix)

    is_ramachandran_file = "ramachandran" in input_filename_stem.lower()
    
    if is_ramachandran_file and data.shape[1] >= 2:
        print(f"Detected Ramachandran data for {args.input}. Generating 2D histogram.")
        plot_ramachandran(data, title, xlabel, ylabel, input_file_path, output_path=output_path, no_show=args.no_show)
    elif data.shape[1] < 2: 
        print(f"Warning: Data in {args.input} has shape {data.shape}, "
              f"which is not suitable for standard time series plotting (requires at least X and one Y column). Skipping.")
        return
    elif args.gaussian:
        print(f"Generating Gaussian analysis plot for {args.input}")
        fig = plot_with_gaussian(data, y_legends, title, xlabel, ylabel, input_filename_stem)
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {output_path}")
        if not args.no_show:
            plt.show()
        plt.close(fig)
    else:
        print(f"Generating standard time series plot for {args.input}")
        plt.figure(figsize=(10, 6))
        for y_col_idx in range(1, data.shape[1]): 
            legend_idx = y_col_idx - 1 
            legend_label = y_legends[legend_idx] if legend_idx < len(y_legends) and y_legends[legend_idx] else f"Series Y{y_col_idx}"
            plt.plot(data[:, 0], data[:, y_col_idx], label=legend_label)
        
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        final_title = f"{title} ({input_filename_stem})" if title != "GROMACS Data" else f"{input_filename_stem}"
        plt.title(final_title)
        
        # Show legend if there are actual legends provided by y_legends or if there are multiple Y series to differentiate
        has_meaningful_legends = any(y_legends) # Check if any legend string is non-empty or not a fallback
        if has_meaningful_legends or (data.shape[1] > 2 and not all(leg.startswith("Y-Col") or leg.startswith("Series Y") for leg in y_legends)):
             plt.legend()
        elif data.shape[1] > 2 : # Multiple series, but maybe all fallbacks, still good to have legend
            plt.legend()


        plt.grid(True, alpha=0.5)
        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {output_path}")
        if not args.no_show:
            plt.show()
        plt.close()

if __name__ == "__main__":
    main()

