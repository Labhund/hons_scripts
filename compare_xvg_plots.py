#!/usr/bin/env python3
"""
Compare Multiple XVG Files - A utility to visualize data from multiple
GROMACS .xvg files on the same plot.

Usage:
    python compare_xvg_plots.py --input-files file1.xvg file2.xvg ... \
                                --output-file combined_plot.png \
                                [--labels label1 label2 ...] \
                                [--plot-type "Gyrate"]
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import re
from pathlib import Path

# Re-using the robust parser from xvg_plotter.py
def parse_xvg_file(filepath):
    """
    Parse an XVG file and extract data and metadata.
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    raw_data = []
    y_legends_dict = {}  # Stores legends like {0: "Potential", 1: "Kinetic En."}
    title = "GROMACS Data" # Default title
    xlabel = "X-axis"     # Default X label
    ylabel = "Y-axis"     # Default Y label

    for line in lines:
        line = line.strip()
        if not line:
            continue

        if line.startswith('@'): # General handler for metadata lines
            if line.startswith('@    title'):
                title = line.split('"')[1] if '"' in line else line.split()[-1]
            elif line.startswith('@    xaxis  label'):
                xlabel = line.split('"')[1] if '"' in line else line.split()[-1]
            elif line.startswith('@    yaxis  label'):
                ylabel = line.split('"')[1] if '"' in line else line.split()[-1]
            else:
                legend_match = re.match(r'@\s+s(\d+)\s+legend\s+"([^"]+)"', line)
                if legend_match:
                    s_index = int(legend_match.group(1)) # GROMACS s-index (s0, s1, ...)
                    legend_text = legend_match.group(2)
                    y_legends_dict[s_index] = legend_text
        elif not line.startswith('#'): # Data line
            parts = line.split()
            numeric_values_for_line = []
            for part in parts:
                try:
                    numeric_values_for_line.append(float(part))
                except ValueError:
                    break # Stop if non-numeric found
            if numeric_values_for_line: # Only add if we extracted some numbers
                raw_data.append(numeric_values_for_line)

    if not raw_data:
        # print(f"Debug: No raw data lines extracted from {filepath}")
        return np.array([]), [], title, xlabel, ylabel

    # Ensure all data rows have the same number of columns as the first data row
    num_cols_first_row = len(raw_data[0])
    consistent_raw_data = [row for row in raw_data if len(row) == num_cols_first_row]
    if len(consistent_raw_data) != len(raw_data):
        print(f"Warning: Inconsistent number of columns in data lines of {filepath}. "
              f"Using {len(consistent_raw_data)} rows out of {len(raw_data)} that matched the first row's column count ({num_cols_first_row}).")
    raw_data = consistent_raw_data

    if not raw_data: # If filtering made it empty
        # print(f"Debug: No consistent data after filtering for {filepath}")
        return np.array([]), [], title, xlabel, ylabel

    data_array = np.array(raw_data)

    final_y_legends = []
    if data_array.ndim == 2 and data_array.shape[1] > 1: # Need at least X and one Y column
        # data_array.shape[1] is total columns. Y-columns start from index 1.
        # s0 corresponds to data_array[:,1], s1 to data_array[:,2] etc.
        # So, for y_legends, index 0 is for s0, index 1 is for s1 etc.
        for y_col_idx_in_data_array in range(1, data_array.shape[1]):
            legend_key_for_dict = y_col_idx_in_data_array - 1 # This is the 's' index (s0, s1, ...)
            fallback_legend_name = f"Y-Col {y_col_idx_in_data_array} (s{legend_key_for_dict})"
            final_y_legends.append(y_legends_dict.get(legend_key_for_dict, fallback_legend_name))
    
    return data_array, final_y_legends, title, xlabel, ylabel

def main():
    parser = argparse.ArgumentParser(description='Compare multiple XVG files on a single plot.')
    parser.add_argument('--input-files', '-i', type=str, nargs='+', required=True,
                        help='List of .xvg files to plot.')
    parser.add_argument('--output-file', '-o', type=str, required=True,
                        help='Path to save the combined plot image.')
    parser.add_argument('--labels', '-l', type=str, nargs='+',
                        help='Optional: List of labels for each input file. Must match number of input files. If not provided, derived from file paths.')
    parser.add_argument('--plot-type', '-pt', type=str, default="Comparison",
                        help='Type of plot (e.g., Gyrate, SASA, RMSF) for title generation if --title is not set explicitly.')
    parser.add_argument('--title', '-t', type=str,
                        help='Custom title for the plot. Overrides default generated title.')
    parser.add_argument('--xlabel', '-xl', type=str, help='Custom X-axis label. Defaults to label from first successfully parsed input file.')
    parser.add_argument('--ylabel', '-yl', type=str, help='Custom Y-axis label. Defaults to label from first successfully parsed input file.')
    parser.add_argument('--y-col-index', '-yc', type=int, default=0,
                        help='0-indexed Y-series to plot from each XVG (e.g., 0 for the s0 series, 1 for s1). Default is 0, which corresponds to the second column in the XVG data block (first Y data series).')
    parser.add_argument('--no-show', action='store_true', help='Do not display the plot, only save it.')

    args = parser.parse_args()

    if args.labels and len(args.labels) != len(args.input_files):
        parser.error("--labels must have the same number of arguments as --input-files.")
        return

    plt.figure(figsize=(12, 7))
    
    plot_xlabel_final = args.xlabel
    plot_ylabel_final = args.ylabel
    plot_title_final = args.title

    first_file_processed_successfully = False
    files_plotted_count = 0

    for idx, filepath_str in enumerate(args.input_files):
        filepath = Path(filepath_str)
        if not filepath.is_file():
            print(f"Warning: File not found: {filepath}. Skipping.")
            continue

        data, file_y_legends, file_title_from_xvg, file_xlabel_from_xvg, file_ylabel_from_xvg = parse_xvg_file(filepath)

        # Determine the actual data column index for Y values in the numpy array
        # data[:, 0] is X. data[:, 1] is Y-series s0 (user's y_col_index 0), data[:, 2] is Y-series s1 (user's y_col_index 1), etc.
        actual_data_y_col_in_array = args.y_col_index + 1

        if data.size == 0 or data.ndim < 2 or data.shape[1] <= actual_data_y_col_in_array:
            print(f"Warning: No suitable data in {filepath} for Y-series index {args.y_col_index} (expected data column {actual_data_y_col_in_array}). File has {data.shape[1]} data columns. Skipping.")
            continue
        
        if not first_file_processed_successfully: # Use metadata from the first successfully parsed file
            if not plot_xlabel_final: plot_xlabel_final = file_xlabel_from_xvg
            if not plot_ylabel_final:
                # If the XVG has a specific legend for the chosen Y-column, prefer that for the Y-axis label
                if args.y_col_index < len(file_y_legends) and "Y-Col" not in file_y_legends[args.y_col_index]:
                    plot_ylabel_final = file_y_legends[args.y_col_index]
                else: # Fallback to the general Y-axis label from the file
                    plot_ylabel_final = file_ylabel_from_xvg

            if not plot_title_final:
                 plot_title_final = f"{args.plot_type} Comparison" # Default title based on plot type
            first_file_processed_successfully = True

        current_label = None
        if args.labels:
            current_label = args.labels[idx]
        else:
            # Derive label from path: e.g., variability_testing/g1_rep1/analysis_output/gyrate.xvg -> g1_rep1
            try:
                # Assumes structure like .../replicate_folder_name/analysis_output/file.xvg
                current_label = filepath.parent.parent.name
            except IndexError: # Fallback if path structure is different
                current_label = filepath.stem 

        # Append specific Y-series legend to the replicate label if it's informative and y_col_index > 0 (or if there are multiple y series)
        if args.y_col_index < len(file_y_legends):
            specific_y_legend_from_file = file_y_legends[args.y_col_index]
            if specific_y_legend_from_file and "Y-Col" not in specific_y_legend_from_file:
                 # Only append if the legend is not generic and provides extra info
                if current_label != specific_y_legend_from_file: # Avoid "Rep1 (Rep1)"
                    current_label = f"{current_label} ({specific_y_legend_from_file})"
            elif not current_label: # If auto-label failed and no specific legend
                 current_label = f"File {idx+1}, Y-series {args.y_col_index}"


        plt.plot(data[:, 0], data[:, actual_data_y_col_in_array], label=current_label, alpha=0.8)
        files_plotted_count +=1

    if files_plotted_count == 0:
        print("Error: No data could be plotted from any of the input files provided.")
        # Optionally create an empty plot or just exit
        if args.output_file and not args.no_show:
            # Create a blank plot to indicate failure if desired, or just don't save.
            plt.title("No Data Plotted")
            plt.xlabel("X")
            plt.ylabel("Y")
            plt.savefig(args.output_file, dpi=300)
            print(f"Empty plot saved to {args.output_file} as no data was plotted.")
        if not args.no_show:
            if files_plotted_count > 0 : # Only show if something was plotted
                plt.show()
        plt.close() # Close the figure regardless
        return


    plt.xlabel(plot_xlabel_final if plot_xlabel_final else "X-axis")
    plt.ylabel(plot_ylabel_final if plot_ylabel_final else "Y-axis")
    plt.title(plot_title_final if plot_title_final else "Comparison Plot")
    
    if files_plotted_count > 0 :
        plt.legend(loc='best', fontsize='small')
    plt.grid(True, alpha=0.4)
    plt.tight_layout()

    if args.output_file:
        plt.savefig(args.output_file, dpi=300)
        print(f"Combined plot saved to {args.output_file}")
    
    if not args.no_show:
        plt.show()
    
    plt.close() # Ensure figure is closed

if __name__ == "__main__":
    main()

