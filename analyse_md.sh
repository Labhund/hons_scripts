#!/bin/bash

# Parse arguments
defnm=""
files=()

while [[ $# -gt 0 ]]; do
    case $1 in
        --defnm)
            defnm="$2"
            shift 2
            ;;
        *)
            files+=("$1")
            shift
            ;;
    esac
done

# Function for energy analysis
analyze_energy() {
    local edr_file=$1
    local term=$2
    local output_xvg="energy_${term// /_}.xvg"
    echo "$term" | gmx_mpi energy -f $edr_file -o $output_xvg
}

# Function for gyration radius
analyze_gyrate() {
    local xtc_file=$1
    local tpr_file=$2
    local output_xvg="gyrate.xvg"
    echo "1" | gmx_mpi gyrate -f $xtc_file -s $tpr_file -o $output_xvg
}

# Function for RMSD
analyze_rmsd() {
    local xtc_file=$1
    local tpr_file=$2
    local output_xvg="rmsd_backbone.xvg"
    echo "4 4" | gmx_mpi rms -f $xtc_file -s $tpr_file -o $output_xvg
}

# Process XVG files with plotter
process_xvg_files() {
    echo "Processing all XVG files with plotter..."
    
    # Create a plots directory if it doesn't exist
    mkdir -p plots
    
    # Loop through all XVG files in the current directory
    for xvg_file in *.xvg; do
        # Skip if no XVG files found
        [[ -e "$xvg_file" ]] || continue
        
        echo "Plotting $xvg_file..."
        
        # Generate output filename in plots directory
        output_file="plots/${xvg_file%.xvg}.png"
        
        # Run the plotter with appropriate options
        python xvg_plotter.py --input "$xvg_file" --gaussian --no-show --output "$output_file"
    done
    
    echo "All plots saved to the 'plots' directory"
}

# Energy terms to analyze
terms=(
    "Temperature"
    "Pressure"
    "Density"
    "Potential"
    "Total-Energy"
    "Volume"
    "T-Protein"
    "T-Water_and_ions"
)

# If defnm is provided, construct file names
if [[ -n "$defnm" ]]; then
    edr_file="${defnm}.edr"
    xtc_file="${defnm}.xtc"
    tpr_file="topol.tpr"  # keeping topol.tpr as default
else
    # Check if all required files are provided
    if [ "${#files[@]}" -lt 3 ]; then
        echo "Usage: $0 [--defnm basename] <edr_file> <xtc_file> <tpr_file>"
        echo "Example with --defnm: $0 --defnm md_10_40"
        echo "Example with files: $0 md_10_40.edr md_10_40.xtc topol.tpr"
        exit 1
    fi
    edr_file=${files[1]}
    xtc_file=${files[2]}
    tpr_file=${files[3]}
fi

# Check if files exist
if [[ ! -f $edr_file ]] || [[ ! -f $xtc_file ]] || [[ ! -f $tpr_file ]]; then
    echo "Error: One or more input files not found"
    echo "Looking for:"
    echo "  EDR: $edr_file"
    echo "  XTC: $xtc_file"
    echo "  TPR: $tpr_file"
    exit 1
fi

# Run energy analysis
echo "Running energy analysis..."
for term in $terms; do
    echo "Analyzing $term..."
    analyze_energy $edr_file $term
done

# Run gyration analysis
echo "Running radius of gyration analysis..."
analyze_gyrate $xtc_file $tpr_file

# Run RMSD analysis
echo "Running RMSD analysis..."
analyze_rmsd $xtc_file $tpr_file

# Process all XVG files with the plotter
process_xvg_files

echo "Analysis complete! All plots have been saved to the 'plots' directory."

