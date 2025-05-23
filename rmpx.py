import argparse
import re # For parsing residue numbers from lines that might include names

def parse_residue_numbers_from_file(filepath):
    """
    Parses residue numbers from a file.
    Each line can contain just a number, or a number and other text (e.g., residue name).
    It will attempt to extract the first integer found on each line.
    """
    residue_numbers = []
    try:
        with open(filepath, 'r') as f:
            for line_number, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'): # Skip empty lines or comments
                    continue
                
                found_numbers = re.findall(r'\d+', line)
                
                if found_numbers:
                    try:
                        res_num = int(found_numbers[0])
                        residue_numbers.append(res_num)
                    except ValueError:
                        print(f"Warning: Could not parse an integer from numbers '{found_numbers[0]}' on line {line_number} in {filepath}: {line}")
    except FileNotFoundError:
        print(f"Error: Input residue file not found at {filepath}")
        return None 
    except Exception as e:
        print(f"Error reading input residue file {filepath}: {e}")
        return None 
    return residue_numbers

def extract_residue_info_from_map_data(map_lines, ref_residue_numbers_input, map_filepath_for_error_reporting="map_file"):
    """
    Core logic to extract residue info from map data (list of lines).
    For a given reference residue number (from input), it finds the matching line in map_lines 
    (based on the second column, RefResNum).
    It then outputs the RefResName (first column) and InputResNum (third column) from that map line.
    """
    selected_residues_output = []
    if not ref_residue_numbers_input:
        return []

    ref_residue_numbers_set = set(ref_residue_numbers_input) # These are numbers like 339 from the input file

    for line_number, line in enumerate(map_lines, start=1): # map_lines already excludes header
        parts = line.strip().split('\t')
        if len(parts) == 3:
            ref_res_name_map = parts[0]       # e.g., "ILE" from map file
            try:
                ref_res_num_map = int(parts[1]) # e.g., 339 (RefResNum from map file)
                
                # Check if the RefResNum from the map file is one of the numbers we're looking for
                if ref_res_num_map in ref_residue_numbers_set: 
                    try:
                        # If it matches, we want to output RefResName (parts[0]) and InputResNum (parts[2])
                        input_res_num_map_val = int(parts[2]) # e.g., 357 (InputResNum from map file)
                        selected_residues_output.append([ref_res_name_map, input_res_num_map_val]) # Append ["ILE", 357]
                    except ValueError:
                        # This would happen if the third column (InputResNum) is not a valid integer.
                        print(f"Warning: Skipping line in map file due to non-integer InputResNum '{parts[2]}' for RefResNum {ref_res_num_map} on map file line corresponding to data line #{line_number}: {line.strip()}")
                        continue
            except ValueError:
                # This happens if the second column (RefResNum) in the map file is not an integer.
                # print(f"Warning: Skipping malformed RefResNum on map file line corresponding to data line #{line_number}: {line.strip()}")
                continue
        # else:
            # if line.strip(): 
                # print(f"Warning: Skipping map file line corresponding to data line #{line_number} with unexpected format: {line.strip()}")
    return selected_residues_output

def main():
    parser = argparse.ArgumentParser(
        description="Extracts residue names (RefResName) and target numbers (InputResNum) from a map file based on a list of reference residue numbers (RefResNum). Residue numbers for matching can be provided via an input file (-f) or directly (-n).",
        formatter_class=argparse.RawTextHelpFormatter
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "-f", "--file",
        dest="residue_file",
        help="Path to an input text file containing reference residue numbers (these should match the 2nd column, RefResNum, of the map file). Each line should contain one residue number, optionally with its name (e.g., '339' or 'ILE 339')."
    )
    input_group.add_argument(
        "-n", "--numbers", 
        dest="residue_numbers_direct",
        nargs='+', 
        type=int,
        help="A space-separated list of reference residue numbers (e.g., 339 343 380). Used if -f is not provided."
    )

    parser.add_argument(
        "-m", "--map",
        dest="map_file",
        required=True,
        help="Path to the residue map file (e.g., ecr_resmap.rmp). This file should have a header and be tab-separated with columns: RefResName, RefResNum, InputResNum."
    )
    parser.add_argument(
        "-o", "--output",
        dest="output_file",
        help="Optional: Path to the output file. If not specified, results are printed to the console. Output format is one 'RefResName\\tInputResNum' per line."
    )

    args = parser.parse_args()

    residue_numbers_to_process = []
    if args.residue_file:
        residue_numbers_to_process = parse_residue_numbers_from_file(args.residue_file)
        if residue_numbers_to_process is None: 
            return 
    elif args.residue_numbers_direct: 
        residue_numbers_to_process = args.residue_numbers_direct
    
    if not residue_numbers_to_process: 
        print("Error: No residue numbers provided. Please use -f or -n.")
        parser.print_help()
        return

    map_lines_content = []
    try:
        with open(args.map_file, 'r') as f_map:
            map_lines_content = f_map.readlines()
    except FileNotFoundError:
        print(f"Error: Map file not found at {args.map_file}")
        return
    except Exception as e:
        print(f"Error reading map file {args.map_file}: {e}")
        return

    if not map_lines_content or len(map_lines_content) < 2:
        print(f"Error: Map file {args.map_file} is empty or does not contain enough data (header + data lines).")
        return
        
    header = map_lines_content[0] 
    data_lines = map_lines_content[1:]

    output_array = extract_residue_info_from_map_data(data_lines, residue_numbers_to_process, args.map_file)

    output_lines = []
    if output_array:
        for res_info in output_array: # res_info is now [RefResName_map, InputResNum_map_val]
            output_lines.append(f"{res_info[0]}\t{res_info[1]}")
    else:
        if residue_numbers_to_process and map_lines_content:
             output_lines.append(f"# No matching residues found in '{args.map_file}' for the provided selection.")

    if args.output_file:
        try:
            with open(args.output_file, 'w') as f_out:
                for line in output_lines:
                    f_out.write(line + "\n")
            if output_array : 
                 print(f"Results successfully written to {args.output_file}")
            elif not output_array and output_lines: 
                 print(f"Information written to {args.output_file}")
        except Exception as e:
            print(f"Error writing to output file {args.output_file}: {e}")
    else:
        if not output_lines and not output_array: 
            print(f"No matching residues found in '{args.map_file}' for the provided selection, and no output generated.")
        else:
            for line in output_lines:
                print(line)

if __name__ == "__main__":
    main()

