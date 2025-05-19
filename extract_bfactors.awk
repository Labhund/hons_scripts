# Corrected extract_bfactors.awk
# Extracts C-alpha B-factors for a specific chain from a PDB file.
# Usage: awk -v TARGET_CHAIN="A" -f extract_bfactors.awk your_pdb_file.pdb

BEGIN {
    if (TARGET_CHAIN == "") {
        print "Error: TARGET_CHAIN variable must be set." > "/dev/stderr"
        exit 1
    }
}

/^ATOM/ {
    # Atom name is in columns 13-16 (e.g., " CA ", "  CA")
    atom_name_str = substr($0, 13, 4)
    # Chain identifier is in column 22
    chain_id_str = substr($0, 22, 1)

    # Trim whitespace from atom_name_str for comparison
    gsub(/^[ \t]+|[ \t]+$/, "", atom_name_str)

    if (atom_name_str == "CA" && chain_id_str == TARGET_CHAIN) {
        # Residue sequence number is in columns 23-26
        res_num_str = substr($0, 23, 4)
        # B-factor is in columns 61-66
        bfactor_str = substr($0, 61, 6)

        # Trim whitespace from extracted numerical values
        gsub(/^[ \t]+|[ \t]+$/, "", res_num_str)
        gsub(/^[ \t]+|[ \t]+$/, "", bfactor_str)

        # Ensure residue number and B-factor are valid numbers before printing
        if (res_num_str ~ /^[0-9]+$/ && bfactor_str ~ /^-?[0-9\.]+$/) {
            printf "%s %s\n", res_num_str, bfactor_str
        }
    }
}

