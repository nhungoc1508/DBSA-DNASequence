#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file> <output_file>"
    exit 1
fi

input_file="$1"
output_file="$2"

# Process the input file and write the output to the output file
awk '{
    # Remove leading numbers and spaces from each line
    $1 = ""; 
    # Remove all spaces in the remaining line
    gsub(/ /, "");
    # Append the cleaned line to the result string
    result = result $0
} END {
    # Write chunks of 8000 characters to the output
    for (i = 1; i <= length(result); i += 8000) {
        print substr(result, i, 8000)
    }
}' "$input_file" > "$output_file"

echo "Processing complete. The output is written to $output_file"