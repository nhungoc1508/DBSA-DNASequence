# #!/bin/bash

# # Check if the correct number of arguments is provided
# if [ $# -ne 2 ]; then
#   echo "Usage: $0 <input_file> <output_file>"
#   exit 1
# fi

# # Input and output files from command line arguments
# input_file="$1"
# output_file="$2"

# # Check if the input file exists
# if [ ! -f "$input_file" ]; then
#   echo "Error: File '$input_file' not found."
#   exit 1
# fi

# # Read the file, remove spaces and newlines, and write to the output file
# tr -d ' \n' < "$input_file" > "$output_file"

# echo "Output written to $output_file"

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
    # Concatenate the cleaned line
    printf "%s", $0
}' "$input_file" > "$output_file"

echo "Processing complete. The output is written to $output_file"