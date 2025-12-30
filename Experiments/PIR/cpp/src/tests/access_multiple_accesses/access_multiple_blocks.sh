#!/bin/bash
# Script generated with assistance from AI(Google Gemini) to parse CSVs and check exit codes.

# 1. Check if a CSV file was provided as an argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <filename.csv, which contains all the indices to be fetched seperated by comma>"
    exit 1
fi

CSV_FILE="$1"
EXECUTABLE="/mnt/sumit/pir/Experiments/PIR/cpp/src/build/pir_client"

# 2. Check file existence
if [ ! -f "$CSV_FILE" ]; then
    echo "Error: File '$CSV_FILE' not found."
    exit 1
fi

if [ ! -x "$EXECUTABLE" ]; then
    echo "Error: '$EXECUTABLE' not found or not executable."
    exit 1
fi

echo "Processing $CSV_FILE..."
echo "--------------------------------"

# Initialize failure counter
FAIL_COUNT=0

# 3. Process the loop
# We use '< <(command)' to redirect output into the loop without creating a subshell.
# This ensures FAIL_COUNT updates are remembered.
while read -r raw_input; do

    # Sanitize input: remove spaces
    number=$(echo "$raw_input" | tr -d '[:space:]')

    # Skip empty lines
    if [ -z "$number" ]; then
        continue
    fi

    # Validate integer
    if ! [[ "$number" =~ ^[0-9]+$ ]]; then
        echo "[SKIP] '$number' is not a valid unsigned integer."
        continue
    fi

    # Run executable
    "$EXECUTABLE" "$number"
    
    # Check exit code
    exit_code=$?

    if [ $exit_code -ne 0 ]; then
        echo "[FAIL] Input: $number"
        # Increment the failure counter
        ((FAIL_COUNT++))
    fi
    sleep 1

done < <(tr ',' '\n' < "$CSV_FILE")

echo "--------------------------------"

# 4. Final Summary Check
if [ $FAIL_COUNT -eq 0 ]; then
    echo "SUCCESS: All block are fetched successfully."
    echo "Verify the logs of all three servers. There must not be any error messages."
    exit 0
else
    echo "FAILURE: There were $FAIL_COUNT error(s) found during execution."
    exit 1
fi