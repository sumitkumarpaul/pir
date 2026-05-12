#!/bin/bash

# Define the path to your config file
CONFIG_FILE="config.txt"


# ==========================================
# 1. Define test functions
# ==========================================
access_random_block() {
    echo -e "\nTest A1: Client is accessing a random block:\n"
}

access_specific_block() {
    echo -e "\nTest A2: Client is accessing a specific block:\n"
}

access_all_blocks() {
    echo -e "\nTest A3: Client is accessing all the blocks:\n"

    # This means the server is required to be restarted sqrt_N times
    for (( i=1; i<=sqrt_N; i++ )); do
        #Wait a little more than the server_alpha
        sleep 2
        $EXECUTABLE "process_request"
        j=$(echo "scale=4; $i * $sqrt_N" | bc -l)
        echo -ne "Progress: $j out of $N blocks are tested\r" >&3
    done
    echo "">&3
    echo "All the blocks are tested!">&3
}

# 2. Load the configuration file
if [[ -f "$CONFIG_FILE" ]]; then
    source "$CONFIG_FILE"
else
    echo "Error: Cannot find $CONFIG_FILE" 2>&1
    exit 1
fi

# 2. Check if log_dir was provided
if [[ -z "$log_dir" ]]; then
    echo "Error: log_dir is not set in $CONFIG_FILE" 2>&1
    exit 1
fi

# 3. Create the log directory if it doesn't exist
mkdir -p "$log_dir"

# 4. Generate a timestamp and define the log file path
TIMESTAMP=$(date +"%d%b_%Y_%Hh_%Mm")
LOG_FILE="${log_dir}/test_beta_${TIMESTAMP}.txt"

EXECUTABLE="${bin_dir}/server_beta"

if [ ! -x "$EXECUTABLE" ]; then
    echo "Error: '$EXECUTABLE' not found or not executable."
    exit 1
fi

# 5. Redirect all output (stdout and stderr) to both the terminal and the log file
# Process substitution >(...) feeds the output of the script into tee
# 5.1. SAVE the original terminal output to File Descriptor 3
exec 3>&1

# 5.2. REDIRECT all future standard output and error to tee
exec > >(tee -a "$LOG_FILE") 2>&1

# ==========================================
# From here down, everything is logged and printed automatically
# ==========================================

echo "--- Test started at $(date) ---"
echo "Testing with the loaded variables from the loaded configuration file: $CONFIG_FILE"
echo "Loaded configurations are:"
echo "=============================================================================="
# Use the variables
echo "Number of blocks in the database (N)                : $N"
echo "Size of epoch (sqrt_N)                              : $sqrt_N"
echo "Number of bits in each block (B)                    : $B"
echo "Number of bits in each shuffled database tag (p)    : $p"
echo "Number of bits in each short tag (r)                : $r"
echo "Log file location (log_dir)                         : $log_dir"
echo "Binary file location (bin_dir)                      : $bin_dir"
echo "=============================================================================="
echo "Log file created at: $LOG_FILE"

# ==========================================
# 2. Display a menu to the user
# ==========================================

echo "================================================="
echo "    Please select any of the following tests"
echo "================================================="
echo "A1. Access one random block"
echo "A2. Access one specific block"
echo "A3. Access all the blocks from the remote database"

echo "q. Quit"
echo "======================"

# ==========================================
# 3. Get input and call the right function
# ==========================================

# Read the user's input into the variable 'choice'
read -p "Please select an option (or q): " choice
echo -e "\n"

# Use a case statement to evaluate the input
case $choice in
    A1)
        # Call the show_date function
        access_random_block
        ;;
    A2)
        # Call the show_disk_space function
        # Ask for another piece of input, then pass it to the function
        read -p "Specify the block index you want to access? " index        
        echo -e "\n"
        access_specific_block $index
        ;;
    A3)
        access_all_blocks
        ;;
    q|Q)
        # Handle upper and lowercase Q
        echo "Exiting script. Goodbye!"
        exit 0
        ;;
    *)
        # The *) acts as a catch-all for invalid inputs
        echo "Invalid option. Please run the script again with proper option."
        ;;
esac


echo "--- Test finished at $(date) ---"