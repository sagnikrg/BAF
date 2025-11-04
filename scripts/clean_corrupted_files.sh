#!/bin/bash

# Directory where files are stored (edit if needed)
DATA_DIR="../../data"
# Output log file
L=12
theta=0.25
OUTPUT_FILE="mbldtc_L${L}_theta_${theta}_corrupted_indices.txt"
# Size threshold (in bytes) — 20 MB = 20 * 1024 * 1024
THRESHOLD=$((20 * 1024 * 1024))

# Clear output file
> "$OUTPUT_FILE"

echo "Checking for corrupted files"

for itr in $(seq 1 2000); do
    FILE="${DATA_DIR}/mbldtc_L${L}_theta_${theta}_${itr}.hdf5"
    if [[ -f "$FILE" ]]; then
        FILE_SIZE=$(stat -c%s "$FILE" 2>/dev/null)
        if (( FILE_SIZE < THRESHOLD )); then
            echo "$itr" >> "$OUTPUT_FILE"
            #echo "Deleting $FILE (size: $((FILE_SIZE / 1024 / 1024)) MB)"
            #rm "$FILE"
        fi
    fi
done

echo "Done. Deleted files listed in $OUTPUT_FILE"
