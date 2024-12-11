#!/bin/bash 
# Convert comma-separated string back into an array
IFS=" " read -r -a pair <<< "$pair"

# Iterate through the array elements
for p in "${pair[@]}"; do
    echo "$p"
done
