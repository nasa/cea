#!/bin/bash
# Script to run all RP-1311 sample examples

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "========================================================================"
echo "Running all RP-1311 examples"
echo "========================================================================"
echo ""

examples=(
    "example1.py"
    "example2.py"
    "example3.py"
    "example4.py"
    "example5.py"
    "example6.py"
    "example7.py"
    "example8.py"
    "example9.py"
    "example10.py"
    "example11.py"
    "example12.py"
    "example13.py"
    "example14.py"
)

total=${#examples[@]}
current=0

for example in "${examples[@]}"; do
    current=$((current + 1))
    echo "------------------------------------------------------------------------"
    echo "[$current/$total] Running $example"
    echo "------------------------------------------------------------------------"
    python "$example"
    echo ""
    echo "[$current/$total] Completed $example"
    echo ""
done

echo "========================================================================"
echo "All RP-1311 examples completed successfully!"
echo "========================================================================"
