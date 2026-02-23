#!/bin/bash
# Script to run all derivative convergence examples

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "========================================================================"
echo "Running all derivative convergence examples"
echo "========================================================================"
echo ""

examples=(
    "hp_derivatives.py"
    "sp_derivatives.py"
    "tp_derivatives.py"
    "tv_derivatives.py"
    "sv_derivatives.py"
    "uv_derivatives.py"
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
echo "All derivative examples completed successfully!"
echo "========================================================================"
