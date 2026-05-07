#!/bin/bash
# Benchmark 2D coordinate generation using locally built RDKit
# Usage: ./time_coord_gen.sh <smiles_file>
#        ./time_coord_gen.sh - < smiles.txt

if [ -z "$1" ]; then
    echo "Usage: time_coord_gen.sh <smiles_file>"
    echo "       time_coord_gen.sh - < smiles.txt"
    echo ""
    echo "The SMILES file should contain one SMILES per line."
    exit 1
fi

# Set up environment for local RDKit build
export RDKIT_BUILD="$SCHRODINGER/rdkit"
export PYTHONPATH="$RDKIT_BUILD/lib/python3.11/site-packages:$PYTHONPATH"

# Run the timing script (from ~/bin to avoid importing source rdkit)
$SCHRODINGER/run python3 ~/bin/time_coord_gen.py "$@"
