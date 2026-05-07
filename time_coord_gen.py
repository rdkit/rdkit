#!/usr/bin/env python3
"""
Benchmark 2D coordinate generation using locally built RDKit.

Usage:
    $SCHRODINGER/run python3 time_coord_gen.py <smiles_file>
    $SCHRODINGER/run python3 time_coord_gen.py - < smiles.txt

The SMILES file should contain one SMILES per line.
"""

import sys
import time
from statistics import mean, median, stdev

# Import RDKit from local build
try:
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
except ImportError as e:
    print(f"ERROR: Failed to import RDKit: {e}", file=sys.stderr)
    print("Make sure you're running with: $SCHRODINGER/run python3", file=sys.stderr)
    sys.exit(1)


def generate_2d_coords(smiles, use_templates=True):
    """
    Generate 2D coordinates for a SMILES string.

    Returns (success, time_seconds, error_message)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, 0.0, "Failed to parse SMILES"

    # Time the coordinate generation
    start_time = time.perf_counter()
    try:
        rdDepictor.Compute2DCoords(
            mol,
            True,   # canonOrient
            True,   # clearConfs
            {},     # coordMap
            3,      # nFlipsPerSample
            100,    # nSample
            100,    # sampleSeed
            False,  # permuteDeg4Nodes
            -1.0,   # bondLength
            False,  # forceRDKit
            use_templates  # useRingTemplates
        )
    except Exception as e:
        return False, 0.0, str(e)

    end_time = time.perf_counter()
    elapsed = end_time - start_time

    return True, elapsed, None


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    input_file = sys.argv[1]

    # Read SMILES
    smiles_list = []
    if input_file == "-":
        print("Reading SMILES from stdin...", file=sys.stderr)
        smiles_list = [line.strip() for line in sys.stdin if line.strip()]
    else:
        print(f"Reading SMILES from {input_file}...", file=sys.stderr)
        try:
            with open(input_file, 'r') as f:
                smiles_list = [line.strip() for line in f if line.strip()]
        except FileNotFoundError:
            print(f"ERROR: File not found: {input_file}", file=sys.stderr)
            sys.exit(1)

    if not smiles_list:
        print("ERROR: No SMILES found in input", file=sys.stderr)
        sys.exit(1)

    print(f"Processing {len(smiles_list)} molecules...", file=sys.stderr)
    print()

    # Process each SMILES
    times = []
    success_count = 0
    fail_count = 0

    for i, smiles in enumerate(smiles_list, 1):
        success, elapsed, error = generate_2d_coords(smiles)

        if success:
            success_count += 1
            times.append(elapsed)
            print(f"[{i}/{len(smiles_list)}] ✓ {elapsed*1000:.2f}ms - {smiles[:60]}")
        else:
            fail_count += 1
            print(f"[{i}/{len(smiles_list)}] ✗ FAILED - {smiles[:60]} - {error}")

    # Print statistics
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total molecules:     {len(smiles_list)}")
    print(f"Successful:          {success_count}")
    print(f"Failed:              {fail_count}")

    if times:
        print()
        print(f"Total time:          {sum(times):.3f}s")
        print(f"Mean time:           {mean(times)*1000:.2f}ms")
        print(f"Median time:         {median(times)*1000:.2f}ms")
        if len(times) > 1:
            print(f"Std dev:             {stdev(times)*1000:.2f}ms")
        print(f"Min time:            {min(times)*1000:.2f}ms")
        print(f"Max time:            {max(times)*1000:.2f}ms")
        print(f"Throughput:          {len(times)/sum(times):.1f} molecules/sec")


if __name__ == "__main__":
    main()
