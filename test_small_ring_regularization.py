#!/usr/bin/env python3
"""
Test script for small ring regularization in macrocycles.
Tests that small rings (4, 5, 6 atoms) maintain regular polygon shapes.
"""

import sys
import math
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor

def compute_turn_angle(coords, i):
    """Compute the turn angle at position i given 3 consecutive points."""
    n = len(coords)
    p_prev = coords[i % n]
    p_curr = coords[(i + 1) % n]
    p_next = coords[(i + 2) % n]

    # Compute vectors
    v1_x = p_curr.x - p_prev.x
    v1_y = p_curr.y - p_prev.y
    v2_x = p_next.x - p_curr.x
    v2_y = p_next.y - p_curr.y

    # Compute signed angle using atan2
    angle1 = math.atan2(v1_y, v1_x)
    angle2 = math.atan2(v2_y, v2_x)

    turn_angle = angle2 - angle1

    # Normalize to [-π, π]
    while turn_angle > math.pi:
        turn_angle -= 2 * math.pi
    while turn_angle < -math.pi:
        turn_angle += 2 * math.pi

    return abs(turn_angle)

def test_molecule(smiles, ring_sizes, expected_angles, tolerance_deg=2.0):
    """
    Test a molecule with fused small rings.

    Args:
        smiles: SMILES string of molecule
        ring_sizes: List of ring sizes to check
        expected_angles: List of expected angles in degrees for each ring
        tolerance_deg: Tolerance in degrees (default ±2°)
    """
    print(f"\nTesting: {smiles}")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("  ERROR: Failed to parse SMILES")
        return False

    # Generate 2D coordinates
    rdDepictor.Compute2DCoords(mol)

    # Get conformer
    conf = mol.GetConformer()

    # Extract 2D coordinates
    coords = []
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        coords.append(pos)

    # Find rings
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    tolerance_rad = tolerance_deg * math.pi / 180.0

    all_passed = True

    for ring_idx, ring in enumerate(rings):
        if len(ring) not in ring_sizes:
            continue

        ring_size = len(ring)
        expected_angle_deg = expected_angles[ring_sizes.index(ring_size)]
        expected_angle_rad = expected_angle_deg * math.pi / 180.0

        print(f"  Ring {ring_idx} (size {ring_size}, atoms {ring}):")
        print(f"    Expected angle: {expected_angle_deg:.1f}°")

        # Measure angles at each atom in the ring
        max_deviation = 0.0
        for i in range(len(ring)):
            # Get three consecutive atoms in the ring
            idx_curr = i
            idx_prev = (i - 1) % len(ring)
            idx_next = (i + 1) % len(ring)

            atom_prev = ring[idx_prev]
            atom_curr = ring[idx_curr]
            atom_next = ring[idx_next]

            p_prev = coords[atom_prev]
            p_curr = coords[atom_curr]
            p_next = coords[atom_next]

            # Compute turn angle
            angle_rad = compute_turn_angle(
                [p_prev, p_curr, p_next], 0
            )
            angle_deg = angle_rad * 180.0 / math.pi

            deviation = abs(angle_rad - expected_angle_rad)
            max_deviation = max(max_deviation, deviation)

            status = "✓" if deviation <= tolerance_rad else "✗"
            print(f"    Atom {atom_curr}: {angle_deg:.2f}° (deviation: {deviation * 180 / math.pi:.2f}°) {status}")

        passed = max_deviation <= tolerance_rad
        if passed:
            print(f"    PASSED: Max deviation {max_deviation * 180 / math.pi:.2f}° ≤ {tolerance_deg}°")
        else:
            print(f"    FAILED: Max deviation {max_deviation * 180 / math.pi:.2f}° > {tolerance_deg}°")
            all_passed = False

    return all_passed

def main():
    print("=" * 60)
    print("Small Ring Regularization Test")
    print("=" * 60)

    # Test 1: Hexagon fused to 16-ring (macrocycle)
    # Expected: 120° angles in hexagon
    test1 = test_molecule(
        "C1CCCCCCCCCCCCCCC2CCCCC12",  # 16-ring with fused hexagon (5 atoms)
        [6],  # Check 6-membered rings
        [120.0],  # Expected 120° for hexagon
        tolerance_deg=2.0
    )

    # Test 2: Pentagon fused to 15-ring (macrocycle)
    # Expected: 108° angles in pentagon
    test2 = test_molecule(
        "C1CCCCCCCCCCCCCC2CCCC12",  # 15-ring with fused pentagon (4 atoms)
        [5],  # Check 5-membered rings
        [108.0],  # Expected 108° for pentagon
        tolerance_deg=2.0
    )

    # Test 3: Square fused to 14-ring (macrocycle)
    # Expected: 90° angles in square
    test3 = test_molecule(
        "C1CCCCCCCCCCCCC2CCC12",  # 14-ring with fused square (3 atoms)
        [4],  # Check 4-membered rings
        [90.0],  # Expected 90° for square
        tolerance_deg=2.0
    )

    print("\n" + "=" * 60)
    print("Test Summary")
    print("=" * 60)
    print(f"Test 1 (Hexagon):  {'PASSED' if test1 else 'FAILED'}")
    print(f"Test 2 (Pentagon): {'PASSED' if test2 else 'FAILED'}")
    print(f"Test 3 (Square):   {'PASSED' if test3 else 'FAILED'}")
    print("=" * 60)

    all_passed = test1 and test2 and test3
    sys.exit(0 if all_passed else 1)

if __name__ == "__main__":
    main()
