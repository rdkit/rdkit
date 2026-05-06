#!/usr/bin/env python3
"""
Process molecular templates and fix abnormal bond lengths using simple geometry.

Instead of force field minimization, this adjusts overly long bonds directly.
"""

import sys
import numpy as np
from rdkit import Chem, Geometry


def get_bond_info(mol):
    """Get bond lengths and identify outliers."""
    if mol.GetNumConformers() == 0:
        return [], None, None

    conf = mol.GetConformer()
    bond_info = []

    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()

        pos1 = conf.GetAtomPosition(begin_idx)
        pos2 = conf.GetAtomPosition(end_idx)

        dx = pos1.x - pos2.x
        dy = pos1.y - pos2.y
        length = np.sqrt(dx*dx + dy*dy)

        bond_info.append({
            'bond': bond,
            'begin': begin_idx,
            'end': end_idx,
            'length': length,
            'pos1': (pos1.x, pos1.y),
            'pos2': (pos2.x, pos2.y)
        })

    lengths = [b['length'] for b in bond_info]
    median_length = np.median(lengths)
    max_length = np.max(lengths)

    return bond_info, median_length, max_length


def fix_long_bonds(mol, threshold_factor=1.1, max_iterations=5):
    """
    Fix bonds that are longer than threshold by adjusting atom positions.

    Uses iterative approach: find longest bond, shorten it slightly, repeat.
    """
    for iteration in range(max_iterations):
        bond_info, median_length, max_length = get_bond_info(mol)

        if max_length <= median_length * threshold_factor:
            # All bonds are within acceptable range
            break

        # Find the longest bond
        longest_bond = max(bond_info, key=lambda b: b['length'])

        if longest_bond['length'] <= median_length * threshold_factor:
            break

        # Calculate how much to shorten it
        target_length = median_length * 1.05  # Slightly above median
        scale = target_length / longest_bond['length']

        # Move both atoms towards the midpoint
        begin_idx = longest_bond['begin']
        end_idx = longest_bond['end']

        pos1 = longest_bond['pos1']
        pos2 = longest_bond['pos2']

        # Midpoint
        mid_x = (pos1[0] + pos2[0]) / 2
        mid_y = (pos1[1] + pos2[1]) / 2

        # New positions (scaled from midpoint)
        new_pos1_x = mid_x + (pos1[0] - mid_x) * scale
        new_pos1_y = mid_y + (pos1[1] - mid_y) * scale
        new_pos2_x = mid_x + (pos2[0] - mid_x) * scale
        new_pos2_y = mid_y + (pos2[1] - mid_y) * scale

        # Update conformer
        conf = mol.GetConformer()
        conf.SetAtomPosition(begin_idx, Geometry.Point3D(new_pos1_x, new_pos1_y, 0.0))
        conf.SetAtomPosition(end_idx, Geometry.Point3D(new_pos2_x, new_pos2_y, 0.0))

    return mol


def process_templates(input_file, output_file=None, threshold_factor=1.1, verbose=True):
    """Process templates file and fix structures with abnormal bond lengths."""
    if output_file is None:
        output_file = input_file.replace('.smi', '.min.smi')

    fixed_count = 0
    total_count = 0

    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line_num, line in enumerate(f_in, 1):
            line = line.strip()
            if not line:
                continue

            total_count += 1

            # Parse CXSMILES
            mol = Chem.MolFromSmiles(line)
            if mol is None:
                if verbose:
                    print(f"Line {line_num}: Failed to parse", file=sys.stderr)
                f_out.write(line + '\n')
                continue

            # Get bond info
            bond_info, median_length, max_length = get_bond_info(mol)

            if median_length is None:
                f_out.write(line + '\n')
                continue

            needs_fix = max_length > median_length * threshold_factor

            if needs_fix:
                if verbose:
                    print(f"Line {line_num}: median={median_length:.3f}, max={max_length:.3f} - FIXING")

                try:
                    mol = fix_long_bonds(mol, threshold_factor)

                    # Verify improvement
                    _, new_median, new_max = get_bond_info(mol)
                    if verbose:
                        print(f"         After: median={new_median:.3f}, max={new_max:.3f}")

                    output_cxsmiles = Chem.MolToCXSmiles(mol)
                    f_out.write(output_cxsmiles + '\n')
                    fixed_count += 1

                except Exception as e:
                    if verbose:
                        print(f"Line {line_num}: Fix failed - {e}", file=sys.stderr)
                    f_out.write(line + '\n')
            else:
                if verbose:
                    print(f"Line {line_num}: median={median_length:.3f}, max={max_length:.3f} - OK")
                f_out.write(line + '\n')

    if verbose:
        print(f"\nProcessed {total_count} structures")
        print(f"Fixed {fixed_count} structures")
        print(f"Output written to: {output_file}")

    return fixed_count, total_count


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Fix templates with abnormal bond lengths")
    parser.add_argument("input", help="Input .smi file")
    parser.add_argument("-o", "--output", help="Output .smi file")
    parser.add_argument("-t", "--threshold", type=float, default=1.1,
                        help="Threshold factor above median to trigger fix (default: 1.1)")
    parser.add_argument("-q", "--quiet", action="store_true", help="Quiet mode")

    args = parser.parse_args()

    process_templates(args.input, args.output, args.threshold, verbose=not args.quiet)
