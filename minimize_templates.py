#!/usr/bin/env python3
"""
Process molecular templates and minimize those with abnormal bond lengths.

Reads templates.smi, checks bond lengths, and minimizes structures
where at least one bond is longer than the norm.
"""

import sys
import numpy as np
from rdkit import Chem


def get_bond_lengths(mol):
    """Calculate all bond lengths in a molecule."""
    if mol.GetNumConformers() == 0:
        return []

    conf = mol.GetConformer()
    bond_lengths = []

    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()

        pos1 = conf.GetAtomPosition(begin_idx)
        pos2 = conf.GetAtomPosition(end_idx)

        # Calculate 2D distance
        dx = pos1.x - pos2.x
        dy = pos1.y - pos2.y
        length = np.sqrt(dx*dx + dy*dy)
        bond_lengths.append(length)

    return bond_lengths


def needs_minimization(mol, threshold_factor=1.1):
    """
    Check if molecule needs minimization based on bond lengths.

    Returns (needs_min, norm, max_length) where:
    - needs_min: True if at least one bond is longer than norm * threshold_factor
    - norm: the median bond length (should be ~1.4)
    - max_length: the maximum bond length found
    """
    bond_lengths = get_bond_lengths(mol)

    if not bond_lengths:
        return False, None, None

    # Use median as the norm
    norm = np.median(bond_lengths)
    max_length = np.max(bond_lengths)

    # Check if any bond exceeds the threshold
    needs_min = max_length > norm * threshold_factor

    return needs_min, norm, max_length


def minimize_2d_coords(mol):
    """
    Minimize 2D coordinates using coordgen.

    Args:
        mol: RDKit molecule with conformer

    Returns:
        Minimized molecule
    """
    from rdkit import Geometry
    from rdkit.Chem import rdCoordGen

    if mol.GetNumConformers() == 0:
        raise ValueError("Molecule has no coordinates")

    # Calculate original scale (median bond length)
    bond_lengths = get_bond_lengths(mol)
    original_scale = np.median(bond_lengths) if bond_lengths else 1.5
    target_scale = 1.0  # coordgen works best with bond lengths ~1.0
    scale_factor = target_scale / original_scale

    # Get existing coordinates and scale them down
    conf = mol.GetConformer()
    coord_map = {}
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        coord_map[i] = Geometry.Point2D(pos.x * scale_factor, pos.y * scale_factor)

    # Set up coordgen parameters for minimization
    params = rdCoordGen.CoordGenParams()
    params.SetCoordMap(coord_map)
    params.dbg_useConstrained = True
    params.dbg_useFixed = False

    # Run coordgen
    rdCoordGen.AddCoords(mol, params)

    # Scale coordinates back to original scale
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Geometry.Point3D(pos.x / scale_factor, pos.y / scale_factor, 0.0))

    return mol


def process_templates(input_file, output_file=None, threshold_factor=1.1, verbose=True):
    """
    Process templates file and minimize structures with abnormal bond lengths.

    Args:
        input_file: Path to input .smi file
        output_file: Path to output .smi file (default: input_file with .min.smi)
        threshold_factor: Factor above norm to trigger minimization (default: 1.1)
        verbose: Print progress information
    """
    if output_file is None:
        output_file = input_file.replace('.smi', '.min.smi')

    minimized_count = 0
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

            # Check if minimization is needed
            needs_min, norm, max_length = needs_minimization(mol, threshold_factor)

            if needs_min:
                if verbose:
                    print(f"Line {line_num}: norm={norm:.3f}, max={max_length:.3f} - MINIMIZING")

                try:
                    mol = minimize_2d_coords(mol)
                    output_cxsmiles = Chem.MolToCXSmiles(mol)
                    f_out.write(output_cxsmiles + '\n')
                    minimized_count += 1

                    # Verify improvement
                    _, new_norm, new_max = needs_minimization(mol, threshold_factor)
                    if verbose and new_norm is not None:
                        print(f"         After: norm={new_norm:.3f}, max={new_max:.3f}")

                except Exception as e:
                    if verbose:
                        print(f"Line {line_num}: Minimization failed - {e}", file=sys.stderr)
                    f_out.write(line + '\n')
            else:
                if verbose and norm is not None:
                    print(f"Line {line_num}: norm={norm:.3f}, max={max_length:.3f} - OK")
                f_out.write(line + '\n')

    if verbose:
        print(f"\nProcessed {total_count} structures")
        print(f"Minimized {minimized_count} structures")
        print(f"Output written to: {output_file}")

    return minimized_count, total_count


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Minimize templates with abnormal bond lengths")
    parser.add_argument("input", help="Input .smi file")
    parser.add_argument("-o", "--output", help="Output .smi file")
    parser.add_argument("-t", "--threshold", type=float, default=1.1,
                        help="Threshold factor above norm to trigger minimization (default: 1.1)")
    parser.add_argument("-q", "--quiet", action="store_true", help="Quiet mode")

    args = parser.parse_args()

    process_templates(args.input, args.output, args.threshold, verbose=not args.quiet)
