#!/usr/bin/env python3
"""
Minimize 2D coordinates using RDKit's coordgen minimizer.

Reads CXSMILES with coordinates, runs coordgen minimizer, outputs minimized CXSMILES.
"""

import sys
from rdkit import Chem, Geometry
from rdkit.Chem import rdCoordGen

def minimize_2d_coords(cxsmiles):
    """
    Minimize 2D coordinates using coordgen.

    Args:
        cxsmiles: CXSMILES string with coordinates

    Returns:
        Minimized CXSMILES string
    """
    # Parse CXSMILES
    mol = Chem.MolFromSmiles(cxsmiles)
    if mol is None:
        raise ValueError(f"Failed to parse SMILES: {cxsmiles}")

    # Check if molecule has conformer (coordinates)
    if mol.GetNumConformers() == 0:
        raise ValueError("Input CXSMILES has no coordinates")

    # Get existing coordinates to use as starting point
    conf = mol.GetConformer()
    coord_map = {}
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        coord_map[i] = Geometry.Point2D(pos.x, pos.y)

    # Set up coordgen parameters for minimization
    # Use constrained mode to minimize existing coordinates rather than regenerate
    params = rdCoordGen.CoordGenParams()
    params.SetCoordMap(coord_map)
    params.dbg_useConstrained = True  # Allow coordinates to be optimized
    params.dbg_useFixed = False       # Don't fix them in place

    # Run coordgen - this will minimize the constrained coordinates
    rdCoordGen.AddCoords(mol, params)

    # Convert back to CXSMILES
    cxsmiles_out = Chem.MolToCXSmiles(mol)

    return cxsmiles_out


if __name__ == "__main__":
    # Example usage
    if len(sys.argv) > 1:
        input_cxsmiles = sys.argv[1]
    else:
        # Default example: benzene with non-ideal coordinates
        input_cxsmiles = "c1ccccc1 |(0,0,;1,0,;2,0.5,;2.5,1.5,;2,2.5,;1,3,)|"

    print("Input CXSMILES:")
    print(input_cxsmiles)
    print()

    try:
        output_cxsmiles = minimize_2d_coords(input_cxsmiles)
        print("Minimized CXSMILES:")
        print(output_cxsmiles)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
