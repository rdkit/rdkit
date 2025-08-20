#!/usr/bin/env python3
"""
Debug ETA epsilon values and hydrogen counts for three molecules.
Run on macOS or Linux. To enable C++ ETA internal counters, set:
  OSMO_ETA_DEBUG=1
Recommended (to avoid stale site-packages):
  PYTHONNOUSERSITE=1
Optionally adjust PYTHONPATH to point to your RDKit build, e.g. mac-build/rdkit.
"""

import os
import sys
import platform
from typing import Tuple

try:
    from rdkit import Chem
    from rdkit.Chem import rdOsmordred
except Exception as e:
    print("ERROR: failed to import RDKit or rdOsmordred.\n"
          "Hint: set PYTHONPATH to your build (e.g. mac-build/rdkit) and DYLD_LIBRARY_PATH/LD_LIBRARY_PATH to your libs.")
    raise


def count_h_stats(smiles: str) -> Tuple[int, int, int]:
    """Return (hmol_atoms, hydrogens_in_hmol, excluded_CH_for_eps5)."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError(f"Failed to parse SMILES: {smiles}")
    hmol = Chem.AddHs(m)
    hmol_atoms = hmol.GetNumAtoms()
    hydrogens = 0
    excluded_ch = 0
    for a in hmol.GetAtoms():
        if a.GetAtomicNum() == 1:
            hydrogens += 1
            has_carbon_neighbor = any(n.GetAtomicNum() == 6 for n in a.GetNeighbors())
            if has_carbon_neighbor:
                excluded_ch += 1
    return hmol_atoms, hydrogens, excluded_ch


def show_case(smiles: str) -> None:
    print(f"SMILES: {smiles}")
    # Python-side counts (for quick visibility)
    hmol_atoms, h_count, excl_ch = count_h_stats(smiles)
    print(f"PY_DEBUG  hmol_atoms={hmol_atoms}  hydrogens_in_hmol={h_count}  excluded_CH_for_eps5={excl_ch}")
    # C++-computed ETA values (indices are 1-based in names; 0-based in array)
    m = Chem.MolFromSmiles(smiles)
    vals = rdOsmordred.CalcExtendedTopochemicalAtom(m)
    e1 = vals[31]   # ExtendedTopochemicalAtom_32
    e2 = vals[32]   # _33
    e5 = vals[35]   # _36
    d = vals[39]    # _40 = e2-e5
    print(f"ETA:  e1(32)={e1:.9f}  e2(33)={e2:.9f}  e5(36)={e5:.9f}  e2-e5(40)={d:.9f}")


def main():
    print(f"Platform: {platform.system()} {platform.release()} ({platform.machine()})  Python {platform.python_version()}")
    print(f"OSMO_ETA_DEBUG={'ON' if os.environ.get('OSMO_ETA_DEBUG') else 'OFF'}  PYTHONNOUSERSITE={'ON' if os.environ.get('PYTHONNOUSERSITE') else 'OFF'}")
    print()
    cases = [
        'C1=CC=CC=C1',
        'NC1=CC=C(Cl)C=C1',
        'ClCCCl',
    ]
    for smi in cases:
        show_case(smi)
        print()


if __name__ == "__main__":
    main()