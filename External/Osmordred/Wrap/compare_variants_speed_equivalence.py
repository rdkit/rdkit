#!/usr/bin/env python3
import os
import sys
import time
from typing import List, Tuple

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger

# Silence RDKit logs
RDLogger.DisableLog('rdApp.*')

# Ensure RDKit build paths are available when running directly
rdkit_site = "/Users/guillaume-osmo/Github/rdkit/lib/python3.11/site-packages"
if rdkit_site not in sys.path:
    sys.path.insert(0, rdkit_site)

from rdkit.Chem import rdOsmordred


def load_smiles(csv_path: str) -> List[str]:
    df = pd.read_csv(csv_path)
    return df['smiles'].astype(str).tolist()


def build_mols(smiles_list: List[str]) -> List[Chem.Mol]:
    mols: List[Chem.Mol] = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        mols.append(mol)
    return mols


def time_variant(func, mols: List[Chem.Mol]) -> Tuple[float, List[np.ndarray]]:
    values: List[np.ndarray] = []
    t0 = time.perf_counter()
    for mol in mols:
        if mol is None:
            values.append(np.array([np.nan], dtype=float))
            continue
        try:
            res = func(mol)
            values.append(np.array(res, dtype=float))
        except Exception:
            # Standardize failures as NaNs of length 1 to allow equal_nan logic
            values.append(np.array([np.nan], dtype=float))
    t1 = time.perf_counter()
    return t1 - t0, values


def compare_vectors(a: np.ndarray, b: np.ndarray, tol: float = 1e-6) -> bool:
    if a.shape != b.shape:
        return False
    if a.size == 0 and b.size == 0:
        return True
    return np.allclose(a, b, rtol=tol, atol=tol, equal_nan=True)


def evaluate_pair(name_a: str, name_b: str, mols: List[Chem.Mol], tol: float = 1e-6):
    func_a = getattr(rdOsmordred, name_a)
    func_b = getattr(rdOsmordred, name_b)

    # Warmup
    for mol in mols[:5]:
        try:
            _ = func_a(mol)
            _ = func_b(mol)
        except Exception:
            pass

    t_a, vals_a = time_variant(func_a, mols)
    t_b, vals_b = time_variant(func_b, mols)

    mismatches = 0
    first_mismatches = []
    for i, (va, vb) in enumerate(zip(vals_a, vals_b)):
        if not compare_vectors(va, vb, tol):
            mismatches += 1
            if len(first_mismatches) < 5:
                first_mismatches.append((i, va, vb))

    n = len(mols)
    print(f"Pair: {name_a} vs {name_b}")
    print(f"  A total: {t_a:.3f}s ({t_a/n:.6f} s/mol)")
    print(f"  B total: {t_b:.3f}s ({t_b/n:.6f} s/mol)")
    faster = name_a if t_a < t_b else name_b
    speedup = (max(t_a, t_b) / max(min(t_a, t_b), 1e-12))
    print(f"  Faster: {faster}  (~{speedup:.2f}x vs slower)")
    print(f"  Equal within tol {tol}: {n - mismatches}/{n} ({(n - mismatches)/n*100:.2f}%)")
    if mismatches:
        print("  Examples of mismatches (up to 5):")
        for idx, va, vb in first_mismatches:
            print(f"    idx={idx}, lenA={va.shape}, lenB={vb.shape}")
    print()


def main():
    csv_path = os.path.join(os.path.dirname(__file__), "..", "testosmordred.csv")
    smiles = load_smiles(csv_path)
    mols = build_mols(smiles)

    pairs = [
        ("CalcAdjacencyMatrix", "CalcAdjacencyMatrixEigen"),
        ("CalcDistanceMatrix", "CalcDistanceMatrixEigen"),
        ("CalcDetourMatrix", "CalcDetourMatrixEigen"),
        ("CalcBaryszMatrix", "CalcBaryszMatrixEigen"),
    ]

    print(f"Evaluating {len(pairs)} pairs on {len(mols)} molecules...")
    for a, b in pairs:
        evaluate_pair(a, b, mols, tol=1e-6)


if __name__ == "__main__":
    main()

