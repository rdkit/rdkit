#!/usr/bin/env python3
import time
import sys
import os
import numpy as np
import pandas as pd

site_packages = "/Users/guillaume-osmo/Github/rdkit/lib/python3.11/site-packages"
mac_build = "/Users/guillaume-osmo/Github/rdkit/mac-build"
if site_packages not in sys.path:
    sys.path.insert(0, site_packages)
if mac_build not in sys.path:
    sys.path.insert(0, mac_build)

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from rdkit import Chem
from rdkit.Chem import Osmordred
from rdkit.Chem import rdOsmordred


def benchmark(smiles_list, iterations=1, warmup=5):
    # pre-build molecules once
    mols = [Chem.MolFromSmiles(smi) for smi in smiles_list]

    # warmup
    for m, smi in zip(mols[:warmup], smiles_list[:warmup]):
        _ = Osmordred.CalcOsmordred(smi)
        _ = rdOsmordred.CalcAllFast(m)

    # python path
    t_py = 0.0
    for _ in range(iterations):
        t0 = time.time()
        for smi in smiles_list:
            _ = Osmordred.CalcOsmordred(smi)
        t_py += (time.time() - t0)

    # c++ aggregate
    t_fast = 0.0
    for _ in range(iterations):
        t2 = time.time()
        for m in mols:
            _ = rdOsmordred.CalcAllFast(m)
        t_fast += (time.time() - t2)

    return t_py/iterations, t_fast/iterations


def main():
    ref = pd.read_csv("../testosmordred.csv")
    smiles = ref['smiles'].tolist()
    n = len(smiles)

    print(f"Benchmarking on {n} molecules...")
    t_py, t_fast = benchmark(smiles, iterations=3)
    print(f"Python CalcOsmordred total: {t_py:.3f}s, per mol: {t_py/n:.6f}s")
    print(f"C++ CalcAllFast total: {t_fast:.3f}s, per mol: {t_fast/n:.6f}s")
    speedup = t_py / t_fast if t_fast > 0 else float('inf')
    print(f"Speedup: {speedup:.2f}x")

if __name__ == "__main__":
    main()