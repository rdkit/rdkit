#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd

# Prefer site-packages layout
site_packages = "/Users/guillaume-osmo/Github/rdkit/lib/python3.11/site-packages"
mac_build = "/Users/guillaume-osmo/Github/rdkit/mac-build"
if site_packages not in sys.path:
    sys.path.insert(0, site_packages)
if mac_build not in sys.path:
    sys.path.insert(0, mac_build)

from rdkit import Chem
from rdkit.Chem import Osmordred
from rdkit.Chem import rdOsmordred
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

TOL = 1e-10


def rel_diff(a, b):
    if pd.isna(a) and pd.isna(b):
        return 0.0
    if pd.isna(a) or pd.isna(b):
        return np.inf
    da = float(a); db = float(b)
    abs_diff = abs(da - db)
    if abs_diff <= TOL:
        return 0.0
    if abs(da) <= TOL:
        return np.inf
    return abs_diff / abs(da)


def main():
    csv_path = os.path.join(os.path.dirname(__file__), "../testosmordred.csv")
    ref = pd.read_csv(csv_path, na_values=["-","nan","NaN",""])
    smiles_col = "smiles" if "smiles" in ref.columns else ("SMILES" if "SMILES" in ref.columns else None)
    if smiles_col is None:
        raise RuntimeError("Could not find smiles column in testosmordred.csv")

    smiles_list = ref[smiles_col].tolist()

    # Get names ordering from Python CalcOsmordred
    vals, names = Osmordred.CalcOsmordred(smiles_list[0], names=True)
    n_desc = len(vals)

    mismatches = []
    for idx, smi in enumerate(smiles_list, start=1):
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        py_vals = Osmordred.CalcOsmordred(smi, names=False)
        fast_vals = rdOsmordred.CalcAllFast(m)
        if len(py_vals) != len(fast_vals):
            mismatches.append((idx, smi, "length", len(py_vals), len(fast_vals)))
            continue
        # elementwise compare with tight tolerance and NaN-equality
        diffs = [rel_diff(py_vals[i], fast_vals[i]) for i in range(len(py_vals))]
        bad_idx = [i for i, d in enumerate(diffs) if not (d == 0.0 or d == np.inf and (pd.isna(py_vals[i]) or pd.isna(fast_vals[i])))]
        if bad_idx:
            # collect a compact sample of first 10 differences
            sample = [(i, names[i], float(py_vals[i]) if not pd.isna(py_vals[i]) else np.nan,
                       float(fast_vals[i]) if not pd.isna(fast_vals[i]) else np.nan,
                       diffs[i]) for i in bad_idx[:10]]
            mismatches.append((idx, smi, "values", sample))

    if not mismatches:
        print("OK: CalcAllFast matches CalcOsmordred for all molecules within tolerance")
    else:
        print(f"Found {len(mismatches)} molecules with mismatches")
        for entry in mismatches[:10]:
            if entry[2] == "length":
                _, smi, _, lp, lf = entry
                print(f"Length mismatch for {smi}: py={lp}, fast={lf}")
            else:
                _, smi, _, sample = entry
                print(f"Value mismatches for {smi} (showing up to 10):")
                for i, nm, pv, fv, rd in sample:
                    print(f"  [{i}] {nm}: py={pv} fast={fv} rel={rd}")

if __name__ == "__main__":
    main()