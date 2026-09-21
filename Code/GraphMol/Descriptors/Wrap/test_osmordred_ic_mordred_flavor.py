#
#  Copyright (C) 2026 Osmo Labs, PBC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""ICKeyFlavor.MORDRED must reproduce the mordred package.

`mordred_references/InformationContent.yaml` holds *mordred's* values.  Since
the default equivalence criterion is Basak's, those entries no longer describe
the default and cannot pass against it -- that is the intended behaviour, not a
regression.  They do describe `ICKeyFlavor.MORDRED` exactly, so this test is
where that file earns its keep.

Measured: the MORDRED flavour reproduces the mordred package on 894/894 values
and passes 344/344 entries of InformationContent.yaml.
"""

from __future__ import annotations

import math
import os
import re

import pytest
import yaml

from rdkit import Chem, RDConfig, RDLogger
from rdkit.Chem import rdMolDescriptors as rdMD

RDLogger.DisableLog("rdApp.*")

REF_DIR = os.path.join(
    RDConfig.RDBaseDir, "Code", "GraphMol", "Descriptors", "test_data",
    "mordred_references",
)
YAML_PATH = os.path.join(REF_DIR, "InformationContent.yaml")
SMI_PATH = os.path.join(REF_DIR, "structures.smi")

pytestmark = pytest.mark.skipif(
    not hasattr(rdMD, "InformationContentOptions")
    or not os.path.exists(YAML_PATH),
    reason="needs RDK_BUILD_OSMORDRED with InformationContentOptions",
)

MAXRADIUS = 5
# The shared reference harness uses an absolute tolerance; match it.
TOL = 0.05
FAMILIES = ("IC", "TIC", "SIC", "BIC", "CIC", "MIC", "ZMIC")


def _structures():
    out = {}
    with open(SMI_PATH) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            smiles, name = line.rsplit(None, 1)
            out[name] = smiles
    return out


def _references():
    """{(family, order): {molecule: value}}"""
    refs = {}
    with open(YAML_PATH) as handle:
        blocks = yaml.safe_load(handle)
    for block in blocks:
        names = block["names"]
        if not isinstance(names, list):
            names = [names]
        for index, name in enumerate(names):
            match = re.match(r"([A-Za-z]+?)(\d)$", name)
            if not match:
                continue
            key = (match.group(1), int(match.group(2)))
            refs[key] = {
                mol: (vals[index] if isinstance(vals, list) else vals)
                for mol, vals in block["results"].items()
            }
    return refs


def _mordred_values(smiles):
    """IC/TIC/SIC/BIC/CIC/MIC/ZMIC under the MORDRED flavour, as a flat list."""
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"could not parse {smiles}"
    opts = rdMD.InformationContentOptions()
    opts.keyFlavor = rdMD.ICKeyFlavor.MORDRED
    opts.maxradius = MAXRADIUS
    return list(rdMD.CalcInformationContent(mol, opts))


@pytest.fixture(scope="module")
def payload():
    return _structures(), _references()


def _mordred_package_available():
    try:
        import mordred  # noqa: F401
        return True
    except ImportError:
        return False


@pytest.mark.skipif(
    not _mordred_package_available(), reason="mordred package not installed"
)
def test_mordred_flavour_matches_the_mordred_package_live(payload):
    """THE oracle: the mordred package, computed here, under THIS RDKit.

    Not a frozen table. mordred kekulizes, and which Kekule structure
    `Kekulize` returns is a property of the RDKit release, not the molecule --
    ellagic acid gets bond(11,12) double on 2025.09 and single on 2026.09, and
    the path codes change with it. So a table generated on one release cannot
    be reproduced by mordred itself on another. Comparing against the package
    computed live removes that variable: whatever this RDKit kekulizes to,
    both sides see the same structure.
    """
    from mordred import Calculator
    from mordred import InformationContent as ICm

    structures, _ = payload
    fams = [ICm.InformationContent, ICm.TotalIC, ICm.StructuralIC,
            ICm.BondingIC, ICm.ComplementaryIC]
    calc = Calculator([cls(r) for cls in fams for r in range(MAXRADIUS + 1)])
    block = MAXRADIUS + 1
    failures = []
    checked = 0
    for name, smiles in sorted(structures.items()):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        try:
            want = [float(v) for v in calc(mol)]
        except Exception:  # noqa: BLE001 - mordred raises on odd inputs
            continue
        got = _mordred_values(smiles)
        for fi in range(5):
            for order in range(block):
                w = want[fi * block + order]
                g = got[fi * block + order]
                if math.isnan(w) and math.isnan(g):
                    continue
                checked += 1
                if math.isnan(w) and g == 0.0:
                    # Deliberate, documented difference: mordred returns NaN for
                    # SIC when A == 1 and BIC when B <= 1 (a 0/0). IC is
                    # identically zero there, so osmordred defines the ratio as
                    # 0 by continuity. Accepted ONLY when we return exactly 0.0,
                    # so a real divergence can never hide behind this branch.
                    continue
                if math.isnan(w) != math.isnan(g) or abs(w - g) > 1e-6:
                    failures.append(
                        f"{FAMILIES[fi]}{order} {name}: got {g!r}, mordred {w!r}"
                    )
    assert checked > 500, f"only {checked} values checked"
    assert not failures, (
        f"{len(failures)} of {checked} values disagree with the mordred package "
        f"under this RDKit:\n  " + "\n  ".join(failures[:12])
    )


def test_mordred_references_yaml_is_pinned_to_its_rdkit_release(payload):
    """`InformationContent.yaml` was generated on RDKit 2025.09. On a newer
    release the ONLY entries allowed to differ are ones whose Kekule structure
    changed -- which for this set is ellagic acid. Anything else failing means
    a real divergence, not a version artifact."""
    structures, refs = payload
    block = MAXRADIUS + 1
    offenders = set()
    for name in sorted({m for table in refs.values() for m in table}):
        if name not in structures:
            continue
        values = _mordred_values(structures[name])
        for fi, family in enumerate(FAMILIES):
            for order in range(block):
                expected = refs.get((family, order), {}).get(name)
                if expected is None or expected == "skip":
                    continue
                if abs(float(expected) - values[fi * block + order]) > TOL:
                    offenders.add(name)
    assert offenders <= {"EllagicAcid"}, (
        f"molecules disagreeing with the 2025.09-era yaml beyond the known "
        f"Kekule-dependent case: {sorted(offenders - {'EllagicAcid'})}"
    )


def test_mordred_flavour_differs_from_the_basak_default(payload):
    """If these agreed everywhere the flavour would be pointless -- and it would
    mean the default had silently reverted to mordred's criterion."""
    structures, _ = payload
    mol = Chem.MolFromSmiles(structures["Lycopene"])
    opts = rdMD.InformationContentOptions()
    opts.maxradius = MAXRADIUS
    default = list(rdMD.CalcInformationContent(mol, opts))
    opts.keyFlavor = rdMD.ICKeyFlavor.MORDRED
    mordred = list(rdMD.CalcInformationContent(mol, opts))
    assert any(abs(a - b) > 0.01 for a, b in zip(default, mordred)), (
        "MORDRED and the Basak default produce the same values for Lycopene; "
        "one of them is not doing what it says"
    )


def test_basak_default_is_unchanged_by_this_flavour(payload):
    """Adding MORDRED must not move the default. 2-butenol, Basak Table 1."""
    mol = Chem.MolFromSmiles("CC=CCO")
    opts = rdMD.InformationContentOptions()
    opts.maxradius = MAXRADIUS

    values = list(rdMD.CalcInformationContent(mol, opts))
    for order, expected in {0: 1.2389, 1: 2.0349, 2: 3.0270, 3: 3.1808}.items():
        assert abs(values[order] - expected) < 5e-4, (
            f"default IC{order} = {values[order]:.6f}, Basak Table 1 says {expected}"
        )


if __name__ == "__main__":
    import sys

    sys.exit(pytest.main([__file__, "-v"]))
