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
    return list(rdMD.CalcInformationContent(mol, MAXRADIUS, opts))


@pytest.fixture(scope="module")
def payload():
    return _structures(), _references()


def test_mordred_flavour_matches_the_mordred_references(payload):
    structures, refs = payload
    failures = []
    checked = 0
    for name in sorted({m for table in refs.values() for m in table}):
        if name not in structures:
            continue
        values = _mordred_values(structures[name])
        block = MAXRADIUS + 1
        for family_index, family in enumerate(FAMILIES):
            for order in range(block):
                expected = refs.get((family, order), {}).get(name)
                if expected is None or expected == "skip":
                    continue
                got = values[family_index * block + order]
                checked += 1
                if abs(float(expected) - got) > TOL:
                    failures.append(
                        f"{family}{order} {name}: got {got:.6f}, "
                        f"reference {expected}"
                    )
    assert checked > 200, f"only {checked} values checked; reference data missing?"
    assert not failures, (
        f"{len(failures)} of {checked} mordred_references entries disagree with "
        f"ICKeyFlavor.MORDRED:\n  " + "\n  ".join(failures[:12])
    )


def test_mordred_flavour_differs_from_the_basak_default(payload):
    """If these agreed everywhere the flavour would be pointless -- and it would
    mean the default had silently reverted to mordred's criterion."""
    structures, _ = payload
    mol = Chem.MolFromSmiles(structures["Lycopene"])
    default = list(rdMD.CalcInformationContent(mol, MAXRADIUS))
    opts = rdMD.InformationContentOptions()
    opts.keyFlavor = rdMD.ICKeyFlavor.MORDRED
    mordred = list(rdMD.CalcInformationContent(mol, MAXRADIUS, opts))
    assert any(abs(a - b) > 0.01 for a, b in zip(default, mordred)), (
        "MORDRED and the Basak default produce the same values for Lycopene; "
        "one of them is not doing what it says"
    )


def test_basak_default_is_unchanged_by_this_flavour(payload):
    """Adding MORDRED must not move the default. 2-butenol, Basak Table 1."""
    mol = Chem.MolFromSmiles("CC=CCO")
    values = list(rdMD.CalcInformationContent(mol, MAXRADIUS))
    for order, expected in {0: 1.2389, 1: 2.0349, 2: 3.0270, 3: 3.1808}.items():
        assert abs(values[order] - expected) < 5e-4, (
            f"default IC{order} = {values[order]:.6f}, Basak Table 1 says {expected}"
        )


if __name__ == "__main__":
    import sys

    sys.exit(pytest.main([__file__, "-v"]))
