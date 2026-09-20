#
#  Copyright (C) 2026 Osmo Labs, PBC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""Expected IC/TIC/SIC/BIC/CIC tables for both ICKeyFlavor options.

Two CSVs, over the 31 molecules of `mordred_references/structures.smi`, with
IC/TIC/SIC/BIC/CIC at r=0..5.  They do NOT have the same standing, and it
matters which is which:

`ic_expected_mordred.csv` is an **oracle**.  It is the output of the mordred
package itself, so `ICKeyFlavor.MORDRED` failing against it means this
implementation is wrong, not that the file is stale.

`ic_expected_basak.csv` is a **regression baseline**, not an oracle.  It is the
output of this implementation's default flavour.  It catches a change to the
default; it cannot tell you the default is right.  What independently validates
BASAK lives elsewhere:

  * `basak_polly_ic_413.csv` -- 413 molecules from POLLY, Basak's own software
  * Basak 1983 Table 1 -- 2-butenol IC0..IC3, reproduced exactly
  * Basak 1983 Table 2 -- ten aliphatic alcohols, IC0 and CIC1

So: if the BASAK test fails, first ask whether the default *should* have moved.
If the MORDRED test fails, the implementation is wrong.
"""

from __future__ import annotations

import csv
import math
import os

import pytest

from rdkit import Chem, RDConfig, RDLogger
from rdkit.Chem import rdMolDescriptors as rdMD

RDLogger.DisableLog("rdApp.*")

DATA_DIR = os.path.join(
    RDConfig.RDBaseDir, "Code", "GraphMol", "Descriptors", "test_data"
)
BASAK_CSV = os.path.join(DATA_DIR, "ic_expected_basak.csv")
MORDRED_CSV = os.path.join(DATA_DIR, "ic_expected_mordred.csv")

pytestmark = pytest.mark.skipif(
    not hasattr(rdMD, "InformationContentOptions") or not os.path.exists(BASAK_CSV),
    reason="needs RDK_BUILD_OSMORDRED with InformationContentOptions",
)

FAMILIES = ("IC", "TIC", "SIC", "BIC", "CIC")
MAXRADIUS = 5
# The tables carry 6 decimals; allow for the last of them.
TOL = 2e-5


def _table(path):
    with open(path) as handle:
        return list(csv.DictReader(handle))


def _computed(smiles, flavor):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    opts = rdMD.InformationContentOptions()
    opts.keyFlavor = flavor
    values = list(rdMD.CalcInformationContent(mol, MAXRADIUS, opts))
    block = MAXRADIUS + 1
    return {
        f"{fam}{order}": values[index * block + order]
        for index, fam in enumerate(FAMILIES)
        for order in range(block)
    }


def _compare(path, flavor):
    """Returns (failures, checked)."""
    failures = []
    checked = 0
    for row in _table(path):
        got = _computed(row["smiles"], flavor)
        if got is None:
            continue
        for column, expected in row.items():
            if column in ("name", "smiles"):
                continue
            want = float(expected)
            have = got[column]
            # NaN is a legitimate value here (single atoms, disconnected salts);
            # it must match NaN, and must never silently pass against a number.
            if math.isnan(want) or math.isnan(have):
                checked += 1
                if math.isnan(want) != math.isnan(have):
                    failures.append(
                        f"{column} {row['name']}: got {have!r}, expected {want!r}"
                    )
                continue
            checked += 1
            if abs(want - have) > TOL:
                failures.append(
                    f"{column} {row['name']}: got {have:.6f}, expected {want:.6f}"
                )
    return failures, checked


def test_mordred_flavor_matches_the_mordred_table():
    """ORACLE. The table is the mordred package's own output."""
    failures, checked = _compare(MORDRED_CSV, rdMD.ICKeyFlavor.MORDRED)
    assert checked > 700, f"only {checked} values checked"
    assert not failures, (
        f"{len(failures)} of {checked} values disagree with the mordred package:\n  "
        + "\n  ".join(failures[:12])
    )


def test_basak_flavor_matches_its_baseline():
    """REGRESSION BASELINE. A failure here means the default moved -- decide
    whether it should have, then regenerate the table deliberately."""
    failures, checked = _compare(BASAK_CSV, rdMD.ICKeyFlavor.BASAK)
    assert checked > 700, f"only {checked} values checked"
    assert not failures, (
        f"{len(failures)} of {checked} values differ from the BASAK baseline. "
        f"If the default was changed on purpose, regenerate the table; otherwise "
        f"this is a regression:\n  " + "\n  ".join(failures[:12])
    )


def test_the_two_tables_actually_differ():
    """If these agreed the options would be pointless, and it would mean one of
    them is not doing what it says."""
    basak = {r["name"]: r for r in _table(BASAK_CSV)}
    mordred = {r["name"]: r for r in _table(MORDRED_CSV)}
    differing = 0
    for name in set(basak) & set(mordred):
        for column in basak[name]:
            if column in ("name", "smiles"):
                continue
            a, b = float(basak[name][column]), float(mordred[name][column])
            if math.isnan(a) or math.isnan(b):
                continue
            if abs(a - b) > 1e-4:
                differing += 1
    assert differing > 100, (
        f"only {differing} values differ between the two tables; the flavours "
        f"should disagree substantially"
    )
