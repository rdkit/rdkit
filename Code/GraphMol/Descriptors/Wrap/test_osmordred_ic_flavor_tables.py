#
#  Copyright (C) 2026 Osmo Labs, PBC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""Expected-value tables for the ICKeyFlavor options.

Three CSVs.  They do NOT have the same standing, and it matters which is which:

`ic_expected_basak_polly.csv` is the **oracle for BASAK**.  Its expected values
are POLLY's -- the program Basak's own group wrote -- for the 319 of 411 POLLY
molecules where the compiled default agreed at every order (r=0..5, tol 0.0015)
when the table was generated (rdkit 2026.09.1pre).  It pins the part of the
default we know is right against Basak himself, so a failure here is a
regression against Basak, not against ourselves.  The remaining 92 molecules --
where the default does not yet match POLLY -- are deliberately absent; the full
comparison lives in test_osmordred_ic_basak_polly.py.  BASAK does not kekulize,
so this table does not need regenerating across RDKit releases.

`ic_expected_mordred.csv` is the **oracle for MORDRED**: the output of the
mordred package itself, over the 31 `structures.smi` molecules, IC/TIC/SIC/BIC/
CIC at r=0..5.  mordred kekulizes, and which Kekule structure RDKit returns is a
property of the release, so this table is generated under rdkit 2026.09.1pre and
MUST be regenerated per release.  Its NaNs (degenerate SIC/BIC denominators)
are mordred's; osmordred returns 0 there by a documented convention.

`ic_expected_basak.csv` is a **regression baseline**, not an oracle: the
default's own output over the same 31 molecules.  It catches a change to the
default; it cannot tell you the default is right.

So: BASAK failing against the POLLY oracle means we broke something Basak got
right.  BASAK failing against its baseline means the default moved -- decide
whether it should have.  MORDRED failing against its table means the
implementation is wrong.
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
BASAK_POLLY_CSV = os.path.join(DATA_DIR, "ic_expected_basak_polly.csv")

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
    opts.maxradius = MAXRADIUS
    values = list(rdMD.CalcInformationContent(mol, opts))
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
                if math.isnan(want) and have == 0.0:
                    # mordred's NaN for a degenerate SIC/BIC denominator (single
                    # atom, or <= 1 bond) is a 0/0 where IC is identically zero;
                    # osmordred defines that ratio as 0 by continuity. Accepted
                    # only when we return exactly 0.0 -- never for another value.
                    continue
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


@pytest.mark.skipif(
    not os.path.exists(BASAK_POLLY_CSV), reason="POLLY-validated table missing"
)
def test_basak_flavor_matches_the_polly_oracle():
    """ORACLE. Expected values are POLLY's own (Basak's software), on the POLLY
    molecules where the default agreed at every order when generated. POLLY
    prints 3-4 decimals, hence the tolerance."""
    failures = []
    checked = 0
    for row in _table(BASAK_POLLY_CSV):
        mol = Chem.MolFromSmiles(row["smiles"])
        if mol is None:
            continue
        opts = rdMD.InformationContentOptions()
        opts.maxradius = MAXRADIUS
        values = list(rdMD.CalcInformationContent(mol, opts))
        for order in range(MAXRADIUS + 1):
            want = float(row[f"IC{order}"])
            got = values[order]
            checked += 1
            if abs(want - got) > 0.0015:
                failures.append(
                    f"IC{order} POLLY#{row['polly_no']}: got {got:.4f}, POLLY {want:.4f}"
                )
    assert checked >= 1800, f"only {checked} values checked; table truncated?"
    assert not failures, (
        f"{len(failures)} of {checked} values disagree with POLLY on molecules "
        f"the default previously matched -- a regression against Basak:\n  "
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

if __name__ == "__main__":
    import sys

    sys.exit(pytest.main([__file__, "-v"]))
