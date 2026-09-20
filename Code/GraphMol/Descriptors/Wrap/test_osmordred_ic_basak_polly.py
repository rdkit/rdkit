#
#  Copyright (C) 2026 Osmo Labs, PBC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""Validate InformationContent against POLLY -- Basak's own software.

`test_data/basak_polly_ic_413.csv` holds IC0..IC5 for 413 molecules as computed
by POLLY, the program Basak's group wrote and used.  It is the supplementary
data of:

    Majumdar, Basak et al. (2019), "Finding Needles in a Haystack: Determining
    Key Molecular Descriptors...", Molecular Informatics, minf201800164.

This matters because `mordred_references/InformationContent.yaml` encodes
*Mordred's* criterion, and Mordred is not a faithful Basak implementation.
Measured on this set (411 parseable molecules x 6 orders = 2466 values):

    osmordred (neighbour degree excluded from the key)   88.0%
    osmordred (neighbour degree included, pre-revert)    66.1%
    mordred                                              49.8%

so Mordred disagrees with Basak's own software on half of all values.

Agreement is not yet 100%.  The residual is concentrated in molecules whose
RDKit graph breaks a symmetry the molecule actually has:

    plain molecules                                      92.3%
    tautomer-ambiguous                                   88.2%
    resonance-asymmetric (nitro, carboxylate, sulfonate) 33.3%

For example RDKit writes a nitro group as [N+](=O)[O-], so the two oxygens land
in different canonical ranks even though they are equivalent; all 8 nitro
molecules in this set disagree with POLLY at one or more orders.  That is a
representation problem, not an algorithm problem, and it is tracked separately.

The thresholds below are therefore regression guards at the level currently
achieved, not statements that the implementation is finished.  Raise them when
the representation issues are fixed; never lower them silently.
"""

from __future__ import annotations

import csv
import math
import os
from collections import defaultdict

import pytest

from rdkit import Chem, RDConfig, RDLogger
from rdkit.Chem import rdMolDescriptors as rdMD

RDLogger.DisableLog("rdApp.*")

DATA = os.path.join(
    RDConfig.RDBaseDir, "Code", "GraphMol", "Descriptors", "test_data",
    "basak_polly_ic_413.csv",
)

pytestmark = pytest.mark.skipif(
    not hasattr(rdMD, "CalcInformationContent") or not os.path.exists(DATA),
    reason="needs RDK_BUILD_OSMORDRED and the POLLY reference data",
)

MAXRADIUS = 5
# POLLY values are printed to 3-4 decimals.
TOL = 0.0015

# Regression floors, per order, from the measured agreement.  IC1 is the order
# the key change moved most (80.3% here versus 16.1% with neighbour degree in
# the key), so it is the load-bearing guard.
MIN_AGREEMENT = {0: 0.90, 1: 0.75, 2: 0.78, 3: 0.82, 4: 0.85, 5: 0.87}
MIN_OVERALL = 0.85

RESONANCE_SMARTS = (
    "[N+](=O)[O-]",          # nitro
    "[CX3](=O)[O-]",         # carboxylate
    "[SX4](=O)(=O)[O-]",     # sulfonate
)


def _rows():
    with open(DATA) as handle:
        return list(csv.DictReader(handle))


def _agreement():
    """{order: (hits, total)} plus the same split by resonance asymmetry."""
    per = defaultdict(lambda: [0, 0])
    split = defaultdict(lambda: [0, 0])
    queries = [Chem.MolFromSmarts(s) for s in RESONANCE_SMARTS]
    for row in _rows():
        mol = Chem.MolFromSmiles(row["SMILES"])
        if mol is None:
            continue
        values = list(rdMD.CalcInformationContent(mol, MAXRADIUS))
        resonant = any(q is not None and mol.HasSubstructMatch(q) for q in queries)
        key = "resonance-asymmetric" if resonant else "plain"
        for order in range(MAXRADIUS + 1):
            ok = abs(values[order] - float(row[f"IC{order}"])) < TOL
            per[order][1] += 1
            split[key][1] += 1
            if ok:
                per[order][0] += 1
                split[key][0] += 1
    return per, split


@pytest.fixture(scope="module")
def agreement():
    return _agreement()


def test_reference_data_is_intact(agreement):
    per, _ = agreement
    assert per[0][1] >= 400, f"expected ~411 usable molecules, got {per[0][1]}"


@pytest.mark.parametrize("order", range(MAXRADIUS + 1))
def test_polly_agreement_per_order(agreement, order):
    per, _ = agreement
    hits, total = per[order]
    frac = hits / total
    assert frac >= MIN_AGREEMENT[order], (
        f"IC{order} agreement with POLLY fell to {frac:.1%} "
        f"({hits}/{total}); floor is {MIN_AGREEMENT[order]:.0%}. The usual cause "
        f"is a change to the equivalence key -- check generateKey."
    )


def test_polly_agreement_overall(agreement):
    per, _ = agreement
    hits = sum(h for h, _ in per.values())
    total = sum(t for _, t in per.values())
    frac = hits / total
    assert frac >= MIN_OVERALL, (
        f"overall POLLY agreement fell to {frac:.1%} ({hits}/{total}); "
        f"floor is {MIN_OVERALL:.0%}"
    )


def test_plain_molecules_agree_better_than_resonance_asymmetric(agreement):
    """Documents where the residual disagreement lives.

    RDKit writes nitro as [N+](=O)[O-] and carboxylate as C(=O)[O-], which makes
    two chemically equivalent oxygens inequivalent in the graph.  POLLY does not
    see them that way.  This is a representation difference; the test pins the
    gap so it is noticed if it widens or silently closes.
    """
    _, split = agreement
    plain = split["plain"][0] / split["plain"][1]
    assert plain >= 0.88, f"plain-molecule agreement fell to {plain:.1%}"
    if split["resonance-asymmetric"][1]:
        res = split["resonance-asymmetric"][0] / split["resonance-asymmetric"][1]
        assert res <= plain, (
            "resonance-asymmetric molecules now agree better than plain ones; "
            "the representation issue may have been fixed -- update this test"
        )


# ---------------------------------------------------------------------------
# Options.  Skipped on builds that predate InformationContentOptions.
# ---------------------------------------------------------------------------

has_options = hasattr(rdMD, "InformationContentOptions")
needs_options = pytest.mark.skipif(
    not has_options, reason="build predates InformationContentOptions"
)

# 2-butenol, from Basak's own Table 1 (Roy/Basak/Harriss/Magnuson 1983).
BUTENOL = "CC=CCO"
BASAK_TABLE1 = {0: 1.2389, 1: 2.0349, 2: 3.0270, 3: 3.1808}


@needs_options
def test_default_options_are_basak():
    """The default must reproduce Basak's published worked example."""
    mol = Chem.MolFromSmiles(BUTENOL)
    values = list(rdMD.CalcInformationContent(mol, MAXRADIUS))
    for order, expected in BASAK_TABLE1.items():
        assert abs(values[order] - expected) < 5e-4, (
            f"default options give IC{order} = {values[order]:.6f}, "
            f"Basak Table 1 says {expected}"
        )


@needs_options
def test_extended_flavour_differs_at_order_one():
    """EXTENDED is the osmordred v3 key; it lands on mordred's answer at r=1."""
    mol = Chem.MolFromSmiles(BUTENOL)
    opts = rdMD.InformationContentOptions()
    opts.keyFlavor = rdMD.ICKeyFlavor.EXTENDED
    extended = list(rdMD.CalcInformationContent(mol, MAXRADIUS, opts))
    assert abs(extended[1] - 2.4997) < 5e-4, (
        f"EXTENDED IC1 = {extended[1]:.6f}, expected mordred's 2.4997"
    )
    assert abs(extended[1] - BASAK_TABLE1[1]) > 0.4, "EXTENDED should not be Basak"


@needs_options
def test_options_round_trip():
    opts = rdMD.InformationContentOptions()
    assert opts.keyFlavor == rdMD.ICKeyFlavor.BASAK
    assert opts.aromaticHandling == rdMD.ICAromaticHandling.DISTINCT
    assert opts.vertexLabel == rdMD.ICVertexLabel.DEGREE
    assert opts.equalizeDelocalizedBonds is False
    opts.equalizeDelocalizedBonds = True
    assert opts.equalizeDelocalizedBonds is True


@needs_options
def test_delocalized_equalisation_merges_nitro_oxygens():
    """Nitro is written [N+](=O)[O-], so its two equivalent oxygens differ by
    bond order.  Equalising the delocalised bonds should lower IC (fewer
    classes) rather than raise it."""
    mol = Chem.MolFromSmiles("CC[N+](=O)[O-]")
    plain = list(rdMD.CalcInformationContent(mol, MAXRADIUS))
    opts = rdMD.InformationContentOptions()
    opts.equalizeDelocalizedBonds = True
    merged = list(rdMD.CalcInformationContent(mol, MAXRADIUS, opts))
    assert merged[1] <= plain[1] + 1e-9, (
        "equalising delocalised bonds should not increase IC1: "
        f"{plain[1]:.6f} -> {merged[1]:.6f}"
    )


@needs_options
def test_negative_radius_raises():
    """Regression: the radius used to reach initializeMatrixAndSP unguarded
    and write out of bounds, segfaulting the interpreter."""
    with pytest.raises(Exception):
        rdMD.CalcInformationContent(Chem.MolFromSmiles("c1ccccc1"), -1)


if __name__ == "__main__":
    import sys

    sys.exit(pytest.main([__file__, "-v"]))
