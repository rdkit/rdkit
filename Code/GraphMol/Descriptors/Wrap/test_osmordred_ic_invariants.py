#
#  Copyright (C) 2026 Osmo Labs, PBC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""Oracle-free invariants for the Osmordred InformationContent descriptors.

Every assertion here follows from the *definition* of an information-content
descriptor, so none of them depends on choosing Mordred, Basak or any other
implementation as the reference.  That matters: the mordred_references values
for IC/TIC/SIC/BIC/CIC pin one particular equivalence criterion, so they stop
being a valid oracle the moment the underlying algorithm is changed.  These
invariants do not.

Background
----------
IC_r is the Shannon entropy of the partition of the (explicit-H) atoms into
equivalence classes determined by their radius-r environment::

    IC_r = -sum_i (n_i/A) * log2(n_i/A)

and the other four families are *defined* from it (Basak et al.)::

    TIC_r = A * IC_r
    SIC_r = IC_r / log2(A)
    CIC_r = log2(A) - IC_r
    BIC_r = IC_r / log2(B)          # B = bond count, weighted by bond order

The load-bearing invariant is the **orbit ceiling**.  Any legitimate
equivalence criterion is defined on the molecular graph, hence is invariant
under graph automorphism, hence two atoms in the same automorphism orbit must
land in the same class at *every* radius.  Therefore::

    IC_r <= IC_orbits   for all r

where IC_orbits is the entropy of the true automorphism-orbit partition, which
RDKit gives us directly via CanonicalRankAtoms(breakTies=False).  A value above
that ceiling proves the partition is not a graph invariant, regardless of which
IC algorithm is intended.
"""

from __future__ import annotations

import math
from collections import Counter

import pytest

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as rdMD

pytestmark = pytest.mark.skipif(
    not hasattr(rdMD, "CalcInformationContent"),
    reason="RDKit built without RDK_BUILD_OSMORDRED",
)

MAXRADIUS = 5

# CalcInformationContent returns 7 blocks of (MAXRADIUS + 1) values, in this
# order (see ShannonEntropies in OsmordredMatrixAutocorrEStateFragments.cpp).
FAMILIES = ("IC", "TIC", "SIC", "BIC", "CIC", "MIC", "ZMIC")

# The 13 connected molecules of test_data/mordred_references/structures.smi that
# carry InformationContent reference values.
MOLECULES = {
    "Hexane": "CCCCCC",
    "Benzene": "c1ccccc1",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Cyanidin": "C1=CC(=C(C=C1C2=C(C=C3C(=CC(=CC3=[O+]2)O)O)O)O)O",
    "Lycopene": (
        "CC(=CCC/C(=C/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/CCC=C(C)C)"
        "\\C)\\C)\\C)/C)/C)/C)C"
    ),
    "Epicatechin": "C1[C@H]([C@H](OC2=CC(=CC(=C21)O)O)C3=CC(=C(C=C3)O)O)O",
    "Limonene": "CC1=CCC(CC1)C(=C)C",
    "Allicin": "C=CCSS(=O)CC=C",
    "Glutathione": "C(CC(=O)N[C@@H](CS)C(=O)NCC(=O)O)[C@@H](C(=O)O)N",
    "Digoxin": (
        "C[C@@H]1[C@H]([C@H](C[C@@H](O1)O[C@@H]2[C@H](O[C@H](C[C@@H]2O)O"
        "[C@@H]3[C@H](O[C@H](C[C@@H]3O)O[C@H]4CC[C@]5([C@@H](C4)CC[C@@H]6"
        "[C@@H]5C[C@H]([C@]7([C@@]6(CC[C@@H]7C8=CC(=O)OC8)O)C)O)C)C)C)O)O"
    ),
    "Capsaicin": "CC(C)/C=C/CCCCC(=O)NCC1=CC(=C(C=C1)O)OC",
    "EllagicAcid": "C1=C2C3=C(C(=C1O)O)OC(=O)C4=CC(=C(C(=C43)OC2=O)O)O",
    "Astaxanthin": (
        "CC1=C(C(C[C@@H](C1=O)O)(C)C)/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/"
        "C2=C(C(=O)[C@H](CC2(C)C)O)C)\\C)\\C)/C)/C"
    ),
}

# Automorphism-orbit ceiling, computed from RDKit alone: (A, n_orbits, IC_orbits).
# test_orbit_ceiling_table_is_current re-derives these so the table cannot rot.
ORBIT_CEILING = {
    "Hexane": (20, 6, 2.446439),
    "Benzene": (12, 2, 1.000000),
    "Caffeine": (24, 18, 3.990602),
    "Cyanidin": (32, 32, 5.000000),
    "Lycopene": (96, 34, 4.839662),
    "Epicatechin": (35, 34, 5.072140),
    "Limonene": (26, 18, 4.026987),
    "Allicin": (19, 15, 3.826875),
    "Glutathione": (37, 32, 4.939183),
    "Digoxin": (119, 97, 6.493352),
    "Capsaicin": (49, 36, 4.956247),
    "EllagicAcid": (28, 14, 3.807355),
    "Astaxanthin": (96, 35, 4.881328),
}

# Molecules/radii exempted from the ceiling check.  Empty, and it should stay
# that way: a breach is a defect under any IC definition, not a convention
# difference.  The one former entry -- EllagicAcid IC5 at 3.878783 against a
# ceiling of 3.807355 -- was the M-row overflow fixed in #43.
#
# Note these are skipped via an imperative pytest.xfail(), which does NOT run
# the assertion, so an entry here silently removes coverage rather than
# reporting an unexpected pass.  Add one only with a tracking note.
KNOWN_CEILING_VIOLATIONS: set[tuple[str, int]] = set()


def _mol(name):
    mol = Chem.MolFromSmiles(MOLECULES[name])
    assert mol is not None, f"{name} failed to parse"
    return mol


def _ic(mol, maxradius=MAXRADIUS):
    """Return {family: [value at r=0..maxradius]}."""
    values = list(rdMD.CalcInformationContent(mol, maxradius))
    block = maxradius + 1
    assert len(values) == len(FAMILIES) * block, (
        f"expected {len(FAMILIES) * block} values, got {len(values)}"
    )
    return {
        fam: values[i * block:(i + 1) * block] for i, fam in enumerate(FAMILIES)
    }


def _orbit_entropy(mol):
    """(A, n_orbits, IC of the true automorphism-orbit partition)."""
    molh = Chem.AddHs(mol)
    natoms = molh.GetNumAtoms()
    ranks = list(Chem.CanonicalRankAtoms(molh, breakTies=False))
    sizes = Counter(ranks).values()
    entropy = -sum((n / natoms) * math.log2(n / natoms) for n in sizes)
    return natoms, len(set(ranks)), entropy


@pytest.mark.parametrize("name", sorted(MOLECULES))
def test_orbit_ceiling_table_is_current(name):
    """The hard-coded ceiling must still match what RDKit derives."""
    natoms, norbits, entropy = _orbit_entropy(_mol(name))
    exp_a, exp_n, exp_ic = ORBIT_CEILING[name]
    assert natoms == exp_a
    assert norbits == exp_n
    assert entropy == pytest.approx(exp_ic, abs=1e-6)


@pytest.mark.parametrize("name", sorted(MOLECULES))
@pytest.mark.parametrize("radius", range(MAXRADIUS + 1))
def test_ic_never_exceeds_orbit_ceiling(name, radius):
    """IC_r <= IC_orbits, for any IC definition whatsoever.

    Exceeding it proves the partition is not automorphism-invariant: it has
    separated two atoms that no graph-defined criterion is allowed to separate.
    """
    if (name, radius) in KNOWN_CEILING_VIOLATIONS:
        pytest.xfail(f"known open defect: {name} IC{radius} exceeds orbit ceiling")
    mol = _mol(name)
    # Derive the ceiling at full precision.  ORBIT_CEILING is documentation and
    # is only rot-checked (to 1e-6) by test_orbit_ceiling_table_is_current --
    # comparing against its rounded value would false-positive on every
    # molecule whose IC has legitimately converged onto the ceiling.
    ceiling = _orbit_entropy(mol)[2]
    got = _ic(mol)["IC"][radius]
    assert got <= ceiling + 1e-9, (
        f"{name} IC{radius} = {got!r} exceeds the automorphism-orbit ceiling "
        f"{ceiling!r}; the radius-{radius} partition splits atoms that are "
        f"genuinely equivalent"
    )


@pytest.mark.parametrize("name", sorted(MOLECULES))
def test_ic_is_monotonic_in_radius(name):
    """Each radius subdivides the previous partition, so entropy cannot drop."""
    ic = _ic(_mol(name))["IC"]
    for radius in range(1, MAXRADIUS + 1):
        assert ic[radius] >= ic[radius - 1] - 1e-9, (
            f"{name}: IC{radius} = {ic[radius]!r} < IC{radius - 1} = "
            f"{ic[radius - 1]!r}; refinement merged classes"
        )


@pytest.mark.parametrize("name", sorted(MOLECULES))
def test_ic_is_invariant_under_atom_renumbering(name):
    """A descriptor of a graph cannot depend on how the atoms are numbered.

    Guards a property no IC algorithm may break.  Note this currently PASSES:
    the refinement is order-independent, so the orbit-ceiling violation above is
    NOT input-order dependence -- the r=5 key simply distinguishes two atoms
    that r=4 had correctly merged.  Kept so a future rework cannot regress it.
    """
    import random

    mol = _mol(name)
    reference = _ic(mol)["IC"]
    rng = random.Random(0xC0FFEE)
    for trial in range(5):
        order = list(range(mol.GetNumAtoms()))
        rng.shuffle(order)
        shuffled = _ic(Chem.RenumberAtoms(mol, order))["IC"]
        for radius in range(MAXRADIUS + 1):
            assert shuffled[radius] == pytest.approx(reference[radius], abs=1e-9), (
                f"{name} trial {trial}: IC{radius} changed under renumbering, "
                f"{reference[radius]!r} -> {shuffled[radius]!r}"
            )


@pytest.mark.parametrize("name", sorted(MOLECULES))
def test_derived_families_follow_their_definitions(name):
    """TIC/SIC/CIC are defined from IC, so they must track it exactly."""
    mol = _mol(name)
    fam = _ic(mol)
    natoms = Chem.AddHs(mol).GetNumAtoms()
    log2a = math.log2(natoms)
    for radius in range(MAXRADIUS + 1):
        ic = fam["IC"][radius]
        assert fam["TIC"][radius] == pytest.approx(natoms * ic, abs=1e-9), (
            f"{name}: TIC{radius} != A * IC{radius}"
        )
        assert fam["SIC"][radius] == pytest.approx(ic / log2a, abs=1e-9), (
            f"{name}: SIC{radius} != IC{radius} / log2(A)"
        )
        assert fam["CIC"][radius] == pytest.approx(log2a - ic, abs=1e-9), (
            f"{name}: CIC{radius} != log2(A) - IC{radius}"
        )


@pytest.mark.skip(
    reason="CalcInformationContent(mol, -1) SEGFAULTS: calcInformationContent "
    "takes a signed radius and reaches initializeMatrixAndSP with no guard, "
    "creating zero-length SP rows and then writing SP[i][0]. Verified to crash "
    "the interpreter, so it cannot be xfailed -- it would take the suite with "
    "it. Un-skip once the radius is validated."
)
def test_negative_radius_is_rejected():
    """A negative radius must raise, not write out of bounds."""
    with pytest.raises(Exception):
        rdMD.CalcInformationContent(_mol("Benzene"), -1)


if __name__ == "__main__":
    import sys

    sys.exit(pytest.main([__file__, "-v"]))
