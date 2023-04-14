# Copyright 2017, D. E. Shaw Research. All rights reserved.
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import collections

from rdkit import Chem
from rdkit.Chem import AllChem


def GetHeterocycleReactionSmarts():
  """
    Return a list of the individual patterns for mutating individual
    atoms in aromatic rings to generate a new aromatic ring. The
    elements are collections.namedtuple objects with the following
    fields:

    SMARTS - the left side of the reaction SMARTS pattern, matches the atom to mutate
    CONVERT_FROM - the element type being converted: c, n, o, s
    CONVERT_TO - the right side of the reaction SMARTS pattern, there can be multiple destination types separated with a comma, these should map to multiple actual reaction objects
    EXAMPLE - an example aromatic ring system that SMARTS should match against, used in test cases
    NEGATIVE_EXAMPLE - an example aromatic ring system that SMART should NOT match against, used in test cases
    DESCRIPTION - a human readable description of the SMARTS pattern matching
    """
  HeteroAtomReaction = collections.namedtuple(
    'HeteroAtomReaction',
    ['SMARTS', 'CONVERT_FROM', 'CONVERT_TO', 'EXAMPLE', 'NEGATIVE_EXAMPLE', 'DESCRIPTION'])
  return [
    HeteroAtomReaction(
      "[c;x2;D2;r5:1]([n;x2;r5:2])[n;x2;r5:3]",
      "c",
      "[n:1]([n:2])[n:3]",
      "c1nc[nH]c1,c1ncn([*])c1,c1nncn([*])1",
      "",
      "aromatic carbon in 5 membered ring, bonded to 2 nitrogens",
    ),
    HeteroAtomReaction(
      "[c;x2;D2;r5:1]([n;h1;!+;!-;x2;r5;!$(*=*);!$(*[*;r5]=*):2])[n;x2;r5;!$(*=*);!$(*[*;r5]=*);!$(*[*;r5;x3]):3]",
      "c",
      "[o:1]([n;H0:2])[n;H0:3],[s:1]([n;H0:2])[n;H0:3]",
      "c1nc[nH]c1,[*]c1nc[nH]c1",
      "[*]c1ncno1,[*]n1[nH]cnc1=O,[*]c1nc2[nH]cnn2c(=O)c1[*],c1nc2cscc2[nH]1,c1c[nH]c2c1[nH]cn2,c1ccc2c(c1)[nH]cn2",
      "aromatic carbon in 5 membered ring, bonded to 1 nitrogen not bonded to another atom that is in 2 rings, bonded to 1 nitrogen with H and no charge",
    ),
    HeteroAtomReaction(
      "[c;x2;D2;r5:1]([n;h1;!+;!-;x2;r5;!$(*=*);!$(*[*;r5]=*):2])[n;x2;r5;!$(*=*);!$(*[*;r5]=*);!$(*[n;r5;x3]);$(*[c;x3;$(*[*;r6])]):3] ",
      "c",
      "[o:1]([n;H0:2])[n;H0:3],[s:1]([n;H0:2])[n;H0:3]",
      "c1ccc2c(c1)[nH]cn2",
      "[*]c1ncno1,[*]n1[nH]cnc1=O,[*]c1nc2[nH]cnn2c(=O)c1[*],c1nc2cscc2[nH]1,c1c[nH]c2c1[nH]cn2",
      "aromatic carbon in 5 membered ring, bonded to 1 nitrogen bonded to a bridgehead carbon in 6-member ring, bonded to 1 nitrogen with H and no charge",
    ),
    HeteroAtomReaction(
      "[c;x2;D2;r5:1]([c;x2;r5:2])[o,s;x2;r5:3]",
      "c",
      "[n:1]([c:2])[*:3]",
      "[*]c1ccoc1,c1ccoc1[*],n1ccoc1[*]",
      "",
      "aromatic carbon in 5 membered ring, bonded to 1 carbon and 1 oxygen or sulfur",
    ),
    HeteroAtomReaction(
      "[c;x2;D2;r5:1]([c;r5:2])[n;r5;h0;!-;!+:3]",
      "c",
      "[n:1]([c:2])[n:3]",
      "c1nc[nH]c1[*],n1cc[nH]c1[*],n1cc[nH][n+]1[*],[*]n1cc2c(c1)CCCC2,[*]N1Cc2cn([*])cc2C1,[*]c1c2c(cn1[*])CCCC2,[*]c1cn2c(c1[*])CCC2",
      "",
      "aromatic carbon in 5 membered ring, bonded to 1 carbon, bonded to 1 nitrogen with no H no charge",
    ),
    HeteroAtomReaction(
      "[c;x2;D2;r5:1]([c;x2;r5:2])[n;x2;r5;h1;+:3]",
      "c",
      "[n:1]([c:2])[n:3]",
      "c1c([*])sc[nH+]1,c1csc([*])[nH+]1",
      "",
      "aromatic carbon in 5 membered ring, bonded to 1 carbon, bonded to 1 nitrogen with H and charge",
    ),
    HeteroAtomReaction(
      "[c;x2;D2;r5:1]([c:2])[n;x2;r5;h1;!+;!-:3]",
      "c",
      "[n:1]([c:2])[n:3]",
      "[*]n1[nH]ccc1=O,[*]n1cc[nH]c1=O,[*]c1nc2[nH]ccn2c(=O)c1[*],[*]c1ccc[nH]1,[*]c1ncc[nH]1,[*][n+]1ccc[nH]1,[*][n+]1cc[nH]c1,c1ccc2c(c1)cc[nH]2",
      "[*]c1cn2ccnc2nc1[*]",
      "aromatic carbon in 5 membered ring, bonded to 1 carbon, bonded to 1 nitrogen with H and no charge",
    ),
    HeteroAtomReaction(
      "[c;x2;D2;r5:1]([c;x2;r5;!$(*=*);!$(*[*;r5]=*);!$(*[*;r5;x3]):2])[n;x2;r5;!$(*=*);!$(*[*;r5]=*);!$(*[*;r5;x3]);h1;!+;!-:3]",
      "c",
      "[o:1]([c:2])[n;H0:3],[s:1]([c:2])[n;H0:3]",
      "[*]c1ccc[nH]1,[*]c1ncc[nH]1,[*][n+]1ccc[nH]1,[*][n+]1cc[nH]c1",
      "c1ccc2c(c1)cc[nH]2,[*]n1[nH]ccc1=O,[*]n1cc[nH]c1=O,[*]c1nc2[nH]ccn2c(=O)c1[*],[*]c1cn2ccnc2nc1[*]",
      "aromatic carbon in 5 membered ring, bonded to 1 carbon, bonded to 1 nitrogen with H and no charge, no exocyclic double bonds or fused rings",
    ),
    HeteroAtomReaction(
      "[c;x2;D2;r5:1]([c;x2;r5:2])[n;h0;+;x2;r5:3]",
      "c",
      "[n:1]([c:2])[n:3]",
      "CC[n+]1ccco1,[*][n+]1ccco1,[*][n+]1ccc[nH]1",
      "",
      "aromatic carbon in 5 membered ring, bonded to 1 carbon, bonded to 1 nitrogen with no H and with charge",
    ),
    HeteroAtomReaction(
      "[c;x2;D2;r5:1]([c;x2;r5:2])[n;h0;+;x2;r5:3][n;x2;r5;h1:4]",
      "c",
      "[o:1]([c:2])[n:3][n;H0:4],[s:1]([c:2])[n:3][n;H0:4]",
      "[*][n+]1ccc[nH]1",
      "CC[n+]1ccco1,[*][n+]1ccco1",
      "aromatic carbon in 5 membered ring, bonded to 1 carbon, bonded to 1 nitrogen with no H and with charge, that nitrogen is then bonded to a nitrogen with a hydrogen",
    ),
    HeteroAtomReaction(
      "[c;x2;D2;r5:1]([c;x2;r5:2])[c;x2;r5:3][o,s;x2;r5:4]",
      "c",
      "[n:1]([c:2])[c:3][*:4]",
      "s1cccn1,o1c([*])ccn1",
      "s1cc([*])cn1",
      "aromatic carbon in 5 membered ring, bonded to 2 carbons, oxygen or sulfur in ring",
    ),
    HeteroAtomReaction(
      "[!o;!s;r5:1]1~[c;r5:2]~[c;x2;D2;r5:3]~[c;r5:4]~[!o;!s;r5:5]~1",
      "c",
      "[*:1]1~[c:2]~[n:3]~[c:4]~[*:5]~1",
      "[*]c1cc[nH]c1,[*]n1cccc1,[*]n1ccc2nccc-2c1",
      "",
      "aromatic carbon in 5 membered ring, bonded to 2 carbons, no oxygen or sulfur in ring",
    ),
    HeteroAtomReaction(
      "[n;h1;!+;!-;r5:1]1~[c;!$(*=*);r5:2]~[c;x2;D2;r5:3]~[c;!$(*=*);r5:4]~[!o;!s;r5:5]~1",
      "c",
      "[n;H0:1]1~[c:2]~[o:3]~[c:4]~[*:5]1,[n;H0:1]1~[c:2]~[s:3]~[c:4]~[*:5]~1",
      "c1cc[nH]c1,c1ccc2c(c1)cc[nH]2",
      "[*]n1cccc1,[*]n1cc([*])cc1,[*]n1ccc(=O)[nH]1",
      "aromatic carbon in 5 membered ring, bonded to 2 carbons without carbonyls, no oxygen or sulfur in ring, and a nitrogen with a hydrogen in the ring",
    ),
    HeteroAtomReaction(
      "[c;h1;D2;r6;$(a1ccc[c,n,o][c,n,o]1),$(a1cc[c,n,o][c,n,o]c1),$(a1c[c,n,o][c,n,o]cc1),$(a1[c,n,o][c,n,o]ccc1),$(a1cc[n,o]c[n,o]1),$(a1c[n,o]c[n,o]c1),$(a1[n,o]c[n,o]cc1),$(a1[n,o]ccc[n,o]1),$(a1c[n,o]cc[n,o]1),$(a1[n,o]cc[n,o]c1):1]",
      "c",
      "[n:1]",
      "c1ccccc1,c1cccnc1,c1nccnc1,[*]c1cc(=O)cco1,c1c([*])[nH]c(=O)[nH]c1=O,[*]c1c[nH]c(=O)[nH]c1=O",
      "",
      "aromatic carbon in 6 membered ring, has 1H, ring has less than 3 nitrogens",
    ),
    HeteroAtomReaction(
      "[n;x2;r5;h1;!+;!-;!$(*[o,s]):1]",
      "n",
      "[o:1],[s:1]",
      "c1n[nH]nn1,c1cc[nH]c1",
      "c1ccn([*])c1,[*]c1o[nH]c(=O)c1[*]",
      "aromatic nitrogen in 5 membered ring, 1H, no charge, not next to sulfur or oxygen",
    ),
    HeteroAtomReaction(
      "[n;x2;X2;r5;h0;!+;!-;$(*[a;r5][a;r5]):1]",
      "n",
      "[c:1]",
      "c1cnc[nH]1",
      "c1ccn([*])c1",
      "aromatic nitrogen in 5 membered ring, bonded to 2 heavy atoms, no H, no charge",
    ),
    HeteroAtomReaction(
      "",
      "n",
      "No mutations",
      "",
      "",
      "aromatic nitrogen in 5 membered ring, bonded to 3 heavy atoms, no charge",
    ),
    HeteroAtomReaction(
      "[n;x2;X3;r5;+;$(*[a;r5][a;r5]):1]",
      "n",
      "[c:1]",
      "c1coc[n+]1C",
      "",
      "aromatic nitrogen in 5 membered ring, bonded to 3 heavy atoms, charged",
    ),
    HeteroAtomReaction(
      "[n;x2;X2;r6;h0;!+;!-;$(*[a;r6][a;r6]):1]",
      "n",
      "[c:1]",
      "c1cccnc1,c1nccnc1",
      "",
      "aromatic nitrogen in 6 membered ring, bonded to 2 heavy atoms, no hydrogens, no charge",
    ),
    HeteroAtomReaction(
      "",
      "n",
      "No mutations",
      "",
      "",
      "aromatic nitrogen in 6 membered ring, bonded to 2 heavy atoms, 1H, no charge",
    ),
    HeteroAtomReaction(
      "[n;x2;r6;h1;+;$(*[a;r6][a;r6]):1]",
      "n",
      "[c:1]",
      "c1ccc2c(c1)ccc[nH+]2",
      "",
      "aromatic nitrogen in 6 membered ring, bonded to 2 heavy atoms, 1H, charged",
    ),
    HeteroAtomReaction(
      "[n;x2;r6;X3;+;$(*[a;r6][a;r6]):1]",
      "n",
      "[c:1]",
      "[*][n+]1cccc2c1cccc2",
      "",
      "aromatic nitrogen in 6 membered ring, bonded to 3 heavy atoms, charged",
    ),
    HeteroAtomReaction(
      "",
      "n",
      "No mutations",
      "",
      "",
      "aromatic nitrogen in 6 membered ring, bonded to 3 heavy atoms, no charge",
    ),
    HeteroAtomReaction(
      "[s;x2;r5;!+;!-;$(*([c;r5])[c;r5]):1]",
      "s",
      "[n;H1:1],[o:1]",
      "c1ccsc1,c1ncsc1",
      "",
      "aromatic sulfur in 5 membered ring, bonded to 2 carbons",
    ),
    HeteroAtomReaction(
      "[o;x2;r5;!+;!-;$(*([c;r5])[c;r5]):1]",
      "o",
      "[n;H1:1],[s:1]",
      "c1ccoc1,c1ncoc1",
      "",
      "aromatic oxygen in 5 membered ring, bonded to 2 carbons",
    ),
    HeteroAtomReaction(
      "[o,s;x2;r5;!+;!-:1]([c;x2;r5;!$(*=*):2])[n;x2;r5;h0;D2;!+;!-:3]",
      "o,s",
      "[c:1]([c:2])[n;H1:3]",
      "c1cnoc1,c1csnc1,c1c([*])snc1",
      "[*]n1occc1=O,[*][n+]1ccco1,[*]n1c(=O)[nH]sc1=O,[*]c1nsc(=O)o1",
      "aromatic oxygen or sulfur in 5 membered ring, bonded to 1 carbon, 1 nitrogen with no charge with no hydrogen and only bonded to 2 atoms",
    ),
    HeteroAtomReaction(
      "[o,s;x2;r5;!+;!-:1]([c;x2;r5:2])[n;x2;r5;+:3]",
      "o,s",
      "[c:1]([c:2])[n;+0:3]",
      "[*][n+]1ccco1",
      "",
      "aromatic oxygen or sulfur in 5 membered ring, bonded to 1 carbon, 1 nitrogen in ring with a charge",
    ),
    HeteroAtomReaction(
      "[s;x2;r5;!+;!-;$(*([c;r5])[n;h1;r5]):1]",
      "s",
      "[n;H1:1],[o:1]",
      "[*]c1s[nH]c(=O)c1[*]",
      "",
      "aromatic sulfur in 5 membered ring, bonded to 1 carbon, 1 nitrogen with a hydrogen",
    ),
    HeteroAtomReaction(
      "[o;x2;r5;!+;!-;$(*([c;r5])[n;h1;r5]):1]",
      "o",
      "[n;H1:1],[s:1]",
      "[*]c1o[nH]c(=O)c1[*]",
      "",
      "aromatic oxygen in 5 membered ring, bonded to 1 carbon, 1 nitrogen with a hydrogen",
    ),
    HeteroAtomReaction(
      "[s;x2;r5;!+;!-:1]([c;x2;r5:2])[n;x2;r5;h0;D3;!+;!-:3]",
      "s",
      "[n;H1:1]([c:2])[n:3],[o:1]([c:2])[n:3]",
      "[*]n1sc(=O)n([*])c1=O",
      "",
      "aromatic sulfur in 5 membered ring, bonded to 1 carbon, 1 nitrogen with no charge bonded to 3 atoms ",
    ),
    HeteroAtomReaction(
      "[o;x2;r5;!+;!-:1]([c;x2;r5:2])[n;x2;r5;h0;D3;!+;!-:3]",
      "o",
      "[n;H1:1]([c:2])[n:3],[s:1]([c:2])[n:3]",
      "[*]n1oc(=O)n([*])c1=O",
      "",
      "aromatic oxygen in 5 membered ring, bonded to 1 carbon, 1 nitrogen with no charge bonded to 3 atoms ",
    ),
    HeteroAtomReaction(
      "[o;r6:1]",
      "o",
      "[n;H1:1]",
      "[*]c1cccc(=O)o1",
      "",
      "aromatic oxygen in 6 membered ring",
    )
  ]


def _ParseReactions():
  for row in GetHeterocycleReactionSmarts():
    smarts = row.SMARTS
    if not smarts:
      continue

    for product in row.CONVERT_TO.split(','):
      reaction = smarts + '>>' + product
      yield AllChem.ReactionFromSmarts(reaction)


REACTION_CACHE = None


def GetHeterocycleReactions():
  """
    Return RDKit ChemicalReaction objects of the reaction SMARTS
    returned from GetHeterocyleReactionSmarts.
    """
  global REACTION_CACHE
  if REACTION_CACHE is None:
    REACTION_CACHE = list(_ParseReactions())
  return REACTION_CACHE


def EnumerateHeterocycles(inputmol, depth=None):
  """
    Enumerate possible relevant heterocycles on the given input
    molecule.

    >>> from rdkit.Chem.EnumerateHeterocycles import EnumerateHeterocycles
    >>> from rdkit import Chem
    >>> for smi in sorted(Chem.MolToSmiles(m) for m in EnumerateHeterocycles(Chem.MolFromSmiles('c1ccccc1'))):
    ...     print(smi)
    c1ccccc1
    c1ccncc1
    c1ccnnc1
    c1cnccn1
    c1cncnc1
    c1cnncn1
    c1cnnnc1
    c1ncncn1

    The algorithm works by mutating only one atom at a time. The depth
    parameter can be used to control the level of this recursion. For
    example, only enumerating aromatic rings that are one atom different:

    >>> for smi in sorted(Chem.MolToSmiles(m) for m in EnumerateHeterocycles(Chem.MolFromSmiles('n1ccccc1'), depth=1)):
    ...     print(smi)
    c1ccccc1
    c1ccnnc1
    c1cnccn1
    c1cncnc1
    """
  starting_points = [(inputmol, 0)]
  seen = set()
  while starting_points:
    curmol, curdepth = starting_points.pop(0)
    if depth is not None and curdepth >= depth:
      return

    for rxn in GetHeterocycleReactions():
      for newmol in rxn.RunReactants((curmol, )):
        newmol = newmol[0]
        Chem.SanitizeMol(newmol)
        newmol_smiles = Chem.MolToSmiles(newmol, isomericSmiles=True)
        if newmol_smiles in seen:
          continue

        starting_points.append((newmol, curdepth + 1))
        seen.add(newmol_smiles)
        yield newmol
