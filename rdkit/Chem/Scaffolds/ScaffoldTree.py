#
# Copyright (C) 2019 Greg Landrum and T5 Informatics GmbH
#  All Rights Reserved
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

from rdkit import Chem
from rdkit.Chem import rdChemReactions


class ScaffoldTreeParams(object):
  bondBreakers = ('[!#0;R:1]-!@[!#0:2]>>[*:1]-[#0].[#0]-[*:2]', )
  includeGenericScaffolds = True
  includeGenericBondScaffolds = False
  keepOnlyFirstFragment = True
  includeScaffoldsWithoutAttachments = True
  pruneFirst = True
  flattenIsotopes = True
  flattenChirality = True
  flattenKeepLargest = True


def _addReactionsToParams(params):
  rxns = []
  for sma in params.bondBreakers:
    rxn = rdChemReactions.ReactionFromSmarts(sma)
    if rxn:
      rxns.append(rxn)
  params._breakers = tuple(rxns)


def _updateMolProps(m):
  ''' at the moment this just calls SanitizeMol.
  it's here in case we want to change that in the future
  '''
  Chem.SanitizeMol(m)


def getMolFragments(mol, params):
  """
    >>> ps = ScaffoldTreeParams()
    >>> m = Chem.MolFromSmiles('c1ccccc1CC1NC(=O)CCC1')
    >>> frags = getMolFragments(m,ps)
    >>> len(frags)
    2

    The results are 2-tuples with SMILES for the parent and then the fragment as a molecule:
    >>> frags[0]
    ('O=C1CCCC(Cc2ccccc2)N1', <rdkit.Chem.rdchem.Mol object at 0x...>)
    >>> sorted((x,Chem.MolToSmiles(y)) for x,y in frags)
    [('O=C1CCCC(Cc2ccccc2)N1', '*C1CCCC(=O)N1'), ('O=C1CCCC(Cc2ccccc2)N1', '*c1ccccc1')]

    Here's what the actual results look like:

    Setting keepOnlyFirstFragment results in us getting all the fragments with linkers:
    >>> ps.keepOnlyFirstFragment = False
    >>> frags = getMolFragments(m,ps)
    >>> len(frags)
    8
    >>> sorted((x,Chem.MolToSmiles(y)) for x,y in frags)
    [('*CC1CCCC(=O)N1', '*C*'), ('*CC1CCCC(=O)N1', '*C1CCCC(=O)N1'), ('*Cc1ccccc1', '*C*'), 
    ('*Cc1ccccc1', '*c1ccccc1'), ('O=C1CCCC(Cc2ccccc2)N1', '*C1CCCC(=O)N1'), ('O=C1CCCC(Cc2ccccc2)N1', 
    '*CC1CCCC(=O)N1'), ('O=C1CCCC(Cc2ccccc2)N1', '*Cc1ccccc1'), ('O=C1CCCC(Cc2ccccc2)N1', '*c1ccccc1')]

  """
  if not hasattr(params, '_breakers'):
    _addReactionsToParams(params)
  res = []
  stack = [mol]
  while stack:
    wmol = stack.pop(0)
    parent_smi = Chem.MolToSmiles(wmol)
    for rxn in params._breakers:
      ps = rxn.RunReactants((wmol, ))
      for p in ps:
        _updateMolProps(p[0])
        stack.append(p[0])
        res.append((parent_smi, p[0]))
        if not params.keepOnlyFirstFragment:
          _updateMolProps(p[1])
          stack.append(p[1])
          res.append((parent_smi, p[1]))
  return res


def removeAttachmentPoints(mol):
  """
  
  >>> m = Chem.MolFromSmiles('*c1ccc(*)cc1')
  >>> m.GetNumAtoms()
  8
  >>> sm = removeAttachmentPoints(m)
  >>> sm.GetNumAtoms()
  6
  >>> Chem.MolToSmiles(sm)
  'c1ccccc1'
  
  """
  res = Chem.RWMol(mol)
  aids = sorted(
    (x.GetIdx() for x in mol.GetAtoms() if x.GetAtomicNum() == 0 and x.GetDegree() == 1),
    reverse=True)
  for aid in aids:
    res.RemoveAtom(aid)
  return Chem.Mol(res)


def makeScaffoldGeneric(mol, doAtoms=True, doBonds=False):
  """

    >>> m = Chem.MolFromSmiles('*c1ncc(C(=O)O)cc1')
    >>> gm = makeScaffoldGeneric(m)
    >>> Chem.SanitizeMol(gm)
    rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
    >>> Chem.MolToSmiles(gm)
    '**(=*)*1:*:*:*(*):*:*:1'
    >>> gm2 = makeScaffoldGeneric(m,doBonds=True)
    >>> Chem.SanitizeMol(gm2)
    rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
    >>> Chem.MolToSmiles(gm2)
    '**1***(*(*)*)**1'
    >>> gm3 = makeScaffoldGeneric(m,doAtoms=False,doBonds=True)
    >>> Chem.SanitizeMol(gm3)
    rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
    >>> Chem.MolToSmiles(gm3)
    '*C1CCC(C(O)O)CN1'

    The original molecule is not affected:
    >>> Chem.MolToSmiles(m)
    '*c1ccc(C(=O)O)cn1'
    
  """
  res = Chem.Mol(mol)
  if doAtoms:
    for at in res.GetAtoms():
      at.SetAtomicNum(0)
  if doBonds:
    for bond in res.GetBonds():
      bond.SetBondType(Chem.BondType.SINGLE)
  return res


def pruneMol(mol, params):
  """

  >>> ps = ScaffoldTreeParams()
  >>> m = Chem.MolFromSmiles('O=C(O)C1C(=O)CC1')
  >>> Chem.MolToSmiles(pruneMol(m,ps))
  'O=C1CCC1'

  """
  mol = Chem.MurckoDecompose(mol)
  mol.UpdatePropertyCache()
  Chem.FastFindRings(mol)
  return mol

  # atomsToPrune = Chem.MolFromSmarts(params.atomsToPrune)
  # res = Chem.RWMol(mol)
  # remove = sorted(res.GetSubstructMatches(atomsToPrune), reverse=True)
  # while remove:
  #   for aid in remove:
  #     res.RemoveAtom(aid[0])
  #   Chem.FastFindRings(res)
  #   remove = sorted(res.GetSubstructMatches(atomsToPrune), reverse=True)
  # return res


def flattenMol(mol, params):
  """

  >>> m = Chem.MolFromSmiles('Cl.[13CH3][C@H](F)/C=C/C')
  >>> ps = ScaffoldTreeParams()
  >>> Chem.MolToSmiles(flattenMol(m,ps))
  'CC=CC(C)F'

  >>> ps = ScaffoldTreeParams()
  >>> ps.flattenIsotopes=False
  >>> Chem.MolToSmiles(flattenMol(m,ps))
  'CC=CC([13CH3])F'

  >>> ps = ScaffoldTreeParams()
  >>> ps.flattenChirality=False
  >>> Chem.MolToSmiles(flattenMol(m,ps))
  'C/C=C/[C@H](C)F'

  >>> ps = ScaffoldTreeParams()
  >>> ps.flattenKeepLargest=False
  >>> Chem.MolToSmiles(flattenMol(m,ps))
  'CC=CC(C)F.Cl'

  """
  if params.flattenKeepLargest:
    frags = sorted(((x.GetNumAtoms(), x) for x in Chem.GetMolFrags(mol, asMols=True)), reverse=True)
    mol = frags[0][1]
  else:
    mol = Chem.Mol(mol)
  for atom in mol.GetAtoms():
    if params.flattenIsotopes:
      atom.SetIsotope(0)
    if params.flattenChirality:
      if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        if atom.GetNoImplicit() and atom.GetAtomicNum() in (6, 7, 15, 16):
          atom.SetNoImplicit(False)
          atom.SetNumExplicitHs(0)
  for bond in mol.GetBonds():
    if params.flattenChirality:
      bond.SetBondDir(Chem.BondDir.NONE)
      bond.SetStereo(Chem.BondStereo.STEREONONE)
  return mol


from collections import namedtuple
NetworkEdge = namedtuple('NetworkEdge', ('start', 'end', 'type'))


class EdgeTypes:
  FragmentEdge = 1
  GenericEdge = 2
  GenericBondEdge = 3
  RemoveAttachmentEdge = 4
  InitializeEdge = 5


def generateScaffoldNetwork(inMol, params=None):
  """

  Run with default settings. The result is a 2-tuple with nodes and edges:
  >>> m = Chem.MolFromSmiles('c1ccccc1CC1NC(=O)CCC1')
  >>> nodes,edges = generateScaffoldNetwork(m)

  nodes is a list of canonical SMILES describing the nodes:
  >>> len(nodes)
  9
  >>> nodes
  ['O=C1CCCC(Cc2ccccc2)N1', '*c1ccccc1', '**1:*:*:*:*:*:1', '*1:*:*:*:*:*:1', 
   'c1ccccc1', '*C1CCCC(=O)N1', '**1****(=*)*1', '*1*****1', 'O=C1CCCCN1']

  edges is a list of NetworkEdge objects with indices into the nodes list for
  the start and end of each edge and the edge type as the last element. Edge
  types are defined in the EdgeType class

  >>> len(edges)
  8
  >>> sorted(x for x in edges if x.type==EdgeTypes.FragmentEdge)          
  [NetworkEdge(start=0, end=1, type=1), 
   NetworkEdge(start=0, end=5, type=1)]
  >>> sorted(x for x in edges if x.type==EdgeTypes.GenericEdge)           
  [NetworkEdge(start=1, end=2, type=2), 
   NetworkEdge(start=5, end=6, type=2)]
  >>> sorted(x for x in edges if x.type==EdgeTypes.RemoveAttachmentEdge)  
  [NetworkEdge(start=1, end=4, type=4), 
   NetworkEdge(start=2, end=3, type=4), 
   NetworkEdge(start=5, end=8, type=4), 
   NetworkEdge(start=6, end=7, type=4)]


  Here's a more complex example, flucloxacillin:
  >>> m = Chem.MolFromSmiles('Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12')
  >>> params = ScaffoldTreeParams()
  >>> params.includeGenericScaffolds = False
  >>> params.includeScaffoldsWithoutAttachments = False
  >>> nodes,edges = generateScaffoldNetwork(m,params)
  >>> nodes
  ['Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12', 
   'O=C(NC1C(=O)N2CCSC12)c1conc1-c1ccccc1', '*c1nocc1C(=O)NC1C(=O)N2CCSC12', 
   '*c1ccccc1', '*c1conc1-c1ccccc1', '*C1C(=O)N2CCSC12', '*c1conc1*']
  >>> len(edges)
  9
  >>> edges
  [NetworkEdge(start=0, end=1, type=5), NetworkEdge(start=1, end=2, type=1), 
   NetworkEdge(start=1, end=3, type=1), NetworkEdge(start=1, end=4, type=1), 
   NetworkEdge(start=1, end=5, type=1), NetworkEdge(start=2, end=6, type=1), 
   NetworkEdge(start=2, end=5, type=1), NetworkEdge(start=4, end=6, type=1), 
   NetworkEdge(start=4, end=3, type=1)]


  """
  if params is None:
    params = ScaffoldTreeParams()

  nodes = []
  edges = []

  ismi = Chem.MolToSmiles(inMol)
  mol = flattenMol(inMol, params)

  if params.pruneFirst:
    mol = pruneMol(mol, params)
  smi = Chem.MolToSmiles(mol)

  if smi != ismi:
    nodes.append(ismi)
    nodes.append(smi)
    edges.append(NetworkEdge(0, 1, EdgeTypes.InitializeEdge))

  frags = getMolFragments(mol, params)
  for smi, fmol in frags:
    if smi not in nodes:
      nodes.append(smi)
    fsmi = Chem.MolToSmiles(fmol)
    if fsmi not in nodes:
      nodes.append(fsmi)
    edges.append(NetworkEdge(nodes.index(smi), nodes.index(fsmi), EdgeTypes.FragmentEdge))

    if params.includeGenericScaffolds:
      gmol = makeScaffoldGeneric(fmol, doAtoms=True, doBonds=False)
      gsmi = Chem.MolToSmiles(gmol)
      if gsmi not in nodes:
        nodes.append(gsmi)
      edges.append(NetworkEdge(nodes.index(fsmi), nodes.index(gsmi), EdgeTypes.GenericEdge))

      if params.includeScaffoldsWithoutAttachments:
        asmi = Chem.MolToSmiles(removeAttachmentPoints(gmol))
        if asmi not in nodes:
          nodes.append(asmi)
        edges.append(
          NetworkEdge(nodes.index(gsmi), nodes.index(asmi), EdgeTypes.RemoveAttachmentEdge))

      if params.includeGenericBondScaffolds:
        gbmol = makeScaffoldGeneric(fmol, doAtoms=True, doBonds=True)
        gbsmi = Chem.MolToSmiles(gbmol)
        if gbsmi not in nodes:
          nodes.append(gbsmi)
        edges.append(NetworkEdge(nodes.index(fsmi), nodes.index(gbsmi), EdgeTypes.GenericBondEdge))

        if params.includeScaffoldsWithoutAttachments:
          asmi = Chem.MolToSmiles(removeAttachmentPoints(gbmol))
          if asmi not in nodes:
            nodes.append(asmi)
          edges.append(
            NetworkEdge(nodes.index(gbsmi), nodes.index(asmi), EdgeTypes.RemoveAttachmentEdge))

    if params.includeScaffoldsWithoutAttachments:
      asmi = Chem.MolToSmiles(removeAttachmentPoints(fmol))
      if asmi not in nodes:
        nodes.append(asmi)
      edges.append(NetworkEdge(nodes.index(fsmi), nodes.index(asmi),
                               EdgeTypes.RemoveAttachmentEdge))

  return nodes, edges


# ------------------------------------
#
#  doctest boilerplate
#
def _runDoctests(verbose=None):  # pragma: nocover
  import sys
  import doctest
  failed, _ = doctest.testmod(optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE,
                              verbose=verbose)
  sys.exit(failed)


if __name__ == '__main__':  # pragma: nocover
  _runDoctests()
