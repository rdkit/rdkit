# $Id: MurckoScaffold.py 3672 2010-06-14 17:10:00Z landrgr1 $
#
# Created by Peter Gedeck, September 2008
#
"""
  Generation of Murcko scaffolds from a molecule
"""

from rdkit import Chem
from rdkit.Chem import AllChem

murckoTransforms = [
  AllChem.ReactionFromSmarts('[*:1]-[!#1;D1]>>[*:1][H]'),
  AllChem.ReactionFromSmarts('[*:1]-[!#1;D2]#[AD1]>>[*:1][H]'),
  AllChem.ReactionFromSmarts('[*:1]-[!#1;D2]=[AD1]>>[*:1][H]'),
  AllChem.ReactionFromSmarts('[*:1]-[!#1;D3](=[AD1])=[AD1]>>[*:1][H]')
]


def MakeScaffoldGeneric(mol):
  """ Makes a Murcko scaffold generic (i.e. all atom types->C and all bonds ->single

  >>> Chem.MolToSmiles(MakeScaffoldGeneric(Chem.MolFromSmiles('c1ccccc1')))
  'C1CCCCC1'
  >>> Chem.MolToSmiles(MakeScaffoldGeneric(Chem.MolFromSmiles('c1ncccc1')))
  'C1CCCCC1'

  The following were associated with sf.net issue 246
  >>> Chem.MolToSmiles(MakeScaffoldGeneric(Chem.MolFromSmiles('c1[nH]ccc1')))
  'C1CCCC1'
  >>> Chem.MolToSmiles(MakeScaffoldGeneric(Chem.MolFromSmiles('C1[NH2+]C1')))
  'C1CC1'
  >>> Chem.MolToSmiles(MakeScaffoldGeneric(Chem.MolFromSmiles('C1[C@](Cl)(F)O1')))
  'CC1(C)CC1'

  """
  res = Chem.Mol(mol)
  for atom in res.GetAtoms():
    if atom.GetAtomicNum() != 1:
      atom.SetAtomicNum(6)
    atom.SetIsAromatic(False)
    atom.SetIsotope(0)
    atom.SetFormalCharge(0)
    atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    atom.SetNoImplicit(0)
    atom.SetNumExplicitHs(0)
  for bond in res.GetBonds():
    bond.SetBondType(Chem.BondType.SINGLE)
    bond.SetIsAromatic(False)
  return Chem.RemoveHs(res)


murckoPatts = [
  '[!#1;D3;$([D3]-[!#1])](=[AD1])=[AD1]', '[!#1;D2;$([D2]-[!#1])]=,#[AD1]',
  '[!#1;D1;$([D1]-[!#1;!n])]'
]
murckoQ = '[' + ','.join(['$(%s)' % x for x in murckoPatts]) + ']'
murckoQ = Chem.MolFromSmarts(murckoQ)
murckoPatts = [Chem.MolFromSmarts(x) for x in murckoPatts]
aromaticNTransform = AllChem.ReactionFromSmarts('[n:1]-[D1]>>[nH:1]')


def GetScaffoldForMol(mol):
  """ Return molecule object containing scaffold of mol

  >>> m = Chem.MolFromSmiles('Cc1ccccc1')
  >>> GetScaffoldForMol(m)
  <rdkit.Chem.rdchem.Mol object at 0x...>
  >>> Chem.MolToSmiles(GetScaffoldForMol(m))
  'c1ccccc1'

  >>> m = Chem.MolFromSmiles('Cc1cc(Oc2nccc(CCC)c2)ccc1')
  >>> Chem.MolToSmiles(GetScaffoldForMol(m))
  'c1ccc(Oc2ccccn2)cc1'

  """
  if 1:
    res = Chem.MurckoDecompose(mol)
    res.ClearComputedProps()
    res.UpdatePropertyCache()
    Chem.GetSymmSSSR(res)
  else:
    # This code cannot be reached
    res = _pyGetScaffoldForMol(mol)
  return res


def _pyGetScaffoldForMol(mol):
  while mol.HasSubstructMatch(murckoQ):
    for patt in murckoPatts:
      mol = Chem.DeleteSubstructs(mol, patt)
  for atom in mol.GetAtoms():
    if atom.GetAtomicNum() == 6:
      if atom.GetNoImplicit() and atom.GetExplicitValence() < 4:
        atom.SetNoImplicit(False)
  h = Chem.MolFromSmiles('[H]')
  mol = Chem.ReplaceSubstructs(mol, Chem.MolFromSmarts('[D1;$([D1]-n)]'), h, True)[0]
  mol = Chem.RemoveHs(mol)
  return mol


def MurckoScaffoldSmiles(smiles=None, mol=None, includeChirality=False):
  """ Returns MurckScaffold Smiles from smiles

  >>> MurckoScaffoldSmiles('Cc1cc(Oc2nccc(CCC)c2)ccc1')
  'c1ccc(Oc2ccccn2)cc1'

  >>> MurckoScaffoldSmiles(mol=Chem.MolFromSmiles('Cc1cc(Oc2nccc(CCC)c2)ccc1'))
  'c1ccc(Oc2ccccn2)cc1'

  """
  if smiles:
    mol = Chem.MolFromSmiles(smiles)
  if mol is None:
    raise ValueError('No molecule provided')
  scaffold = GetScaffoldForMol(mol)
  if not scaffold:
    return None
  return Chem.MolToSmiles(scaffold, includeChirality)


def MurckoScaffoldSmilesFromSmiles(smiles, includeChirality=False):
  """ Returns MurckScaffold Smiles from smiles

  >>> MurckoScaffoldSmilesFromSmiles('Cc1cc(Oc2nccc(CCC)c2)ccc1')
  'c1ccc(Oc2ccccn2)cc1'

  """
  return MurckoScaffoldSmiles(smiles=smiles, includeChirality=includeChirality)


# ------------------------------------
#
#  doctest boilerplate
#
def _runDoctests(verbose=None):  # pragma: nocover
  import doctest
  import sys
  failed, _ = doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=verbose)
  sys.exit(failed)


if __name__ == '__main__':  # pragma: nocover
  _runDoctests()
