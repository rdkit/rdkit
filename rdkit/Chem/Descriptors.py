#
# Copyright (C) 2001-2017 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from collections import abc  # this won't work in python2, but we don't support that any more

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as _rdMolDescriptors
from rdkit.Chem import rdPartialCharges, rdMolDescriptors
import rdkit.Chem.ChemUtils.DescriptorUtilities as _du
from rdkit.Chem.EState.EState import (MaxEStateIndex, MinEStateIndex, MaxAbsEStateIndex,
                                      MinAbsEStateIndex)
from rdkit.Chem.QED import qed


def _isCallable(thing):
  return isinstance(thing, abc.Callable) or \
              hasattr(thing, '__call__')


_descList = []


def _setupDescriptors(namespace):
  global _descList, descList
  from rdkit.Chem import GraphDescriptors, MolSurf, Lipinski, Fragments, Crippen, Descriptors3D
  from rdkit.Chem.EState import EState_VSA
  _descList.clear()

  mods = [GraphDescriptors, MolSurf, EState_VSA, Lipinski, Crippen, Fragments]

  otherMods = [Chem]

  for nm, thing in tuple(namespace.items()):
    if nm[0] != '_' and _isCallable(thing):
      _descList.append((nm, thing))

  others = []
  for mod in otherMods:
    tmp = dir(mod)
    for name in tmp:
      if name[0] != '_':
        thing = getattr(mod, name)
        if _isCallable(thing):
          others.append(name)

  for mod in mods:
    tmp = dir(mod)

    for name in tmp:
      if name[0] != '_' and name[-1] != '_' and name not in others:
        # filter out python reference implementations:
        if name[:2] == 'py' and name[2:] in tmp:
          continue
        if name == 'print_function':
          continue
        thing = getattr(mod, name)
        if _isCallable(thing):
          namespace[name] = thing
          _descList.append((name, thing))
  descList = _descList


MolWt = lambda *x, **y: _rdMolDescriptors._CalcMolWt(*x, **y)
MolWt.version = _rdMolDescriptors._CalcMolWt_version
MolWt.__doc__ = """The average molecular weight of the molecule

  >>> MolWt(Chem.MolFromSmiles('CC'))
  30.07
  >>> MolWt(Chem.MolFromSmiles('[NH4+].[Cl-]'))
  53.49...

"""


def HeavyAtomMolWt(x):
  return MolWt(x, True)


HeavyAtomMolWt.__doc__ = """The average molecular weight of the molecule ignoring hydrogens

  >>> HeavyAtomMolWt(Chem.MolFromSmiles('CC'))
  24.02...
  >>> HeavyAtomMolWt(Chem.MolFromSmiles('[NH4+].[Cl-]'))
  49.46

"""
HeavyAtomMolWt.version = "1.0.0"

ExactMolWt = lambda *x, **y: _rdMolDescriptors.CalcExactMolWt(*x, **y)
ExactMolWt.version = _rdMolDescriptors._CalcExactMolWt_version
ExactMolWt.__doc__ = """The exact molecular weight of the molecule

  >>> ExactMolWt(Chem.MolFromSmiles('CC'))
  30.04...
  >>> ExactMolWt(Chem.MolFromSmiles('[13CH3]C'))
  31.05...

"""


def NumValenceElectrons(mol):
  """ The number of valence electrons the molecule has

    >>> NumValenceElectrons(Chem.MolFromSmiles('CC'))
    14
    >>> NumValenceElectrons(Chem.MolFromSmiles('C(=O)O'))
    18
    >>> NumValenceElectrons(Chem.MolFromSmiles('C(=O)[O-]'))
    18
    >>> NumValenceElectrons(Chem.MolFromSmiles('C(=O)'))
    12

    """
  tbl = Chem.GetPeriodicTable()
  return sum(
    tbl.GetNOuterElecs(atom.GetAtomicNum()) - atom.GetFormalCharge() + atom.GetTotalNumHs()
    for atom in mol.GetAtoms())


NumValenceElectrons.version = "1.1.0"


def NumRadicalElectrons(mol):
  """ The number of radical electrons the molecule has
      (says nothing about spin state)

    >>> NumRadicalElectrons(Chem.MolFromSmiles('CC'))
    0
    >>> NumRadicalElectrons(Chem.MolFromSmiles('C[CH3]'))
    0
    >>> NumRadicalElectrons(Chem.MolFromSmiles('C[CH2]'))
    1
    >>> NumRadicalElectrons(Chem.MolFromSmiles('C[CH]'))
    2
    >>> NumRadicalElectrons(Chem.MolFromSmiles('C[C]'))
    3

    """
  return sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())


NumRadicalElectrons.version = "1.1.0"


def _ChargeDescriptors(mol, force=False):
  if not force and hasattr(mol, '_chargeDescriptors'):
    return mol._chargeDescriptors
  chgs = rdPartialCharges.ComputeGasteigerCharges(mol)
  minChg = 500.
  maxChg = -500.
  for at in mol.GetAtoms():
    chg = float(at.GetProp('_GasteigerCharge'))
    minChg = min(chg, minChg)
    maxChg = max(chg, maxChg)
  res = (minChg, maxChg)
  mol._chargeDescriptors = res
  return res


def MaxPartialCharge(mol, force=False):
  _, res = _ChargeDescriptors(mol, force)
  return res


MaxPartialCharge.version = "1.0.0"


def MinPartialCharge(mol, force=False):
  res, _ = _ChargeDescriptors(mol, force)
  return res


MinPartialCharge.version = "1.0.0"


def MaxAbsPartialCharge(mol, force=False):
  v1, v2 = _ChargeDescriptors(mol, force)
  return max(abs(v1), abs(v2))


MaxAbsPartialCharge.version = "1.0.0"


def MinAbsPartialCharge(mol, force=False):
  v1, v2 = _ChargeDescriptors(mol, force)
  return min(abs(v1), abs(v2))


MinAbsPartialCharge.version = "1.0.0"


def _FingerprintDensity(mol, func, *args, **kwargs):
  fp = func(*((mol, ) + args), **kwargs)
  if hasattr(fp, 'GetNumOnBits'):
    val = fp.GetNumOnBits()
  else:
    val = len(fp.GetNonzeroElements())
  num_heavy_atoms = mol.GetNumHeavyAtoms()
  if num_heavy_atoms == 0:
    return 0
  return float(val) / num_heavy_atoms


def FpDensityMorgan1(x):
  return _FingerprintDensity(x, _rdMolDescriptors.GetMorganFingerprint, 1)


def FpDensityMorgan2(x):
  return _FingerprintDensity(x, _rdMolDescriptors.GetMorganFingerprint, 2)


def FpDensityMorgan3(x):
  return _FingerprintDensity(x, _rdMolDescriptors.GetMorganFingerprint, 3)


_du.setDescriptorVersion('1.0.0')(FpDensityMorgan1)
_du.setDescriptorVersion('1.0.0')(FpDensityMorgan2)
_du.setDescriptorVersion('1.0.0')(FpDensityMorgan3)

if hasattr(rdMolDescriptors, 'BCUT2D'):
  names = [
    "BCUT2D_%s" % s
    for s in ('MWHI', "MWLOW", "CHGHI", "CHGLO", "LOGPHI", "LOGPLOW", "MRHI", "MRLOW")
  ]
  _du.VectorDescriptorWrapper(_rdMolDescriptors.BCUT2D, names=names, version="1.0.0",
                              namespace=locals())

_setupDescriptors(locals())

if hasattr(rdMolDescriptors, 'CalcAUTOCORR2D'):
  names = ["AUTOCORR2D_%s" % str(i + 1) for i in range(192)]
  autocorr = _du.VectorDescriptorWrapper(_rdMolDescriptors.CalcAUTOCORR2D, names=names,
                                         version="1.0.0", namespace=locals())

  def setupAUTOCorrDescriptors():
    """Adds AUTOCORR descriptors to the default descriptor lists"""
    _setupDescriptors(namespace=autocorr.namespace)


class PropertyFunctor(rdMolDescriptors.PythonPropertyFunctor):
  """Creates a python based property function that can be added to the
    global property list.  To use, subclass this class and override the
    __call__ method.  Then create an instance and add it to the
    registry.  The __call__ method should return a numeric value.

    Example:

      class NumAtoms(Descriptors.PropertyFunctor):
        def __init__(self):
          Descriptors.PropertyFunctor.__init__(self, "NumAtoms", "1.0.0")
        def __call__(self, mol):
          return mol.GetNumAtoms()

      numAtoms = NumAtoms()
      rdMolDescriptors.Properties.RegisterProperty(numAtoms)
    """

  def __init__(self, name, version):
    rdMolDescriptors.PythonPropertyFunctor.__init__(self, self, name, version)

  def __call__(self, mol):
    raise NotImplementedError("Please implement the __call__ method")


# ------------------------------------
#
#  doctest boilerplate
#
def _runDoctests(verbose=None):  # pragma: nocover
  import sys
  import doctest
  failed, _ = doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=verbose)
  sys.exit(failed)


if __name__ == '__main__':  # pragma: nocover
  _runDoctests()
