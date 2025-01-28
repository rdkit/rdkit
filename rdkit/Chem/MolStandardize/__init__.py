# -*- coding: utf-8 -*-
"""
Molecule Validation and Standardization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a C++ reimplementation and exapansion of Matt Swain's MolVS tool

"""
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize


def ReorderTautomers(molecule):
  """Returns the list of the molecule's tautomers
    so that the canonical one as determined by the canonical
    scoring system in TautomerCanonicalizer appears first.

    :param molecule: An RDKit Molecule object.
    :return: A list of Molecule objects.
    """
  enumerator = rdMolStandardize.TautomerEnumerator()
  canon = enumerator.Canonicalize(molecule)
  csmi = Chem.MolToSmiles(canon)
  res = [canon]
  tauts = enumerator.Enumerate(molecule)
  smis = [Chem.MolToSmiles(x) for x in tauts]
  stpl = sorted((x, y) for x, y in zip(smis, tauts) if x != csmi)
  res += [y for x, y in stpl]
  return res
