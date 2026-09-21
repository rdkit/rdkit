#  Copyright (C) 2026  Steven Kearnes and other RDKit contributors
#         All Rights Reserved
"""Tests the signatures of the nanobind wrappers."""

import ast
import json
import re
import subprocess
import sys
import unittest

from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import (rdChemReactions, rdEnumerateStereoisomers, rdFMCS, rdMolDescriptors,
                        rdSubstructLibrary, rdSynthonSpaceSearch)
from rdkit.Chem.MolStandardize import rdMolStandardize


def signature_annotations(signature):
  """Returns a {parameter: annotation} dict for one nanobind signature string.

  Raises:
    SyntaxError: If the signature names a numpy.ndarray, which nanobind writes with keyword
      arguments inside the brackets.
  """
  args = parse_signature(signature).args
  return {arg.arg: ast.unparse(arg.annotation)
          for arg in args.args + args.kwonlyargs if arg.annotation}


def parameter_annotations(func):
  """Returns a {parameter: annotation} dict for each overload of a nanobind function."""
  return [signature_annotations(signature) for signature, _, _ in func.__nb_signature__]


def return_annotations(func):
  """Returns the return annotation of each overload of a nanobind function."""
  return [ast.unparse(parse_signature(signature).returns)
          for signature, _, _ in func.__nb_signature__]


def parse_signature(signature):
  """Returns the ast.FunctionDef for one nanobind signature string."""
  # nanobind writes each default value as a \N placeholder, which is not valid Python, and leaves
  # the name out of a property's getter and setter.
  signature = re.sub(r'\\\d+', '...', re.sub(r'^def \(', 'def property(', signature))
  return ast.parse(signature + ': ...').body[0]


class TestSignatures(unittest.TestCase):

  def testParameterAnnotations(self):
    cases = [
      (Chem.MolFromSmiles, 0, 'SMILES', 'str | bytes'),
      (Chem.MolFromMolBlock, 0, 'molBlock', 'str | bytes'),
      (Chem.MolFromMolFile, 0, 'molFileName', 'str | os.PathLike'),
      (Chem.rdchem.CreateStereoGroup, 0, 'atomIds', 'collections.abc.Iterable[int]'),
      (Chem.SubstructMatchParameters.setExtraAtomCheckFunc, 0, 'func',
       'collections.abc.Callable[[rdkit.Chem.rdchem.Atom, rdkit.Chem.rdchem.Atom], bool]'),
      (Chem.AddMetadataToPNGString, 0, 'metadata', 'dict[str, str]'),
      (DataStructs.BulkTanimotoSimilarity, 1, 'bvList',
       'collections.abc.Iterable[rdkit.DataStructs.cDataStructs.ExplicitBitVect]'),
      (Chem.MolAddRecursiveQueries, 0, 'queries', 'dict[str, rdkit.Chem.rdchem.Mol]'),
      (rdSynthonSpaceSearch.SynthonSpace.SubstructureSearchIncremental, 0, 'callback',
       'collections.abc.Callable[[list[rdkit.Chem.rdchem.Mol]], bool | None]'),
      (rdSynthonSpaceSearch.ShapeBuildParams.setUserConformerGenerator, 0, 'func',
       'collections.abc.Callable[[str, int], rdkit.Chem.rdchem.Mol | None]'),
      (rdMolStandardize.TautomerEnumerator.Canonicalize, 1, 'scoreFunc',
       'collections.abc.Callable[[rdkit.Chem.rdchem.Mol], int]'),
      (rdMolStandardize.CleanupInPlace, 1, 'mols',
       'collections.abc.Iterable[rdkit.Chem.rdchem.Mol]'),
      (Chem.SetDoubleBondNeighborDirections, 0, 'conf', 'rdkit.Chem.rdchem.Conformer | None'),
      (rdChemReactions.EnumerateLibrary.__init__, 1, 'reagents',
       'collections.abc.Iterable[collections.abc.Iterable[rdkit.Chem.rdchem.Mol]]'),
      # A list the wrapper fills by index has Any elements, so a caller can preallocate it with
      # placeholders of any type.
      (rdMolDescriptors.CalcHallKierAlpha, 0, 'atomContribs', 'list[typing.Any] | None'),
      (rdMolDescriptors._CalcCrippenContribs, 0, 'atomTypeLabels', 'list[typing.Any] | None'),
    ]
    for func, overload, parameter, expected in cases:
      with self.subTest(func=func.__name__, parameter=parameter):
        self.assertEqual(parameter_annotations(func)[overload][parameter], expected)

  def testNullableReturnAnnotations(self):
    mol = 'rdkit.Chem.rdchem.Mol'
    cases = [
      (Chem.MolFromSmiles, mol + ' | None'),
      (Chem.MolFromMolBlock, mol + ' | None'),
      (Chem.AtomFromSmarts, 'rdkit.Chem.rdchem.Atom | None'),
      (Chem.SDMolSupplier.__next__, mol + ' | None'),
      (Chem.Mol.GetBondBetweenAtoms, 'rdkit.Chem.rdchem.Bond | None'),
      (Chem.ReplaceCore, mol + ' | None'),
      (rdChemReactions.ReactionFromRxnFile,
       'rdkit.Chem.rdChemReactions.ChemicalReaction | None'),
      (rdFMCS.MCSResult.queryMol.fget, mol + ' | None'),
      (rdSubstructLibrary.SubstructLibrary.GetMol, mol + ' | None'),
      (rdEnumerateStereoisomers.StereoisomerEnumerator.next, mol + ' | None'),
      # A function that never returns None keeps a plain annotation.
      (Chem.AddHs, mol),
    ]
    for func, expected in cases:
      with self.subTest(func=func.__qualname__):
        self.assertEqual(set(return_annotations(func)), {expected})

  def testNullableReturnsReturnNone(self):
    with rdBase.BlockLogs():
      self.assertIsNone(Chem.MolFromSmiles('C1CC'))
      self.assertIsNone(Chem.AtomFromSmarts('['))
    mol = Chem.MolFromSmiles('CCO')
    self.assertIsNone(mol.GetBondBetweenAtoms(0, 2))
    self.assertIsNone(Chem.ReplaceCore(mol, Chem.MolFromSmiles('c1ccccc1')))
    self.assertIsNone(rdFMCS.FindMCS([Chem.MolFromSmiles('C'), Chem.MolFromSmiles('N')]).queryMol)

  def testSignaturesNameTypesFromOtherModules(self):
    # A fresh interpreter loads only the modules these two import themselves.
    code = ('import json\n'
            'from rdkit.Chem import rdMolProcessing, rdSynthonSpaceSearch\n'
            'functions = [rdMolProcessing.GetFingerprintsForMolsInFile,\n'
            '             rdSynthonSpaceSearch.SynthonSpace.FingerprintSearch,\n'
            '             rdSynthonSpaceSearch.SynthonSpace.RascalSearch]\n'
            'print(json.dumps([[sig for sig, _, _ in f.__nb_signature__] for f in functions]))\n')
    output = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True,
                            check=True).stdout
    processing, fingerprint, rascal = json.loads(output)
    for signature in processing + fingerprint + rascal:
      # An unregistered type is rendered as its C++ name, such as RDKit::RascalMCES::RascalOptions.
      self.assertNotIn('::', signature)
    self.assertEqual(
      signature_annotations(processing[1])['generator'],
      'rdkit.Chem.rdFingerprintGenerator.FingerprintGenerator64 | None')
    self.assertEqual(
      signature_annotations(fingerprint[0])['fingerprintGenerator'],
      'rdkit.Chem.rdFingerprintGenerator.FingerprintGenerator64')
    self.assertEqual(
      signature_annotations(rascal[0])['rascalOptions'], 'rdkit.Chem.rdRascalMCES.RascalOptions')


if __name__ == '__main__':
  unittest.main()
