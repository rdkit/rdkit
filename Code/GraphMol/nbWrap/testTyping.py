#  Copyright (C) 2026  Steven Kearnes and other RDKit contributors
#         All Rights Reserved
"""Tests the nanobind signatures and the type stubs installed beside the modules."""

import ast
import concurrent.futures
import difflib
import importlib.machinery
import json
import os
import pathlib
import re
import subprocess
import sys
import tempfile
import unittest

import nanobind

import rdkit
from rdkit import Chem, DataStructs, RDConfig, rdBase
from rdkit.Chem import (rdChemReactions, rdEnumerateStereoisomers, rdFMCS, rdSubstructLibrary,
                        rdSynthonSpaceSearch)
from rdkit.Chem.MolStandardize import rdMolStandardize

PATTERNS = pathlib.Path(RDConfig.RDBaseDir, 'Code', 'RDBoost', 'nanobind_stub_patterns.txt')

# The committed stubs were generated with this version of nanobind. Other versions can format a
# stub differently, so the stubs are compared with freshly generated ones only under this one.
STUBGEN_VERSION = '2.15.0'

REGENERATE = ('Regenerate the stubs by installing RDKit and building the nanobind_stubs target, '
              'then commit them.')

# Modules with no committed stub, because the build the stubs come from cannot build them.
WITHOUT_STUBS = frozenset(['rdkit.Chem.Draw.rdMolDraw2DQt'])  # needs Qt


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


def compiled_modules():
  """Returns a {module name: path} dict for the compiled modules in the installed rdkit package."""
  package = pathlib.Path(rdkit.__file__).parent
  suffixes = sorted(importlib.machinery.EXTENSION_SUFFIXES, key=len, reverse=True)
  modules = {}
  for path in package.rglob('*'):
    suffix = next((s for s in suffixes if path.name.endswith(s)), None)
    if suffix:
      parts = path.relative_to(package.parent).parent.parts
      modules['.'.join(parts + (path.name[:-len(suffix)], ))] = path
  return modules


def installed_stub(name, path):
  """Returns the path of the stub installed beside a compiled module."""
  return path.with_name(name.rsplit('.', 1)[-1] + '.pyi')


def generate_stub(name, directory):
  """Generates the stub for a module with the options the nanobind_stubs target uses.

  Returns:
    The text of the generated stub.
  """
  output = pathlib.Path(directory, name + '.pyi')
  subprocess.run([
    sys.executable, '-m', 'nanobind.stubgen', '-q', '-p',
    str(PATTERNS), '-m', name, '-o',
    str(output)
  ], check=True, capture_output=True, text=True)
  return output.read_text()


def unresolved_names(stub):
  """Returns the type names that a stub's annotations quote because stubgen could not resolve them.

  nanobind names a type by its C++ name, such as RDKit::MrvWriterParams, until the module that
  registers it is imported, and stubgen quotes a name it cannot import.
  """
  annotations = []
  for node in ast.walk(ast.parse(stub)):
    if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
      args = node.args
      annotations += [arg.annotation for arg in args.posonlyargs + args.args + args.kwonlyargs]
      annotations += [args.vararg and args.vararg.annotation,
                      args.kwarg and args.kwarg.annotation, node.returns]
    elif isinstance(node, ast.AnnAssign):
      annotations.append(node.annotation)
  names = set()
  pending = [a for a in annotations if a is not None]
  while pending:
    node = pending.pop()
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
      names.add(node.value)
    elif not isinstance(node, ast.Call):
      # A call holds metadata, such as the dict(shape=..., order='C') of a NumPy array.
      pending.extend(ast.iter_child_nodes(node))
  return sorted(names)


class TestStubs(unittest.TestCase):

  @classmethod
  def setUpClass(cls):
    cls.modules = compiled_modules()

  def testEveryModuleHasAStub(self):
    self.assertTrue(self.modules)
    missing = sorted(name for name, path in self.modules.items()
                     if name not in WITHOUT_STUBS and not installed_stub(name, path).exists())
    self.assertEqual(missing, [], REGENERATE)

  def testStubsAreValidPython(self):
    for name, path in sorted(self.modules.items()):
      stub = installed_stub(name, path)
      if not stub.exists():
        continue
      with self.subTest(module=name):
        text = stub.read_text()
        ast.parse(text, filename=str(stub))
        self.assertEqual(unresolved_names(text), [])

  @unittest.skipIf(nanobind.__version__ != STUBGEN_VERSION,
                   f'the stubs were generated with nanobind {STUBGEN_VERSION}')
  def testStubsAreCurrent(self):
    names = sorted(name for name, path in self.modules.items()
                   if installed_stub(name, path).exists())
    with tempfile.TemporaryDirectory() as directory, \
        concurrent.futures.ThreadPoolExecutor(os.cpu_count()) as pool:
      generated = dict(zip(names, pool.map(lambda name: generate_stub(name, directory), names)))
    for name in names:
      with self.subTest(module=name):
        installed = installed_stub(name, self.modules[name]).read_text()
        if installed == generated[name]:
          continue
        diff = [
          line for line in difflib.unified_diff(installed.splitlines(),
                                                generated[name].splitlines(), 'installed',
                                                'generated', lineterm='')
        ]
        if not any(line.startswith('+') and not line.startswith('+++') for line in diff):
          # The committed stubs come from a build with every optional feature, so a build that
          # leaves one out registers fewer names than its stub covers.
          self.skipTest(f'{name} is built without features the committed stub covers')
        self.fail('\n'.join(diff[:60] + [REGENERATE]))


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
