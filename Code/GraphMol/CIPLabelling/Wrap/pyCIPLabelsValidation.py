# coding: utf-8
import os
import sys
from rdkit import Chem

EXPECTED_LABELS_OVERRIDES = {
  # p/m labels are filtered out on input, so they are not
  # included here

  # allene-likes
  'VS063': ['26R', '2S'],
  'VS118': [],
  'VS135': [],
  'VS154': ['6E', '7E'],
  'VS164': ['21R'],
  'VS188': [],
  'VS189': [],
  'VS190': [],
  'VS214': [],
  'VS231': [],

  # Chiralities not flagged by RDKit
  'VS132': [],
  'VS207': ['10r', '2s', '5s', '7r'],
  'VS226': ['12r', '15r', '19s', '22s', '28s', '2r', '31s', '5r'],
  'VS246': [],
  'VS247': [],
  'VS248': [],

  # This one is worrying: the missing chirality '2R'
  # causes '7S', '15S' to flip into '7r', '15s'
  'VS279': ['12s', '15s', '4r', '7r'],

  # Chiralities not flagged by RDKit (Rule 6)
  'VS280': [],
  'VS281': [],
  'VS282': [],
  'VS283': [],
  'VS284': [],
  'VS285': [],
  'VS286': [],
  'VS288': [],
  'VS289': [],
  'VS290': [],
  'VS291': [],
  'VS292': [],
  'VS293': ['11Z', '12Z', '3Z', '4Z', '7Z', '8Z'],
  'VS294': ['10S', '2S', '6R', '9R'],
  'VS295': ['10R', '2S', '6R', '9S'],
  'VS296': [],
  'VS297': [],
  'VS298': [],
  'VS299': [],
  'VS300': [],
}


def split_label_string(labels):
  """
    Splits the string of expected labels and filters out helical
    labels, which are not supported by rdkit
    """
  return {label for label in labels.split() if label[-1] not in 'mpMP'}


def supplier(fname):
  """
    Read the smiles, name and expected labels from the input file.
    If we know of an override for the expected labels, we do the
    override here.
    """
  with open(fname) as f:
    for line in f:
      smiles, name, expected, _ = line.split('\t', 3)

      try:
        expected = set(EXPECTED_LABELS_OVERRIDES[name])
      except KeyError:
        expected = split_label_string(expected)

      mol = Chem.MolFromSmiles(smiles)

      yield mol, name, expected


def getLabels(mol):
  """
    Calculates and extracts the CIP labels for the mol
    """

  Chem.rdCIPLabelling.AssignCIPLabels(mol)

  labels = set()
  for atom in mol.GetAtoms():
    try:
      label = atom.GetProp('_CIPCode')
    except KeyError:
      continue
    else:
      atom_idx = atom.GetIdx() + 1
      labels.add(f'{atom_idx}{label}')

  for bond in mol.GetBonds():
    try:
      label = bond.GetProp('_CIPCode')
    except KeyError:
      continue
    else:
      begin_idx = bond.GetBeginAtomIdx() + 1
      end_idx = bond.GetEndAtomIdx() + 1
      labels.add(f'{begin_idx}{label}')
      labels.add(f'{end_idx}{label}')

  return labels


if __name__ == '__main__':
  # The structures used for the validation come from
  # https://github.com/CIPValidationSuite/ValidationSuite
  # at commit 28d0fe05073905e74a1ba5e06b3bd6298686f6af
  fname = 'compounds.smi'
  fpath = os.path.join(os.environ['RDBASE'], 'Code', 'GraphMol', 'test_data', fname)

  failed = 0
  for mol, name, expected in supplier(fpath):
    actual = getLabels(mol)

    if actual == expected:
      print(f'{name}: PASSED')
    else:
      print(f'{name}: FAILED')
      print(f'    Expected: {sorted(expected)}')
      print(f'    Actual:   {sorted(actual)}')
      failed += 1
  print(f'Check finished: {0} molecules failed.')

  sys.exit(failed > 0)
