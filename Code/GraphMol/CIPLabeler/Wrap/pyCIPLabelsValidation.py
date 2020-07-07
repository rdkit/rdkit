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
  'VS231': [],

  # Chiralities not flagged by RDKit
  'VS132': [],
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

      mol = Chem.MolFromSmiles(smiles, sanitize=False)
      Chem.SanitizeMol(mol)
      Chem.SetBondStereoFromDirections(mol)

      yield mol, name, expected


def getLabels(mol):
  """
    Calculates and extracts the CIP labels for the mol
    """

  Chem.rdCIPLabeler.AssignCIPLabels(mol)

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

  failed = []
  for mol, name, expected in supplier(fpath):
    actual = getLabels(mol)

    if actual == expected:
      print(f'{name}: PASSED')
    else:
      print(f'{name}: FAILED')
      print(f'    Expected: {sorted(expected)}')
      print(f'    Actual:   {sorted(actual)}')
      failed.append(name)
  print(f'Check finished: {len(failed)} molecules failed.')
  print(f'Failed: {", ".join(failed)}')

  sys.exit(len(failed) > 0)
