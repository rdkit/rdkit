import builtins
import importlib.util
import math
import unittest
from unittest import mock

import numpy as np

from rdkit import Chem
from rdkit.Chem import CoordinationDepict


def _xy(mol):
  conf = mol.GetConformer()
  return np.array([
    (conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y)
    for i in range(mol.GetNumAtoms())
  ])


class TestCoordinationDepict(unittest.TestCase):

  def test_module_imports_without_coordgen(self):
    real_import = builtins.__import__

    def import_without_coordgen(name, globals=None, locals=None, fromlist=(), level=0):
      if name == "rdkit.Chem" and "rdCoordGen" in fromlist:
        raise ImportError("CoordGen support disabled for test")
      return real_import(name, globals, locals, fromlist, level)

    spec = importlib.util.spec_from_file_location(
      "CoordinationDepictWithoutCoordGen", CoordinationDepict.__file__)
    module = importlib.util.module_from_spec(spec)
    with mock.patch("builtins.__import__", side_effect=import_without_coordgen):
      spec.loader.exec_module(module)

    self.assertIsNone(module.rdCoordGen)

  @unittest.skipIf(CoordinationDepict.rdCoordGen is None,
                   "CoordGen support not available")
  def test_octahedral_donors_are_distributed_around_metal(self):
    mol = Chem.MolFromSmiles(
      "N->[Co+3](<-N)(<-N)(<-N)(<-N)<-N", sanitize=False)
    self.assertIsNotNone(mol)

    conf_id = CoordinationDepict.Compute2DCoordinationCoords(mol)

    self.assertEqual(conf_id, 0)
    metal = next(atom for atom in mol.GetAtoms() if atom.GetSymbol() == "Co")
    coords = _xy(mol)
    vectors = np.array([
      coords[bond.GetOtherAtomIdx(metal.GetIdx())] - coords[metal.GetIdx()]
      for bond in metal.GetBonds()
    ])
    angles = np.sort(np.arctan2(vectors[:, 1], vectors[:, 0]))
    gaps = np.diff(np.r_[angles, angles[0] + 2.0 * math.pi])
    self.assertGreater(np.min(gaps), math.radians(35.0))

  @unittest.skipIf(CoordinationDepict.rdCoordGen is None,
                   "CoordGen support not available")
  def test_chelating_ligands_form_a_regular_coordination_sphere(self):
    smiles = (
      "Clc1cc2nc3c4ccc[n]5->[Ru+2]67(<-[n]8cccc(c3nc2cc1Cl)c8c45)"
      "(<-[n]1cccc2c3nc4cc(Cl)c(Cl)cc4nc3c3ccc[n]->6c3c21)"
      "<-[n]1cccc2c3nc4cc(Cl)c(Cl)cc4nc3c3ccc[n]->7c3c21")
    mol = Chem.MolFromSmiles(smiles)
    self.assertIsNotNone(mol)
    CoordinationDepict.Compute2DCoordinationCoords(mol)

    metal = next(atom for atom in mol.GetAtoms() if atom.GetSymbol() == "Ru")
    coords = _xy(mol)
    vectors = np.array([
      coords[bond.GetOtherAtomIdx(metal.GetIdx())] - coords[metal.GetIdx()]
      for bond in metal.GetBonds()
    ])
    distances = np.linalg.norm(vectors, axis=1)
    angles = np.sort(np.arctan2(vectors[:, 1], vectors[:, 0]))
    gaps = np.diff(np.r_[angles, angles[0] + 2.0 * math.pi])
    self.assertLess(np.ptp(distances), 1.0e-6)
    self.assertGreater(np.min(gaps), math.radians(20.0))

  @unittest.skipIf(CoordinationDepict.rdCoordGen is None,
                   "CoordGen support not available")
  def test_haptic_groups_use_existing_rdkit_representation(self):
    mol = Chem.MolFromSmiles(
      "C12->[Fe+2]3456789(<-C1=C->3[CH-]->4C->5=2)"
      "<-C1=C->6[CH-]->7C->8=C->91")
    self.assertIsNotNone(mol)

    prepared = CoordinationDepict.PrepareCoordinationMolForDrawing(mol)

    dative_bonds = [
      bond for bond in prepared.GetBonds()
      if bond.GetBondType() == Chem.BondType.DATIVE
    ]
    dummy_atoms = [atom for atom in prepared.GetAtoms()
                   if atom.GetAtomicNum() == 0]
    self.assertEqual(len(dative_bonds), 2)
    self.assertEqual(len(dummy_atoms), 2)
    for bond in dative_bonds:
      self.assertTrue(bond.HasProp("_MolFileBondEndPts"))

  @unittest.skipIf(CoordinationDepict.rdCoordGen is None,
                   "CoordGen support not available")
  def test_adjacent_donors_remain_separate_coordination_sites(self):
    mol = Chem.MolFromSmiles("N1N->[Co+2]<-1", sanitize=False)
    CoordinationDepict.Compute2DCoordinationCoords(mol)

    metal = mol.GetAtomWithIdx(2)
    coords = _xy(mol)
    vectors = np.array([
      coords[bond.GetOtherAtomIdx(metal.GetIdx())] - coords[metal.GetIdx()]
      for bond in metal.GetBonds()
    ])
    self.assertGreater(np.linalg.norm(vectors[0] - vectors[1]), 0.5)

  @unittest.skipIf(CoordinationDepict.rdCoordGen is None,
                   "CoordGen support not available")
  def test_non_coordination_molecule_falls_back_to_coordgen(self):
    mol = Chem.MolFromSmiles("c1ccccc1")
    conf_id = CoordinationDepict.Compute2DCoordinationCoords(mol)
    self.assertEqual(conf_id, 0)
    self.assertEqual(mol.GetNumConformers(), 1)

  def test_invalid_arguments(self):
    with self.assertRaises(ValueError):
      CoordinationDepict.Compute2DCoordinationCoords(None)
    mol = Chem.MolFromSmiles("N->[Pt+2]<-N", sanitize=False)
    with self.assertRaises(ValueError):
      CoordinationDepict.Compute2DCoordinationCoords(mol, metalBondLength=0)

  def test_coordgen_required_for_valid_molecule(self):
    mol = Chem.MolFromSmiles("c1ccccc1")
    with mock.patch.object(CoordinationDepict, "rdCoordGen", None):
      with self.assertRaisesRegex(ImportError, "CoordGen support is required"):
        CoordinationDepict.Compute2DCoordinationCoords(mol)


if __name__ == "__main__":
  unittest.main()
