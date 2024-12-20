import unittest
import numpy
from rdkit.Chem.AllChem import get_rotation_angles_to_align, RotateGroup
from rdkit.Chem import AllChem
from rdkit.Chem import *


class TestRotationAngles(unittest.TestCase):

 def test_zero_rotation(self):
  """Test case where the point is already aligned along the x-axis."""
  point = numpy.array([1.0, 0.0, 0.0])
  theta_1, theta_2 = get_rotation_angles_to_align(point)
  self.assertAlmostEqual(theta_1, 0.0, msg=f"Expected 0.0 for theta_1, got {theta_1}")
  self.assertAlmostEqual(theta_2, 0.0, msg=f"Expected 0.0 for theta_2, got {theta_2}")

 def test_point_with_positive_angles(self):
  """Test case where the point requires rotation to align with the x-axis."""
  point = numpy.array([1.0, 1.0, 1.0])
  theta_1, theta_2 = get_rotation_angles_to_align(point)
  self.assertGreaterEqual(theta_1, 0, msg=f"Expected positive theta_1, got {theta_1}")
  self.assertGreaterEqual(theta_2, 0, msg=f"Expected positive theta_2, got {theta_2}")

 def test_negative_angles(self):
  """Test case where the angles result in negative rotations."""
  point = numpy.array([1.0, -1.0, 1.0])
  theta_1, theta_2 = get_rotation_angles_to_align(point)
  self.assertGreaterEqual(theta_1, 0, msg=f"Expected positive theta_1, got {theta_1}")
  self.assertLessEqual(theta_2, 0, msg=f"Expected negative theta_2, got {theta_2}")

 def test_large_angles(self):
  """Test case where the angles exceed 360 degrees."""
  point = numpy.array([1.0, 1.0, 1.0])
  theta_1, theta_2 = get_rotation_angles_to_align(point)
  self.assertLessEqual(theta_1, 2 * numpy.pi, msg=f"Expected theta_1 <= 2*pi, got {theta_1}")
  self.assertLessEqual(theta_2, 2 * numpy.pi, msg=f"Expected theta_2 <= 2*pi, got {theta_2}")

 def test_edge_case_on_axis(self):
  """Test case where the point lies exactly on one of the axes (x-axis)."""
  point = numpy.array([1.0, 0.0, 0.0])
  theta_1, theta_2 = get_rotation_angles_to_align(point)
  self.assertAlmostEqual(theta_1, 0.0, msg=f"Expected theta_1 to be 0, got {theta_1}")
  self.assertAlmostEqual(theta_2, 0.0, msg=f"Expected theta_2 to be 0, got {theta_2}")


class TestParameterizedCases(unittest.TestCase):

 def test_parameterized_cases(self):
  """Test multiple cases with different inputs."""
  test_cases = [
   (numpy.array([1.0, 0.0, 0.0]), 0.0, 0.0),  # Point already aligned
   (numpy.array([1.0, -1.0, 1.0]), 0.7853981633974483, -0.9553166181245093),  # Negative point
   (numpy.array([1.0, 1.0, 0.0]), 0.0, numpy.arctan(1)),  # Point in the x-y plane
  ]
  for point, expected_theta_1, expected_theta_2 in test_cases:
   with self.subTest(point=point):
    theta_1, theta_2 = get_rotation_angles_to_align(point)
    self.assertAlmostEqual(theta_1, expected_theta_1, msg=f"Expected {expected_theta_1} for theta_1, got {theta_1}")
    self.assertAlmostEqual(theta_2, expected_theta_2, msg=f"Expected {expected_theta_2} for theta_2, got {theta_2}")


class TestMolecule(unittest.TestCase):

 @staticmethod
 def create_molecule():
  """Create a molecule from SMILES and add hydrogen atoms."""
  smiles = "CC(C)C1=CC=CC=C1"
  mol = MolFromSmiles(smiles)
  mol = AddHs(mol)
  AllChem.EmbedMolecule(mol)
  AllChem.MMFFOptimizeMolecule(mol)
  return mol

 def test_rotation_of_bond(self):
  """Test that rotating a single bond (1-3) changes rotated atoms group positions."""
  molecule = self.create_molecule()

  conformer_before = molecule.GetConformer()
  atom_positions_before = [conformer_before.GetAtomPosition(i) for i in range(molecule.GetNumAtoms())]

  RotateGroup(molecule, 1, 3, 60)

  conformer_after = molecule.GetConformer()
  atom_positions_after = [conformer_after.GetAtomPosition(i) for i in range(molecule.GetNumAtoms())]

  self.assertTrue(numpy.allclose(atom_positions_before[1], atom_positions_after[1]), "Atom 1 moved")
  self.assertTrue(numpy.allclose(atom_positions_before[3], atom_positions_after[3]), "Atom 3 moved")

  self.assertFalse(numpy.allclose(atom_positions_before[0], atom_positions_after[0]), "Atom 0 did not move")
  self.assertFalse(numpy.allclose(atom_positions_before[2], atom_positions_after[2]), "Atom 2 did not move")


if __name__ == "__main__":
 unittest.main()
