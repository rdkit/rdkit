from rdkit import RDConfig
import os, sys, math
import unittest
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Geometry import rdGeometry as geom
from rdkit.Chem.MolStandardize import rdMolStandardize


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def test1Cleanup(self):
    mol = Chem.MolFromSmiles("CCC(=O)O[Na]")
    nmol = rdMolStandardize.Cleanup(mol)
    self.assertEqual(Chem.MolToSmiles(nmol), "CCC(=O)[O-].[Na+]")

  def test2StandardizeSmiles(self):
    self.assertEqual(rdMolStandardize.StandardizeSmiles("CCC(=O)O[Na]"), "CCC(=O)[O-].[Na+]")


if __name__ == "__main__":
  unittest.main()
