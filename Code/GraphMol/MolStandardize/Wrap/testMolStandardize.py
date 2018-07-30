from rdkit import RDConfig
import os, sys, math
import unittest
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Geometry import rdGeometry as geom
from rdkit.Chem.MolStandardize import rdMolStandardize, Metal, Charge


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def test1Cleanup(self):
    mol = Chem.MolFromSmiles("CCC(=O)O[Na]")
    nmol = rdMolStandardize.Cleanup(mol)
    self.assertEqual(Chem.MolToSmiles(nmol), "CCC(=O)[O-].[Na+]")

  def test2StandardizeSmiles(self):
    self.assertEqual(rdMolStandardize.StandardizeSmiles("CCC(=O)O[Na]"), "CCC(=O)[O-].[Na+]")

  def test3FragmentParent(self):
    mol = Chem.MolFromSmiles("[Na]OC(=O)c1ccccc1")
    nmol = rdMolStandardize.FragmentParent(mol)
    self.assertEqual(Chem.MolToSmiles(nmol), "O=C([O-])c1ccccc1")

  def test4ChargeParent(self):
    mol = Chem.MolFromSmiles("C[NH+](C)(C).[Cl-]")
    nmol = rdMolStandardize.ChargeParent(mol)
    self.assertEqual(Chem.MolToSmiles(nmol), "CN(C)C")

  def test4Normalize(self):
    mol = Chem.MolFromSmiles("C[N+](C)=C\C=C\[O-]")
    nmol = rdMolStandardize.Normalize(mol)
    self.assertEqual(Chem.MolToSmiles(nmol), "CN(C)C=CC=O")

  def test4Reionize(self):
    mol = Chem.MolFromSmiles("C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O")
    nmol = rdMolStandardize.Reionize(mol)
    self.assertEqual(Chem.MolToSmiles(nmol), "O=S(O)c1ccc(S(=O)(=O)[O-])cc1")

  def test5Metal(self):
    mol = Chem.MolFromSmiles("C1(CCCCC1)[Zn]Br")
    md = Metal.MetalDisconnector()
    nm = md.Disconnect(mol)
#    Metal.MetalDisconnector.Disconnect(mol)
    self.assertEqual(Chem.MolToSmiles(nm), "[Br-].[CH-]1CCCCC1.[Zn+2]")

    # test user defined metal_nof
    md.SetMetalNof(Chem.MolFromSmarts("[Li,K,Rb,Cs,Fr,Be,Mg,Ca,Sr,Ba,Ra,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Al,Ga,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi]~[N,O,F]"))
    mol2 = Chem.MolFromSmiles("CCC(=O)O[Na]")
    nm2 = md.Disconnect(mol2)
    self.assertEqual(Chem.MolToSmiles(nm2), "CCC(=O)O[Na]")

    def test6Charge(self):
        print("Hello world")

if __name__ == "__main__":
  unittest.main()
#  print(Charge.CHARGE_CORRECTIONS)
