from __future__ import print_function
import os, sys
import unittest
import copy
from rdkit import DataStructs, RDConfig
from rdkit.Geometry import rdGeometry as geom
from rdkit.Chem import rdMIF, AllChem
import math
import numpy as np
import pickle


def feq(v1, v2, tol=1.0e-4):
  return abs(v1 - v2) < tol


class testfunctor():

  def __call__(self, x, y, z):
    return 1.0

  def __call__(self, pt):
    return 1.0


class TestCase(unittest.TestCase):

  def SetUp(self):
    pass

  def test1ConstructGrid(self):
    mol = AllChem.MolFromMolFile(
      os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/MolInteractionFields/Wrap/testData/HCl.mol'),
      removeHs=False)

    grd = rdMIF.ConstructGrid(mol, margin=5.0, spacing=0.5)

    bond = mol.GetConformer().GetAtomPosition(1) - mol.GetConformer().GetAtomPosition(0)
    self.assertTrue(feq(grd.GetSpacing(), 0.5))
    self.assertTrue(grd.GetNumX() == int((abs(bond.x) + 10.0) / 0.5 + 0.5))
    self.assertTrue(grd.GetNumY() == int((abs(bond.y) + 10.0) / 0.5 + 0.5))
    self.assertTrue(grd.GetNumZ() == int((abs(bond.z) + 10.0) / 0.5 + 0.5))

  def test2CubeFiles(self):
    mol = AllChem.MolFromMolFile(
      os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/MolInteractionFields/Wrap/testData/HCl.mol'),
      removeHs=False)

    grd = geom.UniformRealValueGrid3D(5.0, 5.0, 5.0, 1.0, geom.Point3D(0.0, 0.0, 0.0))
    for i in range(grd.GetSize()):
      grd.SetVal(i, float(i / 10.0))

    rdMIF.WriteToCubeFile(
      grd,
      os.path.join(RDConfig.RDBaseDir,
                   'Code/GraphMol/MolInteractionFields/Wrap/testData/test3.cube'),
      mol)

    grd2, mol2 = rdMIF.ReadFromCubeFile(
      os.path.join(RDConfig.RDBaseDir,
                   'Code/GraphMol/MolInteractionFields/Wrap/testData/test3.cube'))

    self.assertTrue(grd.GetSize() == grd2.GetSize())

    for i in range(grd.GetSize()):
      self.assertTrue(feq(grd2.GetVal(i), float(i / 10.0)))

    self.assertTrue(grd.CompareGrids(grd2))
    self.assertTrue(mol.GetNumAtoms() == mol2.GetNumAtoms())

    for i in range(mol.GetNumAtoms()):
      self.assertTrue(mol.GetAtomWithIdx(i).GetAtomicNum() == mol2.GetAtomWithIdx(i).GetAtomicNum())
      self.assertTrue(
        feq(mol.GetConformer().GetAtomPosition(i).x,
            mol2.GetConformer().GetAtomPosition(i).x))
      self.assertTrue(
        feq(mol.GetConformer().GetAtomPosition(i).y,
            mol2.GetConformer().GetAtomPosition(i).y))
      self.assertTrue(
        feq(mol.GetConformer().GetAtomPosition(i).z,
            mol2.GetConformer().GetAtomPosition(i).z))

    rdMIF.WriteToCubeFile(
      grd,
      os.path.join(RDConfig.RDBaseDir,
                   'Code/GraphMol/MolInteractionFields/Wrap/testData/test4.cube'))

    grd3, mol3 = rdMIF.ReadFromCubeFile(
      os.path.join(RDConfig.RDBaseDir,
                   'Code/GraphMol/MolInteractionFields/Wrap/testData/test4.cube'))

    self.assertTrue(grd.GetSize() == grd3.GetSize())

    for i in range(grd.GetSize()):
      self.assertTrue(feq(grd3.GetVal(i), float(i / 10.0)))

    self.assertTrue(grd.CompareGrids(grd3))
    self.assertIsNone(mol3)

  def test3Coulomb(self):
    mol = AllChem.MolFromMolFile(
      os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/MolInteractionFields/Wrap/testData/HCl.mol'),
      removeHs=False)
    AllChem.ComputeGasteigerCharges(mol)

    conf = mol.GetConformer(0)
    charges = [a.GetDoubleProp("_GasteigerCharge") for a in mol.GetAtoms()]
    pos = conf.GetPositions()

    grd = rdMIF.ConstructGrid(mol)

    grd2 = copy.deepcopy(grd)

    coul = rdMIF.Coulomb(mol)

    rdMIF.CalculateDescriptors(grd, coul)
    coul1 = rdMIF.Coulomb(charges, pos)
    chargesMismatchingLen = charges[:-1]
    self.assertRaises(ValueError, rdMIF.Coulomb, chargesMismatchingLen, pos)
    rdMIF.CalculateDescriptors(grd2, coul1)
    self.assertTrue(grd.CompareParams(grd2))
    self.assertTrue(grd.CompareVectors(grd2))
    self.assertTrue(grd.CompareGrids(grd2))
    self.assertTrue(feq(coul(0.0, 0.0, 0.0, 1000), 0.0))
    self.assertTrue(coul(2.0, 0.0, 0.0, 1000) < 0)
    self.assertTrue(coul(-2.0, 0.0, 0.0, 1000) > 0)

    rdMIF.CalculateDescriptors(grd, rdMIF.Coulomb(mol, absVal=True))

    for i in range(grd.GetSize()):
      self.assertTrue(grd.GetVal(i) <= 0.0)

    coul1 = rdMIF.Coulomb(mol, probeCharge=-1.0)
    self.assertTrue(coul1(-2.0, 0.0, 0.0, 1000) < 0)
    self.assertTrue(coul1(2.0, 0.0, 0.0, 1000) > 0)

    coul2 = rdMIF.Coulomb(mol, probeCharge=-.5)
    self.assertTrue(coul1(-2.0, 0.0, 0.0, 1000) < coul2(-2.0, 0.0, 0.0, 1000))

    coul3 = rdMIF.Coulomb(mol, confId=0, probeCharge=1.0, absVal=False,
                          chargeKey="_GasteigerCharge", softcoreParam=0.01, cutoff=1.0)
    self.assertTrue(coul3(0.0, 0.0, 0.0, 1000) > coul3(0.1, 0.0, 0.0, 1000))
    self.assertTrue(coul3(0.66, 0.0, 0.0, 1000) > coul3(0.68, 0.0, 0.0, 1000))
    self.assertTrue(coul3(0.70, 0.0, 0.0, 1000) > coul3(0.68, 0.0, 0.0, 1000))

  def test4CoulombDielectric(self):
    mol = AllChem.MolFromMolFile(
      os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/MolInteractionFields/Wrap/testData/HCl.mol'),
      removeHs=False)
    AllChem.ComputeGasteigerCharges(mol)

    conf = mol.GetConformer(0)
    charges = [a.GetDoubleProp("_GasteigerCharge") for a in mol.GetAtoms()]
    pos = conf.GetPositions()

    grd = rdMIF.ConstructGrid(mol, confId=0)
    grd2 = rdMIF.ConstructGrid(mol, confId=0)

    couldiele = rdMIF.CoulombDielectric(mol, confId=0)
    couldiele1 = rdMIF.CoulombDielectric(charges, pos)
    chargesMismatchingLen = charges[:-1]
    self.assertRaises(ValueError, rdMIF.CoulombDielectric, chargesMismatchingLen, pos)

    rdMIF.CalculateDescriptors(grd, couldiele)
    rdMIF.CalculateDescriptors(grd2, couldiele1)

    self.assertTrue(grd.CompareGrids(grd2))
    self.assertTrue(feq(couldiele(0.0, 0.0, 0.0, 1000), 0.0))
    self.assertTrue(couldiele(2.0, 0.0, 0.0, 1000) < 0)
    self.assertTrue(couldiele(-2.0, 0.0, 0.0, 1000) > 0)

    rdMIF.CalculateDescriptors(grd, rdMIF.CoulombDielectric(mol, absVal=True))

    for i in range(grd.GetSize()):
      self.assertTrue(grd.GetVal(i) <= 0.0)

    couldiele1 = rdMIF.CoulombDielectric(mol, probeCharge=-1.0)
    self.assertTrue(couldiele1(-2.0, 0.0, 0.0, 1000) < 0)
    self.assertTrue(couldiele1(2.0, 0.0, 0.0, 1000) > 0)

    couldiele2 = rdMIF.CoulombDielectric(mol, probeCharge=-.5)
    self.assertTrue(couldiele1(-2.0, 0.0, 0.0, 1000) < couldiele2(-2.0, 0.0, 0.0, 1000))

    couldiele3 = rdMIF.CoulombDielectric(mol, confId=0, probeCharge=1.0, absVal=False,
                                         chargeKey="_GasteigerCharge", softcoreParam=0.01,
                                         cutoff=1.0)
    self.assertTrue(couldiele3(0.0, 0.0, 0.0, 1000) > couldiele3(0.1, 0.0, 0.0, 1000))
    self.assertTrue(couldiele3(0.66, 0.0, 0.0, 1000) > couldiele3(0.68, 0.0, 0.0, 1000))
    self.assertTrue(couldiele3(0.70, 0.0, 0.0, 1000) > couldiele3(0.68, 0.0, 0.0, 1000))

    mol = AllChem.MolFromMolFile(
      os.path.join(RDConfig.RDBaseDir,
                   'Code/GraphMol/MolInteractionFields/Wrap/testData/glucose.mol'), removeHs=False)
    AllChem.ComputeGasteigerCharges(mol)

    couldiele4 = rdMIF.CoulombDielectric(mol, confId=0, probeCharge=1.0, absVal=False,
                                         chargeKey="_GasteigerCharge", softcoreParam=0.01,
                                         cutoff=1.0, epsilon=80.0, xi=4.0)
    couldiele5 = rdMIF.CoulombDielectric(mol, confId=0, probeCharge=1.0, absVal=False,
                                         chargeKey="_GasteigerCharge", softcoreParam=0.01,
                                         cutoff=1.0, epsilon=200.0, xi=4.0)
    couldiele6 = rdMIF.CoulombDielectric(mol, confId=0, probeCharge=1.0, absVal=False,
                                         chargeKey="_GasteigerCharge", softcoreParam=0.01,
                                         cutoff=1.0, epsilon=80.0, xi=10.0)

    self.assertTrue(couldiele5(-1.0, 0.0, 0.0, 1000) < couldiele4(-1.0, 0.0, 0.0, 1000))
    self.assertTrue(couldiele6(-1.0, 0.0, 0.0, 1000) < couldiele4(-1.0, 0.0, 0.0, 1000))

  def test5VdWaals(self):
    mol = AllChem.MolFromMolFile(
      os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/MolInteractionFields/Wrap/testData/HCN.mol'),
      removeHs=False)
    vdw = rdMIF.MMFFVdWaals(mol, confId=0, probeAtomType=6, scaling=False, cutoff=1.0)

    self.assertTrue(vdw(-5.0, 0, 0, 1000) < 0)
    self.assertTrue(vdw(-1.68, 0, 0, 1000) > vdw(-5.0, 0, 0, 1000))
    self.assertTrue(vdw(-5.0, 0, 0, 1000) < vdw(-10.0, 0, 0, 1000))

    mol2 = AllChem.MolFromMolFile(
      os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/MolInteractionFields/Wrap/testData/h2o.mol'),
      removeHs=False)
    vdw = rdMIF.MMFFVdWaals(mol2, scaling=False)
    vdw2 = rdMIF.MMFFVdWaals(mol2, scaling=True)

    self.assertTrue(abs(vdw2(-3.0, 0, 0, 1000) - vdw(-3.0, 0, 0, 1000)) > 0.0001)

    vdw3 = rdMIF.UFFVdWaals(mol, confId=0, probeAtomType="O_3", cutoff=1.0)

    self.assertTrue(vdw3(-5.0, 0, 0, 1000) < 0)
    self.assertTrue(vdw3(-1.68, 0, 0, 1000) > vdw3(-5.0, 0, 0, 1000))
    self.assertTrue(vdw3(-5.0, 0, 0, 1000) < vdw3(-10.0, 0, 0, 1000))

  def test6HBond(self):
    mol = AllChem.MolFromMolFile(
      os.path.join(RDConfig.RDBaseDir,
                   'Code/GraphMol/MolInteractionFields/Wrap/testData/ethane.mol'), removeHs=False)
    grd = rdMIF.ConstructGrid(mol, margin=5.0, spacing=0.5)
    for i in range(grd.GetSize()):
      grd.SetVal(i, 1.0)
    grd1 = copy.deepcopy(grd)

    hbonddes = rdMIF.HBond(mol)
    rdMIF.CalculateDescriptors(grd, hbonddes)
    self.assertTrue(not grd.CompareGrids(grd1))
    self.assertTrue(
      abs(int((grd.GetOccupancyVect() - grd1.GetOccupancyVect()).GetTotalVal())) == grd.GetSize())

    hbonddes = rdMIF.HBond(mol, probeAtomType='O')
    for i in range(grd.GetSize()):
      grd.SetVal(i, 1.0)
    rdMIF.CalculateDescriptors(grd, hbonddes)
    self.assertTrue(not grd.CompareGrids(grd1))
    self.assertTrue(
      abs(int((grd.GetOccupancyVect() - grd1.GetOccupancyVect()).GetTotalVal())) == grd.GetSize())

    mol = AllChem.MolFromMolFile(
      os.path.join(RDConfig.RDBaseDir,
                   'Code/GraphMol/MolInteractionFields/Wrap/testData/aceticacid.mol'),
      removeHs=False)
    grd = rdMIF.ConstructGrid(mol, margin=5.0, spacing=0.5)

    hbonddes = rdMIF.HBond(mol, probeAtomType="OH")
    rdMIF.CalculateDescriptors(grd, hbonddes)

    hbonddes1 = rdMIF.HBond(mol, probeAtomType="O", fixed=True)
    rdMIF.CalculateDescriptors(grd, hbonddes1)

    hbonddes2 = rdMIF.HBond(mol, probeAtomType="O", fixed=False)
    self.assertTrue(hbonddes1(4.0, 0.0, 1.0, 1000) > hbonddes2(4.0, 0.0, 1.0, 1000))

    hbonddes3 = rdMIF.HBond(mol, probeAtomType="NH")
    self.assertTrue(hbonddes(2.0, 2.0, 1.0, 1000) < hbonddes3(2.0, 2.0, 1.0, 1000))

    hbonddes4 = rdMIF.HBond(mol, probeAtomType="N")
    self.assertTrue(hbonddes1(3.0, 0.0, 0.0, 1000) < hbonddes4(3.0, 0.0, 0.0, 1000))

  def test7Hydrophilic(self):
    mol = AllChem.MolFromMolFile(
      os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/MolInteractionFields/Wrap/testData/h2o.mol'),
      removeHs=False)
    hydro = rdMIF.Hydrophilic(mol)
    hbondOH = rdMIF.HBond(mol, probeAtomType="OH")
    hbondO = rdMIF.HBond(mol, probeAtomType="O")

    pt = geom.Point3D(0.0, 0.0, 0.0)
    hyd = hydro(pt.x, pt.y, pt.z, 1000)
    hOH = hbondOH(pt.x, pt.y, pt.z, 1000)
    hO = hbondO(pt.x, pt.y, pt.z, 1000)
    self.assertTrue(feq(min(hOH, hO), hyd))

    pt = geom.Point3D(1.0, 1.5, 2.0)
    hyd = hydro(pt.x, pt.y, pt.z, 1000)
    hOH = hbondOH(pt.x, pt.y, pt.z, 1000)
    hO = hbondO(pt.x, pt.y, pt.z, 1000)
    self.assertTrue(feq(min(hOH, hO), hyd))

    pt = geom.Point3D(2.0, 1.5, -3.0)
    hyd = hydro(pt.x, pt.y, pt.z, 1000)
    hOH = hbondOH(pt.x, pt.y, pt.z, 1000)
    hO = hbondO(pt.x, pt.y, pt.z, 1000)
    self.assertTrue(feq(min(hOH, hO), hyd))

    pt = geom.Point3D(-2.5, 0.5, 3.0)
    hyd = hydro(pt.x, pt.y, pt.z, 1000)
    hOH = hbondOH(pt.x, pt.y, pt.z, 1000)
    hO = hbondO(pt.x, pt.y, pt.z, 1000)
    self.assertTrue(feq(min(hOH, hO), hyd))

    pt = geom.Point3D(10.0, 1.5, 1.0)
    hyd = hydro(pt.x, pt.y, pt.z, 1000)
    hOH = hbondOH(pt.x, pt.y, pt.z, 1000)
    hO = hbondO(pt.x, pt.y, pt.z, 1000)
    self.assertTrue(feq(min(hOH, hO), hyd))

    pt = geom.Point3D(6.0, -5.0, 0.0)
    hyd = hydro(pt.x, pt.y, pt.z, 1000)
    hOH = hbondOH(pt.x, pt.y, pt.z, 1000)
    hO = hbondO(pt.x, pt.y, pt.z, 1000)
    self.assertTrue(feq(min(hOH, hO), hyd))

    pt = geom.Point3D(-3.0, -3.0, 7.0)
    hyd = hydro(pt.x, pt.y, pt.z, 1000)
    hOH = hbondOH(pt.x, pt.y, pt.z, 1000)
    hO = hbondO(pt.x, pt.y, pt.z, 1000)
    self.assertTrue(feq(min(hOH, hO), hyd))

    pt = geom.Point3D(1.0, 0.0, 0.0)
    hyd = hydro(pt.x, pt.y, pt.z, 1000)
    hOH = hbondOH(pt.x, pt.y, pt.z, 1000)
    hO = hbondO(pt.x, pt.y, pt.z, 1000)
    self.assertTrue(feq(min(hOH, hO), hyd))

    pt = geom.Point3D(0.0, 2.0, 2.0)
    hyd = hydro(pt.x, pt.y, pt.z, 1000)
    hOH = hbondOH(pt.x, pt.y, pt.z, 1000)
    hO = hbondO(pt.x, pt.y, pt.z, 1000)
    self.assertTrue(feq(min(hOH, hO), hyd))

    pt = geom.Point3D(2.0, -2.0, 0.0)
    hyd = hydro(pt.x, pt.y, pt.z, 1000)
    hOH = hbondOH(pt.x, pt.y, pt.z, 1000)
    hO = hbondO(pt.x, pt.y, pt.z, 1000)
    self.assertTrue(feq(min(hOH, hO), hyd))

    pt = geom.Point3D(2.0, -2.0, -3.0)
    hyd = hydro(pt.x, pt.y, pt.z, 1000)
    hOH = hbondOH(pt.x, pt.y, pt.z, 1000)
    hO = hbondO(pt.x, pt.y, pt.z, 1000)
    self.assertTrue(feq(min(hOH, hO), hyd))


if __name__ == '__main__':
  print("Testing MIF wrapper")
  unittest.main()
