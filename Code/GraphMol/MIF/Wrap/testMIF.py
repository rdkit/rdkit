from __future__ import print_function
import os,sys
import unittest
import copy
from rdkit import DataStructs, RDConfig
from rdkit.Geometry import rdGeometry as geom
from rdkit.Chem import rdMIF, AllChem
import math
import numpy as np

def feq(v1, v2, tol=1.0e-4):
    return abs(v1-v2) < tol

class testfunctor():
    def __call__ (self, x, y, z):
  	return 1.0
    def __call__(self, pt):
	return 1.0

class TestCase(unittest.TestCase):
    def SetUp(self):
        pass

    def test1ConstructGrid(self):
        mol = AllChem.MolFromMolFile(os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/MIF/Wrap/testData/HCl.mol'), removeHs=False)

        grd = rdMIF.ConstructGrid(mol, margin=5.0, spacing=0.5)

        bond = mol.GetConformer().GetAtomPosition(1)-mol.GetConformer().GetAtomPosition(0)
        self.failUnless(feq(grd.GetSpacing (), 0.5))
        self.failUnless(grd.GetNumX () == int((abs(bond.x) + 10.0) / 0.5 + 0.5))
        self.failUnless(grd.GetNumY () == int((abs(bond.y) + 10.0) / 0.5 + 0.5))
        self.failUnless(grd.GetNumZ () == int((abs(bond.z) + 10.0) / 0.5 + 0.5))
        
        
    def test2CubeFiles(self):
        mol = AllChem.MolFromMolFile(os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/MIF/Wrap/testData/HCl.mol'), removeHs=False)

        grd = geom.UniformRealValueGrid3D(5.0, 5.0, 5.0, 1.0, geom.Point3D(0.0, 0.0, 0.0))        
        for i in range(grd.GetSize()):
            grd.SetVal (i, float(i / 10.0))

        rdMIF.WriteToCubeFile(grd, mol, os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/MIF/Wrap/testData/test3.cube'))

        grd2 = geom.UniformRealValueGrid3D(1.0,1.0,1.0,1.0)
        mol2 = rdMIF.ReadFromCubeFile (grd2, os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/MIF/Wrap/testData/test3.cube'))

        self.failUnless(grd.GetSize() == grd2.GetSize())

        for i in range(grd.GetSize()):
            self.failUnless(feq(grd2.GetVal(i), float(i/10.0)))

        self.failUnless(grd.CompareGrids(grd2))
        self.failUnless(mol.GetNumAtoms() == mol2.GetNumAtoms())

        for i in range(mol.GetNumAtoms()):
            self.failUnless(mol.GetAtomWithIdx(i).GetAtomicNum() == mol2.GetAtomWithIdx(i).GetAtomicNum())
            self.failUnless(feq(mol.GetConformer().GetAtomPosition(i).x, mol2.GetConformer().GetAtomPosition(i).x))
            self.failUnless(feq(mol.GetConformer().GetAtomPosition(i).y, mol2.GetConformer().GetAtomPosition(i).y))
            self.failUnless(feq(mol.GetConformer().GetAtomPosition(i).z, mol2.GetConformer().GetAtomPosition(i).z))


    def test3Coulomb(self):
        mol = AllChem.MolFromMolFile(os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/MIF/Wrap/testData/HCl.mol'), removeHs=False)
        AllChem.ComputeGasteigerCharges(mol)

        charges = []
        pos = []
        conf = mol.GetConformer(0)

        for i in range (mol.GetNumAtoms()):
            charges.append(float(mol.GetAtomWithIdx(i).GetProp("_GasteigerCharge")))
            pos.append(conf.GetAtomPosition(i))

        grd  = rdMIF.ConstructGrid(mol)
        grd2 = copy.deepcopy(grd)

        coul = rdMIF.Coulomb(mol)

        rdMIF.CalculateDescriptors(grd, coul)
        coul1 = rdMIF.Coulomb_(charges, pos)
        pt = geom.Point3D(3.0, 3.0, 3.0)
        rdMIF.CalculateDescriptors(grd2, coul1)
        self.failUnless(grd.CompareParams(grd2))
        self.failUnless(grd.CompareVectors(grd2))
        self.failUnless(grd.CompareGrids(grd2))
        self.failUnless(feq(coul(geom.Point3D(0.0, 0.0, 0.0)), 0.0))
        self.failUnless(coul(geom.Point3D(2.0, 0.0, 0.0)) < 0)
        self.failUnless(coul(geom.Point3D(-2.0, 0.0, 0.0)) > 0)

        rdMIF.CalculateDescriptors(grd, rdMIF.Coulomb(mol, absVal=True))

        for i in range(grd.GetSize()):
            self.failUnless(grd.GetVal(i) <= 0.0)

        coul1 = rdMIF.Coulomb(mol, probeCharge=-1.0)
        self.failUnless(coul1(geom.Point3D(-2.0, 0.0, 0.0)) < 0)
        self.failUnless(coul1(geom.Point3D(2.0, 0.0, 0.0)) > 0)

        coul2 = rdMIF.Coulomb(mol, probeCharge= -.5)
        self.failUnless(coul1(geom.Point3D(-2.0, 0.0, 0.0)) < coul2(geom.Point3D(-2.0, 0.0, 0.0)))

        coul3 = rdMIF.Coulomb(mol, confId=0, probeCharge=1.0, absVal=False, chargeKey="_GasteigerCharge", softcoreParam=0.01, cutoffDist=1.0)
        self.failUnless(coul3 (geom.Point3D (0.0, 0.0, 0.0)) > coul3 (geom.Point3D (0.1, 0.0, 0.0)))
        self.failUnless(coul3 (geom.Point3D (0.66, 0.0, 0.0)) > coul3 (geom.Point3D (0.68, 0.0, 0.0)))
        self.failUnless(coul3 (geom.Point3D (0.70, 0.0, 0.0)) > coul3 (geom.Point3D (0.68, 0.0, 0.0)))
        

    def test4CoulombDielectric(self):
        mol = AllChem.MolFromMolFile(os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/MIF/Wrap/testData/HCl.mol'), removeHs=False)
        AllChem.ComputeGasteigerCharges(mol)

        charges = []
        pos = []
        conf = mol.GetConformer(0)

        for i in range (mol.GetNumAtoms()):
            charges.append(float(mol.GetAtomWithIdx(i).GetProp("_GasteigerCharge")))
            pos.append(conf.GetAtomPosition(i))

        grd  = rdMIF.ConstructGrid(mol, confId=0)
        grd2 = rdMIF.ConstructGrid(mol, confId=0)

        couldiele = rdMIF.CoulombDielectric(mol, confId=0)
        couldiele1 =  rdMIF.CoulombDielectric_(charges, pos)
        pt = geom.Point3D(3.0, 3.0, 3.0)

        rdMIF.CalculateDescriptors(grd, couldiele)
        rdMIF.CalculateDescriptors(grd2, rdMIF.CoulombDielectric_(charges, pos))

        self.failUnless(grd.CompareGrids(grd2))
        self.failUnless(feq(couldiele(geom.Point3D(0.0, 0.0, 0.0)), 0.0))
        self.failUnless(couldiele(geom.Point3D(2.0, 0.0, 0.0)) < 0)
        self.failUnless(couldiele(geom.Point3D(-2.0, 0.0, 0.0)) > 0)

        rdMIF.CalculateDescriptors(grd, rdMIF.CoulombDielectric(mol, absVal=True))

        for i in range(grd.GetSize()):
            self.failUnless(grd.GetVal(i) <= 0.0)

        couldiele1 = rdMIF.CoulombDielectric(mol, probeCharge=-1.0)
        self.failUnless(couldiele1(geom.Point3D(-2.0, 0.0, 0.0)) < 0)
        self.failUnless(couldiele1(geom.Point3D(2.0, 0.0, 0.0)) > 0)

        couldiele2 = rdMIF.CoulombDielectric(mol, probeCharge= -.5)
        self.failUnless(couldiele1(geom.Point3D(-2.0, 0.0, 0.0)) < couldiele2(geom.Point3D(-2.0, 0.0, 0.0)))

        couldiele3 = rdMIF.CoulombDielectric(mol, confId=0, probeCharge=1.0, absVal=False, chargeKey="_GasteigerCharge", softcoreParam=0.01, cutoffDist=1.0)
        self.failUnless(couldiele3 (geom.Point3D (0.0, 0.0, 0.0)) > couldiele3 (geom.Point3D (0.1, 0.0, 0.0)))
        self.failUnless(couldiele3 (geom.Point3D (0.66, 0.0, 0.0)) > couldiele3 (geom.Point3D (0.68, 0.0, 0.0)))
        self.failUnless(couldiele3 (geom.Point3D (0.70, 0.0, 0.0)) > couldiele3 (geom.Point3D (0.68, 0.0, 0.0)))


        mol = AllChem.MolFromMolFile(os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/MIF/Wrap/testData/glucose.mol'), removeHs=False)
        AllChem.ComputeGasteigerCharges(mol);

        couldiele4 = rdMIF.CoulombDielectric(mol, confId=0, probeCharge=1.0, absVal=False, chargeKey="_GasteigerCharge", softcoreParam=0.01, cutoffDist=1.0, epsilon=80.0, xi=4.0)
        couldiele5 = rdMIF.CoulombDielectric(mol, confId=0, probeCharge=1.0, absVal=False, chargeKey="_GasteigerCharge", softcoreParam=0.01, cutoffDist=1.0, epsilon=200.0, xi=4.0)
        couldiele6 = rdMIF.CoulombDielectric(mol, confId=0, probeCharge=1.0, absVal=False, chargeKey="_GasteigerCharge", softcoreParam=0.01, cutoffDist=1.0, epsilon=80.0, xi=10.0)
        
        self.failUnless(couldiele5(geom.Point3D(-1.0, 0.0, 0.0)) < couldiele4(geom.Point3D(-1.0, 0.0, 0.0)))
        self.failUnless(couldiele6(geom.Point3D(-1.0, 0.0, 0.0)) < couldiele4(geom.Point3D(-1.0, 0.0, 0.0)))

    def test5VdWaals(self):
        mol = AllChem.MolFromMolFile(os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/MIF/Wrap/testData/HCN.mol'), removeHs=False)
        vdw = rdMIF.ConstructVdWaalsMMFF(mol, confId=0, probeType=6, scaling=False, cutoffDist=1.0)
        
        self.failUnless(vdw(geom.Point3D(-5.0, 0, 0)) < 0)
        self.failUnless(vdw(geom.Point3D(-1.68, 0, 0)) > vdw (geom.Point3D (-5.0, 0, 0)))
        self.failUnless(vdw(geom.Point3D(-5.0, 0, 0)) < vdw (geom.Point3D (-10.0, 0, 0)))
        
        mol2 = AllChem.MolFromMolFile(os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/MIF/Wrap/testData/h2o.mol'), removeHs=False)
        vdw  = rdMIF.ConstructVdWaalsMMFF(mol2, scaling=False)
        vdw2 = rdMIF.ConstructVdWaalsMMFF(mol2, scaling=True)
        
        self.failUnless(abs(vdw2(geom.Point3D(-3.0, 0, 0)) - vdw(geom.Point3D(-3.0, 0, 0))) > 0.0001)

        vdw3 = rdMIF.ConstructVdWaalsUFF(mol, confId=0, probeType="O_3", cutoffDist=1.0)

        self.failUnless(vdw3(geom.Point3D (-5.0, 0, 0)) < 0)
        self.failUnless(vdw3(geom.Point3D(-1.68, 0, 0)) > vdw3(geom.Point3D(-5.0, 0, 0)))
        self.failUnless(vdw3(geom.Point3D(-5.0, 0, 0)) < vdw3(geom.Point3D(-10.0, 0, 0)))

#   std::string names[] = { "acetone", "aceticacid", "phenol", "phenolate",
#       "serine", "threonine", "ethanol", "diethylether", "h2o", "ammonia",
#       "ethylamine", "imine", "acetonitrile", "histidine", "phenylamine",
#       "methylammonium", "fluoromethane", "chloromethane", "bromomethane",
#       "glycine", "glyphe", "glysergly", "glythrgly", "glucose" };
#   for ( unsigned int i = 0; i < 24; i++ ){
#     mol = *MolFileToMol (path + names[i] + ".mol", true, false);
#     vdw = constructVdWaalsMMFF(mol);
#     self.failUnless(vdw(geom.Point3D(0.0,0.0,0.0)), "VdWMMFF: crashed with " + names[i]);
#     vdw = constructVdWaalsUFF(mol);
#     self.failUnless(vdw(geom.Point3D(0.0,0.0,0.0)), "VdWUFF: crashed with " + names[i]);

    def test6HBond(self):
        mol = AllChem.MolFromMolFile(os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/MIF/Wrap/testData/ethane.mol'), removeHs=False)
        grd = rdMIF.ConstructGrid (mol, margin=5.0, spacing=0.5)
        for i in range(grd.GetSize()):
            grd.SetVal (i, 1.0)
        grd1=copy.deepcopy(grd)

        hbonddes = rdMIF.HBond(mol)
        rdMIF.CalculateDescriptors(grd, hbonddes)
        self.failUnless( not grd.CompareGrids(grd1))
        self.failUnless(abs(int((grd.GetOccupancyVect() - grd1.GetOccupancyVect()).GetTotalVal())) == grd.GetSize())

        hbonddes = rdMIF.HBond(mol, probeType='O')
        for i in range(grd.GetSize()):
            grd.SetVal (i, 1.0)
        rdMIF.CalculateDescriptors(grd, hbonddes)
        self.failUnless( not grd.CompareGrids(grd1))
        self.failUnless(abs(int((grd.GetOccupancyVect() - grd1.GetOccupancyVect()).GetTotalVal())) == grd.GetSize())


        mol = AllChem.MolFromMolFile(os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/MIF/Wrap/testData/aceticacid.mol'), removeHs=False)
        grd = rdMIF.ConstructGrid (mol, margin=5.0, spacing=0.5)

        hbonddes = rdMIF.HBond(mol, probeType="OH")
        rdMIF.CalculateDescriptors(grd, hbonddes)

        hbonddes1 = rdMIF.HBond(mol, probeType="O", fixed=True)
        rdMIF.CalculateDescriptors(grd, hbonddes1)

        hbonddes2 = rdMIF.HBond(mol, probeType="O", fixed=False)
        self.failUnless(hbonddes1(geom.Point3D(4.0, 0.0, 1.0)) > hbonddes2(geom.Point3D(4.0, 0.0, 1.0)))

        hbonddes3 = rdMIF.HBond(mol, probeType="NH")
        self.failUnless(hbonddes(geom.Point3D(2.0, 2.0, 1.0)) < hbonddes3(geom.Point3D(2.0, 2.0, 1.0)))

        hbonddes4 = rdMIF.HBond(mol, probeType="N")
        self.failUnless(hbonddes1(geom.Point3D(3.0, 0.0,0.0)) < hbonddes4(geom.Point3D(3.0,0.0,0.0)))
        
#   mol = *MolFileToMol (path + "acetone.mol", true, false);     //Acetone
#   grd = *constructGrid (mol, 0, 5.0, 1);
#   hbonddes = HBond (mol, 0, "OH", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 1, "");
#   calculateDescriptors<HBond> (grd, hbonddes);
#   hbonddes = HBond (mol, 0, "O", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 0, "");
#   calculateDescriptors<HBond> (grd, hbonddes);

#   mol = *MolFileToMol (path + "diethylether.mol", true, false);                //Et2O
#   grd = *constructGrid (mol, 0, 5.0, 1);
#   hbonddes = HBond (mol, 0, "OH", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 1, "");
#   calculateDescriptors<HBond> (grd, hbonddes);
#   hbonddes = HBond (mol, 0, "O", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 0, "");
#   calculateDescriptors<HBond> (grd, hbonddes);

#   mol = *MolFileToMol (path + "h2o.mol", true, false);         //H2O
#   grd = *constructGrid (mol, 0, 5.0, 1);
#   hbonddes = HBond (mol, 0, "OH", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 1, "");
#   calculateDescriptors<HBond> (grd, hbonddes);
#   hbonddes = HBond (mol, 0, "O", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 2, "");
#   calculateDescriptors<HBond> (grd, hbonddes);

#   mol = *MolFileToMol (path + "ammonia.mol", true, false);             //ammonia
#   grd = *constructGrid (mol, 0, 5.0, 1);
#   hbonddes = HBond (mol, 0, "OH", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 1, "");
#   calculateDescriptors<HBond> (grd, hbonddes);
#   hbonddes = HBond (mol, 0, "O", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 3, "");
#   calculateDescriptors<HBond> (grd, hbonddes);


#   mol = *MolFileToMol (path + "imine.mol", true, false);               //imine
#   grd = *constructGrid (mol, 0, 5.0, 1);
#   hbonddes = HBond (mol, 0, "OH", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 1, "");
#   calculateDescriptors<HBond> (grd, hbonddes);
#   hbonddes = HBond (mol, 0, "O", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 0, "");
#   calculateDescriptors<HBond> (grd, hbonddes);

#   mol = *MolFileToMol (path + "methylammonium.mol", true, false);//methylammonium
#   grd = *constructGrid (mol, 0, 5.0, 1);
#   hbonddes = HBond (mol, 0, "OH", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 0, "");
#   calculateDescriptors<HBond> (grd, hbonddes);
#   hbonddes = HBond (mol, 0, "O", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 3, "");
#   calculateDescriptors<HBond> (grd, hbonddes);

#   mol = *MolFileToMol (path + "chloromethane.mol", true, false);       //Chloromethane
#   grd = *constructGrid (mol, 0, 5.0, 1);
#   hbonddes = HBond (mol, 0, "OH", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 1, "");
#   calculateDescriptors<HBond> (grd, hbonddes);
#   hbonddes = HBond (mol, 0, "O", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 0, "");
#   calculateDescriptors<HBond> (grd, hbonddes);

#   mol = *MolFileToMol (path + "phosphonate.mol", true, false);         //Phosphonate
#   grd = *constructGrid (mol, 0, 5.0, 1);
#   hbonddes = HBond (mol, 0, "OH", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 3, "");
#   calculateDescriptors<HBond> (grd, hbonddes);
#   hbonddes = HBond (mol, 0, "O", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 0, "");
#   calculateDescriptors<HBond> (grd, hbonddes);

#   mol = *MolFileToMol (path + "phosphatediester.mol", true, false);//Phosphatediester
#   grd = *constructGrid (mol, 0, 5.0, 1);
#   hbonddes = HBond (mol, 0, "OH", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 4, "");
#   calculateDescriptors<HBond> (grd, hbonddes);
#   hbonddes = HBond (mol, 0, "O", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 0, "");
#   calculateDescriptors<HBond> (grd, hbonddes);

#   mol = *MolFileToMol (path + "hydrogenphosphatediester.mol", true, false);//Hydrogenphosphatediester
#   grd = *constructGrid (mol, 0, 5.0, 1);
#   hbonddes = HBond (mol, 0, "OH", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 4, "");
#   calculateDescriptors<HBond> (grd, hbonddes);
#   hbonddes = HBond (mol, 0, "O", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 1, "");
#   calculateDescriptors<HBond> (grd, hbonddes);

#   mol = *MolFileToMol (path + "mustardgas.mol", true, false);          //mustard gas
#   grd = *constructGrid (mol, 0, 5.0, 1);
#   hbonddes = HBond (mol, 0, "OH", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 2, "");
#   calculateDescriptors<HBond> (grd, hbonddes);
#   hbonddes = HBond (mol, 0, "O", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 0, "");
#   calculateDescriptors<HBond> (grd, hbonddes);

#   mol = *MolFileToMol (path + "alicin.mol", true, false);              //Alicin
#   grd = *constructGrid (mol, 0, 5.0, 1);
#   hbonddes = HBond (mol, 0, "OH", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 1, "");
#   calculateDescriptors<HBond> (grd, hbonddes);
#   hbonddes = HBond (mol, 0, "O", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 0, "");
#   calculateDescriptors<HBond> (grd, hbonddes);

#   mol = *MolFileToMol (path + "sulfanilamide.mol", true, false);       //Sulfanilamide
#   grd = *constructGrid (mol, 0, 5.0, 1);
#   hbonddes = HBond (mol, 0, "OH", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 3, "");
#   calculateDescriptors<HBond> (grd, hbonddes);
#   hbonddes = HBond (mol, 0, "O", true, 0.001);
#   self.failUnless(hbonddes.getNumInteractions () == 4, "");
#   calculateDescriptors<HBond> (grd, hbonddes);
# }

    def test7Hydrophilic(self):
       mol = AllChem.MolFromMolFile(os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/MIF/Wrap/testData/h2o.mol'), removeHs=False)
       hydro = rdMIF.Hydrophilic(mol)
       hbondOH = rdMIF.HBond(mol, probeType="OH")
       hbondO = rdMIF.HBond(mol, probeType="O")

       pt = geom.Point3D(0.0,0.0,0.0)
       hyd = hydro(pt)
       hOH = hbondOH(pt)
       hO  = hbondO(pt)
       self.failUnless(feq(min(hOH, hO), hyd))

       pt = geom.Point3D(1.0, 1.5, 2.0)
       hyd = hydro(pt)
       hOH = hbondOH(pt)
       hO  = hbondO(pt)
       self.failUnless(feq(min(hOH, hO), hyd))
       
       pt = geom.Point3D(2.0, 1.5, -3.0)
       hyd = hydro(pt)
       hOH = hbondOH(pt)
       hO  = hbondO(pt)
       self.failUnless(feq(min(hOH, hO), hyd))
       
       pt = geom.Point3D(-2.5, 0.5, 3.0)
       hyd = hydro(pt)
       hOH = hbondOH(pt)
       hO  = hbondO(pt)
       self.failUnless(feq(min(hOH, hO), hyd))
       
       pt = geom.Point3D(10.0, 1.5, 1.0)
       hyd = hydro(pt)
       hOH = hbondOH(pt)
       hO  = hbondO(pt)
       self.failUnless(feq(min(hOH, hO), hyd))
       
       pt = geom.Point3D(6.0, -5.0, 0.0)
       hyd = hydro(pt)
       hOH = hbondOH(pt)
       hO  = hbondO(pt)
       self.failUnless(feq(min(hOH, hO), hyd))

       pt = geom.Point3D(-3.0, -3.0, 7.0)
       hyd = hydro(pt)
       hOH = hbondOH(pt)
       hO  = hbondO(pt)
       self.failUnless(feq(min(hOH, hO), hyd))
       
       pt = geom.Point3D(1.0, 0.0, 0.0)
       hyd = hydro(pt)
       hOH = hbondOH(pt)
       hO  = hbondO(pt)
       self.failUnless(feq(min(hOH, hO), hyd))

       pt = geom.Point3D(0.0, 2.0, 2.0)
       hyd = hydro(pt)
       hOH = hbondOH(pt)
       hO  = hbondO(pt)
       self.failUnless(feq(min(hOH, hO), hyd))

       pt = geom.Point3D(2.0, -2.0, 0.0)
       hyd = hydro(pt)
       hOH = hbondOH(pt)
       hO  = hbondO(pt)
       self.failUnless(feq(min(hOH, hO), hyd))

       pt = geom.Point3D(2.0, -2.0, -3.0)
       hyd = hydro(pt)
       hOH = hbondOH(pt)
       hO  = hbondO(pt)
       self.failUnless(feq(min(hOH, hO), hyd))


if __name__=='__main__':
    print("Testing MIF wrapper")
    unittest.main()
