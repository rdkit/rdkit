# $Id$
#
#  Copyright (C) 2004-2006 Rational Discovery LLC
#
#     @@  All Rights Reserved  @@
#
from rdkit import RDConfig
import os,sys
import unittest
import math
from rdkit import Chem
from rdkit.Chem import rdMolAlign,rdDistGeom,ChemicalForceFields

def lstFeq(l1, l2, tol=1.e-4):
  if (len(list(l1)) != len(list(l2))):
    return 0
  for i in range(len(list(l1))):
    if not feq(l1[i], l2[i], tol):
      return 0
  return 1

def feq(v1,v2,tol2=1e-4):
  return abs(v1-v2)<=tol2

class TestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test1Basic(self):
        file1 = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                             'MolAlign', 'test_data', '1oir.mol')
        file2 = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                             'MolAlign', 'test_data', '1oir_conf.mol')

        mol1 = Chem.MolFromMolFile(file1)
        mol2 = Chem.MolFromMolFile(file2)

        rmsd = rdMolAlign.AlignMol(mol2, mol1)
        self.failUnless(feq(rmsd, 0.6578))

        file3 = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                             'MolAlign', 'test_data', '1oir_trans.mol')
        mol3 = Chem.MolFromMolFile(file3)
        conf2 = mol2.GetConformer()
        conf3 = mol3.GetConformer()

        for i in range(mol2.GetNumAtoms()):
            self.failUnless(lstFeq(conf2.GetAtomPosition(i), conf3.GetAtomPosition(i)))

        rmsd, trans = rdMolAlign.GetAlignmentTransform(mol2, mol1)
        self.failUnlessAlmostEqual(rmsd, 0.6579,4)

    def test2AtomMap(self) :
        atomMap = ((18,27), (13,23), (21,14), (24,7), (9,19), (16,30))
        file1 = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                             'MolAlign', 'test_data', '1oir.mol')
        file2 = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                             'MolAlign', 'test_data', '1oir_conf.mol')

        mol1 = Chem.MolFromMolFile(file1)
        mol2 = Chem.MolFromMolFile(file2)
        rmsd = rdMolAlign.AlignMol(mol2, mol1, 0, 0, atomMap)
        self.failUnlessAlmostEqual(rmsd, 0.8525,4)

    def test3Weights(self):
        atomMap = ((18,27), (13,23), (21,14), (24,7), (9,19), (16,30))
        file1 = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                             'MolAlign', 'test_data', '1oir.mol')
        file2 = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                             'MolAlign', 'test_data', '1oir_conf.mol')

        mol1 = Chem.MolFromMolFile(file1)
        mol2 = Chem.MolFromMolFile(file2)
        wts = (1.0, 1.0, 1.0, 1.0, 1.0, 2.0)
        rmsd = rdMolAlign.AlignMol(mol2, mol1, 0, 0, atomMap, wts)
        self.failUnlessAlmostEqual(rmsd, 0.9513,4)

    def test4AlignConfs(self):
      mol = Chem.MolFromSmiles('C1CC1CNc(n2)nc(C)cc2Nc(cc34)ccc3[nH]nc4')
      
      cids = rdDistGeom.EmbedMultipleConfs(mol,10,30,100)
      writer = Chem.SDWriter('mol_899.sdf')
    
      for cid in cids:
        print 'cid:',repr(cid)
        ff = ChemicalForceFields.UFFGetMoleculeForceField(mol, confId=cid)
        ff.Initialize()
        more = 1
        while more :
          more = ff.Minimize()
        # FIX: this should not be necessary but somehow more comes out to be 0
        # even with the structure still being crappy
        ff.Minimize() 
      aids = [12, 13, 14, 15, 16, 17, 18]
      rdMolAlign.AlignMolConformers(mol, aids)

      # now test that the atom location of these atom are consistent
      confs = mol.GetConformers()
      for aid in aids:
        mpos = 0
        for i,conf in enumerate(confs):
          if (i == 0):
            mpos = list(conf.GetAtomPosition(aid))
            continue
          else :
            pos = list(conf.GetAtomPosition(aid))
            
            self.failUnless(lstFeq(mpos, pos, .5))

    def test5O3A(self):
      sdf = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                         'MolAlign', 'test_data', 'ref_e2.sdf')
      # alignedSdf = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
      #                           'MolAlign', 'test_data', 'ref_e2_pyO3A.sdf')
      molS = Chem.SDMolSupplier(sdf, True, False)
      # molW = Chem.SDWriter(alignedSdf)
      refNum = 48
      refMol = molS[refNum]
      cumScore = 0.0
      cumMsd = 0.0
      refPyMP = ChemicalForceFields.MMFFGetMoleculeProperties(refMol)
      for prbMol in molS:
        prbPyMP = ChemicalForceFields.MMFFGetMoleculeProperties(prbMol)
        pyO3A = rdMolAlign.GetO3A(prbMol, refMol, prbPyMP, refPyMP)
        cumScore += pyO3A.Score()
        rmsd = pyO3A.Align()
        cumMsd += rmsd * rmsd
        # molW.write(prbMol)
      cumMsd /= len(molS)
      self.failUnlessAlmostEqual(cumScore,6772,0)
      self.failUnlessAlmostEqual(math.sqrt(cumMsd),.385,3)

    def test6O3A(self):
      " now test where the mmff parameters are generated on call "
      sdf = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                         'MolAlign', 'test_data', 'ref_e2.sdf')
      molS = Chem.SDMolSupplier(sdf, True, False)
      refNum = 48
      refMol = molS[refNum]
      cumScore = 0.0
      cumMsd = 0.0
      for prbMol in molS:
        pyO3A = rdMolAlign.GetO3A(prbMol, refMol)
        cumScore += pyO3A.Score()
        rmsd = pyO3A.Align()
        cumMsd += rmsd * rmsd
      cumMsd /= len(molS)
      self.failUnlessAlmostEqual(cumScore,6772,0)
      self.failUnlessAlmostEqual(math.sqrt(cumMsd),.385,3)
          
if __name__ == '__main__':
    print "Testing MolAlign Wrappers"
    unittest.main()
