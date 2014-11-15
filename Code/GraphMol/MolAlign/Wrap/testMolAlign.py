# $Id$
#
#  Copyright (C) 2004-2006 Rational Discovery LLC
#
#     @@  All Rights Reserved  @@
#
from __future__ import print_function
from rdkit import RDConfig
import os,sys,copy
import unittest
import math
from rdkit import Chem
from rdkit.Chem import rdMolAlign,rdMolTransforms,rdMolDescriptors,rdDistGeom,ChemicalForceFields

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
        self.assertTrue(feq(rmsd, 0.6578))

        file3 = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                             'MolAlign', 'test_data', '1oir_trans.mol')
        mol3 = Chem.MolFromMolFile(file3)
        conf2 = mol2.GetConformer()
        conf3 = mol3.GetConformer()

        for i in range(mol2.GetNumAtoms()):
            self.assertTrue(lstFeq(conf2.GetAtomPosition(i), conf3.GetAtomPosition(i)))

        rmsd, trans = rdMolAlign.GetAlignmentTransform(mol2, mol1)
        self.assertAlmostEqual(rmsd, 0.6579,4)

    def test2AtomMap(self) :
        atomMap = ((18,27), (13,23), (21,14), (24,7), (9,19), (16,30))
        file1 = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                             'MolAlign', 'test_data', '1oir.mol')
        file2 = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                             'MolAlign', 'test_data', '1oir_conf.mol')

        mol1 = Chem.MolFromMolFile(file1)
        mol2 = Chem.MolFromMolFile(file2)
        rmsd = rdMolAlign.AlignMol(mol2, mol1, 0, 0, atomMap)
        self.assertAlmostEqual(rmsd, 0.8525,4)

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
        self.assertAlmostEqual(rmsd, 0.9513,4)

    def test4AlignConfs(self):
      mol = Chem.MolFromSmiles('C1CC1CNc(n2)nc(C)cc2Nc(cc34)ccc3[nH]nc4')
      
      cids = rdDistGeom.EmbedMultipleConfs(mol,10,30,100)
      #writer = Chem.SDWriter('mol_899.sdf')
    
      for cid in cids:
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
            
            self.assertTrue(lstFeq(mpos, pos, .5))

      # now test that we can get a list of RMS values
      rmsvals = []
      rdMolAlign.AlignMolConformers(mol, aids, RMSlist=rmsvals)
      self.assertTrue((len(rmsvals)==mol.GetNumConformers()-1))

      # make sure something sensible happens if we provide a stupid
      # argument:
      rmsvals = 4
      self.assertRaises(AttributeError,rdMolAlign.AlignMolConformers,mol, atomIds=aids, RMSlist=rmsvals)

    def test5MMFFO3A(self):
      sdf = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                         'MolAlign', 'test_data', 'ref_e2.sdf')
      # alignedSdf = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
      #                           'MolAlign', 'test_data', 'ref_e2_pyMMFFO3A.sdf')
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
      self.assertAlmostEqual(cumScore,6942,0)
      self.assertAlmostEqual(math.sqrt(cumMsd),.345,3)

    def test6MMFFO3A(self):
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
      self.assertAlmostEqual(cumScore,6942,0)
      self.assertAlmostEqual(math.sqrt(cumMsd),.345,3)

    def test7MMFFO3A(self):
      " make sure we generate an error if parameters are missing (github issue 158) "

      m1 = Chem.MolFromSmiles('c1ccccc1Cl')
      rdDistGeom.EmbedMolecule(m1)
      m2 = Chem.MolFromSmiles('c1ccccc1B(O)O')
      rdDistGeom.EmbedMolecule(m1)

      self.assertRaises(ValueError,lambda :rdMolAlign.GetO3A(m1, m2))
      self.assertRaises(ValueError,lambda :rdMolAlign.GetO3A(m2, m1))

    def test8MMFFO3A(self):
      " test MMFFO3A with constraints "

      #we superimpose two identical coplanar 4-phenylpyridines:
      #1) the usual way
      #2) forcing the pyridine nitrogen to match with the para
      #   carbon of the phenyl ring
      m = Chem.MolFromSmiles('n1ccc(cc1)-c1ccccc1')
      m1 = Chem.AddHs(m)
      rdDistGeom.EmbedMolecule(m1)
      mp = ChemicalForceFields.MMFFGetMoleculeProperties(m1)
      ff = ChemicalForceFields.MMFFGetMoleculeForceField(m1, mp)
      ff.Minimize()
      sub1 = m1.GetSubstructMatch(Chem.MolFromSmarts('nccc-cccc'))
      nIdx = sub1[0]
      cIdx = sub1[-1]
      dihe = sub1[2:6]
      rdMolTransforms.SetDihedralDeg(m1.GetConformer(),
        dihe[0], dihe[1], dihe[2], dihe[3], 0)
      m2 = copy.copy(m1)
      rdMolAlign.RandomTransform(m2)
      m3 = copy.copy(m2)
      pyO3A = rdMolAlign.GetO3A(m2, m1)
      pyO3A.Align()
      d = m2.GetConformer().GetAtomPosition(cIdx). \
        Distance(m1.GetConformer().GetAtomPosition(cIdx))
      self.assertAlmostEqual(d, 0, 0)
      pyO3A = rdMolAlign.GetO3A(m3, m1, constraintMap = [[cIdx, nIdx]])
      pyO3A.Align()
      d = m3.GetConformer().GetAtomPosition(cIdx). \
        Distance(m1.GetConformer().GetAtomPosition(cIdx))
      self.assertAlmostEqual(d, 7, 0)
      #alignedSdf = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
      #                          'MolAlign', 'test_data',
      #                          '4-phenylpyridines_MMFFO3A.sdf')
      #sdW = Chem.SDWriter(alignedSdf)
      #sdW.write(m1)
      #sdW.write(m2)
      #sdW.write(m3)
      #sdW.close()

    def test9MMFFO3A(self):
      " test MMFFO3A with variable weight constraints followed by local-only optimization "

      sdf = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                         'MolAlign', 'test_data', 'ref_e2.sdf')
      # alignedSdf = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
      #                           'MolAlign', 'test_data', 'localonly.sdf')
      molS = Chem.SDMolSupplier(sdf, True, False)
      refNum = 23
      prbNum = 32
      refMol = molS[refNum]
      prbMol = molS[prbNum]
      refPyMP = ChemicalForceFields.MMFFGetMoleculeProperties(refMol)
      prbPyMP = ChemicalForceFields.MMFFGetMoleculeProperties(prbMol)
      refSIdx = refMol.GetSubstructMatch(Chem.MolFromSmarts('S'))[0]
      prbOIdx = prbMol.GetSubstructMatch(Chem.MolFromSmarts('O'))[0]
      # molW = Chem.SDWriter(alignedSdf)
      # molW.write(refMol)
      weights = [10.0, 100.0]
      distOS = [3.2, 0.3]
      for i in [0, 1]:
        pyO3A = rdMolAlign.GetO3A(prbMol, refMol,
          prbPyMP, refPyMP, constraintMap = [[prbOIdx, refSIdx]],
          constraintWeights = [weights[i]])
        pyO3A.Align()
        # molW.write(prbMol)
        pyO3A = rdMolAlign.GetO3A(prbMol, refMol,
          prbPyMP, refPyMP, options = 4)
        pyO3A.Align()
        # molW.write(prbMol)
        d = prbMol.GetConformer().GetAtomPosition(prbOIdx). \
          Distance(refMol.GetConformer().GetAtomPosition(refSIdx))
        self.assertAlmostEqual(d, distOS[i], 1)
      # molW.close()

    def test10CrippenO3A(self):
      sdf = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                         'MolAlign', 'test_data', 'ref_e2.sdf')
      alignedSdf = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                                'MolAlign', 'test_data', 'ref_e2_pyCrippenO3A.sdf')
      molS = Chem.SDMolSupplier(sdf, True, False)
      molW = Chem.SDWriter(alignedSdf)
      refNum = 48
      refMol = molS[refNum]
      cumScore = 0.0
      cumMsd = 0.0
      refList = rdMolDescriptors._CalcCrippenContribs(refMol, True)
      for prbMol in molS:
        prbList = rdMolDescriptors._CalcCrippenContribs(prbMol, True)
        pyO3A = rdMolAlign.GetCrippenO3A(prbMol, refMol, prbList, refList)
        cumScore += pyO3A.Score()
        rmsd = pyO3A.Align()
        cumMsd += rmsd * rmsd
        molW.write(prbMol)
      cumMsd /= len(molS)
      self.assertAlmostEqual(cumScore,4918,0)
      self.assertAlmostEqual(math.sqrt(cumMsd),.304,3)

    def test11CrippenO3A(self):
      " now test where the Crippen parameters are generated on call "
      sdf = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                         'MolAlign', 'test_data', 'ref_e2.sdf')
      molS = Chem.SDMolSupplier(sdf, True, False)
      refNum = 48
      refMol = molS[refNum]
      cumScore = 0.0
      cumMsd = 0.0
      for prbMol in molS:
        pyO3A = rdMolAlign.GetCrippenO3A(prbMol, refMol)
        cumScore += pyO3A.Score()
        rmsd = pyO3A.Trans()[0]
        cumMsd += rmsd * rmsd
      cumMsd /= len(molS)
      self.assertAlmostEqual(cumScore,4918,0)
      self.assertAlmostEqual(math.sqrt(cumMsd),.304,3)

    def test12CrippenO3A(self):
      " test CrippenO3A with constraints "

      #we superimpose two identical coplanar 4-phenylpyridines:
      #1) the usual way
      #2) forcing the pyridine nitrogen to match with the para
      #   carbon of the phenyl ring
      m = Chem.MolFromSmiles('n1ccc(cc1)-c1ccccc1')
      m1 = Chem.AddHs(m)
      rdDistGeom.EmbedMolecule(m1)
      mp = ChemicalForceFields.MMFFGetMoleculeProperties(m1)
      ff = ChemicalForceFields.MMFFGetMoleculeForceField(m1, mp)
      ff.Minimize()
      sub1 = m1.GetSubstructMatch(Chem.MolFromSmarts('nccc-cccc'))
      nIdx = sub1[0]
      cIdx = sub1[-1]
      dihe = sub1[2:6]
      rdMolTransforms.SetDihedralDeg(m1.GetConformer(),
        dihe[0], dihe[1], dihe[2], dihe[3], 0)
      m2 = copy.copy(m1)
      rdMolAlign.RandomTransform(m2)
      m3 = copy.copy(m2)
      pyO3A = rdMolAlign.GetCrippenO3A(m2, m1)
      pyO3A.Align()
      d = m2.GetConformer().GetAtomPosition(cIdx). \
        Distance(m1.GetConformer().GetAtomPosition(cIdx))
      self.assertAlmostEqual(d, 0, 0)
      pyO3A = rdMolAlign.GetCrippenO3A(m3, m1, constraintMap = [[cIdx, nIdx]])
      pyO3A.Align()
      d = m3.GetConformer().GetAtomPosition(cIdx). \
        Distance(m1.GetConformer().GetAtomPosition(cIdx))
      self.assertAlmostEqual(d, 7, 0)
      #alignedSdf = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
      #                          'MolAlign', 'test_data',
      #                          '4-phenylpyridines_CrippenO3A.sdf')
      #sdW = Chem.SDWriter(alignedSdf)
      #sdW.write(m1)
      #sdW.write(m2)
      #sdW.write(m3)
      #sdW.close()

    def test13CrippenO3A(self):
      " test CrippenO3A with variable weight constraints followed by local-only optimization "

      sdf = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                         'MolAlign', 'test_data', 'ref_e2.sdf')
      # alignedSdf = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
      #                           'MolAlign', 'test_data', 'localonly.sdf')
      molS = Chem.SDMolSupplier(sdf, True, False)
      refNum = 23
      prbNum = 32
      refMol = molS[refNum]
      prbMol = molS[prbNum]
      refPyMP = ChemicalForceFields.MMFFGetMoleculeProperties(refMol)
      prbPyMP = ChemicalForceFields.MMFFGetMoleculeProperties(prbMol)
      refSIdx = refMol.GetSubstructMatch(Chem.MolFromSmarts('S'))[0]
      prbOIdx = prbMol.GetSubstructMatch(Chem.MolFromSmarts('O'))[0]
      # molW = Chem.SDWriter(alignedSdf)
      # molW.write(refMol)
      weights = [0.1, 100.0]
      distOS = [2.7, 0.4]
      for i in [0, 1]:
        pyO3A = rdMolAlign.GetCrippenO3A(prbMol, refMol,
          constraintMap = [[prbOIdx, refSIdx]],
          constraintWeights = [weights[i]])
        pyO3A.Align()
        # molW.write(prbMol)
        pyO3A = rdMolAlign.GetCrippenO3A(prbMol, refMol, options = 4)
        pyO3A.Align()
        # molW.write(prbMol)
        d = prbMol.GetConformer().GetAtomPosition(prbOIdx). \
          Distance(refMol.GetConformer().GetAtomPosition(refSIdx))
        self.assertAlmostEqual(d, distOS[i], 1)
      # molW.close()

    def test14Github385(self):
      """ test github issue 385:
        O3A code generating incorrect results for multiconformer molecules
      """
      def _multiConfFromSmiles(smiles, nConfs=10, maxIters=500):
          """Adds hydrogens to molecule and optimises a chosen number of conformers.  Returns the optimised RDKit mol."""
          idea = Chem.MolFromSmiles(smiles)
          idea = Chem.AddHs(idea)
          confs = rdDistGeom.EmbedMultipleConfs(idea, nConfs)

          for conf in confs:
              opt = ChemicalForceFields.MMFFOptimizeMolecule(idea, confId=conf, maxIters=maxIters)
          return idea
      def _confsToAlignedMolsList(multiConfMol):
          """Input is a multiconformer RDKit mol.  Output is an aligned set of conformers as a list of RDKit mols."""
          rdMolAlign.AlignMolConformers(multiConfMol)
          ms = []
          cids = [x.GetId() for x in multiConfMol.GetConformers()]
          for cid in cids:
              newmol = Chem.MolToMolBlock(multiConfMol, confId=cid)
              newmol = Chem.MolFromMolBlock(newmol, removeHs=False)
              ms.append(newmol)
          return ms      
      reference = Chem.MolFromSmiles("c1ccccc1N2CCC(NS(=O)(=O)C(F)(F)F)CC2")
      reference = Chem.AddHs(reference)
      rdDistGeom.EmbedMolecule(reference)
      idea1 = _multiConfFromSmiles("c1ccccc1C2CCCCC2", 10)      
      
      idea1_mols = _confsToAlignedMolsList(idea1)
      cids = [x.GetId() for x in idea1.GetConformers()]

      refParams = ChemicalForceFields.MMFFGetMoleculeProperties(reference)
      prbParams = ChemicalForceFields.MMFFGetMoleculeProperties(idea1)

      for i in range(len(cids)):
          o3a1 = rdMolAlign.GetO3A(idea1_mols[i],reference,prbParams,refParams)
          score1 = o3a1.Score()

          o3a2 = rdMolAlign.GetO3A(idea1,reference,prbParams,refParams,prbCid=cids[i])
          score2 = o3a2.Score()
          self.assertAlmostEqual(score1,score2,3)


if __name__ == '__main__':
    print("Testing MolAlign Wrappers")
    unittest.main()
