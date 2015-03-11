# $Id$
#
#  Copyright (C) 2004-2008  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from __future__ import print_function
from rdkit import RDConfig
import unittest,sys,os
import io
from rdkit.six import PY3
from rdkit.six.moves import cPickle
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures,rdDistGeom
import EmbedLib
import gzip
from rdkit import DistanceGeometry as DG
from rdkit import Geometry
import Pharmacophore
import numpy

def feq(n1,n2,tol=1e-5):
  return abs(n1-n2)<=tol

class TestCase(unittest.TestCase):
  def setUp(self):
    self.dataDir = os.path.join(RDConfig.RDCodeDir,'Chem/Pharm3D/test_data')
    self.fdefBlock = \
                   """DefineFeature HAcceptor1 [N,O;H0]
                      Family HBondAcceptor
                      Weights 1.0
                   EndFeature
                   DefineFeature HDonor1 [N,O;!H0]
                      Family HBondDonor
                      Weights 1.0
                   EndFeature
                   DefineFeature Aromatic1 c1ccccc1
                      Family Aromatic
                      Weights 1.,1.,1.,1.,1.,1.
                   EndFeature\n"""

    self.featFactory = ChemicalFeatures.BuildFeatureFactoryFromString(self.fdefBlock)
    self.feats = [ChemicalFeatures.FreeChemicalFeature('HBondAcceptor', 'HAcceptor1',
                                           Geometry.Point3D(0.0, 0.0, 0.0)),
                  ChemicalFeatures.FreeChemicalFeature('HBondDonor', 'HDonor1',
                                           Geometry.Point3D(2.65, 0.0, 0.0)),
                  ChemicalFeatures.FreeChemicalFeature('Aromatic', 'Aromatic1',
                                           Geometry.Point3D(5.12, 0.908, 0.0)),
                  ]
    self.pcophore=Pharmacophore.Pharmacophore(self.feats)
    self.pcophore.setLowerBound(0,1, 2.0)
    self.pcophore.setUpperBound(0,1, 3.3)
    self.pcophore.setLowerBound(0,2, 5.0)
    self.pcophore.setUpperBound(0,2, 5.4)
    self.pcophore.setLowerBound(1,2, 2.6)
    self.pcophore.setUpperBound(1,2, 3.0)

  def _matchMol(self,tpl,pcophore,featFactory,downSample):
    name,molPkl,boundsMat = tpl
    mol = Chem.Mol(molPkl)
    matched,matches = EmbedLib.MatchPharmacophoreToMol(mol,featFactory,pcophore)
    if matched:
      r = EmbedLib.MatchPharmacophore(matches,boundsMat,pcophore,
                                      useDownsampling=downSample)
      if r[0]:
        return 0
      else:
        return 1
    else:
      return 0
    
  def test1SearchFullMat(self):
    inF = gzip.open(os.path.join(self.dataDir,'cdk2-syn-clip100.pkl.gz'),'rb')
    #outF = gzip.open(os.path.join(self.dataDir,'cdk2-syn-clip100.pkl.new.gz'),'wb+')
    nDone = 0
    nHits = 0
    while 1:
      try:
        tpl = cPickle.load(inF, encoding='latin1')
        if PY3:
          tpl = tpl[0], tpl[1].encode('latin1'), tpl[2]
        #tpl=tpl[0],tpl[1],numpy.array(tpl[2])
        #cPickle.dump(tpl,outF)
      except:
        break
      if self._matchMol(tpl,self.pcophore,self.featFactory,0):
        nHits+=1
      nDone += 1
    self.assertEqual(nDone,100)
    #print 'nHits:',nHits
    self.assertEqual(nHits,47)
    
  def test2SearchDownsample(self):
    inF = gzip.open(os.path.join(self.dataDir,'cdk2-syn-clip100.pkl.gz'),'rb')
    nDone = 0
    nHits = 0
    hits = []
    while 1:
      try:
        tpl = cPickle.load(inF, encoding='latin1')
        if PY3:
          tpl = tpl[0], tpl[1].encode('latin1'), tpl[2]
      except:
        break
      if self._matchMol(tpl,self.pcophore, self.featFactory,1):
        nHits+=1
      nDone += 1
    self.assertEqual(nDone,100)
    #print 'nHits:',nHits
    self.assertEqual(nHits,47)
    
  def test3Embed(self):
    testResults={
      'mol_197':(218.80,35.75,110.33,11.58,109.66,11.09,90.35,2.95,0.00),
      'mol_223':(259.19,6.27,134.13,1.12,134.06,1.12,85.74,0.61,0.00),
      'mol_269':(204.51,7.89,103.89,1.20,102.66,1.20,88.07,1.21,6.00),
      }
    inF = gzip.open(os.path.join(self.dataDir,'cdk2-syn-clip100.pkl.gz'),'rb')
    nDone = 0
    nHits = 0
    while 1:
      try:
        name,molPkl,boundsMat = cPickle.load(inF, encoding='latin1')
        if PY3:
          molPkl = bytes(molPkl, encoding='latin1')
      except:
        break

      nDone += 1

      mol = Chem.Mol(molPkl)
      nboundsMat = rdDistGeom.GetMoleculeBoundsMatrix(mol)
      DG.DoTriangleSmoothing(nboundsMat)
      matched,matches = EmbedLib.MatchPharmacophoreToMol(mol,self.featFactory,
                                                         self.pcophore)
      if matched:
        failed,bm,match,stats = EmbedLib.MatchPharmacophore(matches,nboundsMat,
                                                            self.pcophore,
                                                            useDownsampling=1)
        if not failed:
          nHits += 1
          
          if name in testResults:
            stats = EmbedLib.EmbedOne(mol,name,match,self.pcophore,count=10,
                                      silent=1,randomSeed=23)
            tgt = testResults[name]
            self.assertEqual(len(tgt),len(stats))
            print(name)
            print(','.join(['%.2f'%x for x in stats]))
            # we'll use different tolerances for the different values:
            self.assertTrue(feq(tgt[0],stats[0],5.0),(tgt[0],stats[0]))
            for i in range(2,len(tgt)):
              self.assertTrue(feq(tgt[i],stats[i],5.0),(tgt[i],stats[i]))
        
    self.assertEqual(nDone,100)
    #print 'nHits:',nHits
    self.assertEqual(nHits,50)
    
  def test4Search(self):
    featFactory = ChemicalFeatures.BuildFeatureFactory(os.path.join(self.dataDir,
                                                        'BaseFeatures.fdef'))

    activeFeats = [ChemicalFeatures.FreeChemicalFeature('Acceptor',
                                            Geometry.Point3D(0.0, 0.0, 0.0)),
                   ChemicalFeatures.FreeChemicalFeature('Donor',
                                            Geometry.Point3D(0.0, 0.0, 0.0)),
                   ChemicalFeatures.FreeChemicalFeature('Aromatic',
                                            Geometry.Point3D(0.0, 0.0, 0.0))]
    pcophore= Pharmacophore.Pharmacophore(activeFeats)
    pcophore.setLowerBound(0,1,2.251)
    pcophore.setUpperBound(0,1,2.451)
    pcophore.setUpperBound2D(0,1,3)

    pcophore.setLowerBound(0,2,4.970)
    pcophore.setUpperBound(0,2,5.170)
    pcophore.setUpperBound2D(0,2,6)

    pcophore.setLowerBound(1,2,2.681)
    pcophore.setUpperBound(1,2,2.881)
    pcophore.setUpperBound2D(1,2,6)

    inF = gzip.open(os.path.join(self.dataDir,'cdk2-syn-clip100.pkl.gz'),'rb')
    nDone = 0
    nMatches = 0
    nHits = 0

    while 1:
      try:
        name,molPkl,boundsMat = cPickle.load(inF, encoding='latin1')
        if PY3:
          molPkl = bytes(molPkl, encoding='latin1')
      except:
        break

      nDone += 1

      mol = Chem.Mol(molPkl)
      boundsMat = rdDistGeom.GetMoleculeBoundsMatrix(mol)
      DG.DoTriangleSmoothing(boundsMat)
    
      canMatch,matches = EmbedLib.MatchPharmacophoreToMol(mol,featFactory,
                                                          pcophore)
      if canMatch:
        nMatches+=1
        r = EmbedLib.MatchPharmacophore(matches,boundsMat,pcophore,
                                        useDownsampling=True,use2DLimits=True,
                                        mol=mol)
        failed,bm,match,details = r
        if not failed:
          nHits+=1

    self.assertEqual(nDone,100)
    self.assertEqual(nMatches,93)
    #print 'nhits:',nHits
    self.assertEqual(nHits,67)

  def testIssue268(self):
    from rdkit import RDLogger
    #RDLogger.EnableLog('rdApp.debug')
    featFactory = ChemicalFeatures.BuildFeatureFactory(os.path.join(self.dataDir,
                                                        'Issue268.fdef'))
    m1 = Chem.MolFromMolFile(os.path.join(self.dataDir,
                                          'Issue268_Mol1.mol'))
    m2 = Chem.MolFromMolFile(os.path.join(self.dataDir,
                                          'Issue268_Mol2.mol'))
    with open(os.path.join(self.dataDir,
                           'Issue268_Pcop.pkl'),'r') as inTF:
      buf = inTF.read().replace('\r\n', '\n').encode('utf-8')
      inTF.close()
    with io.BytesIO(buf) as inF:
      pcop = cPickle.load(inF, encoding='latin1')
    #pcop._boundsMat=numpy.array(pcop._boundsMat)
    #pcop._boundsMat2D=numpy.array(pcop._boundsMat2D)
    #cPickle.dump(pcop,file(os.path.join(self.dataDir,
    #                                    'Issue268_Pcop.new.pkl'),'wb+'))
    match,mList1 = EmbedLib.MatchFeatsToMol(m1,featFactory,pcop.getFeatures())
    match,mList2 = EmbedLib.MatchFeatsToMol(m2,featFactory,pcop.getFeatures())
    b1 = rdDistGeom.GetMoleculeBoundsMatrix(m1)
    b2 = rdDistGeom.GetMoleculeBoundsMatrix(m2)

    self.assertEqual(len(EmbedLib.MatchPharmacophore(mList1,b1,pcop)[2]),4)
    self.assertEqual(len(EmbedLib.MatchPharmacophore(mList2,b2,pcop)[2]),4)


    self.assertEqual(len(EmbedLib.MatchPharmacophore(mList1,b1,pcop,
                                                    mol=m1,use2DLimits=True)[2]),4)
    self.assertEqual(len(EmbedLib.MatchPharmacophore(mList2,b2,pcop,
                                                    mol=m2,use2DLimits=True)[2]),4)

    from rdkit import DistanceGeometry as DG
    self.assertTrue(DG.DoTriangleSmoothing(b1))
    self.assertTrue(DG.DoTriangleSmoothing(b2))

    self.assertEqual(len(EmbedLib.MatchPharmacophore(mList1,b1,pcop)[2]),4)
    self.assertEqual(len(EmbedLib.MatchPharmacophore(mList2,b2,pcop)[2]),4)
    
    self.assertEqual(len(EmbedLib.MatchPharmacophore(mList1,b1,pcop,
                                                    mol=m1,use2DLimits=True)[2]),4)
    self.assertEqual(len(EmbedLib.MatchPharmacophore(mList2,b2,pcop,
                                                    mol=m2,use2DLimits=True)[2]),4)

    
 
    
if __name__ == '__main__':
  unittest.main()

