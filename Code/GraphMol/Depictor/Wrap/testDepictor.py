## Automatically adapted for numpy.oldnumeric Jun 27, 2008 by -c

#
#  $Id$
#
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit import Geometry
from rdkit import RDConfig
import unittest
import os,sys
import cPickle as pickle
from rdkit.Chem.ChemUtils import AlignDepict
import numpy.oldnumeric as Numeric

def feq(v1,v2,tol2=1e-4):
  return abs(v1-v2)<=tol2

def ptEq(pt1, pt2, tol=1e-4):
  return feq(pt1.x,pt2.x,tol) and feq(pt1.y,pt2.y,tol) and feq(pt1.z,pt2.z,tol)

def getDistMat(mol):
    conf = mol.GetConformer()
    nat = mol.GetNumAtoms()
    nl = nat*(nat-1)/2
    res = Numeric.zeros(nl, Numeric.Float)

    for i in range(1,nat):
        pi = conf.GetAtomPosition(i)
        id = i*(i-1)/2
        for j in range(i):
            pj = conf.GetAtomPosition(j)
            pj -= pi
            res[id + j] = pj.Length()
                
    return res
  
def compareCoords(m, molFile):
  mo = Chem.MolFromMolFile(molFile)
  co = mo.GetConformer()

  ci = m.GetConformer()
  nat = m.GetNumAtoms()
  if (nat != mo.GetNumAtoms()):
    return 0

  for i in range(nat) :
      pos = ci.GetAtomPosition(i)
      opos = co.GetAtomPosition(i)
      if not ptEq(pos, opos):
        return 0
  return 1

def compareWithOld(smilesFile, sdFile) :
  smiSup = Chem.SmilesMolSupplier(smilesFile, ",", 0, -1)
  sdsup = Chem.SDMolSupplier(sdFile)
  im = 0
  for mol in smiSup :
    omol = sdsup[im]
    rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()
    oconf = omol.GetConformer()
    nat = mol.GetNumAtoms()
    for i in range(nat) :
      pos = conf.GetAtomPosition(i)
      opos = oconf.GetAtomPosition(i)
      if not ptEq(pos, opos):
        print >>sys.stderr,Chem.MolToMolBlock(omol)
        print >>sys.stderr,'> <Failed>\n%d\n'%i
        print >>sys.stderr,"$$$$"
        print >>sys.stderr,Chem.MolToMolBlock(mol)
        print >>sys.stderr,'> <Failed>\n%d\n'%i
        print >>sys.stderr,"$$$$"
        return 0
      
    im += 1
  return 1
        
class TestCase(unittest.TestCase) :
    def setUp(self):
        pass

    def _test0First200(self):
        # this test is disabled because it's not particularly useful and
        # causes problems every time anything changes.
        fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','Depictor',
                                            'test_data','first_200.tpsa.csv')
        #smiSup = Chem.SmilesMolSupplier(fileN, ",", 0, -1)

        ofile = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','Depictor',
                             'test_data', 'first_200.python.sdf')
        self.failUnless(compareWithOld(fileN, ofile))

    def test1CisTrans(self) :
        fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','Depictor',
                             'test_data', "cis_trans_cases.csv")
        ofile = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','Depictor',
                             'test_data', 'cis_trans_cases.sdf')
        
        self.failUnless(compareWithOld(fileN, ofile))
        
    def test2Coords(self) :
        m1 = Chem.MolFromSmiles('C1CCC1CC')
        coordMap = {0:Geometry.Point2D(0,0),
                    1:Geometry.Point2D(1.5,0),
                    2:Geometry.Point2D(1.5,1.5),
                    3:Geometry.Point2D(0,1.5)
                    }
        rdDepictor.Compute2DCoords(m1,coordMap=coordMap)
        conf = m1.GetConformer(0)
        for i in range(4):
            self.failUnless(ptEq(conf.GetAtomPosition(i),Geometry.Point3D(coordMap[i].x,
                                                                          coordMap[i].y,
                                                                          0.0)))

        m1 = Chem.MolFromSmiles('CCC')
        try:
            rdDepictor.Compute2DCoords(m1,coordMap=coordMap)
            ok = 0
        except ValueError:
            ok=1
        self.failUnless(ok)

    def test3IssueSF1526844(self):
      t = Chem.MolFromSmiles('c1nc(N)ccc1')
      rdDepictor.Compute2DCoords(t)
      
      m2 = Chem.MolFromSmiles('c1nc(NC=O)ccc1')
      AlignDepict.AlignDepict(m2,t)
      expected = [Geometry.Point3D(1.5, 0.0, 0.0),
                  Geometry.Point3D(0.75, -1.299, 0.0),
                  Geometry.Point3D(-0.75, -1.299, 0.0),
                  Geometry.Point3D(-1.5, -2.5981, 0.0),
                  Geometry.Point3D(-3.0, -2.5981, 0.0),
                  Geometry.Point3D(-3.75, -3.8971, 0.0),
                  Geometry.Point3D(-1.5, 0.0, 0.0),
                  Geometry.Point3D(-0.75, 1.2990, 0.0),
                  Geometry.Point3D(0.75, 1.2990, 0.0)]

      nat = m2.GetNumAtoms()
      conf = m2.GetConformer()
      for i in range(nat) :
        pos = conf.GetAtomPosition(i)
        self.failUnless(ptEq(pos, expected[i], 0.001))

    def test4SamplingSpread(self):
      mol= Chem.MolFromMolFile('../test_data/7UPJ_xtal.mol')

      # default mode
      rdDepictor.Compute2DCoords(mol)
      self.failUnless(compareCoords(mol, '../test_data/7UPJ_default.mol'))

      # spread the structure as much as possible by sampling
      rdDepictor.Compute2DCoords(mol, nFlipsPerSample=3, nSample=100,
                                 sampleSeed=100, permuteDeg4Nodes=1)
      self.failUnless(compareCoords(mol, '../test_data/7UPJ_spread.mol'))
      
    def test5SamplingMimic3D(self):
      mol = Chem.MolFromMolFile('../test_data/7UPJ_xtal.mol')
      dmat3D = getDistMat(mol)

      # now mimic the coordinate with a very small weight
      rdDepictor.Compute2DCoordsMimicDistmat(mol, dmat3D, weightDistMat=0.001)
      self.failUnless(compareCoords(mol, '../test_data/7UPJ_mimic3D_1.mol'))
      
      # now mimic the coordinate with a very small weight
      rdDepictor.Compute2DCoordsMimicDistmat(mol, dmat3D, weightDistMat=0.003)
      self.failUnless(compareCoords(mol, '../test_data/7UPJ_mimic3D_2.mol'))

      #mb = Chem.MolToMolBlock(mol)
      #ofile = open('../test_data/7UPJ_mimic3D_2.mol', 'w')
      #ofile.write(mb)
      #ofile.close()
      
if __name__ == '__main__':
  unittest.main()
