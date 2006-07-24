#
#  $Id$
#
import Chem
from Chem import rdDepictor
import Geometry
import RDConfig
import unittest
import os,sys
import cPickle as pickle
from Chem.ChemUtils import AlignDepict

def feq(v1,v2,tol2=1e-4):
  return abs(v1-v2)<=tol2

def ptEq(pt1, pt2, tol=1e-4):
  return feq(pt1.x,pt2.x,tol) and feq(pt1.y,pt2.y,tol) and feq(pt1.z,pt2.z,tol)


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

    def test0First200(self):
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
        
if __name__ == '__main__':
  unittest.main()
