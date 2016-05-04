# $Id$
#
#  Copyright (C) 2004  Rational Discovery LLC
#         All Rights Reserved
#
from rdkit import RDConfig
import os,sys
import unittest
from rdkit import Chem
from rdkit.Chem import ChemicalForceFields
from rdkit import Geometry
from rdkit.Geometry import Point3D

def feq(v1,v2,tol2=1e-4):
  return abs(v1-v2)<=tol2

def ptEq(pt1, pt2, tol=1e-4):
  return feq(pt1.x,pt2.x,tol) and feq(pt1.y,pt2.y,tol) and feq(pt1.z,pt2.z,tol)

def addConf(mol):
  conf = Chem.Conformer(mol.GetNumAtoms())
  for i in range(mol.GetNumAtoms()):
    conf.SetAtomPosition(i,(0.,0.,0.))
  mol.AddConformer(conf)
  mb = Chem.MolToMolBlock(mol)
  mb = Chem.MolToMolBlock(mol)
  
class TestCase(unittest.TestCase):
  def setUp(self):
    pass

  def test0Conformers(self) :
    """Test the conformer data structure"""
    mol = Chem.MolFromSmiles("CC")
    conf = Chem.Conformer(2)
    conf.SetAtomPosition(0, (-0.5, 0.0, 0.0))
    conf.SetAtomPosition(1, (1.0, 0.0, 0.0))
    conf.SetId(0)
    cid = mol.AddConformer(conf)

    self.assertTrue(cid == 0)
    
    conf2 = mol.GetConformer(0)
    self.assertTrue(conf2.GetId() == cid)
    pt1 = conf2.GetAtomPosition(0)
    self.assertTrue(ptEq(pt1, Point3D(-0.5, 0.0, 0.0)))
    
    pt2 = conf2.GetAtomPosition(1)
    self.assertTrue(ptEq(pt2, Point3D(1.0, 0.0, 0.0)))
    
    #changing conf should not change conf2 - related to issue 217
    conf.SetAtomPosition(1, Point3D(2.0, 0.0, 0.0))
    pt2 = conf2.GetAtomPosition(1)
    self.assertTrue(feq(pt2[0], 1.0))

    conf = Chem.Conformer(2)
    conf.SetAtomPosition(0, Point3D(-0.5, 0.0, 0.0))
    conf.SetAtomPosition(1, Point3D(1.0, 0.0, 0.0))
    conf.SetId(2)

    cid = mol.AddConformer(conf, 0)
    self.assertTrue(cid == 2)
    conf3 = mol.GetConformer(2)
    
  def test0AddHds(self) :
    mol = Chem.MolFromSmiles("CC")
    conf = Chem.Conformer(1)
    conf.SetAtomPosition(0, Point3D(-0.5, 0.0, 0.0))
    conf.SetAtomPosition(1, Point3D(1.0, 0.0, 0.0))
    cid = mol.AddConformer(conf)

    conf2 = mol.GetConformer()
    self.assertTrue(conf2.GetNumAtoms() == 2)

    nmol = Chem.AddHs(mol, 0,1)
    conf3 = nmol.GetConformer()
    self.assertTrue(conf3.GetNumAtoms() == 8)
    self.assertTrue(conf2.GetNumAtoms() == 2)

    targetCoords = [[-0.5, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [-0.8667, 0.0, 1.03709],
                    [-0.8667, 0.8981, -0.5185],
                    [-0.8667, -0.8981, -0.5185],
                    [1.3667, 0.0, -1.0371],
                    [1.36667, 0.8981, 0.5185],
                    [1.36667, -0.8981, 0.5185]]

    for i in range(8) :
      pt = conf3.GetAtomPosition(i)
      self.assertTrue(ptEq(pt, Point3D(*tuple(targetCoords[i]))))

  def test2Issue217(self) :
    smi = 'c1ccccc1'
    m = Chem.MolFromSmiles(smi)
    addConf(m)
    self.assertTrue(m.GetNumConformers()==1);
    mb2 = Chem.MolToMolBlock(m)

  def test3Exceptions(self) :
    smi = 'c1ccccc1'
    m = Chem.MolFromSmiles(smi)
    addConf(m)
    self.assertTrue(m.GetNumConformers()==1)
    self.assertRaises(ValueError,lambda:m.GetConformer(2))

  def test4ConfTuple(self):
    smi = 'c1ccccc1'
    m = Chem.MolFromSmiles(smi)
    for i in range(10):
      addConf(m)

    confs = m.GetConformers()
    self.assertTrue(len(confs) == 10)

    for conf in confs:
      for i in range(6):
        pt = conf.GetAtomPosition(i)
        self.assertTrue(ptEq(pt, Point3D(0.0, 0.0, 0.0)))

    m.RemoveAllConformers()
    self.assertTrue(m.GetNumConformers() == 0)
    confs = m.GetConformers()
    self.assertTrue(confs == ())

  def testAddConformersFromTrajectory(self):
    molBlock = \
      '\n' \
      '     RDKit          3D\n' \
      '\n' \
      ' 71 74  0  0  0  0  0  0  0  0999 V2000\n' \
      '    8.2543    3.1901   -0.3005 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    7.4558    1.9712    0.0938 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    7.3934    1.0441   -0.9483 O   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    6.6660   -0.0533   -0.4641 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    5.1928    0.2346   -0.4609 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    4.3713   -0.9410   -0.5770 N   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    3.1852   -1.0034   -1.2291 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    2.2914    0.1276   -1.6316 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    0.9308   -0.4468   -1.9908 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    0.1417   -0.7821   -0.7545 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -0.1848    0.3695    0.0456 N   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -1.5661    0.7686   -0.0745 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -2.4768   -0.0640    0.8206 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -3.8874    0.1143    0.3941 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -4.6333   -0.9984    0.0264 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -6.0127   -0.9516   -0.0400 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -6.7062    0.1599    0.3963 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -8.0408    0.4828   -0.1977 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -7.7914    1.1180   -1.5591 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -8.7622    1.4403    0.7265 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -8.8409   -0.7397   -0.4395 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -8.9121   -1.6637    0.4258 O   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -9.7414   -0.7636   -1.5059 O   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -5.9736    1.2357    0.8565 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -4.5843    1.2252    0.8530 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    0.6263    1.4884   -0.3942 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    2.0541    1.0258   -0.4230 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    2.9225   -2.3317   -1.2963 N   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    3.6061   -2.9745   -0.3180 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    3.3554   -4.1536    0.3735 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    3.7653   -4.2712    1.6948 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    4.8254   -3.4613    2.0796 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    5.1978   -2.3436    1.3419 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    4.5694   -2.0799    0.1305 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    9.3138    3.1372    0.0031 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    7.8117    4.0754    0.1798 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    8.2358    3.3535   -1.4074 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    6.4027    2.2146    0.3634 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    7.9270    1.5444    1.0040 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    7.0677   -0.2415    0.5615 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    6.9530   -0.9105   -1.1025 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    4.9578    0.7259    0.5137 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    4.9985    0.9430   -1.3033 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    2.7171    0.7264   -2.4494 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    0.3994    0.2339   -2.6810 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    1.1342   -1.4171   -2.5076 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -0.7632   -1.3370   -1.0391 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    0.7845   -1.4394   -0.1311 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    0.0125    0.1989    1.0673 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -1.6672    1.8215    0.2925 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -1.8705    0.7271   -1.1337 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -2.3045    0.3159    1.8590 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -2.1980   -1.1367    0.7635 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -4.1513   -1.9468   -0.2114 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -6.6138   -1.7460   -0.4718 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -7.0727    0.4399   -2.0858 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -7.3144    2.1076   -1.4482 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -8.7609    1.1720   -2.1135 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -8.3137    2.4504    0.5729 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -8.6170    1.0817    1.7580 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -9.8244    1.4444    0.4200 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -6.4629    2.0541    1.3719 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '   -4.0445    2.0563    1.3058 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    0.3329    1.8224   -1.3991 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    0.4920    2.3164    0.3160 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    2.2025    0.3766    0.4766 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    2.7945    1.8369   -0.3969 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    2.4404   -4.6964    0.1303 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    3.3157   -5.0055    2.3587 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    5.4272   -3.7654    2.9380 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '    5.5668   -1.5069    1.9380 H   0  0  0  0  0  0  0  0  0  0  0  0\n' \
      '  1  2  1  0\n' \
      '  2  3  1  0\n' \
      '  3  4  1  0\n' \
      '  4  5  1  0\n' \
      '  5  6  1  0\n' \
      '  6  7  1  0\n' \
      '  7  8  1  0\n' \
      '  8  9  1  0\n' \
      '  9 10  1  0\n' \
      ' 10 11  1  0\n' \
      ' 11 12  1  0\n' \
      ' 12 13  1  0\n' \
      ' 13 14  1  0\n' \
      ' 14 15  2  0\n' \
      ' 15 16  1  0\n' \
      ' 16 17  2  0\n' \
      ' 17 18  1  0\n' \
      ' 18 19  1  0\n' \
      ' 18 20  1  0\n' \
      ' 18 21  1  0\n' \
      ' 21 22  2  0\n' \
      ' 21 23  1  0\n' \
      ' 17 24  1  0\n' \
      ' 24 25  2  0\n' \
      ' 11 26  1  0\n' \
      ' 26 27  1  0\n' \
      '  7 28  2  0\n' \
      ' 28 29  1  0\n' \
      ' 29 30  2  0\n' \
      ' 30 31  1  0\n' \
      ' 31 32  2  0\n' \
      ' 32 33  1  0\n' \
      ' 33 34  2  0\n' \
      ' 34  6  1  0\n' \
      ' 27  8  1  0\n' \
      ' 34 29  1  0\n' \
      ' 25 14  1  0\n' \
      '  1 35  1  0\n' \
      '  1 36  1  0\n' \
      '  1 37  1  0\n' \
      '  2 38  1  0\n' \
      '  2 39  1  0\n' \
      '  4 40  1  0\n' \
      '  4 41  1  0\n' \
      '  5 42  1  0\n' \
      '  5 43  1  0\n' \
      '  8 44  1  0\n' \
      '  9 45  1  0\n' \
      '  9 46  1  0\n' \
      ' 10 47  1  0\n' \
      ' 10 48  1  0\n' \
      ' 11 49  1  0\n' \
      ' 12 50  1  0\n' \
      ' 12 51  1  0\n' \
      ' 13 52  1  0\n' \
      ' 13 53  1  0\n' \
      ' 15 54  1  0\n' \
      ' 16 55  1  0\n' \
      ' 19 56  1  0\n' \
      ' 19 57  1  0\n' \
      ' 19 58  1  0\n' \
      ' 20 59  1  0\n' \
      ' 20 60  1  0\n' \
      ' 20 61  1  0\n' \
      ' 24 62  1  0\n' \
      ' 25 63  1  0\n' \
      ' 26 64  1  0\n' \
      ' 26 65  1  0\n' \
      ' 27 66  1  0\n' \
      ' 27 67  1  0\n' \
      ' 30 68  1  0\n' \
      ' 31 69  1  0\n' \
      ' 32 70  1  0\n' \
      ' 33 71  1  0\n' \
      'M  CHG  2  11   1  23  -1\n' \
      'M  END\n'
    mol = Chem.MolFromMolBlock(molBlock, removeHs = False)
    everySteps = 10
    maxIts = 1000
    gradTol = 0.01
    rdbase = os.environ['RDBASE']
    fName = os.path.join(rdbase, 'Code', 'GraphMol', 'Wrap', 'test_data', 'bilastine_trajectory.sdf')
    w = Chem.SDWriter(fName)
    field = ChemicalForceFields.MMFFGetMoleculeForceField(mol,
      ChemicalForceFields.MMFFGetMoleculeProperties(mol))
    traj = Geometry.Trajectory(3, mol.GetNumAtoms())
    field.MinimizeTrajectory(everySteps, traj, maxIts, gradTol);
    mol.RemoveConformer(0)
    mol.AddConformersFromTrajectory(traj)
    for nConf in range(mol.GetNumConformers()):
      mol.SetProp('ENERGY', '{0:.4f}'.format(traj.GetSnapshot(nConf).GetEnergy()))
      w.write(mol, nConf)
    w.close()

  def testAddConformersFromAmberTrajectory(self):
    mol = Chem.MolFromSmiles('CCC')
    rdbase = os.environ['RDBASE']
    fName = os.path.join(rdbase, 'Code', 'Geometry', 'testData', 'water_coords.trx')
    traj = Geometry.Trajectory(3, mol.GetNumAtoms())
    Geometry.ReadAmberTrajectory(fName, traj)
    self.assertEqual(len(traj), 1)
    for i in range(2):
      mol.AddConformersFromTrajectory(traj)
      self.assertEqual(mol.GetNumConformers(), i + 1)
      self.assertEqual(mol.GetConformer(i).GetNumAtoms(), 3)
      self.assertAlmostEqual(mol.GetConformer(i).GetAtomPosition(0).x, 0.1941767)
      self.assertAlmostEqual(mol.GetConformer(i).GetAtomPosition(2).z, -0.4088006)
    mol.RemoveAllConformers()
    mol.AddConformersFromTrajectory(traj, 0)
    self.assertEqual(mol.GetNumConformers(), 0)
    fName = os.path.join(rdbase, 'Code', 'Geometry', 'testData', 'water_coords2.trx')
    traj = Geometry.Trajectory(3, mol.GetNumAtoms())
    Geometry.ReadAmberTrajectory(fName, traj)
    self.assertEqual(len(traj), 2)
    mol.AddConformersFromTrajectory(traj)
    self.assertEqual(len(traj), 2)
    mol.RemoveAllConformers()
    mol.AddConformersFromTrajectory(traj, 1)
    self.assertEqual(mol.GetNumConformers(), 1)

  def testAddConformersFromGromosTrajectory(self):
    mol = Chem.MolFromSmiles('CCC')
    rdbase = os.environ['RDBASE']
    fName = os.path.join(rdbase, 'Code', 'Geometry', 'testData', 'water_coords.trc')
    traj = Geometry.Trajectory(3, mol.GetNumAtoms())
    Geometry.ReadGromosTrajectory(fName, traj)
    self.assertEqual(len(traj), 1)
    for i in range(2):
      mol.AddConformersFromTrajectory(traj)
      self.assertEqual(mol.GetNumConformers(), i + 1)
      self.assertEqual(mol.GetConformer(i).GetNumAtoms(), 3)
      self.assertAlmostEqual(mol.GetConformer(i).GetAtomPosition(0).x, 1.941767)
      self.assertAlmostEqual(mol.GetConformer(i).GetAtomPosition(2).z, -4.088006)
    mol.RemoveAllConformers()
    mol.AddConformersFromTrajectory(traj, 0)
    self.assertEqual(mol.GetNumConformers(), 0)
    fName = os.path.join(rdbase, 'Code', 'Geometry', 'testData', 'water_coords2.trc')
    traj = Geometry.Trajectory(3, mol.GetNumAtoms())
    Geometry.ReadGromosTrajectory(fName, traj)
    self.assertEqual(len(traj), 2)
    mol.AddConformersFromTrajectory(traj)
    self.assertEqual(len(traj), 2)
    mol.RemoveAllConformers()
    mol.AddConformersFromTrajectory(traj, 1)
    self.assertEqual(mol.GetNumConformers(), 1)

if __name__ == '__main__':
  unittest.main()


