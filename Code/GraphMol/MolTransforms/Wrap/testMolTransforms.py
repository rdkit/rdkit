from rdkit import RDConfig
import os, sys, math
import unittest
import numpy as np
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Geometry import rdGeometry as geom
from rdkit.Chem import rdMolTransforms as rdmt

def feq(v1, v2, tol=1.0e-4):
  return abs(v1 - v2) < tol


def ptEq(pt, tpt, tol=0.0001):
  pt -= tpt
  return feq(pt.Length(), 0.0)


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def test1Canonicalization(self):
    mol = Chem.MolFromSmiles("C")
    conf = Chem.Conformer(1)
    conf.SetAtomPosition(0, (4.0, 5.0, 6.0))
    mol.AddConformer(conf, 1)

    conf = mol.GetConformer()
    pt = rdmt.ComputeCentroid(conf)
    self.assertTrue(ptEq(pt, geom.Point3D(4.0, 5.0, 6.0)))

    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MolTransforms', 'test_data',
                         '1oir.mol')
    m = Chem.MolFromMolFile(fileN)
    cpt = rdmt.ComputeCentroid(m.GetConformer())
    trans = rdmt.ComputeCanonicalTransform(m.GetConformer(), cpt)
    trans2 = rdmt.ComputeCanonicalTransform(m.GetConformer())
    for i in range(4):
      for j in range(4):
        self.assertTrue(feq(trans[i, j], trans2[i, j]))
    rdmt.TransformConformer(m.GetConformer(), trans2)
    m2 = Chem.MolFromMolFile(fileN)
    rdmt.CanonicalizeConformer(m2.GetConformer())
    nats = m.GetNumAtoms()
    cnf1 = m.GetConformer()
    cnf2 = m2.GetConformer()
    for i in range(nats):
      p1 = list(cnf1.GetAtomPosition(i))
      p2 = list(cnf2.GetAtomPosition(i))
      self.assertTrue(feq(p1[0], p2[0]))
      self.assertTrue(feq(p1[1], p2[1]))
      self.assertTrue(feq(p1[2], p2[2]))

    m3 = Chem.MolFromMolFile(fileN)
    rdmt.CanonicalizeMol(m3)
    cnf1 = m.GetConformer()
    cnf2 = m3.GetConformer()
    for i in range(nats):
      p1 = list(cnf1.GetAtomPosition(i))
      p2 = list(cnf2.GetAtomPosition(i))
      self.assertTrue(feq(p1[0], p2[0]))
      self.assertTrue(feq(p1[1], p2[1]))
      self.assertTrue(feq(p1[2], p2[2]))

  def testComputePrincipalAxesAndMoments(self):
    if (not hasattr(rdmt, 'ComputePrincipalAxesAndMoments')):
      return
    molBlock = '''\

     RDKit          3D

  4  3  0  0  0  0  0  0  0  0999 V2000
   -0.1888    1.3224   -0.2048 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0790    0.1671    0.6070 C   0  0  1  0  0  0  0  0  0  0  0  0
   -1.3083   -0.8939   -0.2236 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    1.5761   -0.5956   -0.1786 Br  0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  6
  2  3  1  0
  2  4  1  0
M  END
'''
    axesRef = (
      (-0.9997, -0.0246,  0.0009),
      ( 0.0246, -0.9981,  0.0559),
      ( 0.0004, -0.0559, -0.9984)
    )
    momentsRef = ( 3.4220,  4.7230,  7.1757)
    m = Chem.MolFromMolBlock(molBlock)
    axes, moments = rdmt.ComputePrincipalAxesAndMoments(m.GetConformer())
    self.assertIsNotNone(axes)
    self.assertIsNotNone(moments)
    for y in range(3):
      for x in range(3):
        self.assertAlmostEqual(axes[y][x], axesRef[y][x], 3)
      self.assertAlmostEqual(moments[y], momentsRef[y], 3)
    failed = False
    try:
      axes, moments = rdmt.ComputePrincipalAxesAndMoments(m.GetConformer(), weights = (0.5, 0.5))
    except Exception:
      failed = True
    self.assertTrue(failed)
    axesWeightedRef = (
      (-0.9998, -0.0114, -0.0189),
      (-0.0153,  0.9744,  0.2245),
      ( 0.0158,  0.2247, -0.9743)
    )
    momentsWeightedRef = ( 0.5496,  1.5559,  1.9361)
    axesWeighted, momentsWeighted = rdmt.ComputePrincipalAxesAndMoments(
                                    m.GetConformer(), weights = (0.1, 0.2, 0.3, 0.4))
    self.assertIsNotNone(axesWeighted)
    self.assertIsNotNone(momentsWeighted)
    for y in range(3):
      for x in range(3):
        self.assertAlmostEqual(axesWeighted[y][x], axesWeightedRef[y][x], 3)
      self.assertAlmostEqual(momentsWeighted[y], momentsWeightedRef[y], 3)

  def testComputePrincipalAxesAndMomentsFromGyrationMatrix(self):
    if (not hasattr(rdmt, 'ComputePrincipalAxesAndMomentsFromGyrationMatrix')):
      return
    molBlock = '''\

     RDKit          3D

  4  3  0  0  0  0  0  0  0  0999 V2000
   -0.1888    1.3224   -0.2048 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0790    0.1671    0.6070 C   0  0  1  0  0  0  0  0  0  0  0  0
   -1.3083   -0.8939   -0.2236 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    1.5761   -0.5956   -0.1786 Br  0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  6
  2  3  1  0
  2  4  1  0
M  END
'''
    axesRef = (
      (-0.0009, -0.0246,  0.9997),
      (-0.0559, -0.9981, -0.0246),
      ( 0.9984, -0.0559, -0.0004)
    )
    momentsRef = ( 0.1212,  0.7343,  1.0596)
    m = Chem.MolFromMolBlock(molBlock)
    axes, moments = rdmt.ComputePrincipalAxesAndMomentsFromGyrationMatrix(m.GetConformer())
    self.assertIsNotNone(axes)
    self.assertIsNotNone(moments)
    for y in range(3):
      for x in range(3):
        self.assertAlmostEqual(axes[y][x], axesRef[y][x], 3)
      self.assertAlmostEqual(moments[y], momentsRef[y], 3)
    failed = False
    try:
      axes, moments = rdmt.ComputePrincipalAxesAndMomentsFromGyrationMatrix(m.GetConformer(), weights = (0.5, 0.5))
    except Exception:
      failed = True
    self.assertTrue(failed)
    axesWeightedRef = (
      ( 0.0189, -0.0114,  0.9998),
      (-0.2245,  0.9744,  0.0153),
      ( 0.9743,  0.2247, -0.0158)
    )
    momentsWeightedRef = (  0.0847,  0.4649,  1.4712)
    axesWeighted, momentsWeighted = rdmt.ComputePrincipalAxesAndMomentsFromGyrationMatrix(
                                    m.GetConformer(), weights = (0.1, 0.2, 0.3, 0.4))
    self.assertIsNotNone(axesWeighted)
    self.assertIsNotNone(momentsWeighted)
    for y in range(3):
      for x in range(3):
        self.assertAlmostEqual(abs(axesWeighted[y][x]), abs(axesWeightedRef[y][x]), 3)
      self.assertAlmostEqual(abs(momentsWeighted[y]), abs(momentsWeightedRef[y]), 3)

  def testGetSetBondLength(self):
    file = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MolTransforms', 'test_data',
                        '3-cyclohexylpyridine.mol')

    m = Chem.MolFromMolFile(file, True, False)
    conf = m.GetConformer()
    dist = rdmt.GetBondLength(conf, 0, 19)
    self.assertAlmostEqual(dist, 1.36, 2)
    rdmt.SetBondLength(conf, 0, 19, 2.5)
    dist = rdmt.GetBondLength(conf, 0, 19)
    self.assertAlmostEqual(dist, 2.5, 1)
    rdmt.SetBondLength(conf, 19, 0, 3.0)
    dist = rdmt.GetBondLength(conf, 0, 19)
    self.assertAlmostEqual(dist, 3.0, 1)

  def testGetSetAngle(self):
    file = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MolTransforms', 'test_data',
                        '3-cyclohexylpyridine.mol')

    m = Chem.MolFromMolFile(file, True, False)
    conf = m.GetConformer()
    angle = rdmt.GetAngleDeg(conf, 0, 19, 21)
    self.assertAlmostEqual(angle, 109.7, 1)
    rdmt.SetAngleDeg(conf, 0, 19, 21, 125.0)
    angle = rdmt.GetAngleDeg(conf, 0, 19, 21)
    self.assertAlmostEqual(angle, 125.0, 1)
    rdmt.SetAngleRad(conf, 21, 19, 0, math.pi / 2.)
    angle = rdmt.GetAngleRad(conf, 0, 19, 21)
    self.assertAlmostEqual(angle, math.pi / 2., 1)
    angle = rdmt.GetAngleDeg(conf, 0, 19, 21)
    self.assertAlmostEqual(angle, 90.0, 1)

  def testGetSetDihedral(self):
    file = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MolTransforms', 'test_data',
                        '3-cyclohexylpyridine.mol')

    m = Chem.MolFromMolFile(file, True, False)
    conf = m.GetConformer()
    dihedral = rdmt.GetDihedralDeg(conf, 0, 19, 21, 24)
    self.assertAlmostEqual(dihedral, 176.05, 2)
    rdmt.SetDihedralDeg(conf, 8, 0, 19, 21, 65.0)
    dihedral = rdmt.GetDihedralDeg(conf, 8, 0, 19, 21)
    self.assertAlmostEqual(dihedral, 65.0, 1)
    rdmt.SetDihedralDeg(conf, 8, 0, 19, 21, -130.0)
    dihedral = rdmt.GetDihedralDeg(conf, 8, 0, 19, 21)
    self.assertAlmostEqual(dihedral, -130.0, 1)
    rdmt.SetDihedralRad(conf, 21, 19, 0, 8, -2. / 3. * math.pi)
    dihedral = rdmt.GetDihedralRad(conf, 8, 0, 19, 21)
    self.assertAlmostEqual(dihedral, -2. / 3. * math.pi, 1)
    dihedral = rdmt.GetDihedralDeg(conf, 8, 0, 19, 21)
    self.assertAlmostEqual(dihedral, -120.0, 1)

  def testGetSetDihedralThroughTripleBond(self):
    file = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MolTransforms', 'test_data',
                        'github1262_2.mol')

    m = Chem.MolFromMolFile(file, True, False)
    conf = m.GetConformer()
    rdmt.SetDihedralDeg(conf, 6, 1, 2, 9, 0.0)
    dihedral = rdmt.GetDihedralDeg(conf, 6, 1, 2, 9)
    self.assertAlmostEqual(dihedral, 0.0, 1)
    dist = rdmt.GetBondLength(conf, 6, 9)
    rdmt.SetDihedralDeg(conf, 6, 1, 2, 9, 120.0)
    dihedral = rdmt.GetDihedralDeg(conf, 6, 1, 2, 9)
    self.assertAlmostEqual(dihedral, 120.0, 1)
    dist2 = rdmt.GetBondLength(conf, 6, 7)
    self.assertAlmostEqual(dist, dist2, 1)
    rdmt.SetDihedralDeg(conf, 6, 1, 2, 9, 180.0)
    dihedral = rdmt.GetDihedralDeg(conf, 6, 1, 2, 9)
    self.assertAlmostEqual(dihedral, 180.0, 1)
    dist3 = rdmt.GetBondLength(conf, 6, 9)
    self.assertNotAlmostEqual(dist, dist3, 1)
    exceptionRaised = False
    try:
      rdmt.SetDihedralDeg(conf, 6, 0, 3, 9, 0.0)
    except ValueError:
      exceptionRaised = True
    self.assertTrue(exceptionRaised)

  def testEigen3CanonicalTransformAgainstNumpy(self):
    def canonicalize_conf_rdkit(mol, conf_id=-1):
      mol = Chem.Mol(mol)
      conf = mol.GetConformer(conf_id)
      ctd = rdmt.ComputeCentroid(conf)
      canon_trans = rdmt.ComputeCanonicalTransform(conf, ctd)
      rdmt.TransformConformer(conf, canon_trans)
      return mol

    def canonicalize_conf_numpy(mol, conf_id=-1):
      mol = Chem.Mol(mol)
      conf = mol.GetConformer(conf_id)
      pos = conf.GetPositions()
      ctd = pos.mean(axis=0)
      trans_pos = pos - ctd
      cov_mat = np.cov(trans_pos, bias=1, rowvar=False) * conf.GetNumAtoms()
      eigval, eigvect = np.linalg.eig(cov_mat)
      eigval_sorted = sorted(enumerate(eigval), key=lambda x: x[1], reverse=True)
      eigvect_sorted = [eigvect[:, i] * (1.0 if eigvect[:, i].sum() > 0.0 else -1.0) for i, _ in eigval_sorted]
      canon_trans = np.array([
        [*eigvect_sorted[0], 0.],
        [*eigvect_sorted[1], 0.],
        [*eigvect_sorted[2], 0.],
        [0., 0., 0., 1.],
      ], dtype=np.double)
      canon_trans[0:3, 3] = np.array([*ctd, 1.]).dot(canon_trans)[:3]
      rdmt.TransformConformer(conf, canon_trans)
      return mol

    file = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MolTransforms', 'test_data',
                        'github4302.sdf')
    built_against_eigen3 = hasattr(rdmt, 'ComputePrincipalAxesAndMomentsFromGyrationMatrix')
    with Chem.SDMolSupplier(file) as suppl:
      for mol in suppl:
        mol_orig = Chem.Mol(mol)
        rdkit_canon = canonicalize_conf_rdkit(mol)
        numpy_canon = canonicalize_conf_numpy(mol)
        for i in range(mol.GetNumAtoms()):
          orig_coord = mol_orig.GetConformer().GetAtomPosition(i)
          rdkit_coord = rdkit_canon.GetConformer().GetAtomPosition(i)
          numpy_coord = numpy_canon.GetConformer().GetAtomPosition(i)
          self.assertLess(orig_coord.Distance(numpy_coord), 1.e-4)
          if built_against_eigen3:
            self.assertLess(orig_coord.Distance(rdkit_coord), 1.e-4)


if __name__ == "__main__":
  unittest.main()
