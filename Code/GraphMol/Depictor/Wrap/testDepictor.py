# Automatically adapted for numpy.oldnumeric Jun 27, 2008 by -c

import os
import sys
import tempfile
#
#  $Id: testDepictor.py 2112 2012-07-02 09:47:45Z glandrum $
#
# pylint:disable=E1101,C0111,C0103,R0904
import unittest

import numpy as np

from rdkit import Chem, Geometry, RDConfig
from rdkit.Chem import rdDepictor, rdMolAlign, rdMolTransforms
from rdkit.Chem.ChemUtils import AlignDepict


def feq(v1, v2, tol2=1e-4):
  return abs(v1 - v2) <= tol2


def ptEq(pt1, pt2, tol=1e-4):
  return feq(pt1.x, pt2.x, tol) and feq(pt1.y, pt2.y, tol) and feq(pt1.z, pt2.z, tol)


def getDistMat(mol):
  conf = mol.GetConformer()
  nat = mol.GetNumAtoms()
  nl = nat * (nat - 1) // 2
  res = np.zeros(nl, float)

  for i in range(1, nat):
    pi = conf.GetAtomPosition(i)
    idx = i * (i - 1) // 2
    for j in range(i):
      pj = conf.GetAtomPosition(j)
      pj -= pi
      res[idx + j] = pj.Length()

  return res


def compareCoords(m, molFile):
  mo = Chem.MolFromMolFile(molFile)
  co = mo.GetConformer()

  ci = m.GetConformer()
  nat = m.GetNumAtoms()
  if (nat != mo.GetNumAtoms()):
    return 0

  for i in range(nat):
    pos = ci.GetAtomPosition(i)
    opos = co.GetAtomPosition(i)
    if not ptEq(pos, opos):
      print(Chem.MolToMolBlock(m))
      print(Chem.MolToMolBlock(mo))
      return 0
  return 1


def compareWithOld(smilesFile, sdFile):
  smiSup = Chem.SmilesMolSupplier(smilesFile, ",", 0, -1)
  sdsup = Chem.SDMolSupplier(sdFile)
  im = 0
  for mol in smiSup:
    omol = sdsup[im]
    rdDepictor.Compute2DCoords(mol, canonOrient=False)
    conf = mol.GetConformer()
    oconf = omol.GetConformer()
    nat = mol.GetNumAtoms()
    for i in range(nat):
      pos = conf.GetAtomPosition(i)
      opos = oconf.GetAtomPosition(i)
      if not ptEq(pos, opos):
        print(Chem.MolToMolBlock(omol), file=sys.stderr)
        print('> <Failed>\n%d\n' % i, file=sys.stderr)
        print("$$$$", file=sys.stderr)
        print(Chem.MolToMolBlock(mol), file=sys.stderr)
        print('> <Failed>\n%d\n' % i, file=sys.stderr)
        print("$$$$", file=sys.stderr)
        return 0
    im += 1
  return 1


def stereoCompare(smilesFile):
  smiSup = Chem.SmilesMolSupplier(smilesFile, ",", 0, -1)
  for mol in smiSup:
    rdDepictor.Compute2DCoords(mol, canonOrient=False)
    mb = Chem.MolToMolBlock(mol)
    nmol = Chem.MolFromMolBlock(mb)
    matches = nmol.GetSubstructMatches(mol, False)
    dbnds = [
      x for x in mol.GetBonds()
      if (x.GetBondType() == Chem.BondType.DOUBLE and x.GetStereo() > Chem.BondStereo.STEREOANY)
    ]
    ok = True
    for match in matches:
      for bnd in dbnds:
        obnd = nmol.GetBondBetweenAtoms(match[bnd.GetBeginAtomIdx()], match[bnd.GetEndAtomIdx()])
        assert (obnd.GetBondType() == Chem.BondType.DOUBLE)
      if ok:
        break
    if not ok:
      print(Chem.MolToMolBlock(mol), file=sys.stderr)
      print("$$$$", file=sys.stderr)
      return 0
  return 1


class TestCase(unittest.TestCase):

  def _test0First200(self):
    # this test is disabled because it's not particularly useful and
    # causes problems every time anything changes.
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Depictor', 'test_data',
                         'first_200.tpsa.csv')
    #smiSup = Chem.SmilesMolSupplier(fileN, ",", 0, -1)

    ofile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Depictor', 'test_data',
                         'first_200.python.sdf')
    self.assertTrue(compareWithOld(fileN, ofile))

  def test1CisTrans(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Depictor', 'test_data',
                         "cis_trans_cases.csv")
    self.assertTrue(stereoCompare(fileN))

  def test2Coords(self):
    m1 = Chem.MolFromSmiles('C1CCC1CC')
    coordMap = {
      0: Geometry.Point2D(0, 0),
      1: Geometry.Point2D(1.5, 0),
      2: Geometry.Point2D(1.5, 1.5),
      3: Geometry.Point2D(0, 1.5)
    }
    rdDepictor.Compute2DCoords(m1, coordMap=coordMap)
    conf = m1.GetConformer(0)
    for i in range(4):
      self.assertTrue(
        ptEq(conf.GetAtomPosition(i), Geometry.Point3D(coordMap[i].x, coordMap[i].y, 0.0)))

    m1 = Chem.MolFromSmiles('CCC')
    try:
      rdDepictor.Compute2DCoords(m1, coordMap=coordMap)
      ok = 0
    except ValueError:
      ok = 1
    self.assertTrue(ok)

  def test3IssueSF1526844(self):
    t = Chem.MolFromSmiles('c1nc(N)ccc1')
    rdDepictor.Compute2DCoords(t, canonOrient=False)

    m2 = Chem.MolFromSmiles('c1nc(NC=O)ccc1')
    AlignDepict.AlignDepict(m2, t)
    expected = [
      Geometry.Point3D(1.5, 0.0, 0.0),
      Geometry.Point3D(0.75, -1.299, 0.0),
      Geometry.Point3D(-0.75, -1.299, 0.0),
      Geometry.Point3D(-1.5, -2.5981, 0.0),
      Geometry.Point3D(-3.0, -2.5981, 0.0),
      Geometry.Point3D(-3.75, -3.8971, 0.0),
      Geometry.Point3D(-1.5, 0.0, 0.0),
      Geometry.Point3D(-0.75, 1.2990, 0.0),
      Geometry.Point3D(0.75, 1.2990, 0.0)
    ]

    nat = m2.GetNumAtoms()
    conf = m2.GetConformer()
    for i in range(nat):
      pos = conf.GetAtomPosition(i)
      self.assertTrue(ptEq(pos, expected[i], 0.001))

  def test4SamplingSpread(self):
    # the expected results here were generated with the legacy stereo code,
    # so we need to use that
    origVal = Chem.GetUseLegacyStereoPerception()
    Chem.SetUseLegacyStereoPerception(True)

    mol = Chem.MolFromMolFile(
      os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/Depictor', 'test_data/7UPJ_xtal.mol'))

    # default mode
    rdDepictor.Compute2DCoords(mol, canonOrient=False)
    self.assertTrue(
      compareCoords(
        mol, os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/Depictor',
                          'test_data/7UPJ_default.mol')))

    # spread the structure as much as possible by sampling
    rdDepictor.Compute2DCoords(mol, canonOrient=False, nFlipsPerSample=3, nSample=100,
                               sampleSeed=100, permuteDeg4Nodes=1)
    self.assertTrue(
      compareCoords(
        mol, os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/Depictor',
                          'test_data/7UPJ_spread.mol')))

    Chem.SetUseLegacyStereoPerception(origVal)

  def test5SamplingMimic3D(self):
    # the expected results here were generated with the legacy stereo code,
    # so we need to use that
    origVal = Chem.GetUseLegacyStereoPerception()
    Chem.SetUseLegacyStereoPerception(True)

    mol = Chem.MolFromMolFile(
      os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/Depictor', 'test_data/7UPJ_xtal.mol'))
    dmat3D = getDistMat(mol)

    # now mimic the coordinate with a very small weight
    rdDepictor.Compute2DCoordsMimicDistmat(mol, dmat3D, weightDistMat=0.001)
    self.assertTrue(
      compareCoords(
        mol,
        os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/Depictor', 'test_data/7UPJ_mimic3D_1.mol')))

    # now mimic the coordinate with a very small weight
    rdDepictor.Compute2DCoordsMimicDistmat(mol, dmat3D, weightDistMat=0.003)
    self.assertTrue(
      compareCoords(
        mol,
        os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/Depictor', 'test_data/7UPJ_mimic3D_2.mol')))

    Chem.SetUseLegacyStereoPerception(origVal)
    #mb = Chem.MolToMolBlock(mol)
    #ofile = open('../test_data/7UPJ_mimic3D_2.mol', 'w')
    # ofile.write(mb)
    # ofile.close()

  def test6ChangeBondLength(self):
    m = Chem.MolFromSmiles('CC')
    rdDepictor.Compute2DCoords(m)
    conf = m.GetConformer()
    self.assertAlmostEqual(conf.GetAtomPosition(0).x, -0.750, 3)
    self.assertAlmostEqual(conf.GetAtomPosition(1).x, 0.750, 3)
    rdDepictor.Compute2DCoords(m, bondLength=1.0)
    conf = m.GetConformer()
    self.assertAlmostEqual(conf.GetAtomPosition(0).x, -0.500, 3)
    self.assertAlmostEqual(conf.GetAtomPosition(1).x, 0.500, 3)
    rdDepictor.Compute2DCoords(m)
    conf = m.GetConformer()
    self.assertAlmostEqual(conf.GetAtomPosition(0).x, -0.750, 3)
    self.assertAlmostEqual(conf.GetAtomPosition(1).x, 0.750, 3)

  def testConstrainedCoords(self):
    templ = Chem.MolFromSmiles('c1nccc2n1ccc2')
    rdDepictor.Compute2DCoords(templ)
    m1 = Chem.MolFromSmiles('c1cccc2ncn3cccc3c21')
    rdDepictor.GenerateDepictionMatching2DStructure(m1, templ)
    m2 = Chem.MolFromSmiles('c1cc(Cl)cc2ncn3cccc3c21')
    rdDepictor.Compute2DCoords(m2)
    refPatt1 = Chem.MolFromSmarts('*1****2*1***2')
    rdDepictor.GenerateDepictionMatching2DStructure(m2, templ, -1, refPatt1)
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Depictor', 'test_data',
                         '1XP0_ligand.sdf')

    xp0_lig = Chem.MolFromMolFile(fileN)
    xp0_lig_2d = Chem.Mol(xp0_lig)
    rdDepictor.GenerateDepictionMatching3DStructure(xp0_lig_2d, xp0_lig)
    xp0_ref = Chem.MolFromSmarts('[#6]1~[#7][#6]~[#6]2[#6](=[#8])[#7]~[#6](c3ccccc3)[#7][#7]12')
    rdDepictor.GenerateDepictionMatching3DStructure(xp0_lig_2d, xp0_lig, -1, xp0_ref)

  def testGenerate2DDepictionRefPatternAtomMap(self):
    indazoleMolblock = """
     RDKit          2D

  9 10  0  0  0  0  0  0  0  0999 V2000
   -6.0878    2.4335    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3867    1.6835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3867    0.1833    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.0878   -0.5666    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7887    0.1833    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7887    1.6835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4897   -0.5664    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1906    1.6833    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1906    0.1835    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  8  9  2  0
  6  8  1  0
  7  9  1  0
  7  5  1  0
M  END"""
    indazoleRef = Chem.MolFromMolBlock(indazoleMolblock)
    cycloheptylPyrazole = Chem.MolFromSmiles("c1cc(C2CCCCCC2)[nH]n1")

    # test using refPattern
    refPatt = Chem.MolFromSmarts("a1aan[nH]1")
    rdDepictor.GenerateDepictionMatching2DStructure(cycloheptylPyrazole, indazoleRef,
                                                    refPatt=refPatt)
    self.assertEqual(cycloheptylPyrazole.GetNumConformers(), 1)
    molMatchVect = cycloheptylPyrazole.GetSubstructMatch(refPatt)
    self.assertEqual(len(molMatchVect), refPatt.GetNumAtoms())
    refMatchVect = indazoleRef.GetSubstructMatch(refPatt)
    self.assertEqual(len(refMatchVect), refPatt.GetNumAtoms())
    atomMap = tuple(zip(refMatchVect, molMatchVect))
    msd = 0.0
    for refIdx, molIdx in atomMap:
      msd += (indazoleRef.GetConformer().GetAtomPosition(refIdx) -
              cycloheptylPyrazole.GetConformer().GetAtomPosition(molIdx)).LengthSq()
    msd /= len(molMatchVect)
    self.assertAlmostEqual(msd, 0.0)
    # try with a pattern larger than the reference molecule
    hugePatt = Chem.MolFromSmarts("CCCCCCCCCCCCCCCCCCCCCCCCCCC")
    with self.assertRaises(ValueError):
      rdDepictor.GenerateDepictionMatching2DStructure(cycloheptylPyrazole, indazoleRef,
                                                      refPatt=hugePatt)

    # try with an out of range confId
    with self.assertRaises(ValueError):
      rdDepictor.GenerateDepictionMatching2DStructure(cycloheptylPyrazole, indazoleRef, confId=1,
                                                      refPatt=refPatt)

    # test using atomMap directly
    cycloheptylPyrazole.RemoveAllConformers()
    rdDepictor.GenerateDepictionMatching2DStructure(cycloheptylPyrazole, indazoleRef,
                                                    atomMap=atomMap)
    self.assertEqual(cycloheptylPyrazole.GetNumConformers(), 1)
    msd = 0.0
    for refIdx, molIdx in atomMap:
      msd += (indazoleRef.GetConformer().GetAtomPosition(refIdx) -
              cycloheptylPyrazole.GetConformer().GetAtomPosition(molIdx)).LengthSq()
    msd /= len(atomMap)
    self.assertAlmostEqual(msd, 0.0)

    # try with an atomMap larger than the reference molecule
    atomMapHuge = list(atomMap) + [(0, 0) for i in range(indazoleRef.GetNumAtoms())]
    with self.assertRaises(ValueError):
      rdDepictor.GenerateDepictionMatching2DStructure(cycloheptylPyrazole, indazoleRef,
                                                      atomMap=atomMapHuge)

    # try with an atomMap with out of range indices
    atomMapOutOfRange = list(atomMap) + [(100, 100)]
    with self.assertRaises(ValueError):
      rdDepictor.GenerateDepictionMatching2DStructure(cycloheptylPyrazole, indazoleRef,
                                                      atomMap=atomMapOutOfRange)

    # try with an out of range confId
    with self.assertRaises(ValueError):
      rdDepictor.GenerateDepictionMatching2DStructure(cycloheptylPyrazole, indazoleRef,
                                                      atomMap=atomMap, confId=1)

  def testGenerate2DDepictionAllowRGroups(self):
    templateMolblock = """
     RDKit          2D

  9  9  0  0  0  0  0  0  0  0999 V2000
   -0.8929    1.0942    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8929   -1.9059    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4910    1.0942    0.0000 R1  0  0  0  0  0  0  0  0  0  0  0  0
    1.7051    1.0942    0.0000 R2  0  0  0  0  0  0  0  0  0  0  0  0
   -3.4910   -1.9059    0.0000 R3  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  6  8  1  0
  3  9  1  0
  2  7  1  0
M  RGP  3   7   1   8   2   9   3
M  END"""
    templateRef = Chem.MolFromMolBlock(templateMolblock)
    orthoMeta = Chem.MolFromSmiles("c1ccc(-c2ccc(-c3ccccc3)c(-c3ccccc3)c2)cc1")
    ortho = Chem.MolFromSmiles("c1ccc(-c2ccccc2-c2ccccc2)cc1")
    meta = Chem.MolFromSmiles("c1ccc(-c2cccc(-c3ccccc3)c2)cc1")
    para = Chem.MolFromSmiles("c1ccc(-c2ccc(-c3ccccc3)cc2)cc1")
    biphenyl = Chem.MolFromSmiles("c1ccccc1-c1ccccc1")
    phenyl = Chem.MolFromSmiles("c1ccccc1")

    atomMap = rdDepictor.GenerateDepictionMatching2DStructure(orthoMeta, templateRef)
    self.assertEqual(orthoMeta.GetNumConformers(), 1)

    # test original usage pattern
    for mol in (ortho, meta, biphenyl, phenyl):
      # fails as does not match template
      with self.assertRaises(ValueError):
        rdDepictor.GenerateDepictionMatching2DStructure(mol, templateRef)

      # succeeds with allowRGroups=true
      atomMap = rdDepictor.GenerateDepictionMatching2DStructure(mol, templateRef, allowRGroups=True)
      self.assertEqual(mol.GetNumConformers(), 1)
      msd = 0.0
      for refIdx, molIdx in atomMap:
        msd += (templateRef.GetConformer().GetAtomPosition(refIdx) -
                mol.GetConformer().GetAtomPosition(molIdx)).LengthSq()
      msd /= len(atomMap)
      self.assertAlmostEqual(msd, 0.0)

    for alignOnly in (True, False):
      p = rdDepictor.ConstrainedDepictionParams()
      p.allowRGroups = True
      p.alignOnly = alignOnly
      for mol in (ortho, meta, para, biphenyl, phenyl):
        # fails as does not match template
        with self.assertRaises(ValueError):
          rdDepictor.GenerateDepictionMatching2DStructure(mol, templateRef)

        # succeeds with allowRGroups=true
        atomMap = rdDepictor.GenerateDepictionMatching2DStructure(mol, templateRef, params=p)
        self.assertGreater(len(atomMap), 0)
        self.assertEqual(mol.GetNumConformers(), 1)
        msd = 0.0
        for refIdx, molIdx in atomMap:
          msd += (templateRef.GetConformer().GetAtomPosition(refIdx) -
                  mol.GetConformer().GetAtomPosition(molIdx)).LengthSq()
        msd /= len(atomMap)
        self.assertAlmostEqual(msd, 0.0)
        if not p.alignOnly:
          self.assertEqual(
            atomMap,
            rdDepictor.GenerateDepictionMatching2DStructure(mol, templateRef, allowRGroups=True))

      # test that using a refPattern with R groups and a reference without works
      pyridineRef = Chem.MolFromMolBlock("""
     RDKit          2D

  8  8  0  0  0  0  0  0  0  0999 V2000
    0.0000    1.5469    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3395    0.7734    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3395   -0.7732    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.5469    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3395   -0.7732    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3395    0.7734    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    3.0938    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -3.0938    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  1  7  1  0
  4  8  1  0
M  END""")
      genericRefPatternWithRGroups = Chem.MolFromSmarts("[*:3]a1a([*:1])aa([*:2])aa1")

      for numExpectedMatches, mol in ((8, orthoMeta), (7, ortho), (7, meta), (8, para),
                                      (7, biphenyl), (6, phenyl)):
        atomMap = rdDepictor.GenerateDepictionMatching2DStructure(mol, pyridineRef, -1,
                                                                  genericRefPatternWithRGroups, p)
        self.assertEqual(len(atomMap), numExpectedMatches)
        self.assertEqual(mol.GetNumConformers(), 1)
        msd = 0.0
        for refIdx, molIdx in atomMap:
          msd += (pyridineRef.GetConformer().GetAtomPosition(refIdx) -
                  mol.GetConformer().GetAtomPosition(molIdx)).LengthSq()
        msd /= len(atomMap)
        self.assertLess(msd, 5.e-3 if alignOnly else 1.e-4)

  def testNormalizeStraighten(self):
    noradrenalineMJ = Chem.MolFromMolBlock("""
  MJ201100                      

 12 12  0  0  1  0  0  0  0  0999 V2000
    2.2687    1.0716    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.4437    1.0716    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0312    0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4437   -0.3572    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2062    0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2062   -0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0312   -0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4437   -1.0716    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4437    0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2687    0.3572    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0312    1.0716    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2062    1.0716    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  3  2  1  0  0  0  0
  3  4  1  6  0  0  0
  3  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  7  1  0  0  0  0
  7  8  1  0  0  0  0
  7  9  2  0  0  0  0
  9 10  1  0  0  0  0
  9 11  1  0  0  0  0
 11 12  2  0  0  0  0
  5 12  1  0  0  0  0
M  END)""")

    noradrenalineMJCopy = Chem.Mol(noradrenalineMJ)
    conformer0 = noradrenalineMJCopy.GetConformer(0)
    conformer1 = Chem.Conformer(conformer0)
    noradrenalineMJCopy.AddConformer(conformer1, True)
    conformer1 = noradrenalineMJCopy.GetConformer(1)
    self.assertLess(rdMolAlign.CalcRMS(noradrenalineMJ, noradrenalineMJCopy, 0, 0), 1.e-5)
    self.assertLess(rdMolAlign.CalcRMS(noradrenalineMJ, noradrenalineMJCopy, 0, 1), 1.e-5)
    scalingFactor = rdDepictor.NormalizeDepiction(noradrenalineMJCopy, 1)
    self.assertLess(rdMolAlign.CalcRMS(noradrenalineMJ, noradrenalineMJCopy, 0, 0), 1.e-5)
    self.assertGreater(rdMolAlign.CalcRMS(noradrenalineMJ, noradrenalineMJCopy, 0, 1), 1.e-5)
    self.assertAlmostEqual(scalingFactor, 1.875, 3)
    conformer2 = Chem.Conformer(conformer1)
    noradrenalineMJCopy.AddConformer(conformer2, True)
    conformer2 = noradrenalineMJCopy.GetConformer(2)
    bond10_11Conf0 = conformer0.GetAtomPosition(11) - conformer0.GetAtomPosition(10)
    self.assertAlmostEqual(bond10_11Conf0.x, 0.825, 3)
    self.assertAlmostEqual(bond10_11Conf0.y, 0., 3)
    bond10_11Conf1 = conformer1.GetAtomPosition(11) - conformer1.GetAtomPosition(10)
    self.assertAlmostEqual(bond10_11Conf1.x, 1.513, 3)
    self.assertAlmostEqual(bond10_11Conf1.y, -0.321, 3)
    rdDepictor.StraightenDepiction(noradrenalineMJCopy, 1)
    bond10_11Conf1 = conformer1.GetAtomPosition(11) - conformer1.GetAtomPosition(10)
    self.assertAlmostEqual(bond10_11Conf1.x, 1.340, 3)
    self.assertAlmostEqual(bond10_11Conf1.y, -0.773, 3)
    bond4_11Conf1 = conformer1.GetAtomPosition(11) - conformer1.GetAtomPosition(4)
    self.assertAlmostEqual(bond4_11Conf1.x, 0., 3)
    self.assertAlmostEqual(bond4_11Conf1.y, 1.547, 3)
    rdDepictor.StraightenDepiction(noradrenalineMJCopy, 2, True)
    bond10_11Conf2 = conformer2.GetAtomPosition(11) - conformer2.GetAtomPosition(10)
    self.assertAlmostEqual(bond10_11Conf2.x, 1.547, 3)
    self.assertAlmostEqual(bond10_11Conf2.y, 0.0, 3)
    bond4_11Conf2 = conformer2.GetAtomPosition(11) - conformer2.GetAtomPosition(4)
    self.assertAlmostEqual(bond4_11Conf2.x, -0.773, 3)
    self.assertAlmostEqual(bond4_11Conf2.y, 1.339, 3)

    noradrenalineMJCopy = Chem.Mol(noradrenalineMJ)
    conformer0 = noradrenalineMJCopy.GetConformer(0)
    conformer1 = Chem.Conformer(conformer0)
    noradrenalineMJCopy.AddConformer(conformer1, True)
    conformer1 = noradrenalineMJCopy.GetConformer(1)
    scalingFactor = rdDepictor.NormalizeDepiction(noradrenalineMJCopy, 1, -1)
    self.assertLess(rdMolAlign.CalcRMS(noradrenalineMJ, noradrenalineMJCopy, 0, 0), 1.e-5)
    self.assertGreater(rdMolAlign.CalcRMS(noradrenalineMJ, noradrenalineMJCopy, 0, 1), 1.e-5)
    self.assertAlmostEqual(scalingFactor, 1.875, 3)
    conformer2 = Chem.Conformer(conformer1)
    noradrenalineMJCopy.AddConformer(conformer2, True)
    conformer2 = noradrenalineMJCopy.GetConformer(2)
    bond10_11Conf0 = conformer0.GetAtomPosition(11) - conformer0.GetAtomPosition(10)
    self.assertAlmostEqual(bond10_11Conf0.x, 0.825, 3)
    self.assertAlmostEqual(bond10_11Conf0.y, 0., 3)
    bond10_11Conf1 = conformer1.GetAtomPosition(11) - conformer1.GetAtomPosition(10)
    self.assertAlmostEqual(bond10_11Conf1.x, 0.321, 3)
    self.assertAlmostEqual(bond10_11Conf1.y, 1.513, 3)
    rdDepictor.StraightenDepiction(noradrenalineMJCopy, 1)
    bond10_11Conf1 = conformer1.GetAtomPosition(11) - conformer1.GetAtomPosition(10)
    self.assertAlmostEqual(bond10_11Conf1.x, 0.0, 3)
    self.assertAlmostEqual(bond10_11Conf1.y, 1.547, 3)
    rdDepictor.StraightenDepiction(noradrenalineMJCopy, 2, True)
    bond10_11Conf2 = conformer2.GetAtomPosition(11) - conformer2.GetAtomPosition(10)
    self.assertAlmostEqual(bond10_11Conf2.x, bond10_11Conf1.x, 3)
    self.assertAlmostEqual(bond10_11Conf2.y, bond10_11Conf1.y, 3)

    noradrenalineMJCopy = Chem.Mol(noradrenalineMJ)
    conformer0 = noradrenalineMJCopy.GetConformer(0)
    conformer1 = Chem.Conformer(conformer0)
    noradrenalineMJCopy.AddConformer(conformer1, True)
    conformer1 = noradrenalineMJCopy.GetConformer(1)
    scalingFactor = rdDepictor.NormalizeDepiction(noradrenalineMJCopy, 1, 0, 3.0)
    self.assertLess(rdMolAlign.CalcRMS(noradrenalineMJ, noradrenalineMJCopy, 0, 0), 1.e-5)
    self.assertGreater(rdMolAlign.CalcRMS(noradrenalineMJ, noradrenalineMJCopy, 0, 1), 1.e-5)
    self.assertAlmostEqual(scalingFactor, 3.0, 3)
    conformer2 = Chem.Conformer(conformer1)
    noradrenalineMJCopy.AddConformer(conformer2, True)
    conformer2 = noradrenalineMJCopy.GetConformer(2)
    conformer3 = Chem.Conformer(conformer1)
    noradrenalineMJCopy.AddConformer(conformer3, True)
    conformer3 = noradrenalineMJCopy.GetConformer(3)
    bond10_11Conf0 = conformer0.GetAtomPosition(11) - conformer0.GetAtomPosition(10)
    self.assertAlmostEqual(bond10_11Conf0.x, 0.825, 3)
    self.assertAlmostEqual(bond10_11Conf0.y, 0., 3)
    bond10_11Conf1 = conformer1.GetAtomPosition(11) - conformer1.GetAtomPosition(10)
    self.assertAlmostEqual(bond10_11Conf1.x, 2.475, 3)
    self.assertAlmostEqual(bond10_11Conf1.y, 0., 3)
    rdDepictor.StraightenDepiction(noradrenalineMJCopy, 1)
    bond10_11Conf1 = conformer1.GetAtomPosition(11) - conformer1.GetAtomPosition(10)
    self.assertAlmostEqual(bond10_11Conf1.x, 2.143, 3)
    self.assertAlmostEqual(bond10_11Conf1.y, -1.237, 3)
    bond4_11Conf1 = conformer1.GetAtomPosition(11) - conformer1.GetAtomPosition(4)
    self.assertAlmostEqual(bond4_11Conf1.x, 0., 3)
    self.assertAlmostEqual(bond4_11Conf1.y, 2.475, 3)
    rdDepictor.StraightenDepiction(noradrenalineMJCopy, 2, True)
    bond10_11Conf2 = conformer2.GetAtomPosition(11) - conformer2.GetAtomPosition(10)
    bond10_11Conf3 = conformer3.GetAtomPosition(11) - conformer3.GetAtomPosition(10)
    self.assertAlmostEqual(bond10_11Conf2.x, bond10_11Conf3.x, 3)
    self.assertAlmostEqual(bond10_11Conf2.y, bond10_11Conf3.y, 3)
    bond4_11Conf2 = conformer2.GetAtomPosition(11) - conformer2.GetAtomPosition(4)
    bond4_11Conf3 = conformer3.GetAtomPosition(11) - conformer3.GetAtomPosition(4)
    self.assertAlmostEqual(bond4_11Conf2.x, bond4_11Conf3.x, 3)
    self.assertAlmostEqual(bond4_11Conf2.y, bond4_11Conf3.y, 3)

  @unittest.skipIf(not rdDepictor.IsCoordGenSupportAvailable(), "CoordGen not available, skipping")
  def testUsingCoordGenCtxtMgr(self):
    default_status = rdDepictor.GetPreferCoordGen()

    # This is the default; we shouldn't have changed it
    self.assertEqual(default_status, False)

    with rdDepictor.UsingCoordGen(True):
      current_status = rdDepictor.GetPreferCoordGen()
      self.assertEqual(current_status, True)

    current_status = rdDepictor.GetPreferCoordGen()
    self.assertEqual(current_status, False)

    rdDepictor.SetPreferCoordGen(True)

    with rdDepictor.UsingCoordGen(False):
      current_status = rdDepictor.GetPreferCoordGen()
      self.assertEqual(current_status, False)

    current_status = rdDepictor.GetPreferCoordGen()
    self.assertEqual(current_status, True)

    rdDepictor.SetPreferCoordGen(default_status)

  def molMatchesTemplate(self, mol, template):
    """
    Determines if the shape/layout of the template and mol are the same. It
    is ok if the mol and template are not centered at the same place, or if
    the mol and template have different orientations.
    """
    match = mol.GetSubstructMatch(template)
    if not match or len(match) != template.GetNumAtoms():
      return False

    # get positions of atoms with centroid at origin, it is ok if the
    # template or mol is not centered
    template_match_positions = [
      mol.GetConformer().GetPositions()[mol_at_idx] for mol_at_idx in match
    ]
    template_match_center = sum(template_match_positions) / len(template_match_positions)
    mol_positions = [p - template_match_center for p in mol.GetConformer().GetPositions()]

    template_center = sum(template.GetConformer().GetPositions()) / template.GetNumAtoms()
    template_positions = [p - template_center for p in template.GetConformer().GetPositions()]

    # the mol may match the template but be slightly rotated about the centroid
    # or reflected across the x or y axis
    rotations = [[], [], [], []]
    for template_idx, idx in enumerate(match):
      v1 = mol_positions[idx]

      # no reflection
      v2 = template_positions[template_idx]
      val = round(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), 4)
      rotations[0].append(np.arccos(val))

      # reflect across x-axis
      v2[0] = v2[0] * -1
      val = round(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), 4)
      rotations[1].append(np.arccos(val))

      # reflect across y-axis
      v2[0] = v2[0] * -1
      v2[1] = v2[1] * -1
      val = round(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), 4)
      rotations[2].append(np.arccos(val))

      # reflect across y-axis and x-acis
      v2[0] = v2[0] * -1
      # v2[1] = v2[1] * -1
      val = round(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), 4)
      rotations[3].append(np.arccos(val))

    # if all the rotations are similar, then the shape is the same
    return np.any([np.allclose(r, r[0], atol=.05) for r in rotations])

  def assertMolMatchesCoordMap(self, mol, coord_map):
    for aid, expected_position in coord_map.items():
      actual_position = mol.GetConformer().GetAtomPosition(aid)
      self.assertAlmostEqual(actual_position.x, expected_position.x)
      self.assertAlmostEqual(actual_position.y, actual_position.y)

  def testUseMultipleTemplates(self):
    with rdDepictor.UsingCoordGen(False):
      # templates that will be linked together
      template1 = Chem.MolFromSmiles(
        "C1=CCCC2CCCCC2CCCCC2CCCC(CCCCCCC1)C2 |(-0.04,3.43,;-0.04,1.93,;-1.34,1.18,;-2.64,1.93,;-3.94,1.18,;-5.24,1.93,;-6.54,1.18,;-6.54,-0.32,;-5.24,-1.07,;-3.94,-0.32,;-2.64,-1.07,;-2.64,-2.57,;-1.34,-3.32,;-0.04,-2.57,;1.26,-3.32,;1.26,-4.82,;2.56,-5.56,;3.86,-4.82,;3.86,-3.32,;5.16,-2.57,;5.16,-1.07,;3.86,-0.32,;3.86,1.18,;2.56,1.93,;2.56,3.43,;1.26,4.18,;2.56,-2.57,)|"
      )
      template2 = Chem.MolFromSmiles(
        "C1CCC2C(C1)C1CCN2NN1 |(-2.94,-0.77,;-2.94,0.77,;-1.6,1.54,;-0.27,0.77,;-0.27,-0.77,;-1.6,-1.54,;1.06,-1.54,;2.4,-0.77,;2.4,0.77,;1.06,1.54,;1.33,0.51,;1.33,-0.51,)|"
      )
      template3 = Chem.MolFromSmiles(
        "C1C2CC3CC1CC3C2 |(-7.01,3.13,;-7.71,4.35,;-7.01,5.56,;-5.61,5.56,;-4.91,4.35,;-5.61,3.13,;-4.28,3.57,;-4.28,5.13,;-6.34,4.05,)|"
      )

      # example with 2 templates linked together
      two_linked_templates = Chem.MolFromSmiles(
        "NC(CCC1CCC2C(C1)C1CCN2NN1)CC(=O)CCC1=CCCC2CCCCC2CCCCC2CCCC(CCCCCCC1)C2")
      rdDepictor.Compute2DCoords(two_linked_templates, useRingTemplates=False)
      assert not self.molMatchesTemplate(two_linked_templates, template1)
      assert not self.molMatchesTemplate(two_linked_templates, template2)

      rdDepictor.Compute2DCoords(two_linked_templates, useRingTemplates=True)
      assert self.molMatchesTemplate(two_linked_templates, template1)
      assert self.molMatchesTemplate(two_linked_templates, template2)

      # example with 3 templates linked together
      three_linked_templates = Chem.MolFromSmiles(
        "NC(CCC1CCC2C(C1)C1CC(CCC(=O)CC(N)CC~C3C4CC5CC3CC5C4)N2NN1)CC(=O)CCC1=CCCC2CCCCC2CCCCC2CCCC(CCCCCCC1)C2"
      )
      rdDepictor.Compute2DCoords(three_linked_templates, useRingTemplates=False)
      assert not self.molMatchesTemplate(three_linked_templates, template1)
      assert not self.molMatchesTemplate(three_linked_templates, template2)
      assert not self.molMatchesTemplate(three_linked_templates, template3)

      rdDepictor.Compute2DCoords(three_linked_templates, useRingTemplates=True)
      assert self.molMatchesTemplate(three_linked_templates, template1)
      assert self.molMatchesTemplate(three_linked_templates, template2)
      assert self.molMatchesTemplate(three_linked_templates, template3)

  def testUseTemplateAndCoordMap(self):
    with rdDepictor.UsingCoordGen(False):
      template1 = Chem.MolFromSmiles(
        "C1=CCCC2CCCCC2CCCCC2CCCC(CCCCCCC1)C2 |(-0.04,3.43,;-0.04,1.93,;-1.34,1.18,;-2.64,1.93,;-3.94,1.18,;-5.24,1.93,;-6.54,1.18,;-6.54,-0.32,;-5.24,-1.07,;-3.94,-0.32,;-2.64,-1.07,;-2.64,-2.57,;-1.34,-3.32,;-0.04,-2.57,;1.26,-3.32,;1.26,-4.82,;2.56,-5.56,;3.86,-4.82,;3.86,-3.32,;5.16,-2.57,;5.16,-1.07,;3.86,-0.32,;3.86,1.18,;2.56,1.93,;2.56,3.43,;1.26,4.18,;2.56,-2.57,)|"
      )
      template2 = Chem.MolFromSmiles(
        "C1CCC2C(C1)C1CCN2NN1 |(-2.94,-0.77,;-2.94,0.77,;-1.6,1.54,;-0.27,0.77,;-0.27,-0.77,;-1.6,-1.54,;1.06,-1.54,;2.4,-0.77,;2.4,0.77,;1.06,1.54,;1.33,0.51,;1.33,-0.51,)|"
      )
      two_linked_templates = Chem.MolFromSmiles(
        "NC(CCC1CCC2C(C1)C1CCN2NN1)CC(=O)CCC1=CCCC2CCCCC2CCCCC2CCCC(CCCCCCC1)C2")

      # when a coord map doesn't contain any part of a ring system, ring system
      # templates should still be adhered to
      linker_coord_map = {
        16: Geometry.Point2D(1.5, 0),
        17: Geometry.Point2D(1.5, 1.5),
        19: Geometry.Point2D(0, 1.5)
      }
      rdDepictor.Compute2DCoords(two_linked_templates, coordMap=linker_coord_map,
                                 useRingTemplates=True)
      self.assertMolMatchesCoordMap(two_linked_templates, linker_coord_map)
      assert self.molMatchesTemplate(two_linked_templates, template1)
      assert self.molMatchesTemplate(two_linked_templates, template2)

      # when a coord map contains a partial ring system, ring system templates
      # should not be used because they could be distorted by the user-provided
      # templates
      ring_system_coord_map = {
        31: Geometry.Point2D(1.5, 0),
        32: Geometry.Point2D(1.5, 1.5),
        33: Geometry.Point2D(0, 1.5)
      }
      rdDepictor.Compute2DCoords(two_linked_templates, coordMap=ring_system_coord_map,
                                 useRingTemplates=True)
      self.assertMolMatchesCoordMap(two_linked_templates, ring_system_coord_map)
      # atoms 10, 11, and 13 are in this template so the ring template should not be used
      assert not self.molMatchesTemplate(two_linked_templates, template1)
      assert self.molMatchesTemplate(two_linked_templates, template2)

      # when a coord map contains a single atom, even if it is a part of a ring
      # system, ring system templates should be used and the coord map should be
      # followed
      single_atom_coord_map = {10: Geometry.Point2D(0, 0)}
      rdDepictor.Compute2DCoords(two_linked_templates, coordMap=single_atom_coord_map,
                                 useRingTemplates=True)
      self.assertMolMatchesCoordMap(two_linked_templates, single_atom_coord_map)
      assert self.molMatchesTemplate(two_linked_templates, template1)
      assert self.molMatchesTemplate(two_linked_templates, template2)

  def testSetRingSystemTemplates(self):
    with rdDepictor.UsingCoordGen(False):
      mol = Chem.MolFromSmiles("C1CC2CCOC3OC4CCC(C1)C23OO4")
      default_template = Chem.MolFromSmiles(
        "C1CC2CCOC3OC4CCC(C1)C23OO4 |(3.53,-1.22,;3.53,0.3,;2.21,1.06,;2.21,2.59,;0.89,3.35,;-0.43,2.59,;-0.43,1.06,;-1.9,0.65,;-2.47,-0.76,;-1.71,-2.08,;-0.2,-2.29,;0.89,-1.22,;2.21,-1.99,;0.89,0.3,;0.12,-0.83,;-1.19,-1.25,)|"
      )
      user_provided_template = Chem.MolFromSmiles(
        "C1CC2CCOC3OC4CCC(C1)C23OO4 |(-0.5537,-3.1595,;-1.6057,-2.003,;-1.4262,-0.4072,;-2.9804,0.0271,;-3.5191,1.502,;-2.2028,2.3562,;-0.6818,1.8511,;1.0592,1.4391,;2.6123,1.8366,;3.5191,0.5341,;2.6067,-0.7521,;1.0061,-0.773,;0.7888,-2.3546,;-0.0405,0.5251,;0.4049,2.3,;1.7604,3.1594,)|"
      )

      # default templates are loaded automatically
      rdDepictor.Compute2DCoords(mol, useRingTemplates=True)
      assert self.molMatchesTemplate(mol, default_template)

      # set to user-provided template, this will delete default templates
      fpath = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Depictor', 'test_data',
                           'ring_system_templates.smi')
      rdDepictor.SetRingSystemTemplates(fpath)
      rdDepictor.Compute2DCoords(mol, useRingTemplates=True)
      assert self.molMatchesTemplate(mol, user_provided_template)

      # set back to default ring system templates
      rdDepictor.LoadDefaultRingSystemTemplates()
      rdDepictor.Compute2DCoords(mol, useRingTemplates=True)
      assert self.molMatchesTemplate(mol, default_template)

  def testSetBadRingSystemTemplates(self):
    with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
      tmp_file.write(b"invalidsmiles")
      tmp_file.seek(0)
      with self.assertRaises(ValueError):
        rdDepictor.SetRingSystemTemplates(tmp_file.name)

    with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
      # not a ring system
      tmp_file.write(b"C |(-0.5537,-3.1595)|")
      tmp_file.seek(0)
      with self.assertRaises(ValueError):
        rdDepictor.SetRingSystemTemplates(tmp_file.name)

    with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
      # no coordinates
      tmp_file.write(b"C1CCCCC1")
      tmp_file.seek(0)
      with self.assertRaises(ValueError):
        rdDepictor.SetRingSystemTemplates(tmp_file.name)

    with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
      # bridged ring system
      tmp_file.write(b"c1ccccc1-c1ccccc1")
      tmp_file.seek(0)
      with self.assertRaises(ValueError):
        rdDepictor.SetRingSystemTemplates(tmp_file.name)

    # set back to default ring system templates
    rdDepictor.LoadDefaultRingSystemTemplates()

  def testAddRingSystemTemplates(self):
    with rdDepictor.UsingCoordGen(False):
      # this is what is in the default ring system templates
      mol = Chem.MolFromSmiles("C1CCCN2CCCC(CC1)C2")
      default_template = Chem.MolFromSmiles(
        "C1CCCN2CCCC(CC1)C2 |(2.64,1.53,;3.1,0.11,;2.52,-1.28,;1.17,-1.94,;-0.28,-1.57,;-1.7,-2.03,;-2.82,-1.03,;-2.51,0.44,;-1.08,0.9,;-0.13,2.06,;1.35,2.31,;0.03,-0.1,)|"
      )
      new_template_smi = "C1CCCN2CCCC(CC1)C2 |(3.2585,-0.7733,;3.2585,0.7733,;2.225,1.9568,;0.6347,2.17,;-0.8068,1.4216,;-2.4372,1.3464,;-3.2585,0,;-2.4372,-1.3464,;-0.8068,-1.4216,;0.6347,-2.17,;2.225,-1.9569,;-0.6701,0,)|"
      new_template = Chem.MolFromSmiles(new_template_smi)

      rdDepictor.Compute2DCoords(mol, useRingTemplates=True)
      assert self.molMatchesTemplate(mol, default_template)

      with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
        tmp_file.write(new_template_smi.encode('utf-8'))
        tmp_file.seek(0)
        rdDepictor.AddRingSystemTemplates(tmp_file.name)

      # now the new template should be used
      rdDepictor.Compute2DCoords(mol, useRingTemplates=True)
      assert self.molMatchesTemplate(mol, new_template)

    # set back to default ring system templates
    rdDepictor.LoadDefaultRingSystemTemplates()

  def testGenerateAlignedCoordsAcceptFailure(self):
    template_ref_molblock = """
     RDKit          2D

  9  9  0  0  0  0  0  0  0  0999 V2000
   -0.8929    1.0942    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8929   -1.9059    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4910    1.0942    0.0000 R1  0  0  0  0  0  0  0  0  0  0  0  0
    1.7051    1.0942    0.0000 R2  0  0  0  0  0  0  0  0  0  0  0  0
   -3.4910   -1.9059    0.0000 R3  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  6  8  1  0
  3  9  1  0
  2  7  1  0
M  RGP  3   7   1   8   2   9   3
M  END
"""
    template_ref = Chem.MolFromMolBlock(template_ref_molblock)
    mol_molblock = """
     RDKit          2D

  9  9  0  0  0  0  0  0  0  0999 V2000
   -0.8929    1.0942    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8929   -1.9059    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4910    1.0942    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    1.7051    1.0942    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -3.4910   -1.9059    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  6  8  1  0
  3  9  1  0
  2  7  1  0
M  END
"""
    mol = Chem.MolFromMolBlock(mol_molblock)
    self.assertRaises(
      ValueError, lambda: rdDepictor.GenerateDepictionMatching2DStructure(
        mol, template_ref, -1, None, False, False, True))
    self.assertEqual(Chem.MolToMolBlock(mol), mol_molblock)
    self.assertEqual(
      len(
        rdDepictor.GenerateDepictionMatching2DStructure(mol, template_ref, -1, None, True, False,
                                                        True)), 0)
    self.assertNotEqual(Chem.MolToMolBlock(mol), mol_molblock)
    mol.RemoveAllConformers()
    p = rdDepictor.ConstrainedDepictionParams()
    p.allowRGroups = True
    p.acceptFailure = False
    self.assertEqual(mol.GetNumConformers(), 0)
    self.assertRaises(
      ValueError,
      lambda: rdDepictor.GenerateDepictionMatching2DStructure(mol, template_ref, -1, None, p))
    self.assertEqual(mol.GetNumConformers(), 0)
    mol.RemoveAllConformers()
    p = rdDepictor.ConstrainedDepictionParams()
    p.allowRGroups = True
    p.acceptFailure = True
    self.assertEqual(mol.GetNumConformers(), 0)
    self.assertEqual(
      len(rdDepictor.GenerateDepictionMatching2DStructure(mol, template_ref, -1, None, p)), 0)
    self.assertEqual(mol.GetNumConformers(), 1)

  def testGenerateAlignedCoordsAlignOnly(self):
    template_ref_molblock = """
     RDKit          2D

  6  6  0  0  0  0  0  0  0  0999 V2000
  -13.7477    6.3036    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  -13.7477    4.7567    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -12.6540    3.6628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -13.7477    2.5691    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -14.8414    3.6628    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  -11.1071    3.6628    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  4  5  1  0
  2  5  1  0
  3  6  1  0
M  RGP  2   1   1   6   2
M  END
"""
    template_ref = Chem.MolFromMolBlock(template_ref_molblock)
    mol_molblock = """
     RDKit          2D

 18 22  0  0  0  0  0  0  0  0999 V2000
    4.3922   -1.5699    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9211   -2.0479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5995   -0.5349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3731    0.8046    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8441    1.2825    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0704   -0.0568    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8666    0.7748    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.7736   -0.3197    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7749   -1.8666    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7718   -1.8679    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7731   -0.3208    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8679    0.7718    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0718   -0.0598    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3933   -1.5729    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9222   -2.0509    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6008   -0.5379    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3744    0.8016    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8454    1.2795    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  9 10  1  0
 11 10  1  0
 11  8  1  0
  8  9  1  0
  4  5  1  0
  6  5  1  0
  7  6  1  0
  3  4  1  0
  3  7  1  0
  1  6  1  0
  2  3  1  0
  2  1  1  0
 17 18  1  0
 13 18  1  0
 12 13  1  0
 16 17  1  0
 16 12  1  0
 14 13  1  0
 15 16  1  0
 15 14  1  0
 12 11  1  0
  8  7  1  0
M  END
"""
    mol = Chem.MolFromMolBlock(mol_molblock)
    bondLength11_12 = rdMolTransforms.GetBondLength(mol.GetConformer(), 11, 12)
    bondLength5_6 = rdMolTransforms.GetBondLength(mol.GetConformer(), 5, 6)
    self.assertLess(abs(bondLength11_12 - bondLength5_6), 1.e-4)
    self.assertGreater(bondLength11_12, 2.3)
    #for alignOnly in (False, True):
    for alignOnly in (True, ):
      mol = Chem.MolFromMolBlock(mol_molblock)
      p = rdDepictor.ConstrainedDepictionParams()
      p.allowRGroups = True
      p.alignOnly = alignOnly
      res = rdDepictor.GenerateDepictionMatching2DStructure(mol, template_ref, params=p)
      expectedMolIndices = [11, 10, 7, 8, 9, 6]
      self.assertTrue(
        all(templateRefAtomIdx == i and molAtomIdx == expectedMolIndices[i]
            for i, (templateRefAtomIdx, molAtomIdx) in enumerate(res)))
      self.assertEqual(Chem.MolToSmiles(mol), "C1CC2CCC1N2C1CNC1N1C2CCC1CC2")
      self.assertTrue(
        all((mol.GetConformer().GetAtomPosition(molAtomIdx) -
             template_ref.GetConformer().GetAtomPosition(templateRefAtomIdx)).LengthSq() < 1.e-4
            for templateRefAtomIdx, molAtomIdx in res))
      bondLengthAli11_12 = rdMolTransforms.GetBondLength(mol.GetConformer(), 11, 12)
      bondLengthAli5_6 = rdMolTransforms.GetBondLength(mol.GetConformer(), 5, 6)
      self.assertLess(abs(bondLengthAli11_12 - bondLengthAli5_6), 1.e-4)
      if alignOnly:
        self.assertGreater(bondLengthAli11_12, 2.3)
      else:
        self.assertLess(bondLengthAli11_12, 1.6)

  def testRGroupMatchHeavyHydroNoneCharged(self):
    templateRef = Chem.MolFromMolBlock("""
  MJ201100                      

  7  7  0  0  0  0  0  0  0  0999 V2000
   -0.5804    1.2045    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2948    0.7920    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2948   -0.0330    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5804   -0.4455    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1340   -0.0330    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
    0.1340    0.7920    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8485   -0.4455    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  6  1  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  5  7  1  0  0  0  0
M  RGP  1   7   1
M  END
""")
    p = rdDepictor.ConstrainedDepictionParams()
    p.allowRGroups = True
    for alignOnly in (True, False):
      p.alignOnly = alignOnly
      mol = Chem.MolFromSmiles("Cc1ccccc1")
      self.assertEqual(mol.GetNumAtoms(), 7)
      self.assertEqual(
        len(rdDepictor.GenerateDepictionMatching2DStructure(mol, templateRef, params=p)), 7)
      mol = Chem.MolFromSmiles("c1ccccc1")
      self.assertEqual(mol.GetNumAtoms(), 6)
      self.assertEqual(
        len(rdDepictor.GenerateDepictionMatching2DStructure(mol, templateRef, params=p)), 6)
      smilesParams = Chem.SmilesParserParams()
      smilesParams.removeHs = False
      mol = Chem.MolFromSmiles("[H]c1ccccc1", smilesParams)
      self.assertEqual(mol.GetNumAtoms(), 7)
      self.assertEqual(
        len(rdDepictor.GenerateDepictionMatching2DStructure(mol, templateRef, params=p)), 7)
      mol = Chem.MolFromSmiles("n1ccccc1")
      self.assertEqual(mol.GetNumAtoms(), 6)
      self.assertEqual(
        len(rdDepictor.GenerateDepictionMatching2DStructure(mol, templateRef, params=p)), 6)
      mol = Chem.MolFromSmiles("C[n+]1ccccc1")
      self.assertEqual(mol.GetNumAtoms(), 7)
      self.assertEqual(
        len(rdDepictor.GenerateDepictionMatching2DStructure(mol, templateRef, params=p)), 7)


if __name__ == '__main__':
  unittest.main()
