#
#  Copyright (C) 2003-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import io
import os
import unittest

from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import FragmentCatalog, BuildFragmentCatalog
from rdkit.six.moves import cPickle

doLong = False


class TestCase(unittest.TestCase):

  def setUp(self):
    self.smiList = ["S(SC1=NC2=CC=CC=C2S1)C3=NC4=C(S3)C=CC=C4", "CC1=CC(=O)C=CC1=O",
                    "OC1=C(Cl)C=C(C=C1[N+]([O-])=O)[N+]([O-])=O", "[O-][N+](=O)C1=CNC(=N)S1",
                    "NC1=CC2=C(C=C1)C(=O)C3=C(C=CC=C3)C2=O",
                    "OC(=O)C1=C(C=CC=C1)C2=C3C=CC(=O)C(=C3OC4=C2C=CC(=C4Br)O)Br",
                    "CN(C)C1=C(Cl)C(=O)C2=C(C=CC=C2)C1=O",
                    "CC1=C(C2=C(C=C1)C(=O)C3=CC=CC=C3C2=O)[N+]([O-])=O", "CC(=NO)C(C)=NO"]
    self.smiList2 = ['OCCC',
                     'CCC',
                     'C=CC',
                     'OC=CC',
                     'CC(O)C',
                     'C=C(O)C',
                     'OCCCC',
                     'CC(O)CC',
                     'C=CCC',
                     'CC=CC',
                     'OC=CCC',
                     'CC=C(O)C',
                     'OCC=CC',
                     'C=C(O)CC',
                     'C=CC(O)C',
                     'C=CCCO', ]
    self.list2Acts = [1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1]
    self.list2Obls = [(0, 1, 2), (1, 3), (1, 4, 5), (1, 6, 7), (0, 8), (0, 6, 9), (0, 1, 2, 3, 10),
                      (0, 1, 2, 8, 11), (1, 3, 4, 5, 12), (1, 4, 5, 13), (1, 3, 6, 7, 14),
                      (0, 1, 6, 7, 9, 15)]

    ffile = os.path.join(RDConfig.RDDataDir, 'FunctionalGroups.txt')
    self.catParams = FragmentCatalog.FragCatParams(1, 6, ffile)
    self.fragCat = FragmentCatalog.FragCatalog(self.catParams)
    self.fgen = FragmentCatalog.FragCatGenerator()

  def _fillCat(self, smilList):
    for smi in self.smiList2:
      mol = Chem.MolFromSmiles(smi)
      self.fgen.AddFragsFromMol(mol, self.fragCat)

  def _testBits(self, fragCat):
    fpgen = FragmentCatalog.FragFPGenerator()
    obits = [3, 2, 3, 3, 2, 3, 5, 5, 5, 4, 5, 6]
    obls = self.list2Obls
    suppl = Chem.SmilesMolSupplierFromText('\n'.join(self.smiList2), ',', 0, -1, 0)
    i = 0
    for mol in suppl:
      fp = fpgen.GetFPForMol(mol, fragCat)
      if i < len(obits):
        smi = Chem.MolToSmiles(mol)
        self.assertEqual(fp.GetNumOnBits(), obits[i], msg='%s: %s' % (smi, str(fp.GetOnBits())))
      obl = fp.GetOnBits()
      if i < len(obls):
        self.assertEqual(tuple(obl), obls[i], msg='%s: %s' % (smi, obl))
      i += 1

  def test1CatGen(self):
    self._fillCat(self.smiList2)
    self.assertEqual(self.fragCat.GetNumEntries(), 21)
    self.assertEqual(self.fragCat.GetFPLength(), 21)
    self._testBits(self.fragCat)

  def test2CatStringPickle(self):
    self._fillCat(self.smiList2)

    # test non-binary pickle:
    cat2 = cPickle.loads(cPickle.dumps(self.fragCat))
    self.assertEqual(cat2.GetNumEntries(), 21)
    self.assertEqual(cat2.GetFPLength(), 21)
    self._testBits(cat2)

    # test binary pickle:
    cat2 = cPickle.loads(cPickle.dumps(self.fragCat, 1))
    self.assertEqual(cat2.GetNumEntries(), 21)
    self.assertEqual(cat2.GetFPLength(), 21)
    self._testBits(cat2)

  def test3CatFilePickle(self):
    with open(os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', 'simple_catalog.pkl'),
              'r') as pklTFile:
      buf = pklTFile.read().replace('\r\n', '\n').encode('utf-8')
      pklTFile.close()
    with io.BytesIO(buf) as pklFile:
      cat = cPickle.load(pklFile, encoding='bytes')
    self.assertEqual(cat.GetNumEntries(), 21)
    self.assertEqual(cat.GetFPLength(), 21)
    self._testBits(cat)

  def test4CatGuts(self):
    self._fillCat(self.smiList2)
    self.assertEqual(self.fragCat.GetNumEntries(), 21)
    self.assertEqual(self.fragCat.GetFPLength(), 21)
    #
    # FIX: (Issue 162)
    # bits like 11 and 15 are questionable here because the underlying
    #  fragments are symmetrical, so they can generate one of two
    #  text representations (i.e. there is nothing to distinguish
    #  between 'CC<-O>CC' and 'CCC<-O>C').
    # This ought to eventually be cleaned up.
    descrs = [(0, 'C<-O>C', 1, (34, )),
              (1, 'CC', 1, ()),
              (2, 'C<-O>CC', 2, (34, )),
              (3, 'CCC', 2, ()),
              (4, 'C=C', 1, ()),
              (5, 'C=CC', 2, ()),
              (6, 'C<-O>=C', 1, (34, )),
              (7, 'C<-O>=CC', 2, (34, )),
              (8, 'CC<-O>C', 2, (34, )),
              (9, 'C=C<-O>C', 2, (34, )),
              (10, 'C<-O>CCC', 3, (34, )),
              (11, 'CC<-O>CC', 3, (34, )),
              (12, 'C=CCC', 3, ()),
              (13, 'CC=CC', 3, ()),
              (14, 'C<-O>=CCC', 3, (34, )),
              (15, 'CC=C<-O>C', 3, (34, )),
              (16, 'C=CC<-O>', 2, (34, )), ]
    for i in range(len(descrs)):
      ID, d, order, ids = descrs[i]
      descr = self.fragCat.GetBitDescription(ID)
      self.assertEqual(descr, d, msg='%d: %s != %s' % (ID, descr, d))
      self.assertEqual(self.fragCat.GetBitOrder(ID), order)
      self.assertEqual(
        tuple(self.fragCat.GetBitFuncGroupIds(ID)), ids,
        msg='%d: %s != %s' % (ID, str(self.fragCat.GetBitFuncGroupIds(ID)), str(ids)))

  def test5MoreComplex(self):
    if not doLong:
      raise unittest.SkipTest('Longer running test skipped')
    lastIdx = 0
    ranges = {}
    suppl = Chem.SmilesMolSupplierFromText('\n'.join(self.smiList), ',', 0, -1, 0)
    for i, mol in enumerate(suppl):
      nEnt = self.fgen.AddFragsFromMol(mol, self.fragCat)
      ranges[i] = range(lastIdx, lastIdx + nEnt)
      lastIdx += nEnt
    # now make sure that those bits are contained in the signatures:
    fpgen = FragmentCatalog.FragFPGenerator()
    for i, mol in enumerate(suppl):
      fp = fpgen.GetFPForMol(mol, self.fragCat)
      for bit in ranges[i]:
        self.assertEqual(fp[bit], 1, msg='%s: %s' % (Chem.MolToSmiles(mol), str(bit)))

  def test6Builder(self):
    suppl = Chem.SmilesMolSupplierFromText('\n'.join(self.smiList2), ',', 0, -1, 0)
    cat = BuildFragmentCatalog.BuildCatalog(suppl, minPath=1, reportFreq=20)
    self.assertEqual(cat.GetNumEntries(), 21)
    self.assertEqual(cat.GetFPLength(), 21)
    self._testBits(cat)

  def test7ScoreMolecules(self):
    suppl = Chem.SmilesMolSupplierFromText('\n'.join(self.smiList2), ',', 0, -1, 0)
    cat = BuildFragmentCatalog.BuildCatalog(suppl, minPath=1, reportFreq=20)
    self.assertEqual(cat.GetNumEntries(), 21)
    self.assertEqual(cat.GetFPLength(), 21)

    scores, obls = BuildFragmentCatalog.ScoreMolecules(suppl, cat, acts=self.list2Acts,
                                                       reportFreq=20)
    for i in range(len(self.list2Obls)):
      self.assertEqual(
        tuple(obls[i]), self.list2Obls[i], msg='%d: %s != %s' % (i, str(obls[i]),
                                                                 str(self.list2Obls[i])))

    scores2 = BuildFragmentCatalog.ScoreFromLists(obls, suppl, cat, acts=self.list2Acts,
                                                  reportFreq=20)
    for i in range(len(scores)):
      self.assertTrue((scores[i] == scores2[i]).all(),
                      msg='%d: %s != %s' % (i, str(scores[i]), str(scores2[i])))

  def test8MolRanks(self):
    suppl = Chem.SmilesMolSupplierFromText('\n'.join(self.smiList2), ',', 0, -1, 0)
    cat = BuildFragmentCatalog.BuildCatalog(suppl, minPath=1, reportFreq=20)
    self.assertEqual(cat.GetNumEntries(), 21)
    self.assertEqual(cat.GetFPLength(), 21)

    # new InfoGain ranking:
    bitInfo, _ = BuildFragmentCatalog.CalcGains(suppl, cat, topN=10, acts=self.list2Acts,
                                                reportFreq=20, biasList=(1, ))
    entry = bitInfo[0]
    self.assertEqual(int(entry[0]), 0)
    self.assertEqual(cat.GetBitDescription(int(entry[0])), 'C<-O>C')
    self.assertAlmostEqual(entry[1], 0.4669, delta=1e-4)

    entry = bitInfo[1]
    self.assertIn(int(entry[0]), (2, 6))
    txt = cat.GetBitDescription(int(entry[0]))
    self.assertIn(txt, ('C<-O>CC', 'C<-O>=C'), msg=txt)
    self.assertAlmostEqual(entry[1], 0.1611, delta=1e-4)

    entry = bitInfo[6]
    self.assertEqual(int(entry[0]), 16)
    self.assertEqual(cat.GetBitDescription(int(entry[0])), 'C=CC<-O>')
    self.assertAlmostEqual(entry[1], 0.0560, delta=1e-4)

    # standard InfoGain ranking:
    bitInfo, _ = BuildFragmentCatalog.CalcGains(suppl, cat, topN=10, acts=self.list2Acts,
                                                reportFreq=20)
    entry = bitInfo[0]
    self.assertEqual(int(entry[0]), 0)
    self.assertEqual(cat.GetBitDescription(int(entry[0])), 'C<-O>C')
    self.assertAlmostEqual(entry[1], 0.4669, delta=1e-4)

    entry = bitInfo[1]
    self.assertEqual(int(entry[0]), 5)
    self.assertEqual(cat.GetBitDescription(int(entry[0])), 'C=CC')
    self.assertAlmostEqual(entry[1], 0.2057, delta=1e-4)

  def test9Issue116(self):
    smiList = ['Cc1ccccc1']
    suppl = Chem.SmilesMolSupplierFromText('\n'.join(smiList), ',', 0, -1, 0)
    cat = BuildFragmentCatalog.BuildCatalog(suppl, minPath=2, maxPath=2)
    self.assertEqual(cat.GetFPLength(), 2)
    self.assertEqual(cat.GetBitDescription(0), 'ccC')
    fpgen = FragmentCatalog.FragFPGenerator()
    mol = Chem.MolFromSmiles('Cc1ccccc1')
    fp = fpgen.GetFPForMol(mol, cat)
    self.assertEqual(fp[0], 1)
    self.assertEqual(fp[1], 1)
    mol = Chem.MolFromSmiles('c1ccccc1-c1ccccc1')
    fp = fpgen.GetFPForMol(mol, cat)
    self.assertEqual(fp[0], 0)
    self.assertEqual(fp[1], 1)


if __name__ == '__main__':  # pragma: nocover
  import argparse
  import sys
  parser = argparse.ArgumentParser()
  parser.add_argument('-l', default=False, action='store_true', dest='doLong')
  args = parser.parse_args()
  doLong = args.doLong

  # Remove the -l flag if present so that it isn't interpreted by unittest.main()
  if '-l' in sys.argv:
    sys.argv.remove('-l')
  unittest.main()
