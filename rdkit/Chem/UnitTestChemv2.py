# $Id$
#
#  Copyright (C) 2003-2006  Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""basic unit testing code for the rdkit Boost wrapper

"""
import unittest, os
import pickle
from rdkit import RDConfig
from rdkit import Chem
from rdkit.Chem import AllChem


class TestCase(unittest.TestCase):

  def setUp(self):
    self.bigSmiList = [
      "CC1=CC(=O)C=CC1=O",
      "S(SC1=NC2=CC=CC=C2S1)C3=NC4=C(S3)C=CC=C4",
      "OC1=C(Cl)C=C(C=C1[N+]([O-])=O)[N+]([O-])=O",
      "[O-][N+](=O)C1=CNC(=N)S1",
      "NC1=CC2=C(C=C1)C(=O)C3=C(C=CC=C3)C2=O",
      "OC(=O)C1=C(C=CC=C1)C2=C3C=CC(=O)C(=C3OC4=C2C=CC(=C4Br)O)Br",
      "CN(C)C1=C(Cl)C(=O)C2=C(C=CC=C2)C1=O",
      "CC1=C(C2=C(C=C1)C(=O)C3=CC=CC=C3C2=O)[N+]([O-])=O",
      "CC(=NO)C(C)=NO",
      "C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3",
      "CC(C)(C)C1=C(O)C=C(C(=C1)O)C(C)(C)C",
      "CC1=NN(C(=O)C1)C2=CC=CC=C2",
      "NC1=CC=NC2=C1C=CC(=C2)Cl",
      "CCCCCC[CH]1CCCCN1",
      "O=CC1=C2C=CC=CC2=CC3=C1C=CC=C3",
      "BrN1C(=O)CCC1=O",
      "CCCCCCCCCCCCCCCC1=C(N)C=CC(=C1)O",
      "C(COC1=C(C=CC=C1)C2=CC=CC=C2)OC3=CC=CC=C3C4=CC=CC=C4",
      "CCCCSCC",
      "CC(=O)NC1=NC2=C(C=C1)C(=CC=N2)O",
      "CC1=C2C=CC(=NC2=NC(=C1)O)N",
      "CCOC(=O)C1=CN=C2N=C(N)C=CC2=C1O",
      "CC1=CC(=NC=C1)N=CC2=CC=CC=C2",
      "C[N+](C)(C)CC1=CC=CC=C1",
      "C[N+](C)(C)C(=O)C1=CC=CC=C1",
      "ICCC(C1=CC=CC=C1)(C2=CC=CC=C2)C3=CC=CC=C3",
      "CC1=CC(=C(C[N+](C)(C)C)C(=C1)C)C",
      "C[C](O)(CC(O)=O)C1=CC=C(C=C1)[N+]([O-])=O",
      "CC1=CC=C(C=C1)C(=O)C2=CC=C(Cl)C=C2",
      "ON=CC1=CC=C(O)C=C1",
      "CC1=CC(=C(N)C(=C1)C)C",
      "CC1=CC=C(C=C1)C(=O)C2=CC=C(C=C2)[N+]([O-])=O",
      "CC(O)(C1=CC=CC=C1)C2=CC=CC=C2",
      "ON=CC1=CC(=CC=C1)[N+]([O-])=O",
      "OC1=C2C=CC(=CC2=NC=C1[N+]([O-])=O)Cl",
      "CC1=CC=CC2=NC=C(C)C(=C12)Cl",
      "CCC(CC)([CH](OC(N)=O)C1=CC=CC=C1)C2=CC=CC=C2",
      "ON=C(CC1=CC=CC=C1)[CH](C#N)C2=CC=CC=C2",
      "O[CH](CC1=CC=CC=C1)C2=CC=CC=C2",
      "COC1=CC=C(CC2=CC=C(OC)C=C2)C=C1",
      "CN(C)[CH](C1=CC=CC=C1)C2=C(C)C=CC=C2",
      "COC1=CC(=C(N)C(=C1)[N+]([O-])=O)[N+]([O-])=O",
      "NN=C(C1=CC=CC=C1)C2=CC=CC=C2",
      "COC1=CC=C(C=C1)C=NO",
      "C1=CC=C(C=C1)C(N=C(C2=CC=CC=C2)C3=CC=CC=C3)C4=CC=CC=C4",
      "C1=CC=C(C=C1)N=C(C2=CC=CC=C2)C3=CC=CC=C3",
      "CC1=C(C2=CC=CC=C2)C(=C3C=CC=CC3=N1)O",
      "CCC1=[O+][Cu]2([O+]=C(CC)C1)[O+]=C(CC)CC(=[O+]2)CC",
      "OC(=O)[CH](CC1=CC=CC=C1)C2=CC=CC=C2",
      "CCC1=C(N)C=C(C)N=C1",
    ]

  def test1(self):
    """ basic building stuff """
    m = Chem.MolFromSmiles('COC(=O)O')

    a1 = m.GetAtomWithIdx(1)
    assert a1.GetAtomicNum() == 8
    assert m.GetAtomWithIdx(2).GetAtomicNum() == 6
    b1 = m.GetBondWithIdx(1)
    assert b1.GetBondType() == Chem.BondType.SINGLE
    assert m.GetBondWithIdx(2).GetBondType() == Chem.BondType.DOUBLE
    assert m.GetBondBetweenAtoms(0, 1).GetBondType() == Chem.BondType.SINGLE
    assert m.GetBondBetweenAtoms(2, 3).GetBondType() == Chem.BondType.DOUBLE

  def test2(self):
    """ editing/persistence basics """
    m = Chem.MolFromSmiles('COC(=C)O')

    a1 = m.GetAtomWithIdx(3)
    assert a1.GetAtomicNum() == 6, 'bad atom order'
    a1.SetAtomicNum(7)
    assert a1.GetAtomicNum() == 7, 'bad atom order'
    assert m.GetAtomWithIdx(3).GetAtomicNum() == 7, 'atom order not stored'

  def test3(self):
    """ SMARTS basics """
    m = Chem.MolFromSmiles('COC(=O)O')
    p = Chem.MolFromSmarts('CO')
    assert m.HasSubstructMatch(p)
    p2 = Chem.MolFromSmarts('CS')
    assert not m.HasSubstructMatch(p2)

    #assert p.GetSMARTS()=='CO','bad getsmarts'
    assert p.GetNumAtoms() == 2
    assert p.GetNumBonds() == 1

    assert m.HasSubstructMatch(p)
    matches = m.GetSubstructMatches(p)
    assert len(matches) == 3
    match = matches[0]
    assert len(match) == 2, 'bad match length'
    matches = m.GetSubstructMatches(p, 0)
    assert len(matches) == 3
    match = matches[0]
    assert len(match) == 2, 'bad match length'

    p = Chem.MolFromSmarts('COC')
    assert m.HasSubstructMatch(p)
    matches = m.GetSubstructMatches(p)
    assert len(matches) == 1
    matches = m.GetSubstructMatches(p, 0)
    assert len(matches) == 2

  def test4Pkl2(self):
    """ further pickle tests """
    smis = self.bigSmiList

    for smi in smis:
      m = Chem.MolFromSmiles(smi)
      newM1 = pickle.loads(pickle.dumps(m))
      newM2 = pickle.loads(pickle.dumps(newM1))
      oldSmi = Chem.MolToSmiles(newM1)
      newSmi = Chem.MolToSmiles(newM2)
      assert newM1.GetNumAtoms() == m.GetNumAtoms(), 'num atoms comparison failed'
      assert newM2.GetNumAtoms() == m.GetNumAtoms(), 'num atoms comparison failed'
      assert oldSmi == newSmi, 'string compare failed: %s != %s' % (oldSmi, newSmi)

  def test5Data(self):
    """ testing Get/Set/HasData """
    m = Chem.MolFromSmiles('CCOC')
    try:
      m.SetProp('foo', '3')
    except Exception:
      ok = 0
    else:
      ok = 1
    assert ok

    try:
      v = m.GetProp('foo')
    except Exception:
      ok = 0
    else:
      ok = 1
    assert ok
    assert v == '3'

    try:
      v = m.GetProp('monkey')
    except KeyError:
      ok = 1
    except Exception:
      ok = 0
    else:
      ok = 0
    assert ok

  def testIssue399(self):
    m = Chem.MolFromSmiles('C[C@H]1CO1')
    AllChem.Compute2DCoords(m)
    Chem.WedgeMolBonds(m, m.GetConformer())
    self.assertTrue(m.GetBondWithIdx(0).GetBondDir() == Chem.rdchem.BondDir.BEGINDASH)
    self.assertTrue(m.GetBondWithIdx(1).GetBondDir() == Chem.rdchem.BondDir.NONE)
    self.assertTrue(m.GetBondWithIdx(2).GetBondDir() == Chem.rdchem.BondDir.NONE)
    self.assertTrue(m.GetBondWithIdx(3).GetBondDir() == Chem.rdchem.BondDir.NONE)


if __name__ == '__main__':
  unittest.main()
