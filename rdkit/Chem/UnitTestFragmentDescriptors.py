# $Id$
#
# Copyright (C) 2003-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Unit Test code for Fragment Descriptors

"""
import unittest
from rdkit import Chem
from rdkit.Chem import Fragments


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def _runTest(self, data, fn):
    for smi, tgtVal in data:
      mol = Chem.MolFromSmiles(smi)
      assert mol, "Smiles parsing failed for %s" % (smi)
      count = fn(mol)
      assert count == tgtVal, "bad value (%d != %d) for smiles %s" % (count, tgtVal, smi)

  def test1(self):
    data = [('C=O', 1),
            ('O=CC=O', 2),
            ('O=CC(=O)O', 2),
            ('O=C(O)C(=O)O', 2), ]
    self._runTest(data, Fragments.fr_C_O)

  def test2(self):
    data = [('C=O', 1),
            ('O=CC=O', 2),
            ('O=CC(=O)O', 1),
            ('O=C(O)C(=O)O', 0), ]
    self._runTest(data, Fragments.fr_C_O_noCOO)

  def test3(self):
    data = [('[O-][N+](C1=C([N+]([O-])=O)C=C([N+]([O-])=O)C=C1)=O', 3),
            ('[O-][N+](C1=CC=CC=C1)=O', 1),
            ('[N+]([O-])(=O)C1OC(=CC=1)/C=N/NC(N)=O', 0),
            ('O=C(O)C(=O)O', 0), ]
    self._runTest(data, Fragments.fr_nitro_arom)

  def test4(self):
    data = [('[O-][N+](C1=C([N+]([O-])=O)C=C([N+]([O-])=O)C=C1)=O', 1),
            ('[O-][N+](C1=CC=CC=C1)=O', 1),
            ('[N+]([O-])(=O)C1OC(=CC=1)/C=N/NC(N)=O', 0),
            ('O=C(O)C(=O)O', 0), ]
    self._runTest(data, Fragments.fr_nitro_arom_nonortho)

  def test5(self):
    data = [('C1(N=[N+]=[N-])=C(C(=C(C(=C1F)F)F)F)F', 1), ]
    self._runTest(data, Fragments.fr_azide)

  def test6(self):
    data = [('C1C=CNC=C1', 1),
            ('C1=CCNC=C1', 1),
            ('C1=CCN=CC1', 1),
            ('C1=CC=NCC1', 1),
            ('C1C=CC(=CC=1[N+]([O-])=O)C2C(=C(NC(=C2C(OC)=O)C)C)C(OCCNCC3=CC=CC=C3)=O', 1),
            ('O=C(O)C(=O)O', 0), ]
    self._runTest(data, Fragments.fr_dihydropyridine)

  def test7(self):
    data = [('OC1=C(C(O)=O)C=CC=C1', 0),
            ('OC1=CC=CC=C1', 1),
            ('OC1=C(O)C=CC=C1', 2),
            ('O=C(O)C1=CC(O)=CC=C1', 1),
            ('OC1=C(O)C(C(O)=O)=CC=C1', 1),
            ('OC1=C(C(N)=O)C=CC=C1', 0),
            ('OC1=C(CO)C=CC=C1', 0),
            ('OC1=C(C(OC)=O)C=CC=C1', 1),
            ('C1=CC=NCC1', 0),
            ('C1C=CC(=CC=1[N+]([O-])=O)C2C(=C(NC(=C2C(OC)=O)C)C)C(OCCNCC3=CC=CC=C3)=O', 0),
            ('O=C(O)C(=O)O', 0), ]
    self._runTest(data, Fragments.fr_phenol_noOrthoHbond)

  def test8(self):
    data = [('OC1=C(C(O)=O)C=CC=C1', 0),
            ('CC(C)(O)C', 0),
            ('CC(O)C', 1),
            ('C=C(O)C', 1),  # includes enols
            ('OC1=C(O)C=CC=C1', 0),
            ('OC1=C(CO)C=CC=C1', 1),
            ('C1C=CC(=CC=1[N+]([O-])=O)C2C(=C(NC(=C2C(OC)=O)C)C)C(OCCNCC3=CC=CC=C3)=O', 0)]
    self._runTest(data, Fragments.fr_Al_OH_noTert)

  def test9(self):
    data = [('OC1=C(C(O)=O)C=CC=C1', 0),
            ('C12=CC=CC=C1NCCN=C2', 1),
            ('C1=C(C=C2C(=C1)NC(CN=C2C3=C(C=CC=C3)Cl)=O)[N+]([O-])=O', 1),  # clonazepam
            ('C12=C(C(=NCC(N1C)=O)C3=CC=CC=C3)C=C(C=C2)Cl', 1),  # valium
            ('CCC1=CC2C(=NCC(=O)N(C)C=2S1)C3=CC=CC=C3Cl', 0),  # clotiazepam has a 5-mb ring
            ('CN(C)CCN1C(=O)C2=CC=CC=C2NC3C=C(Cl)C=CC1=3', 0),  #clobenzepam has a third fused ring
            ('C1C=CC(=CC=1[N+]([O-])=O)C2C(=C(NC(=C2C(OC)=O)C)C)C(OCCNCC3=CC=CC=C3)=O', 0)]
    self._runTest(data, Fragments.fr_benzodiazepine)

  def test10(self):
    data = [('OC1=C(C(O)=O)C=CC=C1', 1),
            ('c1ccccc1Oc2ccccc2', 2),  # diphenyl ether has 2 sites
            ('COC1=CC=CC=C1', 1),
            ('CN(C)C2=CC=CC=C2', 1),
            ('O=C(C4=CC=CC=C4)C3=CC=CC=C3', 0),
            ('COC1=CC=C(C)C(OC)=C1C', 0),
            ('CN(C)C1=CC=CC=C1OC', 2),
            ('O=C(C2=CC=CC=C2)CC1=CC=CC=C1', 0),
            ('O=C(NC2=CC=CC=C2)CN(C)C1=CC=CC=C1', 2),
            ('O=C(N(C)C2=CC=CC=C2)CC1=CC=CC=C1', 1),  # methylated amide is included
            ('C12=CC=CC=C1OCCC2', 1),
            ('C12=CC=CC=C1OC=C2', 1),  # benzofuran
            ('C12=CC=CC=C1NC=N2', 2),  # benzimidazole hydroxylated at 5 position
            ('CC1=CC=C(C=CO2)C2=C1', 0)]
    self._runTest(data, Fragments.fr_para_hydroxylation)

  def test11(self):
    data = [('C1C(C=C2[C@](C1)([C@]3([C@@](CC2)([C@]4([C@](CC3)([C@@H](CC4)O)C)[H])[H])[H])C)=O',
             1),  # testosterone 6beta hydroxylation includedbut not dienone
            ('CCC=CCCC', 2),
            ('C=CC', 1),
            ('C=CCO', 0)]
    self._runTest(data, Fragments.fr_allylic_oxid)

  def test12(self):
    data = [('c1ccccc1C', 1), ('c1ccccc1CC', 1), ('c1(C)cccc(C)c1CC', 2), ('c1(N)cccc(O)c1CC', 0),
            ('c1cccc(C)c1CC', 2), ('c1ccccc1CN', 0), ('c1ccccc1CCN', 0), ('c1(C)ccccc1CCN', 1)]
    self._runTest(data, Fragments.fr_aryl_methyl)

  def test13(self):
    data = [('CNC1(CCCCC1=O)C2(=CC=CC=C2Cl)', 1),  # ketamine
            ('C1=C(C(=C(C=C1)C)NC(CN(CC)CC)=O)C', 1),  # lidocaine
            ('c1(C)ccccc1CCN', 0)]
    self._runTest(data, Fragments.fr_Ndealkylation1)

  def test14(self):
    data = [('C12(=C(N=CC=C1[C@@H](O)C3(N4CC(C(C3)CC4)([H])C=C)[H])C=CC(=C2)OC)', 0),  # quinine
            ('N12CCC(CC2)CC1', 0),  # bridged N
            ('N1(CCC)CCCCC1', 1),
            ('CCC1(CCCCN(C)C1)C2=CC=CC(O)=C2', 1),  # meptazinol
            ('C1(=C2C3=C(C=C1)C[C@@H]4[C@]5([C@@]3([C@H]([C@H](CC5)O)O2)CCN4CC6CCC6)O)O',
             1),  # nalbuphine
            ('CC1N=C2CCCCN2C(=O)C=1CCN3CCC(CC3)C4=NOC5C=C(F)C=CC4=5', 1),  # risperidone
            ('N1CCOCC1', 0),  # morpholino group
            ('n1ccccc1', 0)]
    self._runTest(data, Fragments.fr_Ndealkylation2)

  def test15(self):
    data = [
      ('C(COC(NC(C)C)=O)(COC(N)=O)(CCC)C', 1),  # carisoprodol
      ('CN(CC3=CSC(C(C)C)=N3)C(N[C@@H]([C@H](C)C)C(N[C@@H](CC4=CC=CC=C4)C[C@H](O)[C@H](CC2=CC=CC=C2)NC(OCC1=CN=CS1)=O)=O)=O',
       1),  # ritonavir
      ('c1(C)ccccc1CCN', 0)
    ]
    self._runTest(data, Fragments.fr_alkyl_carbamate)

  def test16(self):
    data = [
      (r'C1(CCC(O1)=O)(C)CC/C=C\CC', 1),
      ('CN(CC3=CSC(C(C)C)=N3)C(N[C@@H]([C@H](C)C)C(N[C@@H](CC4=CC=CC=C4)C[C@H](O)[C@H](CC2=CC=CC=C2)NC(OCC1=CN=CS1)=O)=O)=O',
       0),  # ritonavir
      ('c1(C)ccccc1CCN', 0)
    ]
    self._runTest(data, Fragments.fr_lactone)

  def test17(self):
    data = [('O=C(C=CC1=CC=CC=C1)C=CC2=CC=CC=C2', 0),  # a,b-unsat. dienone
            ('c1ccccc1-C(=O)-c2ccccc2', 0),
            ('c1ccccc1-C(=O)-CCO', 1),
            ('CC(=O)NCCC', 0)]
    self._runTest(data, Fragments.fr_ketone_Topliss)

  def test18(self):
    data = [('C1=CC(=CC=C1C(N[C@@H](CCC(O)=O)C(O)=O)=O)N(CC2=NC3C(=NC(=NC=3N=C2)N)N)C',
             2),  # methotrexate
            ('S(NC1=NC=CC=N1)(=O)(=O)C2=CC=C(C=C2)N', 1),  # sulfadiazine
            ('S(NC1=NC=CC=N1)(=O)(=O)C2=CC=C(C=C2)', 0),
            ('c1ccccc1-C(=O)-CCO', 0),
            ('NNC(=O)C1C2CCCCC=2N=C3C=CC=CC=13', 1),
            ('NC(=N)C1C2CCCCC=2N=C3C=CC=CC=13', 1),
            ('NNC1C2CCCCC=2N=C3C=CC=CC=13', 1),
            ('NC1C2CCCCC=2N=C3C=CC=CC=13', 1)  # tacrine
            ]
    self._runTest(data, Fragments.fr_ArN)

  def test19(self):
    data = [('CC(C)(C)NCC(O)C1=CC(O)=CC(O)=C1', 1),  # terbulatine
            (r'OCCN1CCN(CC/C=C2\C3=CC=CC=C3SC4C=CC(Cl)=CC2=4)CC1', 1),  # clopenthixol
            ('c1ccccc1-C(=O)-CCO', 0),
            ('CC(=O)NCCC', 0)]
    self._runTest(data, Fragments.fr_HOCCN)

  def test20(self):
    data = [('c1ccccc1OC', 1),
            ('c1ccccc1Oc2ccccc2', 0),
            ('c1ccccc1OCC', 0), ]
    self._runTest(data, Fragments.fr_methoxy)

  def test21(self):
    data = [(r'C/C(C(C)=O)=N\O', 1),
            ('C(=N/OC(C(=O)O)(C)C)(/C1=CS[N+](=N1)C)C(N[C@@H]2C(N([C@@H]2C)S(O)(=O)=O)=O)=O',
             1),  # aztreonam
            ('c1ccccc1OCC', 0), ]
    self._runTest(data, Fragments.fr_oxime)

  def test22(self):
    data = [('CCCF', 1),
            ('c1ccccc1C(C)(C)Br', 1),
            ('c1ccccc1CC(C)(C)Cl', 1), ]
    self._runTest(data, Fragments.fr_alkyl_halide)


if __name__ == '__main__':
  unittest.main()
