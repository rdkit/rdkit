"""
This file contains tests to be verified for the reorderTautomers() so defined.
The function reorderTautomers() return the list of all tautomers in the order 
in which the canonical Tautomer will always occurs at the first position.
This would help the chemist to identify the canonical tautomer easily by looking 
at the first position.

"""
# twenty molecules/ compounds are selected for testing which are decribed in-line at each test definition.
# 20 molecules/ Compounds are listed as follows in the form of smiles:

#test input smiles =
#['Oc1c(cccc3)c3nc2ccncc12',
#'C=C(O)Cc1ccccc1',
#'OC(C)=C(C)C',
#'Oc1nccc2cc[nH]c(=N)c12',
#'C1(C=CCCC1)=O',
#'C(=C)(O)C1=CC=CC=C1',
#'C1(CCCCC1)=N',
#'C1(=NC=CC=C1)CC',
#'Oc1ccccn1',
#'C(=N)=O',
#'O=Nc1ccc(O)cc1',
#'CC(C)=NO',
#'cc1cc(=O)[nH]c2nccn12',
#'O=C1CC=CO1',
#'n1ccc2c[nH]ccc12',
#'NC1C=CC(=Cc2ccc(O)cc2)C=C1',
#'CN=c1ccnc2nc[nH]n21',
#'C1(=CCCCC1)O',
#'S=C(N)N',
#'CNC(C)=O'
]

# The above tuple enlists the inputted smiles used for testing the reorderTautomers().
# The above compounds in the form of smiles are descripted at each test function().   



import unittest
import reorderTautomers
from rdkit import Chem
from rdkit.Chem.MolStandardize import canonicalize,enumerate 



class TestReorderCanonicalTautomer(unittest.TestCase)
"""20 test function tested for the canonical tautomer to be at the first position"""

      def test_1():

          """ an examplar for 1,5 keto/enol tautomer"""
          # smile('Oc1nccc2cc[nH]c(=N)c12')

          assert reorderTautomers('Oc1nccc2cc[nH]c(=N)c12') == ['Nc1nccc2cc[nH]c(=O)c12','Nc1nccc2ccnc(O)c12','N=c1[nH]ccc2ccnc(O)c12','Nc1[nH]ccc2ccnc(=O)c1-2','N=c1[nH]ccc2cc[nH]c(=O)c12','N=c1nccc2cc[nH]c(O)c1-2']


     def test_2():

      """ Aromatic compound with three membered ring with two nitrogen and single oxygen atoms"""
          #smile('Oc1c(cccc3)c3nc2ccncc12')

          assert reorderTautomers('Oc1c(cccc3)c3nc2ccncc12') == ['O=c2c1ccccc1nc3cc[nH]cc23','O=c2c1ccccc1[nH]c3ccnccc23','Oc2c1ccccc1nc3ccncc23']


      def test_3():


         """an examplar 1-phenyl-2-propanone enol/keto"""
         #smile('C=C(O)Cc1ccccc1')

          assert reorderTautomers('C=C(O)Cc1ccccc1') == ['CC(=O)Cc1ccccc1','C=C(O)Cc1ccccc1','CC(O)=Cc1ccccc1']


      def test_4():
         """a simple illustration for keto/ enol tautomer"""
         #a compound with only setreocenter
         #smile('OC(C)=C(C)C')

          assert reorderTautomers('OC(C)=C(C)C') == ['CC(=O)C(C)C','C=C(O)C(C)C','CC(C)=C(C)O']
  
      def test_5():
          
          """a family of six-membered ring containg one double bond and a single oxygen atom"""
         # examplar for 1,5 keto/enol tautomer
         #ssmile('C1(C=CCCC1)=O')

          assert reorderTautomers('C1(C=CCCC1)=O') == ['O=C1CC=CCC1','O=C1C=CCCC1','OC1=CCC=CC1','OC1=CC=CCC1','OC1=CCCC=C1']

      def test_6():

          """an example for acetophenone tautomer"""
          # a keto/ enol tautomer examplar
          #smile('C(=C)(O)C1=CC=CC=C1')

          assert reorderTautomers('C(=C)(O)C1=CC=CC=C1') == ['CC(=O)c1ccccc1','C=C(O)c1ccccc1']
  
      def test_7():
          
          """an examplar for pyridine tautomer"""
          # aliphatic imine tautomer
          #smile('C1(CCCCC1)=N')

          assert reorderTautomers('C1(CCCCC1)=N') == ['N=C1CCCCC1','NC1=CCCCC1']

      def test_8():
           
          """A pyridine compund with a side branch for tautomerization"""
          #imine tautomer
          #smile('C1(=NC=CC=C1)CC')

          assert reorderTautomers('C1(=NC=CC=C1)CC') == ['CC=C1C=CCC=N1','CC=C1C=CC=CN1','CCc1ccccn1']

      def test_9():

          # an example for testing 1,5 heteroatom H shift
          # smile('Oc1ccccn1')

          assert reorderTautomers('Oc1ccccn1') == ['Oc1ccccn1','O=c1cccc[nH]1']

      def test_10():

              """an examplar for cyano/iso-cyanic acid tautomer"""
              # smile('C(=N)=O')

          assert reorderTautomers('C(=N)=O') == ['N=c=O','N#CO']

      def test_11():

          """an examplar for phenolic tautomer"""
          """oxim/nitroso tautomer via phenol"""
          # smile('O=Nc1ccc(O)cc1')

          assert reorderTautomers('O=Nc1ccc(O)cc1') == ['O=NC1C=CC(=O)C=C1','O=C1C=CC(=NO)C=C1','O=Nc1ccc(O)cc1']

      def test_12():

          """oxim nitroso tautomer"""
          # smile('CC(C)=NO')

          assert reorderTautomers('CC(C)=NO') == ['CC(C)NO','CC(C)=NO','C=C(C)NO']

      def test_13():

          # furanone tautomer
          # smile('O=C1CC=CO1')
          """shows two tautomer containing two oxygen"""

          assert reorderTautomers('O=C1CC=CO1') == ['O=C1CC=CO1','Oc1cccO1']

      def test_14():

          # aromatic heteroatom tautomer
           """1,5 aromatic heteroatom H shift"""
          #smile('cc1cc(=O)[nH]c2nccn12')

          assert reorderTautomers('cc1cc(=O)[nH]c2nccn12') == ['O=c1ccn2cc[nH]c2n1','Oc1ccn2ccnc2n1','O=c1ccn2ccnc2[nH]1']

      def test_15():

          """An examplar for heterocyclic tautomer"""
          #smile('n1ccc2c[nH]ccc12')

          assert reorderTautomers('n1ccc2c[nH]ccc12') == ['c1cc2[nH]ccc2cn1','c1cc2c[nH]ccc-2n1']

      def test_16():
          
          """aromatic heteroatom tautomer"""
          # an examplar for 1,11 aromatic heteroatom H shift
          # smile('NC1C=CC(=Cc2ccc(O)cc2)C=C1')

          assert reorderTautomers('NC1C=CC(=Cc2ccc(O)cc2)C=C1') == ['N=C1C=CC(C=C2C=CC(=O)C=C2)C=C1','Nc1ccc(C=C2C=CC(=O)C=C2)cc1','N=C1C=CC(=CC2C=CC(=O)C=C2)C=C1','N=C1C=CC(Cc2ccc(O)cc2)C=C1']

      def test_17():

          """An examplar for Aromtaic tautomer"""
          # 1,9 aromatic heteroatom H shift
          # smile('CN=c1ccnc2nc[nH]n21')
          assert reorderTautomers('CN=c1ccnc2nc[nH]n21') == ['CN=c1cc[nH]c2ncnn12','CN=c1ccnc2nc[nH]n12','CN=c1ccnc2[nH]cnn12','CNc1ccnc2ncnn12']

      def test_18():

          """ an examplar for 1,3 keto/enol tautomer"""
          #smile('C1(=CCCCC1)O')
          # a simple test for the keto/ enol tautomer

          assert reorderTautomers('C1(=CCCCC1)O') == ['O=C1CCCCC1','OC1=CCCCC1']

      def test_19():

          """an examplar for thiocyno tautomer"""
              #smile('S=C(N)N')
              # an examplar for 1,3 heteroatom H shift

          assert reorderTautomers('S=C(N)N') == ['N=C(N)S','NC(N)=S']

      def test_20():

              """an examplar for oxim nitroso tautomer"""
              #smile('CNC(C)=O')

          assert reorderTautomers('CNC(C)=O') == ['CNC(C)=O','CN=C(C)O','C=C(O)NC']
      







           
