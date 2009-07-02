#
#  Copyright (C) 2003  Greg Landrum and Rational Discovery LLC
#

"""unit testing code for the SimilarityScreeners

"""
import unittest
from rdkit import Chem
from rdkit.Chem.Fingerprints import SimilarityScreener
from rdkit import DataStructs

def feq(v1,v2,tol=1e-4):
  return abs(v1-v2)<=tol


class TestCase(unittest.TestCase):
  def setUp(self):
    #print '\n%s: '%self.shortDescription(),
    pass

  def test1(self):
    """ TopN screener
    """
    smis = ['C1CCCCC1','C1OCCCC1','C1NCCCC1','c1ccccc1','C1C(C)CCCC1','C1C(C)C(C)CCC1']
    suppl = Chem.SmilesMolSupplierFromText('\n'.join(smis),
                                           delimiter=",",
                                           smilesColumn=0,
                                           nameColumn=-1,
                                           titleLine=0)
    metric = DataStructs.TanimotoSimilarity
    fingerprinter = lambda x:Chem.RDKFingerprint(x,minPath=2,maxPath=7,fpSize=2048)
    probe = fingerprinter(Chem.MolFromSmiles('C1OCCCC1'))

    screener = SimilarityScreener.TopNScreener(3,probe=probe,metric=metric,
                                               fingerprinter=fingerprinter,
                                               dataSource=suppl)
    matches1 = [x for x in screener]
    assert len(matches1)==3
    matches2 = [x for x in screener]
    assert len(matches2)==3
    assert matches1==matches2

  def test2(self):
    """ threshold screener
    """
    smis = ['C1CCCCC1','C1OCCCC1','C1NCCCC1','c1ccccc1','C1C(C)CCCC1','C1C(C)C(C)CCC1']
    suppl = Chem.SmilesMolSupplierFromText('\n'.join(smis),
                                           delimiter=",",
                                           smilesColumn=0,
                                           nameColumn=-1,
                                           titleLine=0)

    metric = DataStructs.TanimotoSimilarity
    fingerprinter = lambda x:Chem.RDKFingerprint(x,minPath=2,maxPath=7,fpSize=2048)
    probe = fingerprinter(Chem.MolFromSmiles('C1OCCCC1'))

    screener = SimilarityScreener.ThresholdScreener(0.09,probe=probe,metric=metric,
                                                    fingerprinter=fingerprinter,
                                                    dataSource=suppl)
    matches1 = [x[0] for x in screener]
    assert len(matches1)==5
    matches2 = [x[0] for x in screener]
    assert len(matches2)==5
    assert matches1==matches2

  def test3(self):
    """ threshold screener, including folding
    """
    smis = ['C1CCCCC1','C1OCCCC1','C1NCCCC1','c1ccccc1','C1C(C)CCCC1','C1C(C)C(C)CCC1']
    suppl = Chem.SmilesMolSupplierFromText('\n'.join(smis),
                                           delimiter=",",
                                           smilesColumn=0,
                                           nameColumn=-1,
                                           titleLine=0)


    metric = DataStructs.TanimotoSimilarity
    fingerprinter = lambda x:Chem.RDKFingerprint(x,minPath=2,maxPath=7,fpSize=2048)
    probe = Chem.RDKFingerprint(Chem.MolFromSmiles('C1OCCCC1'),
                                     minPath=2,maxPath=7,fpSize=4096)

    screener = SimilarityScreener.ThresholdScreener(0.09,probe=probe,metric=metric,
                                                    fingerprinter=fingerprinter,
                                                    dataSource=suppl)
    matches1 = [x[0] for x in screener]
    assert len(matches1)==5
    matches2 = [x[0] for x in screener]
    assert len(matches2)==5
    assert matches1==matches2


      
      
    

if __name__ == '__main__':
  unittest.main()

