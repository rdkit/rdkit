from rdkit import RDConfig
import os, sys
import unittest
from rdkit import Chem
from rdkit.Chem import rdMolHash


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def test1(self):
    m = Chem.MolFromSmiles('C1CCCC(O)C1c1ccnc(OC)c1')
    self.assertEqual(rdMolHash.MolHash(m, rdMolHash.HashFunction.AnonymousGraph),
                     '***1****(*2*****2*)*1')
    self.assertEqual(rdMolHash.MolHash(m, rdMolHash.HashFunction.ElementGraph),
                     'COC1CC(C2CCCCC2O)CCN1')
    self.assertEqual(rdMolHash.MolHash(m, rdMolHash.HashFunction.CanonicalSmiles),
                     'COc1cc(C2CCCCC2O)ccn1')
    self.assertEqual(rdMolHash.MolHash(m, rdMolHash.HashFunction.MurckoScaffold),
                     'c1cc(C2CCCCC2)ccn1')
    self.assertEqual(rdMolHash.MolHash(m, rdMolHash.HashFunction.ExtendedMurcko),
                     '*c1cc(C2CCCCC2*)ccn1')
    self.assertEqual(rdMolHash.MolHash(m, rdMolHash.HashFunction.MolFormula), 'C12H17NO2')
    self.assertEqual(rdMolHash.MolHash(m, rdMolHash.HashFunction.AtomBondCounts), '15,16')
    self.assertEqual(rdMolHash.MolHash(m, rdMolHash.HashFunction.DegreeVector), '0,4,9,2')
    self.assertEqual(rdMolHash.MolHash(m, rdMolHash.HashFunction.Mesomer),
                     'CO[C]1[CH][C](C2CCCCC2O)[CH][CH][N]1_0')
    self.assertEqual(rdMolHash.MolHash(m, rdMolHash.HashFunction.Regioisomer),
                     '*O.*O*.C.C1CCCCC1.c1ccncc1')
    self.assertEqual(rdMolHash.MolHash(m, rdMolHash.HashFunction.NetCharge), '0')
    self.assertEqual(rdMolHash.MolHash(m, rdMolHash.HashFunction.SmallWorldIndexBR), 'B16R2')
    self.assertEqual(rdMolHash.MolHash(m, rdMolHash.HashFunction.SmallWorldIndexBRL), 'B16R2L9')
    self.assertEqual(rdMolHash.MolHash(m, rdMolHash.HashFunction.ArthorSubstructureOrder),
                     '000f001001000c000300005f000000')


if __name__ == "__main__":
  unittest.main()
