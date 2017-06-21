from rdkit import Chem
from rdkit import rdBase
from rdkit import RDConfig
import os

from rdkit.Chem import rdMolDescriptors as rdMD
from rdkit.Chem import AllChem

import time,unittest

def _gen3D(m,is3d,calculator):
    if not is3d:
        m = Chem.AddHs(m)
        ps = AllChem.ETKDG()
        ps.randomSeed = 0xf00d
        AllChem.EmbedMolecule(m,ps)
    return calculator(m)


class TestCase(unittest.TestCase):

  def setUp(self):
    self.dataDir = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
    'Descriptors','test_data')
    self.suppl = Chem.SDMolSupplier(os.path.join(self.dataDir,'PBF_egfr.sdf'),removeHs=False)

  def test1AUTOCORR3D(self):
      with open(os.path.join(self.dataDir,'auto3D_dragon.out')) as refFile:
          for i,m in enumerate(self.suppl):
              if i>10: break
              nm = m.GetProp('_Name')
              inl = refFile.readline()
              split = inl.split('\t')
              self.assertEqual(split[0],nm)
              split.pop(0)
              vs = _gen3D(m,True,rdMD.CalcAUTOCORR3D)
              for rv,nv in zip(split,vs):
                  self.assertAlmostEqual(float(rv),nv,delta=0.05)
  def test1AUTOCORR2D(self):
      with open(os.path.join(self.dataDir,'auto2D.out')) as refFile:
          for i,m in enumerate(self.suppl):
              if i>10: break
              nm = m.GetProp('_Name')
              inl = refFile.readline()
              split = inl.split('\t')
              self.assertEqual(split[0],nm)
              split.pop(0)
              vs = rdMD.CalcAUTOCORR2D(m)
              for rv,nv in zip(split,vs):
                  self.assertAlmostEqual(float(rv),nv,delta=0.05)


if(__name__=='__main__'):
  unittest.main()
