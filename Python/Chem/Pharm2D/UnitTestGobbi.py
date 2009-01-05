# $Id$
#
#  Copyright (C) 2002-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
"""unit testing code for the signatures

"""
import unittest
import Chem
from Chem.Pharm2D import Gobbi_Pharm2D,Generate

class TestCase(unittest.TestCase):
  def setUp(self):
    self.factory = Gobbi_Pharm2D.factory

  def test1Sigs(self):
    probes = [
      ('OCCC=O',{'HA':(1,((0,),(4,))),
                 'HD':(1,((0,),)),
                 'LH':(0,None),
                 'AR':(0,None),
                 'RR':(0,None),
                 'X':(0,None),
                 'BG':(0,None),
                 'AG':(0,None),
                 }
       ),
      ('OCCC(=O)O',{'HA':(1,((0,),(4,))),
                 'HD':(1,((0,),(5,))),
                 'LH':(0,None),
                 'AR':(0,None),
                 'RR':(0,None),
                 'X':(0,None),
                 'BG':(0,None),
                 'AG':(1,((3,),)),
                 }
       ),
      ('CCCN',{'HA':(1,((3,),)),
                 'HD':(1,((3,),)),
                 'LH':(0,None),
                 'AR':(0,None),
                 'RR':(0,None),
                 'X':(0,None),
                 'BG':(1,((3,),)),
                 'AG':(0,None),
                 }
       ),
      ('CCCCC',{'HA':(0,None),
                 'HD':(0,None),
                 'LH':(1,((1,),(3,))),
                 'AR':(0,None),
                 'RR':(0,None),
                 'X':(0,None),
                 'BG':(0,None),
                 'AG':(0,None),
                 }
       ),
      ('CC1CCC1',{'HA':(0,None),
                 'HD':(0,None),
                 'LH':(1,((1,),(3,))),
                 'AR':(0,None),
                 'RR':(1,((1,),)),
                 'X':(0,None),
                 'BG':(0,None),
                 'AG':(0,None),
                 }
       ),
      ('[SiH3]C1CCC1',{'HA':(0,None),
                 'HD':(0,None),
                 'LH':(1,((1,),)),
                 'AR':(0,None),
                 'RR':(1,((1,),)),
                 'X':(1,((0,),)),
                 'BG':(0,None),
                 'AG':(0,None),
                 }
       ),
      ('[SiH3]c1ccccc1',{'HA':(0,None),
                          'HD':(0,None),
                          'LH':(0,None),
                          'AR':(1,((1,),)),
                          'RR':(0,None),
                          'X':(1,((0,),)),
                          'BG':(0,None),
                          'AG':(0,None),
                          }
       ),
      ]
    for smi,d in probes:
      mol = Chem.MolFromSmiles(smi)
      feats=self.factory.featFactory.GetFeaturesForMol(mol)
      for k in d.keys():
        shouldMatch,mapList=d[k]
        feats=self.factory.featFactory.GetFeaturesForMol(mol,includeOnly=k)
        if shouldMatch:
          self.failUnless(feats)
          self.failUnlessEqual(len(feats),len(mapList))
          aids = [(x.GetAtomIds()[0],) for x in feats]
          aids.sort()
          self.failUnlessEqual(tuple(aids),mapList)

  def test2Sigs(self):
    probes = [('O=CCC=O',(149,)),
              ('OCCC=O',(149,156)),
              ('OCCC(=O)O',(22, 29, 149, 154, 156, 184, 28810, 30055)),
              ]
    for smi,tgt in probes:
      sig = Generate.Gen2DFingerprint(Chem.MolFromSmiles(smi),self.factory)
      self.failUnlessEqual(len(sig),39972)
      bs = tuple(sig.GetOnBits())
      self.failUnlessEqual(len(bs),len(tgt))
      self.failUnlessEqual(bs,tgt)

if __name__ == '__main__':
  unittest.main()

