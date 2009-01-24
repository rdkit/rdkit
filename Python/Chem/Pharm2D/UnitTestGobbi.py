# $Id$
#
#  Copyright (C) 2002-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
"""unit testing code for the signatures

"""
import unittest
from pyRDKit import Chem
from pyRDKit.Chem.Pharm2D import Gobbi_Pharm2D,Generate

class TestCase(unittest.TestCase):
  def setUp(self):
    self.factory = Gobbi_Pharm2D.factory

  def testPatts(self):
    patts = Gobbi_Pharm2D.patts
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
                 'AG':(1,((3,4,5),)),
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
      for k in d.keys():
        shouldMatch,mapList=d[k]
        patt = patts[k][0]
        self.failUnless(mol.HasSubstructMatch(patt)==shouldMatch,'bad match (!=%d) for smi %s and patt %s'%(shouldMatch,smi,k))
        if shouldMatch:
          mapL = mol.GetSubstructMatches(patt)
          self.failUnless(mapL==mapList,
                          'bad match (%s!=%s) for smi %s and patt %s'%(str(mapL),
                                                                       str(mapList),
                                                                       smi,k))


  def testSig(self):
    probes = [('O=CCC=O',1)]
    for smi,tgt in probes:
      sig = Generate.Gen2DFingerprint(Chem.MolFromSmiles(smi),self.factory)
      cnt = len(sig.GetOnBits())
      self.failUnless(cnt==tgt,'bad # onbits for smi %s: %d != %d'%(smi,cnt,tgt))
    


if __name__ == '__main__':
  unittest.main()

