# $Id$
#
#  Copyright (C) 2007 Greg Landrum
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from rdkit import RDConfig
import unittest
from rdkit import Chem

from rdkit.Chem import AllChem
from rdkit.Chem.Subshape import SubshapeBuilder,SubshapeObjects,SubshapeAligner
from rdkit.six.moves import cPickle
import copy,sys,os

class TestCase(unittest.TestCase):
  def setUp(self):
    pass
  def test1(self):
    suppl = Chem.SDMolSupplier(os.path.join(RDConfig.RDCodeDir,'Chem','Subshape',
                                        'test_data/5ht3ligs.sdf'))
    builder = SubshapeBuilder.SubshapeBuilder()
    builder.gridDims=(20.,20.,10)
    builder.gridSpacing=0.5
    builder.winRad=4.

    ms = []
    shapes=[]
    for m in suppl:
      m = Chem.AddHs(m,addCoords=True)
      AllChem.CanonicalizeConformer(m.GetConformer())
      ms.append(m)
      shape = builder(m,terminalPtsOnly=True)
      shapes.append(shape)
    
    self.failUnless(len(ms)==4)
    self.failUnless(len(shapes)==4)
    self.failUnless([len(x.skelPts) for x in shapes] == [5,5,5,5])

    refShape = builder.GenerateSubshapeShape(ms[0])
    self.failUnless(len(refShape.skelPts)==15)

    aligner = SubshapeAligner.SubshapeAligner()
    aligner.shapeDistTol=.30

    algStore = []
    for i,s1 in enumerate(shapes):
      if not i or not s1: 
        algStore.append([])
        continue
      m1 = ms[i]
      alignments = aligner.GetSubshapeAlignments(ms[0],refShape,m1,s1,builder)
      algStore.append(alignments)
    self.failUnlessEqual([len(x) for x in algStore],[0,2,39,0])

    algStore = []
    for i,s1 in enumerate(shapes):
      if not i or not s1: 
        algStore.append([])
        continue
      m1 = ms[i]
      alignments = list(aligner(ms[0],refShape,m1,s1,builder))
      algStore.append(alignments)
    self.failUnless([len(x) for x in algStore] == [0,2,39,0])


    
    pruned=[]
    for i,mi in enumerate(ms):
      alignments=algStore[i]
      pruned.append(SubshapeAligner.ClusterAlignments(mi,alignments,builder,
                                                      neighborTol=0.15))
    self.failUnless([len(x) for x in pruned] == [0,2,29,0])


    



if __name__ == '__main__':
  unittest.main()


