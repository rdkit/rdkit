# $Id$
#
# Copyright (C) 2007 by Greg Landrum 
#  All rights reserved
#
import Chem,Geometry
import Numeric
from Numerics import Alignment
from Chem.Subshape import SubshapeObjects

class SubshapeAlignment(object):
  transform=None
  triangleSSD=None
  targetTri=None
  queryTri=None
  alignedConfId=-1
  
def _getAllTriangles(pts,orderedTraversal=False):
  for i in range(len(pts)):
    if orderedTraversal:
      jStart=i+1
    else:
      jStart=0
    for j in range(jStart,len(pts)):
      if j==i:
        continue
      if orderedTraversal:
        kStart=j+1
      else:
        kStart=0
      for k in range(j+1,len(pts)):
        if k==i or k==j:
          continue
        yield (i,j,k)

class SubshapeDistanceMetric(object):
  TANIMOTO=0
  PROTRUDE=1

def TransformMol(mol,tform,confId=-1,newConfId=100):
  """  Applies the transformation to a molecule and sets it up with
  a single conformer

  """
  newConf = Chem.Conformer()
  newConf.SetId(0)
  refConf = mol.GetConformer(confId)
  for i in range(refConf.GetNumAtoms()):
    pos = list(refConf.GetAtomPosition(i))
    pos.append(1.0)
    newPos = Numeric.matrixmultiply(tform,Numeric.array(pos))
    newConf.SetAtomPosition(i,list(newPos)[:3])
    newConf.SetId(newConfId)
  mol.RemoveConformer(newConfId)
  mol.AddConformer(newConf,assignId=False)
  
class SubshapeAligner(object):
  triangleRMSTol=0.5
  distMetric=SubshapeDistanceMetric.PROTRUDE
  shapeDistTol=0.2
  numFeatThresh=3
  dirThresh=2.6
  
  def GetTriangleMatches(self,target,query):
    ssdTol = (self.triangleRMSTol**2)*9
    res = []
    tgtPts = target.skelPts
    queryPts = query.skelPts
    for tgtTri in _getAllTriangles(tgtPts,orderedTraversal=True):
      tgtLocs=[tgtPts[x].location for x in tgtTri]
      for queryTri in _getAllTriangles(queryPts):
        queryLocs=[queryPts[x].location for x in queryTri]
        ssd,tf = Alignment.GetAlignmentTransform(tgtLocs,queryLocs)
        if ssd<=ssdTol:
          alg = SubshapeAlignment()
          alg.transform=tf
          alg.triangleSSD=ssd
          alg.targetTri=tgtTri
          alg.queryTri=queryTri
          res.append(alg)
    return res

  def PruneMatchesUsingFeatures(self,target,query,alignments):
    i = 0
    tgtPts = target.skelPts
    queryPts = query.skelPts
    while i<len(alignments):
      nMatched=0
      alg = alignments[i]
      for j in range(3):
        tgtFeats = tgtPts[alg.targetTri[j]].molFeatures
        qFeats = queryPts[alg.queryTri[j]].molFeatures
        for k,kFeat in enumerate(tgtFeats):
          if kFeat in qFeats:
            nMatched+=1
            break
      if nMatched<self.numFeatThresh:
        del alignments[i]
      else:
        i+=1

  def PruneMatchesUsingDirection(self,target,query,alignments):
    i = 0
    tgtPts = target.skelPts
    queryPts = query.skelPts
    while i<len(alignments):
      removeIt=False
      nMatched=0
      alg = alignments[i]
      prunedTf = alg.transform[:3,:3]
      dot = 0.0
      for j in range(3):
        tgtPt = tgtPts[j]
        queryPt = queryPts[j]
        m1,m2,m3=tgtPt.shapeMoments
        tgtR = m1/(m2+m3)
        m1,m2,m3=queryPt.shapeMoments
        queryR = m1/(m2+m3)

        qv = Numeric.array(queryPt.shapeDirs[0])
        tv = Numeric.array(tgtPt.shapeDirs[0])
        rotV = Numeric.matrixmultiply(prunedTf,qv)
        dot += abs(Numeric.sum(rotV*tv))
        if dot>self.dirThresh:
          removeIt=True
          break
      print dot
      if removeIt:
        del alignments[i]
      else:
        alignments[i].dirMatch=dot
        i+=1
    

  def _addCoarseAndMediumGrids(self,mol,tgt,confId,builder):
    oSpace=builder.gridSpacing
    builder.gridSpacing = oSpace*1.5
    tgt.medGrid = builder.GenerateSubshapeShape(mol,confId,addSkeleton=False)
    builder.gridSpacing = oSpace*2
    tgt.coarseGrid = builder.GenerateSubshapeShape(mol,confId,addSkeleton=False)
    builder.gridSpacing = oSpace

  def _getShapeShapeDistance(self,s1,s2):
    #print s1.grid.GetNumX(),s1.grid.GetNumY(),s1.grid.GetNumZ()
    #print s2.grid.GetNumX(),s2.grid.GetNumY(),s2.grid.GetNumZ()
    #print tuple(s1.grid.GetOffset()),s1.grid.GetSpacing()
    #print tuple(s2.grid.GetOffset()),s2.grid.GetSpacing()
    #print
    if self.distMetric==SubshapeDistanceMetric.PROTRUDE:
      if s1.grid.GetOccupancyVect().GetTotalVal()<s2.grid.GetOccupancyVect().GetTotalVal():
        d = Geometry.ProtrudeDistance(s1.grid,s2.grid)
      else:
        d = Geometry.ProtrudeDistance(s2.grid,s1.grid)
    else:
      d = Geometry.TanimotoDistance(s1.grid,s2.grid)
    return d

  def PruneMatchesUsingShape(self,targetMol,target,queryMol,query,builder,
                             alignments,tgtConf=-1,queryConf=-1):
    if not hasattr(target,'medGrid'):
      self._addCoarseAndMediumGrids(targetMol,target,tgtConf,builder)
      self._addCoarseAndMediumGrids(targetMol,target,tgtConf,builder)
    i=0
    while i < len(alignments):
      removeIt=False
      tConfId=100+i
      alg = alignments[i]
      TransformMol(queryMol,alg.transform,confId=queryConf,newConfId=tConfId)
      alg.alignedConfId=tConfId
      oSpace=builder.gridSpacing
      builder.gridSpacing=oSpace*2
      coarseGrid=builder.GenerateSubshapeShape(queryMol,tConfId,addSkeleton=False)
      d = self._getShapeShapeDistance(coarseGrid,target.coarseGrid)
      if d>self.shapeDistTol:
        print 'coarsePrune'
        removeIt=True
      else:
        builder.gridSpacing=oSpace*1.5
        medGrid=builder.GenerateSubshapeShape(queryMol,tConfId,addSkeleton=False)
        d = self._getShapeShapeDistance(medGrid,target.medGrid)
        if d>self.shapeDistTol:
          print 'medPrune'
          removeIt=True
        else:
          builder.gridSpacing=oSpace
          fineGrid=builder.GenerateSubshapeShape(queryMol,tConfId,addSkeleton=False)
          d = self._getShapeShapeDistance(fineGrid,target)
          if d>self.shapeDistTol:
            print 'finePrune'
            removeIt=True
      builder.gridSpacing=oSpace
      if removeIt:
        del alignments[i]
      else:
        i+=1

  def GetSubshapeAlignments(self,targetMol,target,queryMol,query,builder,tgtConf=-1,queryConf=-1):
    res = self.GetTriangleMatches(target,query)
    self.PruneMatchesUsingFeatures(target,query,res)
    self.PruneMatchesUsingDirection(target,query,res)
    self.PruneMatchesUsingShape(targetMol,target,queryMol,query,builder,res,
                                tgtConf=tgtConf,queryConf=queryConf)
    return res


if __name__=='__main__':
  import cPickle
  tgtMol,tgtShape = cPickle.load(file('target.pkl','rb'))
  queryMol,queryShape = cPickle.load(file('query.pkl','rb'))
  builder = cPickle.load(file('builder.pkl','rb'))
  aligner = SubshapeAligner()
  algs = aligner.GetSubshapeAlignments(tgtMol,tgtShape,queryMol,queryShape,builder)
  print len(algs)

  from Chem.PyMol import MolViewer
  v = MolViewer()
  v.ShowMol(tgtMol,name='Target',showOnly=True)
  v.ShowMol(queryMol,name='Query',showOnly=False)
  SubshapeObjects.DisplaySubshape(v,tgtShape,'target_shape',color=(.8,.2,.2))
  SubshapeObjects.DisplaySubshape(v,queryShape,'query_shape',color=(.2,.2,.8))
  for i,alg in enumerate(algs):
    confId = alg.alignedConfId
    v.ShowMol(queryMol,confId=confId,name='align-%d'%i,showOnly=False)
      
