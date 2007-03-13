# $Id$
#
# Copyright (C) 2007 by Greg Landrum 
#  All rights reserved
#
import Chem,Geometry
from Chem.Subshape import SubshapeObjects
from Chem.Subshape import BuilderUtils
import time,cPickle

#-----------------------------------------------------------------------------
class SubshapeBuilder(object):
  gridDims=(20,15,10)
  gridSpacing=0.5
  winRad=3.0
  nbrCount=7
  terminalPtRadScale=0.75
  fraction=0.25
  stepSize=1.0
  featFactory=None

  def GenerateSubshapeShape(self,cmpd,confId=-1,addSkeleton=True,**kwargs):
    from Chem import AllChem
    shape = SubshapeObjects.ShapeWithSkeleton()
    shape.grid=Geometry.UniformGrid3D(self.gridDims[0],self.gridDims[1],self.gridDims[2],
                                      self.gridSpacing)
    AllChem.EncodeShape(cmpd,shape.grid,ignoreHs=False,confId=confId)
    conf = cmpd.GetConformer(confId)
    if addSkeleton:
      self.GenerateSubshapeSkeleton(shape,conf,kwargs)
    return shape
  
  def GenerateSubshapeSkeleton(self,shape,conf=None,terminalPtsOnly=False,skelFromConf=True):
    if conf and skelFromConf:
      pts = BuilderUtils.FindTerminalPtsFromConformer(conf,self.winRad,self.nbrCount)
    else:
      pts = BuilderUtils.FindTerminalPtsFromShape(shape,self.winRad,self.fraction)
      
    pts = BuilderUtils.ClusterTerminalPts(pts,self.winRad,self.terminalPtRadScale)
    if not terminalPtsOnly:
      pts = BuilderUtils.AppendSkeletonPoints(shape.grid,pts,self.winRad,self.stepSize)
    for i,pt in enumerate(pts):
      BuilderUtils.CalculateDirectionsAtPoint(pt,shape.grid,self.winRad)
    if conf and self.featFactory:
      BuilderUtils.AssignMolFeatsToPoints(pts,conf.GetOwningMol(),self.featFactory,self.winRad)
    shape.skelPts=pts

      
if __name__=='__main__':
  from Chem import AllChem,ChemicalFeatures
  from Chem.PyMol import MolViewer
  #cmpd = Chem.MolFromSmiles('CCCc1cc(C(=O)O)ccc1')
  #cmpd = Chem.AddHs(cmpd)
  if 1:
    cmpd = Chem.MolFromSmiles('C1=CC=C1C#CC1=CC=C1')
    cmpd = Chem.AddHs(cmpd)
    AllChem.EmbedMolecule(cmpd)
    AllChem.UFFOptimizeMolecule(cmpd)
    AllChem.CanonicalizeMol(cmpd)
    print >>file('testmol.mol','w+'),Chem.MolToMolBlock(cmpd)
  else:
    cmpd = Chem.MolFromMolFile('testmol.mol')
  builder=SubshapeBuilder()
  if 1:
    shape=builder.GenerateSubshapeShape(cmpd)
  v = MolViewer()
  if 1:
    import tempfile
    tmpFile = tempfile.mktemp('.grd')
    v.server.deleteAll()
    Geometry.WriteGridToFile(shape.grid,tmpFile)
    time.sleep(1)
    v.ShowMol(cmpd,name='testMol',showOnly=True)
    v.server.loadSurface(tmpFile,'testGrid','',2.5)
  v.server.resetCGO('*')

  cPickle.dump(shape,file('subshape.pkl','w+'))
  for i,pt in enumerate(shape.skelPts):
    v.server.sphere(tuple(pt.location),.5,(1,0,1),'Pt-%d'%i)
    if not hasattr(pt,'shapeDirs'): continue
    momBeg = pt.location-pt.shapeDirs[0]
    momEnd = pt.location+pt.shapeDirs[0]
    v.server.cylinder(tuple(momBeg),tuple(momEnd),.1,(1,0,1),'v-%d'%i)

