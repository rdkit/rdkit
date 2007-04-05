# $Id$
#
# Copyright (C) 2007 by Greg Landrum 
#  All rights reserved
#
import Chem
import Geometry

class SkeletonPoint(object):
  location=None
  shapeMoments=None
  shapeDirs=None
  molFeatures=None
  featmapFeatures=None
  fracVol=0.0
  
  def __init__(self,*args,**kwargs):
    self._initMemberData()
    self.location = kwargs.get('location',None)

  def _initMemberData(self):
    self.shapeMoments=(0.0,)*3
    self.shapeDirs=[None]*3
    self.molFeatures=[]
    self.featmapFeatures=[]

class ShapeWithSkeleton(object):
  grid=None
  skelPts=None

  def __init__(self,*args,**kwargs):
    self._initMemberData()

  def _initMemberData(self):
    self.skelPts=[]


class SubshapeShape(object):
  shapes=None # a list of ShapeWithSkeleton objects at multiple resolutions (low to high)
  featMap=None
  keyFeat=None
  
  def __init__(self,*args,**kwargs):
    self._initMemberData()

  def _initMemberData(self):
    self.shapes=[]

def DisplaySubshapeSkeleton(viewer,shape,name,color=(1,0,1)):
  cgoNm='%s-skeleton'%name
  viewer.server.resetCGO(cgoNm)
  for i,pt in enumerate(shape.skelPts):
    viewer.server.sphere(tuple(pt.location),.5,color,cgoNm)
    if not hasattr(pt,'shapeDirs'): continue
    momBeg = pt.location-pt.shapeDirs[0]
    momEnd = pt.location+pt.shapeDirs[0]
    viewer.server.cylinder(tuple(momBeg),tuple(momEnd),.1,color,cgoNm)
  
def DisplaySubshape(viewer,shape,name,showSkelPts=True,color=(1,0,1)):
  import Geometry
  import os,tempfile
  fName = tempfile.mktemp('.grd')
  Geometry.WriteGridToFile(shape.grid,fName)
  viewer.server.loadSurface(fName,name,'',2.5)
  if showSkelPts:
    DisplaySubshapeSkeleton(viewer,shape,name,color)
  try:
    os.unlink(fName)
  except:
    import time
    time.sleep(.5)
    try:
      os.unlink(fName)
    except:
      pass
    

  
