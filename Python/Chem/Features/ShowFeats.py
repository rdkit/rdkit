# $Id$
#
# Created by Greg Landrum Aug 2006
#
#
_version = "0.1.0"

_usage="""
   ShowFeats [optional args] <filenames>

Optional arguments:
   -x "list":     specify a comma-separated list of feature types to not be excluded from
                  the visualization

   -f "fdef":
   --fdef="fdef": specify the name of the fdef file.

   --sd:
   --sdf:         expect the input to be one or more SD files (otherwise mol files are expected)

   --dirs:        do not calculate and display feature directions
   --heads:       do not put heads on the feature-direction arrows (only cylinders will be drawn)
   
"""

_welcomeMessage="This is ShowFeats version %s"%(_version)


import math

import Geometry
from Chem.Features import FeatDirUtilsRD as FeatDirUtils
_featColors = {
  'Donor':(0,1,1),
  'Acceptor':(1,0,1),
  'NegIonizable':(1,0,0),
  'PosIonizable':(0,0,1),
  'ZnBinder':(1,.5,.5),
  'Aromatic':(1,.8,.2),
  'LumpedHydrophobe':(.5,.25,0),
  'Hydrophobe':(.5,.25,0),
  }

def _getVectNormal(v,tol=1e-4):
  if math.fabs(v.x)>tol:
    res = Geometry.Point3D(v.y,-v.x,0)
  elif math.fabs(v.y)>tol:
    res = Geometry.Point3D(-v.y,v.x,0)
  elif math.fabs(v.z)>tol:
    res = Geometry.Point3D(1,0,0)
  else:
    raise ValueError,'cannot find normal to the null vector'
  res.Normalize()
  return res

_canonArrowhead=None
def _buildCanonArrowhead(headFrac,nSteps,aspect):
  global _canonArrowhead
  startP = RDGeometry.Point3D(0,0,headFrac)
  _canonArrowhead=[startP]

  scale = headFrac*aspect
  baseV = RDGeometry.Point3D(scale,0,0)
  _canonArrowhead.append(baseV)

  twopi = 2*math.pi
  for i in range(1,nSteps):
    v = RDGeometry.Point3D(scale*math.cos(i*twopi),scale*math.sin(i*twopi),0)
    _canonArrowhead.append(v)
  
    
_globalArrowCGO=[]
_globalSphereCGO=[]
# taken from pymol's cgo.py
BEGIN=2
END=3
TRIANGLE_FAN=6
COLOR=6
VERTEX=4
NORMAL=5
SPHERE=7
CYLINDER=9
ALPHA=25

def _cgoArrowhead(viewer,tail,head,radius,color,label,headFrac=0.3,nSteps=10,aspect=.5):
  global _globalArrowCGO
  delta = head-tail
  normal = _getVectNormal(delta)
  delta.Normalize()
  
  dv = head-tail
  dv.Normalize()
  dv *= headFrac
  startP = head

  normal*=headFrac*aspect
  
  cgo = [BEGIN,TRIANGLE_FAN,
         COLOR,color[0],color[1],color[2],
          NORMAL,dv.x,dv.y,dv.z,
         VERTEX,head.x+dv.x,head.y+dv.y,head.z+dv.z]
  base = [BEGIN,TRIANGLE_FAN,
          COLOR,color[0],color[1],color[2],
          NORMAL,-dv.x,-dv.y,-dv.z,
          VERTEX,head.x,head.y,head.z]
  v = startP+normal
  cgo.extend([NORMAL,normal.x,normal.y,normal.z])
  cgo.extend([VERTEX,v.x,v.y,v.z])
  base.extend([VERTEX,v.x,v.y,v.z])
  for i in range(1,nSteps):
    v = FeatDirUtils.ArbAxisRotation(360./nSteps*i,delta,normal)
    cgo.extend([NORMAL,v.x,v.y,v.z])
    v += startP
    cgo.extend([VERTEX,v.x,v.y,v.z])
    base.extend([VERTEX,v.x,v.y,v.z])

  cgo.extend([NORMAL,normal.x,normal.y,normal.z])
  cgo.extend([VERTEX,startP.x+normal.x,startP.y+normal.y,startP.z+normal.z])
  base.extend([VERTEX,startP.x+normal.x,startP.y+normal.y,startP.z+normal.z])
  cgo.append(END)
  base.append(END)
  cgo.extend(base)
  
  #viewer.server.renderCGO(cgo,label)
  _globalArrowCGO.extend(cgo)
  
def ShowArrow(viewer,tail,head,radius,color,label,transparency=0,includeArrowhead=True):
  global _globalArrowCGO
  if transparency:
    _globalArrowCGO.extend([ALPHA,1-transparency])
  else:
    _globalArrowCGO.extend([ALPHA,1])
  _globalArrowCGO.extend([CYLINDER,tail.x,tail.y,tail.z,
                          head.x,head.y,head.z,
                          radius*.10,
                          color[0],color[1],color[2],
                          color[0],color[1],color[2],
                          ])
    
  if includeArrowhead:
    _cgoArrowhead(viewer,tail,head,radius,color,label)
  

def ShowMolFeats(mol,factory,viewer,radius=0.5,confId=-1,showOnly=True,
                 name='',transparency=0.0,colors=None,excludeTypes=[],
                 useFeatDirs=True,featLabel=None,dirLabel=None,includeArrowheads=True):
  global _globalSphereCGO
  if not name:
    if mol.HasProp('_Name'):
      name =  mol.GetProp('_Name')
    else:
      name = 'molecule'
  if not colors:
    colors = _featColors

  viewer.ShowMol(mol,name=name,showOnly=showOnly,confId=confId)

  molFeats=factory.GetFeaturesForMol(mol)
  if not featLabel:
    featLabel='%s-feats'%name
    viewer.server.resetCGO(featLabel)
  if not dirLabel:
    dirLabel=featLabel+"-dirs"
    viewer.server.resetCGO(dirLabel)
  queueIt=hasattr(viewer,'AddSpheres')
  for i,feat in enumerate(molFeats):
    family=feat.GetFamily()
    if family in excludeTypes:
      continue
    pos = feat.GetPos(confId)
    color = colors.get(family,(.5,.5,.5))
    nm = '%s(%d)'%(family,i+1)
    #if transparency:
    #  viewer.server.sphere((pos.x,pos.y,pos.z),radius,color,featLabel,1,1,transparency)
    #else:
    #  viewer.server.sphere((pos.x,pos.y,pos.z),radius,color,featLabel,1,0)

    if transparency:
      _globalSphereCGO.extend([ALPHA,1-transparency])
    else:
      _globalSphereCGO.extend([ALPHA,1])
    _globalSphereCGO.extend([COLOR,color[0],color[1],color[2],
                             SPHERE,pos.x,pos.y,pos.z,
                             radius])

    if useFeatDirs:
      ps = []
      if family=='Aromatic':
        ps,fType = FeatDirUtils.GetAromaticFeatVects(mol.GetConformer(confId),
                                                               feat.GetAtomIds(),pos,
                                                               scale=1.0)
      elif family=='Donor':
        aids = feat.GetAtomIds()
        if len(aids)==1:
          if len(mol.GetAtomWithIdx(aids[0]).GetNeighbors())==1:
            ps,fType = FeatDirUtils.GetDonor1FeatVects(mol.GetConformer(confId),
                                                       aids,scale=1.0)
          elif len(mol.GetAtomWithIdx(aids[0]).GetNeighbors())==2:
            ps,fType = FeatDirUtils.GetDonor2FeatVects(mol.GetConformer(confId),
                                                       aids,scale=1.0)
          elif len(mol.GetAtomWithIdx(aids[0]).GetNeighbors())==3:
            ps,fType = FeatDirUtils.GetDonor3FeatVects(mol.GetConformer(confId),
                                                       aids,scale=1.0)
      elif family=='Acceptor':
        aids = feat.GetAtomIds()
        if len(aids)==1:
          if len(mol.GetAtomWithIdx(aids[0]).GetNeighbors())==1:
            ps,fType = FeatDirUtils.GetAcceptor1FeatVects(mol.GetConformer(confId),
                                                          aids,scale=1.0)
          elif len(mol.GetAtomWithIdx(aids[0]).GetNeighbors())==2:
            ps,fType = FeatDirUtils.GetAcceptor2FeatVects(mol.GetConformer(confId),
                                                       aids,scale=1.0)
          elif len(mol.GetAtomWithIdx(aids[0]).GetNeighbors())==3:
            ps,fType = FeatDirUtils.GetAcceptor3FeatVects(mol.GetConformer(confId),
                                                       aids,scale=1.0)

            

      for tail,head in ps:
        ShowArrow(viewer,tail,head,radius,color,dirLabel,
                  transparency=transparency,includeArrowhead=includeArrowheads)


if __name__=='__main__':
  def Usage():
    print >>sys.stderr,_usage

  import sys,os,getopt
  import RDConfig
  import Chem
  from Chem import AllChem
  from Chem.PyMol import MolViewer

  try:
    args,extras = getopt.getopt(sys.argv[1:],'x:f:',
                                ('sdf','sd','fdef=','heads','dirs'))
  except:
    Usage()
    sys.exit(1)
    
  exclude=[]
  fdef='BaseFeatures.fdef'
  molFormat='mol'
  includeArrowheads=True
  useDirs=True
  for arg,val in args:
    if arg=='-x':
      val = val.replace('[','').replace('(','').replace(']','').replace(')','')
      for ex in val.split(','):
        exclude.append(ex)
    elif arg in ('-f','--fdef'):
      fdef = val
    elif arg in ('--sd','--sdf'):
      molFormat='sdf'
    elif arg in ('--heads',):
      includeArrowheads=False
    elif arg in ('--dirs',):
      useDirs=False
      
  print >>sys.stderr,_welcomeMessage

  try:
    v = MolViewer()
  except:
    print >>sys.stderr,'ERROR: Unable to connect to PyMol server.\nPlease run ~landrgr1/extern/PyMol/launch.sh to start it.'
    sys.exit(1)
  v.DeleteAll()
    

  
  try:
    fdef = file(fdef,'r').read()
  except IOError:
    print >>sys.stderr,'ERROR: Could not open fdef file %s'%fdef
    sys.exit(1)
    
  factory = AllChem.BuildFeatureFactoryFromString(fdef)

  i = 1
  for molN in extras:
    featLabel='%s_Feats'%molN
    v.server.resetCGO(featLabel)
    # this is a big of kludgery to work around what seem to be pymol cgo bug:
    v.server.sphere((0,0,0),.01,(1,0,1),featLabel)
    dirLabel=featLabel+"-dirs"
    if useDirs:
      v.server.resetCGO(dirLabel)
      # this is a big of kludgery to work around what seem to be pymol cgo bug:
      v.server.cylinder((0,0,0),(.01,.01,.01),.01,(1,0,1),dirLabel)

    if molFormat=='mol':
      m = Chem.MolFromMolFile(molN)
      if not m:
        print >>sys.stderr,'Could not parse molecule: %s'%molN
        ms = []
      else:
        ms = [m]
    elif molFormat=='sdf':
      ms = Chem.SDMolSupplier(molN)
    else:
      ms = []

    for m in ms:
      nm = 'Mol_%d'%(i)
      if m.HasProp('_Name'):
        nm += '_'+m.GetProp('_Name')
      ShowMolFeats(m,factory,v,transparency=0.25,excludeTypes=exclude,name=nm,showOnly=False,
                   useFeatDirs=useDirs,
                   featLabel=featLabel,dirLabel=dirLabel,
                   includeArrowheads=includeArrowheads)
      i += 1
    if ms:
      v.server.renderCGO(_globalSphereCGO,featLabel,1)
      if useDirs:
        v.server.renderCGO(_globalArrowCGO,dirLabel,1)
  sys.exit(0)
