# $Id$
#
# Copyright (C) 2007 by Greg Landrum 
#  All rights reserved
#
import Geometry
from Chem.Subshape import SubshapeObjects
import math
import Numeric
import LinearAlgebra


#-----------------------------------------------------------------------------
def ComputeGridIndices(shapeGrid,winRad):
  if getattr(shapeGrid,'_indicesInSphere',None):
    return shapeGrid._indicesInSphere
  gridSpacing = shapeGrid.GetSpacing()
  dX = shapeGrid.GetNumX()
  dY = shapeGrid.GetNumY()
  dZ = shapeGrid.GetNumZ()
  radInGrid = int(winRad/gridSpacing)
  indicesInSphere=[]
  for i in range(-radInGrid,radInGrid+1):
    for j in range(-radInGrid,radInGrid+1):
      for k in range(-radInGrid,radInGrid+1):
        d=int(math.sqrt(i*i + j*j + k*k ))
        if d<=radInGrid:
          idx = (i*dY+j)*dX+k
          indicesInSphere.append(idx)
  shapeGrid._indicesInSphere = indicesInSphere
  return indicesInSphere

#-----------------------------------------------------------------------------
def ComputeShapeGridCentroid(pt,shapeGrid,winRad):
  count=0
  centroid = Geometry.Point3D(0,0,0)
  idxI = shapeGrid.GetGridPointIndex(pt)
  shapeGridVect = shapeGrid.GetOccupancyVect()
  indicesInSphere = ComputeGridIndices(shapeGrid,winRad)
  for idxJ in indicesInSphere:
    idx = idxI+idxJ;
    if idx>=0 and idx<len(shapeGridVect):
      wt = shapeGridVect[idx]
      tPt = shapeGrid.GetGridPointLoc(idx)
      centroid += tPt*wt
      count+=wt
  if not count:
    raise ValueError,'found no weight in sphere'
  centroid /= count
  #print 'csgc:','(%2f,%2f,%2f)'%tuple(pt),'(%2f,%2f,%2f)'%tuple(centroid),count
  return count,centroid
      


#-----------------------------------------------------------------------------
def FindTerminalPtsFromShape(shape,winRad,fraction,maxGridVal=3):
  termPts = []
  shapeGrid=shape.grid
  shapeVect = shapeGrid.GetOccupancyVect()
  nGridPts = len(shapeVect)
  indicesInSphere = ComputeGridIndices(shapeGrid,winRad)

  for i in range(nGridPts):
    if shapeVect[i]<maxGridVal:
      continue
    # compute the volume inside the sphere:
    volInSphere=0
    nPtsHere=0
    for idx in indicesInSphere:
      idx += i
      if idx>=0 and idx<nGridPts:
        volInSphere += shapeVect[idx]
        nPtsHere +=1
    # the shape may be cut off by the edge of the grid, so
    # the actual max volume in the sphere may well be less
    # than the theoretical max:
    maxVolInSphere=maxGridVal*nPtsHere
    fracVol = float(volInSphere)/maxVolInSphere
    if fracVol<fraction:
      # the sphere is towards the edge of the shape, add a skel pt:
      posI = shapeGrid.GetGridPointLoc(i)
      count,centroid = ComputeShapeGridCentroid(posI,shapeGrid,winRad)
      termPts.append(SubshapeObjects.SkeletonPoint(location=centroid))
  return termPts



#-----------------------------------------------------------------------------
def FindTerminalPtsFromConformer(conf,winRad,nbrCount):
  mol = conf.GetOwningMol()
  nAts = conf.GetNumAtoms()
  nbrLists=[[] for x in range(nAts)]
  for i in range(nAts):
    if(mol.GetAtomWithIdx(i).GetAtomicNum()<=1): continue
    pi = conf.GetAtomPosition(i)
    nbrLists[i].append((i,pi))
    for j in range(i+1,nAts):
      if(mol.GetAtomWithIdx(j).GetAtomicNum()<=1): continue
      pj = conf.GetAtomPosition(j)
      dist = pi.Distance(conf.GetAtomPosition(j))
      if dist<winRad:
        nbrLists[i].append((j,pj))
        nbrLists[j].append((i,pi))
  termPts=[]
  for i in range(nAts):
    if not nbrLists[i]: continue
    pos = Geometry.Point3D(0,0,0)
    totWt=0.0
    if len(nbrLists[i])<nbrCount:
      nbrList = nbrLists[i]
      for j in range(0,len(nbrList)):
        nbrJ,posJ=nbrList[j]
        weight = 1.*len(nbrLists[i])/len(nbrLists[nbrJ])
        pos += posJ*weight
        totWt+=weight
      pos /= totWt
      termPts.append(SubshapeObjects.SkeletonPoint(location=pos))
  return termPts

#-----------------------------------------------------------------------------
def ClusterTerminalPts(pts,winRad,scale):
  res = []
  tagged = [(y,x) for x,y in enumerate(pts)]
  while tagged:
    head,headIdx = tagged.pop(0)
    currSet = [head]
    i=0
    while i<len(tagged):
      nbr,nbrIdx=tagged[i]
      if head.location.Distance(nbr.location)<scale*winRad:
        currSet.append(nbr)
        del tagged[i]
      else:
        i+=1
    pt = Geometry.Point3D(0,0,0)
    for o in currSet:
      pt += o.location
    pt /= len(currSet)
    res.append(SubshapeObjects.SkeletonPoint(location=pt))
  return res

   
#-----------------------------------------------------------------------------
def AppendSkeletonPoints(shapeGrid,termPts,winRad,stepDist,maxGridVal=3,
                         maxDistC=15.0,distTol=1.5,symFactor=1.5):
  nTermPts = len(termPts)
  skelPts=[]
  shapeVect = shapeGrid.GetOccupancyVect()
  nGridPts = len(shapeVect)
  # find all possible skeleton points
  print 'generate all possible'
  for i in range(nGridPts):
    if shapeVect[i]<maxGridVal:
      continue
    posI = shapeGrid.GetGridPointLoc(i)
    ok=True
    for pt in termPts:
      dst = posI.Distance(pt.location)
      if dst<stepDist:
        ok=False
        break
    if ok:
      skelPts.append(SubshapeObjects.SkeletonPoint(location=posI))
  # now start removing them
  print 'Compute centroids:',len(skelPts)
  gridBoxVolume=shapeGrid.GetSpacing()**3
  maxVol = 4.0*math.pi/3.0 * winRad**3 * maxGridVal / gridBoxVolume
  i=0
  while i<len(skelPts):
    pt = skelPts[i]
    count,centroid=ComputeShapeGridCentroid(pt.location,shapeGrid,winRad)
    centroidPtDist=centroid.Distance(pt.location)
    if centroidPtDist>maxDistC:
      del skelPts[i]
    else:
      pt.fracVol = float(count)/maxVol
      pt.location.x = centroid.x
      pt.location.y = centroid.y
      pt.location.z = centroid.z
    i+=1

  print 'remove points:',len(skelPts)
  res = termPts+skelPts
  i=0
  while i<len(res):
    p=-1
    mFrac=0.0
    ptI = res[i]

    startJ = max(i+1,nTermPts)
    for j in range(startJ,len(res)):
      ptJ=res[j]
      distC = ptI.location.Distance(ptJ.location)
      if distC<symFactor*stepDist:
        if ptJ.fracVol>mFrac:
          p=j
          mFrac=ptJ.fracVol
    #print i,len(res),p,mFrac
    if p>-1:
      j = startJ
      while j < len(res):
        ptJ=res[j]
        distC = ptI.location.Distance(ptJ.location)
        if distC<symFactor*stepDist and j!=p:
          del res[j]
        else:
          j+=1
    i+=1
  return res

#-----------------------------------------------------------------------------
def CalculateDirectionsAtPoint(pt,shapeGrid,winRad):
  shapeGridVect = shapeGrid.GetOccupancyVect()
  tmp = winRad/shapeGrid.GetSpacing()
  radInGrid=int(tmp)
  radInGrid2=int(tmp*tmp)
  covMat = Numeric.zeros((3,3),Numeric.Float)

  dX = shapeGrid.GetNumX()
  dY = shapeGrid.GetNumY()
  dZ = shapeGrid.GetNumZ()
  idx = shapeGrid.GetGridPointIndex(pt.location)
  idxZ = idx//(dX*dY)
  rem = idx%(dX*dY)
  idxY = rem//dX
  idxX = rem%dX
  totWt=0.0
  for i in range(-radInGrid,radInGrid):
    xi = idxX+i
    for j in range(-radInGrid,radInGrid):
      xj = idxY+j
      for k in range(-radInGrid,radInGrid):
        xk = idxZ+k
        d2 = i*i+j*j+k*k
        if d2>radInGrid2 and int(math.sqrt(d2))>radInGrid:
          continue
        gridIdx = (xk*dY+xj)*dX+xi
        if gridIdx>=0 and gridIdx<len(shapeGridVect):
          wtHere = shapeGridVect[gridIdx]
          totWt += wtHere
          ptInSphere = shapeGrid.GetGridPointLoc(gridIdx)
          ptInSphere -= pt.location
          covMat[0][0]+= wtHere*(ptInSphere.x*ptInSphere.x)
          covMat[0][1]+= wtHere*(ptInSphere.x*ptInSphere.y)
          covMat[0][2]+= wtHere*(ptInSphere.x*ptInSphere.z)
          covMat[1][1]+= wtHere*(ptInSphere.y*ptInSphere.y)
          covMat[1][2]+= wtHere*(ptInSphere.y*ptInSphere.z)
          covMat[2][2]+= wtHere*(ptInSphere.z*ptInSphere.z)
  covMat /= totWt
  covMat[1][0] = covMat[0][1]
  covMat[2][0] = covMat[0][2]
  covMat[2][1] = covMat[1][2]

  eVals,eVects = LinearAlgebra.eigenvectors(covMat)
  sv = [(x,y) for x,y in zip(eVals,eVects)]
  sv.sort(reverse=True)
  eVals = [x for x,y in sv]
  eVects = [y for x,y in sv]
  pt.shapeMoments=tuple(eVals)
  pt.shapeDirs = tuple([Geometry.Point3D(p[0],p[1],p[2]) for p in eVects])

  #print '-------------'
  #print pt.location.x,pt.location.y,pt.location.z
  #for v in covMat:
  # print '  ',v
  #print '---'
  #print eVals
  #for v in eVects:
  #  print '  ',v
  #print '-------------'
  

#-----------------------------------------------------------------------------
def AssignMolFeatsToPoints(pts,mol,featFactory,winRad):
  feats = featFactory.GetFeaturesForMol(mol)
  for i,pt in enumerate(pts):
    for feat in feats:
      if feat.GetPos().Distance(pt.location)<winRad:
        typ = feat.GetFamily()
        if typ not in pt.molFeatures:
          pt.molFeatures.append(typ)
    print i,pt.molFeatures

