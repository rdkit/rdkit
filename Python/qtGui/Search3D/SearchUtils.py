## Automatically adapted for numpy.oldnumeric Sep 23, 2006 by alter_code1.py

# $Id$
#
#  Copyright (C) 2005,2006 Rational Discovery LLC
#   All Rights Reserved
#
from Numeric import *
import sys
from qtGui.qtUtils import logger
import Chem
from Chem import ChemicalForceFields,rdDistGeom,rdShapeHelpers
from Chem.Pharm3D import EmbedLib
EmbedLib.logger = logger
import DistanceGeometry as DG
import numpy.oldnumeric as Numerics.rdAlignment as Aligner
import DataStructs
import cPickle
from qtGui.Search3D import LocalConfig
colors=((1,0,1),(0,.8,.8),(0.2,0.2,1),(1,0,0),(.9,.9,.1),(.1,.7,.1),(.6,.6,.6),
         (.7,.1,.1), (.1,.1,.7), (.3,.3,.3), )

ZEROTOL=1e-3

def GetFeatToPointDistanceRange(mol,feat,pos):
  pos = array(pos)
  atomIds = feat.GetAtomIds()
  if len(atomIds)==1:
    coords = array(feat.GetPos())
    v = pos-coords
    d = sqrt(dot(v,v))
    minD = d
    maxD = d
  else:
    conf = mol.GetConformer()
    minD = 1e8
    maxD = -1e8
    for id in atomIds:
      coords = array(conf.GetAtomPosition(id))
      v = pos-coords
      d = sqrt(dot(v,v))
      minD = min(minD,d)
      maxD = max(maxD,d)
  return minD,maxD

def GetFeatsPerAtom(feats):
  """  Returns a dictionary keyed by atom id, with lists of features as
   the values

  """
  res = {}
  for feat in feats:
    for atomId in feat.GetAtomIds():
      if res.has_key(atomId):
        res[atomId].append(feat)
      else:
        res[atomId] = [feat]
  return res


def GetFeatDistances(mol,feats):
  """  Returns the 3D and 2D distances between a set of features as
    a 2-tuple of lists
    
  """
  dm = Chem.GetDistanceMatrix(mol,False,False)
  res = []
  res2D = []
  nFeats = len(feats)
  for i in range(nFeats):
    featI = feats[i]
    posI = array(featI.GetPos())
    atomsI = list(featI.GetAtomIds())
    for j in range(i+1,nFeats):
      featJ = feats[j]
      posJ = array(featJ.GetPos())
      atomsJ = list(featJ.GetAtomIds())
      delta = posI-posJ
      dist = sqrt(dot(delta,delta))
      dist2D = 0
      for atomI in atomsI:
        for atomJ in atomsJ:
          dist2D = max(dist2D,dm[atomI,atomJ])
      res.append(dist)
      res2D.append(dist2D)
  return res,res2D

def InitializeTargets(supplier,addHs=False,progressCallback=None,maxToGet=-1,
                      pickledBoundsName='BOUNDSMATRIX',res=None,removeHs=True):
  """ Pulls in a set of molecules from a MolSupplier and returns a list
   of (mol,boundsMatrix) 2-tuples.
   
  """
  if res is None:
    res = []
  nDone = 0
  for mol in supplier:
    if not mol: continue
    try:
      Chem.SanitizeMol(mol)
      if removeHs:
        mol = Chem.RemoveHs(mol)
      if addHs:
        mol = Chem.AddHs(mol)
    except ValueError:
      txt = 'problems sanitizing molecule '
      if mol.HasProp('_Name'):
        txt = txt + mol.GetProp('_Name')
        
      logger.debug(txt,exc_info=True)
    else:
      # compute the molecule's distance matrix now (we'll use it later and
      # the molecule itself will cache the value:
      Chem.GetDistanceMatrix(mol,False,False,True)
      if not pickledBoundsName or not mol.HasProp(pickledBoundsName):
        try:
          bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
        except:
          logger.debug('problems building bounds matrix for molecule %s'%Chem.MolToSmiles(mol),
                       exc_info=True)
        else:
          if DG.DoTriangleSmoothing(bm):
            res.append((mol,bm))
          else:
            #print >>sys.stderr,'problems encountered smoothing bounds matrix for molecule %d'%nDone
            logger.info('Problems encountered smoothing bounds matrix for molecule %d\nSMILES: %s'%(nDone,Chem.MolToSmiles(mol)))
          mol.RemoveAllConformers()
      else:
        pkl = mol.GetProp(pickledBoundsName)
        bm = cPickle.loads(pkl)
        res.append((mol,bm))
    if len(res)==maxToGet:
      break
    nDone += 1
    if progressCallback:
      progressCallback(nDone)
  return res


def GetMolFeatMatches(mol,featFactory,pcophore):
  """ Finds all matches for each of the pharmacophore's features in a molecule.

  Returns a list of lists of features

  """
  canMatch,matches = EmbedLib.MatchPharmacophoreToMol(mol,featFactory,pcophore)
  if not canMatch:
    return []
  return matches

def SearchTargets(targets,featFactory,pcophore,progressCallback=None,
                  excludedVolumes=None,useDirs=False,
                  maxToDo=-1,res=None):
  """ Searches through a list of (mol,bounds) tuples for molecules that can
   match a pharmacophore.

   Returns a list(mol,match,bounds) tuples with matches

  """
  nDone = 0
  if res is None:
    res = []
  for mol,bounds in targets:
    matches = GetMolFeatMatches(mol,featFactory,pcophore)
    if matches:
      r = EmbedLib.MatchPharmacophore(matches,bounds,pcophore,
                                      useDownsampling=True,
                                      excludedVolumes=excludedVolumes,
                                      use2DLimits=True,mol=mol,
                                      useDirs=useDirs)
      failed,bm,match,details = r
      #print match
      #print details
      #print '-----'
      if not failed:
        res.append((mol,match,bounds))
    nDone += 1
    if nDone == maxToDo: break
    if progressCallback:
      progressCallback(nDone)
  return res

def GetMolEmbeddings(mol,bounds,atomMatch,pcophore,numDesired=5,
                     maxToTry=20,excludedVolumes=None,randomSeed=-1,
                     useDirs=False,
                     twoStageOptimization=False,
                     progressCallback=None):
  """ Finds a set of optimized embeddings for a molecule to a pharmacophore

  returns a list of 2-tuples: (E2, embedded molecule)

  """
  res = []
  try:
    bm,embeds,bar = EmbedLib.EmbedPharmacophore(mol,atomMatch,pcophore,count=maxToTry,
                                                smoothFirst=False,silent=True,bounds=bounds,
                                                excludedVolumes=excludedVolumes,
                                                randomSeed=randomSeed,
                                                targetNumber=numDesired,
                                                useDirs=useDirs)
  except ValueError:
    embeds = []
  nDone = 0
  for embed in embeds:
    try:
      if twoStageOptimization:
        e1,e2 = EmbedLib.OptimizeMol(embed,bm)
        e3,e2 = EmbedLib.OptimizeMol(embed,bm,atomMatches=atomMatch,
                                     excludedVolumes=excludedVolumes,forceConstant=20.)
      else:
        e1,e2 = EmbedLib.OptimizeMol(embed,bm,atomMatches=atomMatch,
                                     excludedVolumes=excludedVolumes,forceConstant=20.)
    except ValueError,msg:
      logger.info('ValueError: %s'%msg)
    else:
      res.append((e2,embed,excludedVolumes[:]))
    nDone += 1
    if progressCallback: progressCallback(nDone)
  return res  

def GetFeaturesForAtomMatch(mol,atomMatch,pcophore,featFactory):
  """ Given a set of atom indices for a match, returns the features
  that correspond to those indices.
  
  returns a list of MolChemicalFeatures

  """
  matches = GetMolFeatMatches(mol,featFactory,pcophore)

  assert len(atomMatch)==len(matches)
  res = []
  for i in range(len(atomMatch)):
    refIds = atomMatch[i]
    refIds.sort()
    foundOne=False
    for match in matches[i]:
      atomIds = list(match.GetAtomIds())
      atomIds.sort()
      if atomIds==refIds:
        res.append(match)
        foundOne=True
        break
    assert foundOne
  return res 


def TransformMol(mol,tform):
  """  Applies the transformation to a molecule and sets it up with
  a single conformer

  """
  newConf = Chem.Conformer()
  newConf.SetId(0)
  refConf = mol.GetConformer()
  for i in range(refConf.GetNumAtoms()):
    pos = list(refConf.GetAtomPosition(i))
    pos.append(1.0)
    newPos = matrixmultiply(tform,array(pos))
    newConf.SetAtomPosition(i,list(newPos)[:3])
  mol.RemoveAllConformers()
  mol.AddConformer(newConf)

def AlignMatchToReference(mol,probeFeats,ref,refFeats,exVols=[],useDirs=False):
  refPts = [list(x.GetPos()) for x in refFeats]

  if useDirs:
    AddFeatureDirections(mol,probeFeats)
    AddFeatureDirections(ref,refFeats)

  probePts = [list(x.GetPos()) for x in probeFeats]

  refDirs = []
  probeDirs = []
  if useDirs:
    for rFeat,pFeat in zip(refFeats,probeFeats):
      if hasattr(rFeat,'dirPts') and hasattr(pFeat,'dirPts') and \
         len(rFeat.dirPts)==len(pFeat.dirPts):
        refDirs.extend([list(x) for x in rFeat.dirPts])
        probeDirs.extend([list(x) for x in pFeat.dirPts])
      else:
        #print 'not using feat dirs between %s feats'%(rFeat.GetFamily())
        pass
  nRefPts = len(refPts)+len(refDirs)
  nProbePts = len(probePts)+len(probeDirs)

  if nProbePts!=nRefPts:
    raise ValueError,'point array length mismatch'

  probeVols = [list(x.pos) for x in exVols]
  refVols = [list(x.origPos) for x in exVols]

  probeArr = array(probePts+probeDirs+probeVols)
  refArr = array(refPts+refDirs+refVols)

  weights = [8.0]*(nProbePts)+[1.0]*len(exVols)

  #print 'probeArr=',repr(probeArr)
  #print 'refArr=',repr(refArr)
  #print 'weights=',repr(weights)

  ssd,tform = Aligner.GetAlignmentTransform(refArr,probeArr,weights=weights)
  # if the molecule has chiral centers, do not try the reflected alignment:
  if not hasattr(mol,'_chiralCenters'):
    centers=EmbedLib.FindMolChiralCenters(mol)
    mol._chiralCenters=centers
  if not mol._chiralCenters:
    ssd2,tform2 = Aligner.GetAlignmentTransform(refArr,probeArr,weights=weights,
                                              reflect=True)
    if ssd2<ssd:
      tform = tform2
      ssd = ssd2
  rms = sqrt(ssd/sum(weights))
  volPs=[]
  for vol in exVols:
    p = list(vol.pos)
    p.append(1.0)
    p = list(matrixmultiply(tform,array(p)))[:-1]
    volPs.append(p)

  if useDirs:
    for i,pFeat in enumerate(probeDirs):
      rFeat = refDirs[i]
      p = list(pFeat)
      p.append(1.0)
      p = list(matrixmultiply(tform,array(p)))[:-1]
  TransformMol(mol,tform)
  tani = rdShapeHelpers.ShapeTanimotoDist(mol,ref,gridSpacing=0.5)
  shapeScore=1.0-tani

  return rms,volPs,shapeScore,tform
  
def cross(v1,v2):
  res = array([ v1[1]*v2[2] - v1[2]*v2[1],
                -v1[0]*v2[2] + v1[2]*v2[0],
                v1[0]*v2[1] - v1[1]*v2[0]],Float)
  return res


def _AddAromaticFeatDirs(mol,feat,dirLen=EmbedLib.defaultFeatLength,confId=-1):
  if feat.GetFamily()!='Aromatic':
    raise ValueError,'bad feature type'
  conf = mol.GetConformer(confId)

  core = array(feat.GetPos())
  atomIds = feat.GetAtomIds()
  a1 = conf.GetAtomPosition(atomIds[0])
  a2 = conf.GetAtomPosition(atomIds[1])

  v1 = a1-core
  v1 = v1 / sqrt(dot(v1,v1))
  
  v2 = a2-core
  v2 = v2 / sqrt(dot(v2,v2))
  
  dir1 = cross(v1,v2)
  dir1 *= dirLen

  # for the purposes of making alignments smoother, make sure
  # that dir1 points into a "positive" direction:
  if dir1[-1]<-ZEROTOL:
    dir1 *=-1
  elif dir1[-1]<ZEROTOL:
    if dir1[0]<-ZEROTOL:
      dir1*=-1
    elif dir1[0]<ZEROTOL:
      if dir1[1]<-ZEROTOL:
        dir1*=-1
  dir2 = dir1*-1

  feat.dirPts = (core+dir1,core+dir2)
  
def _GetDegree1FeatDirs(conf,atom,dirLen=EmbedLib.defaultFeatLength,nbr=None):
  if nbr is None:
    for nbr in atom.GetNeighbors():
      if nbr.GetAtomicNum()!=1:
        break
  core = array(conf.GetAtomPosition(atom.GetIdx()))
  a1 = array(conf.GetAtomPosition(nbr.GetIdx()))
  v1 = core-a1
  v1 = v1/dot(v1,v1)
  v1 *= dirLen

  return core+v1,

def _AddDonorFeatDirs(mol,feat,dirLen=EmbedLib.defaultFeatLength,confId=-1):
  if feat.GetFamily()!='Donor':
    raise ValueError,'bad feature type'

  # at the moment we only handle single-coordinate points:
  atomIds = feat.GetAtomIds()
  if len(atomIds)!=1:
    return
  atom = mol.GetAtomWithIdx(atomIds[0])
  nbrs = EmbedLib.GetAtomHeavyNeighbors(atom)
  if len(nbrs)!=1:
    #print 'no dir for feat %s with degree %d'%(feat.GetFamily(),atom.GetDegree())
    return
  conf = mol.GetConformer(confId)
  feat.dirPts = _GetDegree1FeatDirs(conf,atom,dirLen=dirLen,nbr=nbrs[0])

def _AddAcceptorFeatDirs(mol,feat,dirLen=EmbedLib.defaultFeatLength,confId=-1):
  if feat.GetFamily()!='Acceptor':
    raise ValueError,'bad feature type'
  conf = mol.GetConformer(confId)
  
  # at the moment we only handle single-coordinate points:
  atomIds = feat.GetAtomIds()
  if len(atomIds)!=1:
    return
  atom = mol.GetAtomWithIdx(atomIds[0])
  nbrs = EmbedLib.GetAtomHeavyNeighbors(atom)
  if len(nbrs)!=1:
    #print 'no dir for feat %s with degree %d'%(feat.GetFamily(),atom.GetDegree())
    return
  conf = mol.GetConformer(confId)
  feat.dirPts = _GetDegree1FeatDirs(conf,atom,dirLen=dirLen,nbr=nbrs[0])


def AddFeatureDirections(mol,feats,confId=-1):
  """ the features and their atoms should have coords already
  """
  for feat in feats:
    if feat.GetFamily()=='Aromatic':
      _AddAromaticFeatDirs(mol,feat)
    #elif feat.GetFamily()=='Donor':
    #  _AddDonorFeatDirs(mol,feat)
    #elif feat.GetFamily()=='Acceptor':
    #  _AddAcceptorFeatDirs(mol,feat)
    else:
      feat.dirPts = []

    
def DiversityPick(hits,numToPick,pickledFpName='SIMILARITY_FP',
                  progressCallback=None):
  import SimDivFilters
  import DataStructs
  from Chem.Fingerprints import FingerprintMols
  if not hits or len(hits)<=numToPick:
    return
  if numToPick <= 0:
    return
  
  nHits = len(hits)
  dMat = []
  params = FingerprintMols.FingerprinterDetails()
  for i in range(1,nHits):
    molI = hits[i][0]
    if not hasattr(molI,'simFp') or not molI.simFp:
      if molI.HasProp(pickledFpName):
        fpI = DataStructs.ExplicitBitVect(molI.GetProp(pickledFpName))
      else:
        fpI = apply(FingerprintMols.FingerprintMol,(molI,),params.__dict__)
      molI.simFp = fpI
    else:
      fpI = molI.simFp
    for j in range(i):
      molJ = hits[j][0]
      if not hasattr(molJ,'simFp') or not molJ.simFp:
        if molJ.HasProp(pickledFpName):
          fpJ = DataStructs.ExplicitBitVect(molJ.GetProp(pickledFpName))
        else:
          fpJ = apply(FingerprintMols.FingerprintMol,(molJ,),params.__dict__)
        molJ.simFp = fpJ
      else:
        fpJ = molJ.simFp
      dMat.append(1.0-DataStructs.FingerprintSimilarity(fpI,fpJ))
      if progressCallback:
        progressCallback(len(dMat))
  dMat = array(dMat)
  picker = SimDivFilters.HierarchicalClusterPicker(SimDivFilters.ClusterMethod.WARD)
  clusters = picker.Pick(dMat,nHits,numToPick)
  logger.debug('clusters: %s'%(str(list(clusters))))
  return [hits[x] for x in clusters]

def GetAlignmentRMSD(refMol,refFeats,probeMol,probeFeats,useDirs=True):
  refPts = [x.GetPos() for x in refFeats]
  probePts = [x.GetPos() for x in probeFeats]
  refDirs = []
  probeDirs = []
  if useDirs:
    AddFeatureDirections(refMol,refFeats)
    AddFeatureDirections(probeMol,probeFeats)
    for rFeat,pFeat in zip(refFeats,probeFeats):
      if hasattr(rFeat,'dirPts') and hasattr(pFeat,'dirPts') and \
         len(rFeat.dirPts)==len(pFeat.dirPts):
        refDirs.extend([x for x in rFeat.dirPts])
        probeDirs.extend([x for x in pFeat.dirPts])
      else:
        pass
  refPts = refPts + refDirs
  probePts = probePts + probeDirs

  if len(probePts)!=len(refPts):
    logger.error('point array length mismatch')
    return

  ssd = 0.0
  for i,probePt in enumerate(probePts):
    refPt = refPts[i]
    d = array(list(probePt))-array(list(refPt))
    ssd += dot(d,d)
  rmsd = sqrt(ssd/len(probePts))
  return rmsd

      
def _MolMolAlignmentRMSD(refConf,probeConf,mapping):
  d2 = 0.0
  nPts = 0
  for probeId,refId in mapping.iteritems():
    refP = refConf.GetAtomPosition(refId)
    probeP = probeConf.GetAtomPosition(probeId)
    refP-=probeP
    d2 += refP.LengthSq()
    nPts +=1
  d2 /= nPts
  return sqrt(d2)

def ImproveMolMolAlignment(refMol,probeMol,mapping,
                           forceConstant=LocalConfig.refineForceConstant,
                           refConfId=-1,probeConfId=-1,maxIts=200,
                           forceTol=1e-4,energyTol=1e-4):
  """ mapping should be a dictionary mapping
    probe atomIds -> reference atomIds
  """
  try:
    ff = ChemicalForceFields.UFFGetMoleculeForceField(probeMol,confId=probeConfId)
  except:
    logger.error('Could not build forcefield for probe molecule',
                 exc_info=True)
    return

  refConf = refMol.GetConformer(refConfId)
  probeConf = probeMol.GetConformer(probeConfId)
  initRMS = _MolMolAlignmentRMSD(refConf,probeConf,mapping)

  refAtomsAdded = {}
  for probeId,refId in mapping.iteritems():
    if not refAtomsAdded.has_key(refId):
      pos = refConf.GetAtomPosition(refId)
      idx =ff.AddExtraPoint(pos.x,pos.y,pos.z)-1
      refAtomsAdded[refId] = idx
    else:
      idx = refAtomsAdded[refId]
    ff.AddDistanceConstraint(probeId,idx,0.0,0.0,forceConstant)
    
  ff.Initialize()
  e1 = ff.CalcEnergy()
  if EmbedLib.isNaN(e1):
    raise ValueError,'bogus energy'
  ff.Minimize(maxIts=maxIts,forceTol=forceTol,energyTol=energyTol)
  e2 = ff.CalcEnergy()
  if EmbedLib.isNaN(e2):
    raise ValueError,'bogus energy'
  finalRMS = _MolMolAlignmentRMSD(refConf,probeConf,mapping)
  
  return (initRMS,e1),(finalRMS,e2)

def GuessMolMolCorrespondence(refMol,probeMol,
                              refConfId=-1,probeConfId=-1,
                              noDupes=False,
                              maxDist=LocalConfig.maxAtomAtomCorrespondenceDistance):
  import bisect
  refConf = refMol.GetConformer(refConfId)
  probeConf = probeMol.GetConformer(probeConfId)

  dist2Tol = maxDist*maxDist
  res = {}

  if not noDupes:
    for probeId in range(probeMol.GetNumAtoms()):
      probePos = probeConf.GetAtomPosition(probeId)
      closest = 1000
      closestId = -1
      for refId in range(refMol.GetNumAtoms()):
        refPos = refConf.GetAtomPosition(refId)
        refPos -= probePos
        d2 = refPos.LengthSq()
        if d2 <= dist2Tol and d2<closest:
          closest = d2
          closestId=refId
      if closestId>-1:
        res[probeId] = closestId
  else:
    # if we need to ensure no dupes, the algorithm is more complex:
    refAtomNbrs = {}
    probeAtomNbrs = {}
    for probeId in range(probeMol.GetNumAtoms()):
      probePos = probeConf.GetAtomPosition(probeId)
      for refId in range(refMol.GetNumAtoms()):
        refPos = refConf.GetAtomPosition(refId)
        refPos -= probePos
        d2 = refPos.LengthSq()
        if d2 <= dist2Tol:
          #print probeId,refId,d2
          refL = refAtomNbrs.get(refId,[])
          probeL = probeAtomNbrs.get(probeId,[])
          bisect.insort(refL,(d2,probeId))
          bisect.insort(probeL,(d2,refId))
          refAtomNbrs[refId]=refL
          probeAtomNbrs[probeId]=probeL
    # refAtomNbrs and probeAtomNbrs now have a list of
    # neighbors for each atomId that is sorted by distance
    # remove the distances:
    for dict in (refAtomNbrs,probeAtomNbrs):
      for k,v in dict.iteritems():
        v = [x[1] for x in v]
        dict[k] = v

    probeIds = probeAtomNbrs.keys()
    changed=True
    while probeIds and changed:
      for probeId in probeIds:
        #print '>>>',probeId,'-',probeIds
        matchId = -1
        seen=False
        for refId in probeAtomNbrs[probeId]:
          nbrL= refAtomNbrs.get(refId,[])
          if not nbrL:
            continue
          #print '\t',refId,nbrL
          if probeId in nbrL:
            seen=True
          if probeId==nbrL[0]:
            # this probe atom is the closest to the reference atom, keep it
            matchId = refId
            del refAtomNbrs[refId]
            seen=True
            break
        if matchId>=0:
          # we found a match.
          res[probeId] = matchId
          probeIds.remove(probeId)
          changed=True
          break
        elif not seen:
          # this atom has no neighbors remaining
          probeIds.remove(probeId)
          changed=True
          break
        elif not len(refAtomNbrs.keys()):
          # nothing more to be found:
          probeIds=[]
          changed=True
          break
        else:
          changed=False
  return res
