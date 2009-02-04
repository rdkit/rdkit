# $Id$
#
#  Copyright (C) 2005  Rational Discovery LLC
#    All Rights Reserved
#
""" for the moment this is using Francois Fleuret's cmim library
 to do the feature selection

 Reference: F. Fleuret "Fast Binary Feature Selection with Conditional
            Mutual Information", J. Machine Learn. Res. 5, 1531-1535
            (2004)
  
"""
from rdkit import RDConfig
from rdkit import DataStructs
import tempfile
import os
import rdFeatSelect

def SelectFeatures(examples,nFeatsToPick,bvCol=1):
  res = rdFeatSelect.selectCMIM(examples,nFeatsToPick)
  if -1 in res:
    res = list(res)
    res = tuple(res[:res.index(-1)])
  return res

def _SelectFeatures(examples,nFeatsToPick,bvCol=1):
  nPts = len(examples)
  nFeats = examples[0][bvCol].GetNumBits()

  exe = os.path.join(RDConfig.RDBaseDir,'External','cmim-1.0','cmim.exe')
  if not os.path.exists(exe):
    raise ValueError,'could not find cmim executable %s'%exe
  
  inFname = tempfile.mktemp('.dat')
  outFname = inFname + '.out'
  inF = open(inFname,'w+')
  print >>inF,nPts,nFeats
  for row in examples:
    print >>inF,row[bvCol].ToBitString()
    print >>inF,row[-1]
  inF.close()
  inF = None

  os.spawnlp(os.P_WAIT,exe,exe,'--nb-features',str(nFeatsToPick),'--train',
            inFname,outFname)

  inD = open(outFname,'r')
  inL = inD.readline()
  nCreated = int(inL)
  inL = inD.readline()
  res = []
  splitL = inL.split(' ')
  for i in range(nFeatsToPick):
    res.append(int(splitL[i]))
  inD.close()
  inD = None
  
  os.unlink(inFname)
  os.unlink(outFname)

  return res
