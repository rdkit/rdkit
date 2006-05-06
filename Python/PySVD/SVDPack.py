#
#  Copyright (C) 2004 Rational Discovery LLC
#    All Rights Reserved
#
import RDConfig
import os
from Numeric import *
import sys

def HBFOut(calc,outF):
  vects = calc._vects
  nRows = len(vects)
  nCols = len(calc._idMap)
  vals = calc._vals
  nVals = len(vals)
  outF.write('% 72s%s\n'%('svdpack input','rd'))
  outF.write('blank\n')
  outF.write('rra %d %d %d 0\n'%(nCols,nRows,nVals))
  outF.write('X X X X\n')
  pos = 1
  for i in range(nRows):
    outF.write('%d '%pos)
    pos += len(vects[i])
    if not (i+1)%20: outF.write('\n')
  outF.write('%d\n'%(nVals+1))

  for i in range(nRows):
    row = vects[i]
    nHere = len(row)
    for entry in row:
      outF.write('%d '%(entry+1))
    outF.write('\n')
  for i in range(nVals):
    outF.write('%d '%vals[i])
    if not (i+1)%20:
      outF.write('\n')

def ReadSVDPackOutputFile(inF):
  import re
  inL = inF.readline()
  while inL and inL.find('(ROWS)')==-1:
    inL = inF.readline()
  if inL.index('(ROWS)')<0:
    raise ValueError('input file contains no (ROWS) line')
  txt = (inL.split('=')[-1]).strip()
  nCols = int(txt)
  inL = inF.readline()
  if not inL:
    raise ValueError('premature EOF hit')
  txt = (inL.split('=')[-1]).strip()
  nRows = int(txt)

  while inL and inL.find('COMPUTED SINGULAR VALUES')==-1:
    inL = inF.readline()
  if inL.index('COMPUTED SINGULAR VALUES')<0:
    raise ValueError('input file contains no singular vals')
  inL = inF.readline()
  inL = inF.readline()
  blankExpr = re.compile(' +')
  vals = []
  inL = inL.strip()
  while inL and inL[0]=='.':
    splitL = blankExpr.split(inL)
    idx = int(splitL[1])
    val = float(splitL[2])
    vals.append(val)
    inL = inF.readline()
    inL = inL.strip()
  return nRows,nCols,vals
  
def ReadSVDPackArrays(inF,nRows,nCols,k):
  import struct
  T = []
  D = []
  for i in range(k):
    d = inF.read(nCols*8)
    if not d:
      raise ValueError('premature EOF hit')
    tmp = struct.unpack('%dd'%nCols,d)
    T.append(tmp)

    d = inF.read(nRows*8)
    if not d:
      raise ValueError('premature EOF hit')
    tmp = struct.unpack('%dd'%nRows,d)
    D.append(tmp)
  T = transpose(array(T))
  D = transpose(array(D))
  return T,D

def DoSVD(calc,k,exe=None,tol=1e-8):
  if exe is None:
    exe = os.path.join(RDConfig.RDBinDir,'sis2-rd')
  parmFilename = "sis2-parms"
  matFilename = "sis2-mat"
  
  nRows = len(calc._vects)
  nCols = len(calc._idMap)

  nExtras = max(2,.1*k)
  if (k+nExtras)>nRows:
    nExtras = nRows-k
  maxIts = max(nRows,nCols)  
  parmText="'AutoGen' %d %d %d %g TRUE\n"%(k,nExtras,maxIts,tol)
  open(parmFilename,'w+').write(parmText)
  matF=open(matFilename,'w+')
  HBFOut(calc,matF)
  matF=None
  res = os.spawnl(os.P_WAIT,exe,exe,parmFilename,matFilename)
  if not res:
    nRows,nCols,vals = ReadSVDPackOutputFile(open('sio2','r'))
    T,D = ReadSVDPackArrays(open('siv2','rb'),nRows,nCols,k)

    # update the lengths:
    while vals[-1]<=tol:
      vals.pop(-1)
    k = len(vals)
    T = T[:,:k]
    D = D[:,:k]
    calc.ForceSingularValues(k,T,D,array(vals))
    
    

if __name__ == '__main__':
  import SVDSimilarity
  m = [[0,1,2],
       [2,3,4,5,6,8],
       [1,3,4,7],
       [0,4,4,7],
       [3,5,6],
       [9],
       [9,10],
       [9,10,11],
       [8,10,11]]
  calc = SVDSimilarity.SimilarityCalculator()
  calc.SetVects(m)
  DoSVD(calc,2)
  print calc.ScorePoint(m[0])

