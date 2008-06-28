# $Id$
#
#  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
import Chem
from rdfragcatalog import *
import sys
import sets

def message(msg,dest=sys.stdout):
  dest.write(msg)


class BitGainsInfo(object):
  id=-1
  description=''
  gain=0.0
  nPerClass=None
  
def ProcessGainsFile(fileName,nToDo=-1,delim=',',haveDescriptions=1):
  inFile = open(fileName,'r')
  nRead = 0
  res = []
  for line in inFile.xreadlines():
    nRead += 1
    splitL = [x.strip() for x in line.split(delim)]
    if nRead != 1 and len(splitL):
      bit = BitGainsInfo()
      bit.id = int(splitL[0])
      col = 1
      if haveDescriptions:
        bit.description = splitL[col]
        col += 1
      bit.gain = float(splitL[col])
      col += 1
      nPerClass = []
      for entry in splitL[col:]:
        nPerClass.append(int(entry))
      bit.nPerClass = nPerClass  
      res.append(bit)
      if len(res)==nToDo:
        break
  return res   

def BuildAdjacencyList(catalog,bits,limitInclusion=1,orderLevels=0):
  adjs = {}
  levels = {}
  bitIds = [bit.id for bit in bits]
  for bitId in bitIds:
    entry = catalog.GetBitEntryId(bitId)
    tmp = []
    order = catalog.GetEntryOrder(entry)
    s = levels.get(order,sets.Set())
    s.add(bitId)
    levels[order] = s
    for down in catalog.GetEntryDownIds(entry):
      id = catalog.GetEntryBitId(down)
      if not limitInclusion or id in bitIds:
        tmp.append(id)
        order = catalog.GetEntryOrder(down)
        s = levels.get(order,sets.Set())
        s.add(id)
        levels[order] = s
    adjs[bitId] = tmp
  if orderLevels:
    # we'll play a little game and sort the indices in each level by
    #  the number of downlinks they have:
    for order in levels.keys():
      ids = levels[order]
      counts = [len(adjs[id]) for id in ids]
      countOrder = argsort(counts)
      l = [ids[x] for x in countOrder]
      l.reverse()
      levels[order] = l
  return adjs,levels

def GetMolsMatchingBit(mols,bit,fps):
  res = []
  if isinstance(bit,BitGainsInfo):
    bitId = bit.id
  else:
    bitId = bit
  for i,mol in enumerate(mols):
    fp = fps[i]
    if fp[bitId]:
      res.append(mol)
  return res

xl = None
def ShowMolsMatchingBit(mols,bit,fps,actName="",wrapper=None,col=1,row=1):
  global xl
  from Excel import Molecules
  if wrapper is not None:
    xl = wrapper
  if xl is None:
    xl = Molecules.ExcelWrapper()
    xl.Workbooks.Add()
    xl.ChemdrawNewWorksheet()
  mols = GetMolsMatchingBit(mols,bit,fps)
  i = row
  for mol in mols:
    xl[i,col] = Chem.MolToSmiles(mol)
    if actName and mol.HasProp(actName):
      xl[i,col+1] = mol.GetProp(actName)
    i+=1
  i = row
  for mol in mols:
    try:
      xl.ChemdrawConvertCellsToMols(i,col)
      xl.ChemdrawShowPictures(i,col)
    except:
      pass
    i+=1
  return mols
        
    

