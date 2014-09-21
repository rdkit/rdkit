# $Id$
#
#  Copyright (C) 2007 Greg Landrum
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import sys
from rdkit import Chem

class FastSDMolSupplier(object):
  """ A wrapper around an SDMolSupplier that precomputes and stores 
      molecular indices (via text processing) to allow quick length 
      calculations and random access.

      NOTE that this class needs to have the entire SD data in memory,
      so it's probably not particularly useful with large files.
  """
  suppl=None
  data=None
  sanitize=True
  def __init__(self,fileN=None,data=None,sanitize=True,removeHs=True):
    if fileN:
      data = open(fileN,'r').read()
    self.sanitize=sanitize
    self.removeHs=removeHs
    if data:
      data = data.replace('\r\n','\n')
      self.init(data)
      
  def init(self,data,recogTxt='$$$$\n'):
    if not data:
      raise ValueError('no data')
    # FIX: it'd be nice to not be caching data locally like this, but it's the easiest
    # way to handle pickle support.
    self.data=data
    self.suppl = Chem.SDMolSupplier()
    self.suppl.SetData(data,sanitize=self.sanitize,removeHs=self.removeHs)

    self._pos = [0]
    p = 0
    while 1:
      try:
        p = data.index(recogTxt,p+1)
        p+=len(recogTxt)
      except:
        break
      else:
        self._pos.append(p)
    self._pos.pop(-1)
    self.suppl._SetStreamIndices(self._pos)
    self._idx=0
    
  def GetItemText(self,idx):
    startOfItem = self._pos[idx]
    if idx+1<len(self._pos):
      endOfItem = self._pos[idx+1]
    else:
      endOfItem = -1
    return self.data[startOfItem:endOfItem]
  
  def reset(self):
    self.suppl.reset()
    self._idx=0

  # ----------------------------------------------------------------
  # support random access and an iterator interface:
  def __iter__(self):
    self.suppl.reset()
    return self
  def next(self):
    self._idx+=1
    return self.suppl.next()
  
  def __len__(self):
    return len(self.suppl)
  def __getitem__(self,idx):
    return self.suppl[idx]
  

