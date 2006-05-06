# $Id$
#
# Copyright (C) 2005 Rational Discovery LLC
#  All Rights Reserved
#
import os,tempfile

class TempFileHandler:
  def __init__(self):
    self.files = []

  def get(self,ext=''):
    fName = tempfile.mktemp(ext)
    self.files.append(fName)
    return fName
  def __del__(self):
    for fileN in self.files:
      try:
        os.unlink(fileN)
      except OSError:
        pass
    

def ModuleVarsToDict(mod):
  """
    >>> import LocalConfig
    >>> LocalConfig.bogusVar1=1
    >>> d = ModuleVarsToDict(LocalConfig)
    >>> d['bogusVar1']
    1
    >>> d['applicationName']==LocalConfig.applicationName
    True
    >>> d['distance3DSlop']==LocalConfig.distance3DSlop
    True

    we don't put callables in the results:
    >>> LocalConfig.tempFunc = lambda x:x
    >>> d = ModuleVarsToDict(LocalConfig)
    >>> d.has_key('tempFunc')
    False

  """
  res = {}
  for name,val in mod.__dict__.iteritems():
    if name[0] != '_' and not callable(val):
      res[name] = val
  return res

def DictToModuleVars(inDict,mod):
  """
    >>> d = {}
    >>> d['bogusVar1']=11
    >>> import LocalConfig
    >>> hasattr(LocalConfig,'bogusVar1')
    False
    >>> DictToModuleVars(d,LocalConfig)
    >>> LocalConfig.bogusVar1
    11

  """
  for name,val in inDict.iteritems():
    setattr(mod,name,val)
  return


#------------------------------------
#
#  doctest boilerplate
#
def _test():
  import doctest,sys
  return doctest.testmod(sys.modules["__main__"])

if __name__ == '__main__':
  import sys
  failed,tried = _test()
  sys.exit(failed)
