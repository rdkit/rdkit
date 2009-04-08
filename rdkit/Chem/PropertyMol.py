# $Id$
#
# Copyright (C) 2007-2009 Greg Landrum
# All Rights Reserved
#
from rdkit import Chem

class PropertyMol(Chem.Mol):
  """ allows rdkit molecules to be pickled with their properties saved.

   >>> import cPickle
   >>> m = Chem.MolFromMolFile('tests/benzene.mol')
   >>> m.GetProp('_Name')
   'benzene.mol'

   by default pickling removes properties:
   >>> m2 = cPickle.loads(cPickle.dumps(m))
   >>> m2.HasProp('_Name')
   0

   Property mols solve this:
   >>> pm = PropertyMol(m)
   >>> pm.GetProp('_Name')
   'benzene.mol'
   >>> pm.SetProp('MyProp','foo')
   >>> pm.HasProp('MyProp')
   1

   >>> pm2 = cPickle.loads(cPickle.dumps(pm))
   >>> Chem.MolToSmiles(pm2)
   'c1ccccc1'
   >>> pm2.GetProp('_Name')
   'benzene.mol'
   >>> pm2.HasProp('MyProp')
   1
   >>> pm2.GetProp('MyProp')
   'foo'
   >>> pm2.HasProp('MissingProp')
   0

   Property mols are a bit more permissive about the types
   of property values:
   >>> pm.SetProp('IntVal',1)

   That wouldn't work with a standard mol:
   >>> m.SetProp('IntVal',1)
   Traceback (most recent call last):
     ...
   ArgumentError: Python argument types in
       Mol.SetProp(Mol, str, int)
   did not match C++ signature:
       SetProp(RDKit::ROMol self, char const* key, std::string val, bool computed=False)

   but the Property mols still convert all values to strings before storing:
   >>> pm.GetProp('IntVal')
   '1'
   
  """
  __getstate_manages_dict__=True
  def __init__(self,mol):
    if not isinstance(mol,Chem.Mol): return
    Chem.Mol.__init__(self,mol.ToBinary())
    self.__propDict={}
    for pn in mol.GetPropNames(includePrivate=True):
      self.__propDict[pn]=mol.GetProp(pn)
  def GetPropNames(self):
    return self.__propDict.keys()
  def GetProp(self,prop):
    return self.__propDict[prop]
  def SetProp(self,prop,val,**kwargs):
    self.__propDict[prop]=str(val)
  def HasProp(self,prop):
    return int(self.__propDict.has_key(prop))
  def __getstate__(self):
    return {'pkl':self.ToBinary(),
            'propD':self.__propDict}
  def __setstate__(self,stateD):
    Chem.Mol.__init__(self,stateD['pkl'])
    self.__propDict=stateD['propD']

    
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
  
