#
# Copyright (C) 2017 Peter Gedeck
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
'''
Collection of utilities to be used with descriptors

'''


def setDescriptorVersion(version='1.0.0'):
  """ Set the version on the descriptor function.

  Use as a decorator """
  def wrapper(func):
    func.version = version
    return func
  return wrapper

class VectorDescriptorWrapper:
    """Wrap a function that returns a vector and make it seem like there
    is one function for each entry.  These functions are added to the global
    namespace with the names provided"""
    def __init__(self, func, names, version):
        self.func = func
        self.names = names
        self.func_key = "__%s"%(func.__name__)
        for i,n in enumerate(names):
            def f(mol, key=self._get_key(i)):
                return self.call_desc(mol, key)
            f.__name__ = n
            f.__qualname__ = n
            f.version = version
            globals()[n] = f

    def _get_key(self, index):
        return "%s%s"%(self.func_key, index)
    
    def call_desc(self, mol, key):
        if mol.HasProp(key):
            return mol.GetDoubleProp(key)
        
        results = self.func(mol)
        _get_key = self._get_key
        for i,v in enumerate(results):
            mol.SetDoubleProp(_get_key(i), v)
        return mol.GetDoubleProp(key)
