#
# Copyright (C) 2007-2010 Greg Landrum
# All Rights Reserved
#
from rdkit import Chem


class PropertyMol(Chem.Mol):
  """ allows rdkit molecules to be pickled with their properties saved.

     >>> import os
     >>> import pickle
     >>> from rdkit import RDConfig
     >>> m = Chem.MolFromMolFile(os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data/benzene.mol'))
     >>> m.GetProp('_Name')
     'benzene.mol'

     by default pickling removes properties:

     >>> m2 = pickle.loads(pickle.dumps(m))
     >>> m2.HasProp('_Name')
     0

     Property mols solve this:

     >>> pm = PropertyMol(m)
     >>> pm.GetProp('_Name')
     'benzene.mol'
     >>> pm.SetProp('MyProp','foo')
     >>> pm.HasProp('MyProp')
     1

     >>> pm2 = pickle.loads(pickle.dumps(pm))
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

     That wouldn't work with a standard mol

     but the Property mols still convert all values to strings before storing:

     >>> pm.GetProp('IntVal')
     '1'

     This is a test for sf.net issue 2880943: make sure properties end up in SD files:

     >>> import tempfile, os
     >>> fn = tempfile.NamedTemporaryFile(suffix='.sdf', delete=False).name
     >>> w = Chem.SDWriter(fn)
     >>> w.write(pm)
     >>> w=None
     >>> with open(fn,'r') as inf:
     ...   txt = inf.read()
     >>> '<IntVal>' in txt
     True
     >>> try:
     ...   os.unlink(fn)
     ... except Exception:
     ...   pass

     The next level of that bug: does writing a *depickled* propertymol
     to an SD file include properties:

     >>> fn = tempfile.NamedTemporaryFile(suffix='.sdf', delete=False).name
     >>> w = Chem.SDWriter(fn)
     >>> pm = pickle.loads(pickle.dumps(pm))
     >>> w.write(pm)
     >>> w=None
     >>> with open(fn,'r') as inf:
     ...   txt = inf.read()
     >>> '<IntVal>' in txt
     True
     >>> try:
     ...   os.unlink(fn)
     ... except Exception:
     ...   pass



    """
  __getstate_manages_dict__ = True

  def __init__(self, mol):
    if not isinstance(mol, Chem.Mol):
      return
    Chem.Mol.__init__(self, mol)
    for pn in mol.GetPropNames(includePrivate=True):
      self.SetProp(pn, mol.GetProp(pn))

  def SetProp(self, nm, val):
    Chem.Mol.SetProp(self, nm, str(val))

  def __getstate__(self):
    pDict = {}
    for pn in self.GetPropNames(includePrivate=True):
      pDict[pn] = self.GetProp(pn)
    return {'pkl': self.ToBinary(), 'propD': pDict}

  def __setstate__(self, stateD):
    Chem.Mol.__init__(self, stateD['pkl'])
    for prop, val in stateD['propD'].items():
      self.SetProp(prop, val)

    # ------------------------------------
    #
    #  doctest boilerplate
    #


def _test():
  import doctest
  import sys
  return doctest.testmod(sys.modules["__main__"], optionflags=doctest.ELLIPSIS)


if __name__ == '__main__':
  import sys
  failed, tried = _test()
  sys.exit(failed)
