#  $Id$
#
#  Copyright (C) 2003 Rational Discovery LLC
#     All Rights Reserved
#
from rdkit import Chem
from rdkit.VLib.Transform import TransformNode


class SmartsRemover(TransformNode):
  """ transforms molecules by removing atoms matching smarts patterns

  Assumptions:

    - inputs are molecules


  Sample Usage:
    >>> smis = ['C1CCC1.C=O','C1CCC1C=O','CCC=O.C=O','NCC=O.C=O.CN']
    >>> mols = [Chem.MolFromSmiles(x) for x in smis]
    >>> from rdkit.VLib.Supply import SupplyNode
    >>> suppl = SupplyNode(contents=mols)
    >>> ms = [x for x in suppl]
    >>> len(ms)
    4

    We can pass in SMARTS strings:
    >>> smas = ['C=O','CN']
    >>> tform = SmartsRemover(patterns=smas)
    >>> tform.AddParent(suppl)
    >>> ms = [x for x in tform]
    >>> len(ms)
    4
    >>> Chem.MolToSmiles(ms[0])
    'C1CCC1'
    >>> Chem.MolToSmiles(ms[1])
    'O=CC1CCC1'
    >>> Chem.MolToSmiles(ms[2])
    'CCC=O'
    >>> Chem.MolToSmiles(ms[3])
    'NCC=O'

    We can also remove pieces of the molecule that are not complete
    fragments:
    >>> tform.Destroy()
    >>> smas = ['C=O','CN']
    >>> smas = [Chem.MolFromSmarts(x) for x in smas]
    >>> tform = SmartsRemover(patterns=smas,wholeFragments=0)
    >>> tform.AddParent(suppl)
    >>> ms = [x for x in tform]
    >>> len(ms)
    4
    >>> Chem.MolToSmiles(ms[0])
    'C1CCC1'
    >>> Chem.MolToSmiles(ms[1])
    'C1CCC1'
    >>> Chem.MolToSmiles(ms[3])
    ''

    Or patterns themselves:
    >>> tform.Destroy()
    >>> smas = ['C=O','CN']
    >>> smas = [Chem.MolFromSmarts(x) for x in smas]
    >>> tform = SmartsRemover(patterns=smas)
    >>> tform.AddParent(suppl)
    >>> ms = [x for x in tform]
    >>> len(ms)
    4
    >>> Chem.MolToSmiles(ms[0])
    'C1CCC1'
    >>> Chem.MolToSmiles(ms[3])
    'NCC=O'


  """

  def __init__(self, patterns=[], wholeFragments=1, **kwargs):
    TransformNode.__init__(self, func=self.transform, **kwargs)
    self._wholeFragments = wholeFragments
    self._initPatterns(patterns)

  def _initPatterns(self, patterns):
    nPatts = len(patterns)
    targets = [None] * nPatts
    for i in range(nPatts):
      p = patterns[i]
      if type(p) in (str, bytes):
        m = Chem.MolFromSmarts(p)
        if not m:
          raise ValueError('bad smarts: %s' % (p))
        p = m
      targets[i] = p
    self._patterns = tuple(targets)

  def transform(self, cmpd):
    # sys.stderr.write('\tTRANSFORM: %s\n'%(Chem.MolToSmiles(cmpd)))
    for patt in self._patterns:
      cmpd = Chem.DeleteSubstructs(cmpd, patt, onlyFrags=self._wholeFragments)
      # sys.stderr.write('\t\tAfter %s: %s\n'%(Chem.MolToSmiles(patt),Chem.MolToSmiles(cmpd)))

    return cmpd


biggerTest = """
>>> smis = ['CCOC','CCO.Cl','CC(=O)[O-].[Na+]','OCC','C[N+](C)(C)C.[Cl-]']
>>> mols = [Chem.MolFromSmiles(x) for x in smis]
>>> from rdkit.VLib.Supply import SupplyNode
>>> suppl = SupplyNode(contents=mols)
>>> ms = [x for x in suppl]
>>> len(ms)
5

#>>> salts = ['[Cl;H1&X1,-]','[Na+]','[O;H2,H1&-,X0&-2]']

>>> salts = ['[Cl;H1&X1,-]','[Na+]','[O;H2,H1&-,X0&-2]']
>>> m = mols[2]
>>> m.GetNumAtoms()
5
>>> patts = [Chem.MolFromSmarts(x) for x in salts]
>>> m2 = Chem.DeleteSubstructs(m,patts[0],1)
>>> m2.GetNumAtoms()
5
>>> m2 = Chem.DeleteSubstructs(m2,patts[1],1)
>>> m2.GetNumAtoms()
4
>>> m2 = Chem.DeleteSubstructs(m2,patts[2],1)
>>> m2.GetNumAtoms()
4

>>> tform = SmartsRemover(patterns=salts)
>>> tform.AddParent(suppl)
>>> ms = [x for x in tform]
>>> len(ms)
5

"""

# ------------------------------------
#
#  doctest boilerplate
#
__test__ = {'bigger': biggerTest}


# ------------------------------------
#
#  doctest boilerplate
#
def _runDoctests(verbose=None):  # pragma: nocover
  import doctest
  import sys
  failed, _ = doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=verbose)
  sys.exit(failed)


if __name__ == '__main__':  # pragma: nocover
  _runDoctests()
