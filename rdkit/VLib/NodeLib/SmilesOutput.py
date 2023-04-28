#  $Id$
#
#  Copyright (C) 2003 Rational Discovery LLC
#     All Rights Reserved
#
from rdkit import Chem
from rdkit.VLib.Output import OutputNode as BaseOutputNode


class OutputNode(BaseOutputNode):
  """ dumps smiles output

    Assumptions:

      - destination supports a write() method

      - inputs (parents) can be stepped through in lockstep


    Usage Example:
      >>> smis = ['C1CCC1','C1CC1','C=O','NCC']
      >>> mols = [Chem.MolFromSmiles(x) for x in smis]
      >>> from rdkit.VLib.Supply import SupplyNode
      >>> suppl = SupplyNode(contents=mols)
      >>> from io import StringIO
      >>> sio = StringIO()
      >>> node = OutputNode(dest=sio,delim=', ')
      >>> node.AddParent(suppl)
      >>> ms = [x for x in node]
      >>> len(ms)
      4
      >>> txt = sio.getvalue()
      >>> repr(txt)
      "'1, C1CCC1\\\\n2, C1CC1\\\\n3, C=O\\\\n4, CCN\\\\n'"

    """

  def __init__(self, dest=None, delim='\t', idField=None, **kwargs):
    BaseOutputNode.__init__(self, dest=dest, strFunc=self.smilesOut)
    self._dest = dest
    self._idField = idField
    self._delim = delim
    self._nDumped = 0

  def reset(self):
    BaseOutputNode.reset(self)
    self._nDumped = 0

  def smilesOut(self, mol):
    self._nDumped += 1
    if isinstance(mol, (tuple, list)):
      args = mol
      mol = args[0]
      if len(args) > 1:
        args = list(args[1:])
      else:
        args = []
    else:
      args = []

    if self._idField and mol.HasProp(self._idField):
      label = mol.GetProp(self._idField)
    else:
      label = str(self._nDumped)
    smi = Chem.MolToSmiles(mol)
    outp = [label, smi] + args
    return '%s\n' % (self._delim.join(outp))


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
