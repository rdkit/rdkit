#
#  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Created by Greg Landrum, Dec 2006
#

import os
import re

from rdkit import Chem
from rdkit import RDConfig
from contextlib import closing


class SaltRemover(object):
  defnFilename = os.path.join(RDConfig.RDDataDir, 'Salts.txt')

  def __init__(self, defnFilename=None, defnData=None):
    if defnFilename:
      self.defnFilename = defnFilename
    self.defnData = defnData
    self.salts = None
    self._initPatterns()

  def _initPatterns(self):
    """

    >>> remover = SaltRemover()
    >>> len(remover.salts)>0
    True

    >>> remover = SaltRemover(defnData="[Cl,Br]")
    >>> len(remover.salts)
    1

    >>> from rdkit import RDLogger
    >>> RDLogger.DisableLog('rdApp.error')
    >>> remover = SaltRemover(defnData="[Cl,fail]")
    Traceback (most recent call last):
      ...
    ValueError: [Cl,fail]

    >>> RDLogger.EnableLog('rdApp.error')
    """
    whitespace = re.compile(r'[\t ]+')
    if self.defnData:
      from rdkit.six.moves import cStringIO as StringIO
      inF = StringIO(self.defnData)
    else:
      inF = open(self.defnFilename, 'r')
    with closing(inF):
      self.salts = []
      for line in inF:
        line = line.strip().split('//')[0]
        if line:
          splitL = whitespace.split(line)
          salt = Chem.MolFromSmarts(splitL[0])
          if salt is None:
            raise ValueError(line)
          self.salts.append(salt)

  def StripMol(self, mol, dontRemoveEverything=False):
    """

    >>> remover = SaltRemover(defnData="[Cl,Br]")
    >>> len(remover.salts)
    1

    >>> mol = Chem.MolFromSmiles('CN(C)C.Cl')
    >>> res = remover.StripMol(mol)
    >>> res is not None
    True
    >>> res.GetNumAtoms()
    4

    Notice that all salts are removed:
    >>> mol = Chem.MolFromSmiles('CN(C)C.Cl.Cl.Br')
    >>> res = remover.StripMol(mol)
    >>> res.GetNumAtoms()
    4

    Matching (e.g. "salt-like") atoms in the molecule are unchanged:
    >>> mol = Chem.MolFromSmiles('CN(Br)Cl')
    >>> res = remover.StripMol(mol)
    >>> res.GetNumAtoms()
    4

    >>> mol = Chem.MolFromSmiles('CN(Br)Cl.Cl')
    >>> res = remover.StripMol(mol)
    >>> res.GetNumAtoms()
    4

    Charged salts are handled reasonably:
    >>> mol = Chem.MolFromSmiles('C[NH+](C)(C).[Cl-]')
    >>> res = remover.StripMol(mol)
    >>> res.GetNumAtoms()
    4


    Watch out for this case (everything removed):
    >>> remover = SaltRemover()
    >>> len(remover.salts)>1
    True
    >>> mol = Chem.MolFromSmiles('CC(=O)O.[Na]')
    >>> res = remover.StripMol(mol)
    >>> res.GetNumAtoms()
    0

    dontRemoveEverything helps with this by leaving the last salt:
    >>> res = remover.StripMol(mol,dontRemoveEverything=True)
    >>> res.GetNumAtoms()
    4

    but in cases where the last salts are the same, it can't choose
    between them, so it returns all of them:
    >>> mol = Chem.MolFromSmiles('Cl.Cl')
    >>> res = remover.StripMol(mol,dontRemoveEverything=True)
    >>> res.GetNumAtoms()
    2

    """

    def _applyPattern(m, salt, notEverything):
      nAts = m.GetNumAtoms()
      if not nAts:
        return m
      res = m

      t = Chem.DeleteSubstructs(res, salt, True)
      if not t or (notEverything and t.GetNumAtoms() == 0):
        return res
      else:
        res = t
      while res.GetNumAtoms() and nAts > res.GetNumAtoms():
        nAts = res.GetNumAtoms()
        t = Chem.DeleteSubstructs(res, salt, True)
        if notEverything and t.GetNumAtoms() == 0:
          break
        else:
          res = t
      return res

    if dontRemoveEverything and len(Chem.GetMolFrags(mol)) <= 1:
      return mol
    modified = False
    natoms = mol.GetNumAtoms()
    for salt in self.salts:
      mol = _applyPattern(mol, salt, dontRemoveEverything)
      if natoms != mol.GetNumAtoms():
        natoms = mol.GetNumAtoms()
        modified = True
        if dontRemoveEverything and len(Chem.GetMolFrags(mol)) <= 1:
          break
    if modified and mol.GetNumAtoms() > 0:
      Chem.SanitizeMol(mol)
    return mol

  def __call__(self, mol, dontRemoveEverything=False):
    """

    >>> remover = SaltRemover(defnData="[Cl,Br]")
    >>> len(remover.salts)
    1

    >>> mol = Chem.MolFromSmiles('CN(C)C.Cl')
    >>> res = remover(mol)
    >>> res is not None
    True
    >>> res.GetNumAtoms()
    4

    """
    return self.StripMol(mol, dontRemoveEverything=dontRemoveEverything)


# ------------------------------------
#
#  doctest boilerplate
#
def _runDoctests(verbose=None):  # pragma: nocover
  import sys
  import doctest
  failed, _ = doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=verbose)
  sys.exit(failed)


if __name__ == '__main__':  # pragma: nocover
  _runDoctests()
