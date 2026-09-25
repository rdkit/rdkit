#
#  Copyright (c) 2011, Novartis Institutes for BioMedical Research Inc.
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

INCHI_AVAILABLE = True

import logging
import re

from rdkit import Chem, Geometry
from rdkit import RDLogger
from rdkit.Chem import rdinchi

logger = RDLogger.logger()

logLevelToLogFunctionLookup = {
  logging.INFO: logger.info,
  logging.DEBUG: logger.debug,
  logging.WARNING: logger.warning,
  logging.CRITICAL: logger.critical,
  logging.ERROR: logger.error
}


class InchiReadWriteError(Exception):
  pass


def _parse_auxinfo_coordinates(auxinfo):
  """Parse the rC: (coordinate) layer from an InChI AuxInfo string.

  Returns (coords_list, is_3d) where coords_list is a list of (x, y, z) tuples
  in original input atom order, or (None, None) if parsing fails or coords are empty.
  """
  if not auxinfo:
    return None, None
  match = re.search(r'/rC:([^/]+)', auxinfo)
  if not match:
    return None, None
  entries = match.group(1).split(';')
  # Remove trailing empty entries from trailing semicolons
  while entries and not entries[-1].strip():
    entries.pop()
  coords = []
  for entry in entries:
    entry = entry.strip()
    if not entry:
      return None, None
    parts = entry.split(',')
    try:
      x = float(parts[0]) if parts[0] else 0.0
      y = float(parts[1]) if len(parts) > 1 and parts[1] else 0.0
      z = float(parts[2]) if len(parts) > 2 and parts[2] else 0.0
    except (ValueError, IndexError):
      return None, None
    coords.append((x, y, z))
  if not coords:
    return None, None
  # All-zero coordinates means no real coords were present
  if all(x == 0.0 and y == 0.0 and z == 0.0 for x, y, z in coords):
    return None, None
  is_3d = any(z != 0.0 for _, _, z in coords)
  return coords, is_3d


def _parse_auxinfo_atom_order(auxinfo):
  """Parse the N: (atom numbering) layer from an InChI AuxInfo string.

  Returns a list of 0-based original atom indices, or None if parsing fails.
  The returned list maps from InChI canonical order to original atom order:
  result[i] is the original atom index for InChI canonical atom i.
  """
  if not auxinfo:
    return None
  match = re.search(r'/N:([^/]+)', auxinfo)
  if not match:
    return None
  # The N: layer contains comma-separated 1-based atom indices
  # possibly with semicolons separating disconnected fragments
  tokens = match.group(1).replace(';', ',').split(',')
  try:
    return [int(t) - 1 for t in tokens]
  except (ValueError, IndexError):
    return None


def _attach_conformer(mol, coords, is_3d):
  """Attach parsed /rC: coordinates to a molecule as a conformer."""
  if coords is not None and len(coords) == mol.GetNumAtoms():
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(is_3d)
    for i, (x, y, z) in enumerate(coords):
      conf.SetAtomPosition(i, Geometry.Point3D(x, y, z))
    mol.AddConformer(conf, assignId=True)
  return mol


def _build_inverse_permutation(atom_order, size):
  """Build the inverse permutation for RenumberAtoms.

  atom_order[inchi_idx] = original_idx. Returns new_order where
  new_order[original_idx] = inchi_idx, or None if any index is out of range.
  """
  new_order = [0] * size
  for inchi_idx, orig_idx in enumerate(atom_order):
    if orig_idx >= size:
      return None
    new_order[orig_idx] = inchi_idx
  return new_order


def MolFromInchiAndAuxInfo(inchi, auxinfo, sanitize=True, removeHs=True, logLevel=None,
                           treatWarningAsError=False):
  """Construct a molecule from an InChI string and its AuxInfo, restoring the
  original atom ordering.

    Keyword arguments:
    sanitize -- set to True to enable sanitization of the molecule. Default is
    True
    removeHs -- set to True to remove Hydrogens from a molecule. This only
    makes sense when sanitization is enabled
    logLevel -- the log level used for logging logs and messages from InChI
    API. set to None to diable the logging completely
    treatWarningAsError -- set to True to raise an exception in case of a
    molecule that generates warning in calling InChI API. The resultant
    molecule and error message are part of the excpetion

    Returns:
    a rdkit.Chem.rdchem.Mol instance with atoms reordered to match the
    original atom ordering encoded in the AuxInfo
    """
  mol = MolFromInchi(inchi, sanitize=sanitize, removeHs=removeHs, logLevel=logLevel,
                     treatWarningAsError=treatWarningAsError)
  if mol is None:
    return None

  # /rC: coordinates are in original input order; attach after reordering.
  coords, is_3d = _parse_auxinfo_coordinates(auxinfo)

  atom_order = _parse_auxinfo_atom_order(auxinfo)
  if atom_order is None:
    return _attach_conformer(mol, coords, is_3d)

  from rdkit.Chem import RenumberAtoms
  num_mol_atoms = mol.GetNumAtoms()

  if removeHs:
    # N: layer lists heavy atoms only; must match molecule size after H removal
    if len(atom_order) != num_mol_atoms:
      return mol
    new_order = _build_inverse_permutation(atom_order, num_mol_atoms)
    if new_order is None:
      return mol
    return _attach_conformer(RenumberAtoms(mol, new_order), coords, is_3d)

  if len(atom_order) > num_mol_atoms:
    return mol

  if len(atom_order) == num_mol_atoms:
    new_order = _build_inverse_permutation(atom_order, num_mol_atoms)
    if new_order is None:
      return mol
    return _attach_conformer(RenumberAtoms(mol, new_order), coords, is_3d)

  # More atoms in molecule than in atom_order (explicit Hs added by InChI).
  # Place heavy atoms in original order, then append Hs.
  num_heavy = len(atom_order)
  new_order = _build_inverse_permutation(atom_order, num_heavy)
  if new_order is None:
    return mol
  new_order.extend([0] * (num_mol_atoms - num_heavy))
  h_slot = num_heavy
  for i in range(num_mol_atoms):
    if i not in atom_order:
      if h_slot >= num_mol_atoms:
        return mol
      new_order[h_slot] = i
      h_slot += 1
  return _attach_conformer(RenumberAtoms(mol, new_order), coords, is_3d)


def MolFromInchi(inchi, sanitize=True, removeHs=True, logLevel=None, treatWarningAsError=False):
  """Construct a molecule from a InChI string

    Keyword arguments:
    sanitize -- set to True to enable sanitization of the molecule. Default is
    True
    removeHs -- set to True to remove Hydrogens from a molecule. This only
    makes sense when sanitization is enabled
    logLevel -- the log level used for logging logs and messages from InChI
    API. set to None to diable the logging completely
    treatWarningAsError -- set to True to raise an exception in case of a
    molecule that generates warning in calling InChI API. The resultant
    molecule  and error message are part of the excpetion

    Returns:
    a rdkit.Chem.rdchem.Mol instance
    """
  try:
    mol, retcode, message, log = rdinchi.InchiToMol(inchi, sanitize, removeHs)
  except ValueError as e:
    logger.error(str(e))
    return None

  if logLevel is not None:
    if logLevel not in logLevelToLogFunctionLookup:
      raise ValueError("Unsupported log level: %d" % logLevel)
    log = logLevelToLogFunctionLookup[logLevel]
    if retcode == 0:
      log(message)

  if retcode != 0:
    if retcode == 1:
      logger.warning(message)
    else:
      logger.error(message)
  if treatWarningAsError and retcode != 0:
    raise InchiReadWriteError(mol, message)
  return mol


def MolToInchiAndAuxInfo(mol, options="", logLevel=None, treatWarningAsError=False):
  """Returns the standard InChI string and InChI auxInfo for a molecule

    Keyword arguments:
    logLevel -- the log level used for logging logs and messages from InChI
    API. set to None to diable the logging completely
    treatWarningAsError -- set to True to raise an exception in case of a
    molecule that generates warning in calling InChI API. The resultant InChI
    string and AuxInfo string as well as the error message are encoded in the
    exception.

    Returns:
    a tuple of the standard InChI string and the auxInfo string returned by
    InChI API, in that order, for the input molecule
    """
  inchi, retcode, message, logs, aux = rdinchi.MolToInchi(mol, options)
  if logLevel is not None:
    if logLevel not in logLevelToLogFunctionLookup:
      raise ValueError("Unsupported log level: %d" % logLevel)
    log = logLevelToLogFunctionLookup[logLevel]
    if retcode == 0:
      log(message)
  if retcode != 0:
    if retcode == 1:
      logger.warning(message)
    else:
      logger.error(message)

  if treatWarningAsError and retcode != 0:
    raise InchiReadWriteError(inchi, aux, message)
  return inchi, aux


def MolBlockToInchiAndAuxInfo(molblock, options="", logLevel=None, treatWarningAsError=False):
  """Returns the standard InChI string and InChI auxInfo for a mol block

    Keyword arguments:
    logLevel -- the log level used for logging logs and messages from InChI
    API. set to None to diable the logging completely
    treatWarningAsError -- set to True to raise an exception in case of a
    molecule that generates warning in calling InChI API. The resultant InChI
    string and AuxInfo string as well as the error message are encoded in the
    exception.

    Returns:
    a tuple of the standard InChI string and the auxInfo string returned by
    InChI API, in that order, for the input molecule
    """
  inchi, retcode, message, logs, aux = rdinchi.MolBlockToInchi(molblock, options)
  if logLevel is not None:
    if logLevel not in logLevelToLogFunctionLookup:
      raise ValueError("Unsupported log level: %d" % logLevel)
    log = logLevelToLogFunctionLookup[logLevel]
    if retcode == 0:
      log(message)
  if retcode != 0:
    if retcode == 1:
      logger.warning(message)
    else:
      logger.error(message)

  if treatWarningAsError and retcode != 0:
    raise InchiReadWriteError(inchi, aux, message)
  return inchi, aux


def MolToInchi(mol, options="", logLevel=None, treatWarningAsError=False):
  """Returns the standard InChI string for a molecule

    Keyword arguments:
    logLevel -- the log level used for logging logs and messages from InChI
    API. set to None to diable the logging completely
    treatWarningAsError -- set to True to raise an exception in case of a
    molecule that generates warning in calling InChI API. The resultant InChI
    string and AuxInfo string as well as the error message are encoded in the
    exception.

    Returns:
    the standard InChI string returned by InChI API for the input molecule
    """
  if options.find('AuxNone') == -1:
    if options:
      options += " /AuxNone"
    else:
      options += "/AuxNone"

  try:
    inchi, aux = MolToInchiAndAuxInfo(mol, options, logLevel=logLevel,
                                      treatWarningAsError=treatWarningAsError)
  except InchiReadWriteError as inst:
    inchi, aux, message = inst.args
    raise InchiReadWriteError(inchi, message)
  return inchi


def MolBlockToInchi(molblock, options="", logLevel=None, treatWarningAsError=False):
  """Returns the standard InChI string for a mol block

    Keyword arguments:
    logLevel -- the log level used for logging logs and messages from InChI
    API. set to None to diable the logging completely
    treatWarningAsError -- set to True to raise an exception in case of a
    molecule that generates warning in calling InChI API. The resultant InChI
    string and AuxInfo string as well as the error message are encoded in the
    exception.

    Returns:
    the standard InChI string returned by InChI API for the input molecule
    """
  if options.find('AuxNone') == -1:
    if options:
      options += " /AuxNone"
    else:
      options += "/AuxNone"

  try:
    inchi, aux = MolBlockToInchiAndAuxInfo(molblock, options, logLevel=logLevel,
                                           treatWarningAsError=treatWarningAsError)
  except InchiReadWriteError as inst:
    inchi, aux, message = inst.args
    raise InchiReadWriteError(inchi, message)
  return inchi


def InchiToInchiKey(inchi):
  """Return the InChI key for the given InChI string. Return None on error"""
  ret = rdinchi.InchiToInchiKey(inchi)
  if ret:
    return ret
  else:
    return None


def MolToInchiKey(mol, options=""):
  """Returns the standard InChI key for a molecule

    Returns:
    the standard InChI key returned by InChI API for the input molecule
    """
  return rdinchi.MolToInchiKey(mol, options)


GetInchiVersion = rdinchi.GetInchiVersion

__all__ = [
  'MolToInchiAndAuxInfo', 'MolToInchi', 'MolBlockToInchiAndAuxInfo', 'MolBlockToInchi',
  'MolFromInchi', 'MolFromInchiAndAuxInfo', 'InchiReadWriteError', 'InchiToInchiKey',
  'MolToInchiKey', 'GetInchiVersion', 'INCHI_AVAILABLE'
]
