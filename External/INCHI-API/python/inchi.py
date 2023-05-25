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


__all__ = [
  'MolToInchiAndAuxInfo', 'MolToInchi', 'MolBlockToInchiAndAuxInfo', 'MolBlockToInchi',
  'MolFromInchi', 'InchiReadWriteError', 'InchiToInchiKey', 'MolToInchiKey', 'INCHI_AVAILABLE'
]
