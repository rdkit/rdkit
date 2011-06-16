# $Id$

# Copyright (C) 2011-2011 Novartis Institutes for BioMedical Research, Inc
#
#  @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

INCHI_AVAILABLE = True

import rdinchi
import logging
from rdkit import RDLogger
logger = RDLogger.logger()

logLevelToLogFunctionLookup = {
        logging.INFO : logger.info,
        logging.DEBUG : logger.debug,
        logging.WARNING : logger.warning,
        logging.CRITICAL : logger.critical,
        logging.ERROR : logger.error
        }

class InchiReadWriteError(Exception):
    pass

def MolFromInchi(inchi, sanitize=True, removeHs=True, logLevel=None,
        treatWarningAsError=False):
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
    except ValueError,e :
        logger.error(str(e))
        return None

    if logLevel is not None:
        if logLogLevel not in logLevelToLogFunctionLookup:
            raise ValueError("Unsupported log level: %d" % logLogLevel)
        log = logLevelToLogFunctionLookup[logLevel]
        log(log)
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

def MolToInchiAndAuxInfo(mol, options="", logLevel=None,
        treatWarningAsError=False):
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
        log(log)
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
    try:
        inchi, aux = MolToInchiAndAuxInfo(mol, options, logLevel=logLevel,
                treatWarningAsError=treatWarningAsError)
    except InchiReadWriteError,inst:
        inchi, aux, message = inst
        raise InchiReadWriteError(inchi, message)
    return inchi

def InchiToInchiKey(inchi):
    """Return the InChI key for the given InChI string. Return None on error"""
    ret = rdinchi.InchiToInchiKey(inchi)
    if ret: return ret
    else: return None

__all__ = ['MolToInchiAndAuxInfo', 'MolToInchi', 'MolFromInchi', 
        'InchiReadWriteError', 'InchiToInchiKey']
