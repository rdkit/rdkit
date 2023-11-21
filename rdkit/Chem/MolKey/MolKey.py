#
#  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
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
# Created by Greg Landrum based on code from Thomas Mueller
# August 2012
#
try:
  from rdkit.Avalon import pyAvalonTools
except ImportError:
  raise ImportError("This code requires the RDKit to be built with AvalonTools support")

import base64
import hashlib
import logging
import os
import re
import tempfile
import uuid
from collections import namedtuple

from rdkit import Chem, RDConfig
from rdkit.Chem.MolKey import InchiInfo


class MolIdentifierException(Exception):
  pass


class BadMoleculeException(Exception):
  pass


MOL_KEY_VERSION = '1'

ERROR_DICT = dict(BAD_MOLECULE=1, ALIAS_CONVERSION_FAILED=2, TRANSFORMED=4, FRAGMENTS_FOUND=8,
                  EITHER_WARNING=16, STEREO_ERROR=32, DUBIOUS_STEREO_REMOVED=64, ATOM_CLASH=128,
                  ATOM_CHECK_FAILED=256, SIZE_CHECK_FAILED=512, RECHARGED=1024,
                  STEREO_FORCED_BAD=2048, STEREO_TRANSFORMED=4096, TEMPLATE_TRANSFORMED=8192,
                  INCHI_COMPUTATION_ERROR=65536, RDKIT_CONVERSION_ERROR=131072,
                  INCHI_READWRITE_ERROR=262144, NULL_MOL=524288)

INCHI_COMPUTATION_ERROR = ERROR_DICT['INCHI_COMPUTATION_ERROR']
RDKIT_CONVERSION_ERROR = ERROR_DICT['RDKIT_CONVERSION_ERROR']
INCHI_READWRITE_ERROR = ERROR_DICT['INCHI_READWRITE_ERROR']
NULL_MOL = ERROR_DICT['NULL_MOL']

BAD_SET = (pyAvalonTools.StruChkResult.bad_set | INCHI_COMPUTATION_ERROR | RDKIT_CONVERSION_ERROR
           | INCHI_READWRITE_ERROR | NULL_MOL)

GET_STEREO_RE = re.compile(r'^InChI=1S(.*?)/(t.*?)/m\d/s1(.*$)')
NULL_SMILES_RE = re.compile(r'^\s*$|^\s*NO_STRUCTURE\s*$', re.IGNORECASE)
PATTERN_NULL_MOL = r'^([\s0]+[1-9]+[\s]+V[\w]*)'

CHIRAL_POS = 12

T_NULL_MOL = (NULL_MOL, '')  # the NULL mol result tuple

stereo_code_dict = {}
stereo_code_dict['DEFAULT'] = 0
stereo_code_dict['S_ACHIR'] = 1
stereo_code_dict['S_ABS'] = 2
stereo_code_dict['S_REL'] = 3
stereo_code_dict['S_PART'] = 4
stereo_code_dict['S_UNKN'] = 5
stereo_code_dict['S_ABS_ACHIR'] = 6
stereo_code_dict['R_ONE'] = 11
stereo_code_dict['R_REL'] = 12
stereo_code_dict['R_OTHER'] = 13
stereo_code_dict['MX_ENANT'] = 21
stereo_code_dict['MX_DIAST'] = 22
stereo_code_dict['MX_SP2'] = 31
stereo_code_dict['MX_DIAST_ABS'] = 32
stereo_code_dict['MX_DIAST_REL'] = 33
stereo_code_dict['OTHER'] = 100
stereo_code_dict['UNDEFINED'] = 200


def _fix_all(pat, sbt, my_string):
  try:
    return re.sub(pat, sbt, my_string)
  except Exception:
    return None


def _fix_line_ends(my_string):
  pat = '\r\n{0,1}'
  sbt = '\n'
  return _fix_all(pat, sbt, my_string)


def _fix_chemdraw_header(my_string):
  pat = '0V2000'
  sbt = 'V2000'
  return _fix_all(pat, sbt, my_string)


def _ctab_has_atoms(ctab_lines):
  ''' look at atom count position (line 4, characters 0:3)
    Return True if the count is > 0, False if 0.
    Throw BadMoleculeException if there are no characters
    at the required position or if they cannot be converted
    to a positive integer
    '''
  try:
    a_count = int(ctab_lines[3][0:3])
    if a_count < 0:
      raise BadMoleculeException('Atom count negative')
    return a_count > 0
  except IndexError:
    raise BadMoleculeException('Invalid molfile format')
  except ValueError:
    raise BadMoleculeException(f'Expected integer')


def _ctab_remove_chiral_flag(ctab_lines):
  ''' read the chiral flag (line 4, characters 12:15)
    and set it to 0. Return True if it was 1, False if 0.
    Throw BadMoleculeException if there are no characters
    at the required position or if they where not 0 or 1
    '''
  str_a_count = ctab_lines[3][12:15]
  try:
    a_count = int(str_a_count)
    if a_count == 0:
      rval = False
    elif a_count == 1:
      rval = True
      orig_line = ctab_lines[3]
      ctab_lines[3] = orig_line[:CHIRAL_POS] + '  0' + orig_line[CHIRAL_POS + 3:]
    else:
      raise BadMoleculeException('Expected chiral flag 0 or 1')
  except IndexError:
    raise BadMoleculeException('Invalid molfile format')
  except ValueError:
    raise BadMoleculeException(f'Expected integer, got {str_a_count}')
  return rval


__initCalled = False


def initStruchk(configDir=None, logFile=None):
  global __initCalled
  if configDir is None:
    configDir = os.path.join(RDConfig.RDDataDir, 'struchk')
  if configDir[-1] != os.path.sep:
    configDir += os.path.sep
  if logFile is None:
    fd = tempfile.NamedTemporaryFile(suffix='.log', delete=False)
    fd.close()
    logFile = fd.name
  struchk_init = '''-tm
-ta {0}checkfgs.trn
-tm
-or
-ca {0}checkfgs.chk
-cc
-cl 3
-cs
-cn 999
-l {1}\n'''.format(configDir, logFile)
  initRes = pyAvalonTools.InitializeCheckMol(struchk_init)
  if initRes:
    raise ValueError(f'bad result from InitializeCheckMol: {initRes}')
  __initCalled = True


def CheckCTAB(ctab, isSmiles=True):
  if not __initCalled:
    initStruchk()
  mol_str = ctab
  if not mol_str:
    raise BadMoleculeException('Unexpected blank or NULL molecule')

  mol_str = _fix_line_ends(mol_str)
  mol_str = _fix_chemdraw_header(mol_str)

  if isSmiles:  # branch for NULL_MOL checks
    if mol_str and NULL_SMILES_RE.match(mol_str):
      return T_NULL_MOL
    return pyAvalonTools.CheckMoleculeString(mol_str, isSmiles)

  # decompose the ctab into lines
  # the line terminator may be \n or \r\n, or even r'\n'
  ctab_lines = mol_str.split('\n')
  if len(ctab_lines) <= 3:
    raise BadMoleculeException('Not enough lines in CTAB')
  _ctab_remove_chiral_flag(ctab_lines)

  if not _ctab_has_atoms(ctab_lines):
    return T_NULL_MOL
  # reassemble the ctab lines into one string.
  mol_str = '\n'.join(ctab_lines)
  return pyAvalonTools.CheckMoleculeString(mol_str, isSmiles)


InchiResult = namedtuple('InchiResult', ['error', 'inchi', 'fixed_ctab'])


def GetInchiForCTAB(ctab):
  """
    >>> from rdkit.Chem.MolKey import MolKey
    >>> from rdkit.Avalon import pyAvalonTools
    >>> res = MolKey.GetInchiForCTAB(pyAvalonTools.Generate2DCoords('c1cn[nH]c1C(Cl)Br',True))
    >>> res.inchi
    'InChI=1/C4H4BrClN2/c5-4(6)3-1-2-7-8-3/h1-2,4H,(H,7,8)/t4?/f/h8H'
    >>> res = MolKey.GetInchiForCTAB(pyAvalonTools.Generate2DCoords('c1c[nH]nc1C(Cl)Br',True))
    >>> res.inchi
    'InChI=1/C4H4BrClN2/c5-4(6)3-1-2-7-8-3/h1-2,4H,(H,7,8)/t4?/f/h7H'
    >>>
    """
  inchi = None
  strucheck_err, fixed_mol = CheckCTAB(ctab, False)
  if strucheck_err & BAD_SET:
    return InchiResult(strucheck_err, None, fixed_mol)

  conversion_err = 0
  try:
    r_mol = Chem.MolFromMolBlock(fixed_mol, sanitize=False)
    if r_mol:
      inchi = Chem.MolToInchi(r_mol, '/FixedH /SUU')
      if not inchi:
        conversion_err = INCHI_COMPUTATION_ERROR
    else:
      conversion_err = RDKIT_CONVERSION_ERROR
  except Chem.InchiReadWriteError:
    conversion_err = INCHI_READWRITE_ERROR
  # keep warnings from strucheck
  return InchiResult(strucheck_err | conversion_err, inchi, fixed_mol)


def _make_racemate_inchi(inchi):
  """ Normalize the stereo information (t-layer) to one selected isomer. """
  # set stereo type = 3 (racemate) for consistency
  # reset inverted flag to m0 - not inverted
  new_stereo = '/m0/s3/'
  stereo_match = GET_STEREO_RE.match(inchi)
  if stereo_match:
    inchi = stereo_match.group(1) + new_stereo + stereo_match.group(2)
  return inchi


def _get_identification_string(err, ctab, inchi, stereo_category=None, extra_stereo=None):
  if err & NULL_MOL:
    return _get_null_mol_identification_string(extra_stereo)
  elif err & BAD_SET:  # bad molecules get special key
    return _get_bad_mol_identification_string(ctab, stereo_category, extra_stereo)

  # make key string
  pieces = []
  if inchi:
    pieces.append(inchi)

  if not stereo_category:
    raise MolIdentifierException('Stereo category may not be left undefined')

  pieces.append(f'ST={stereo_category}')
  if extra_stereo:
    pieces.append(f'XTR={extra_stereo}')
  return '/'.join(pieces)


def _get_null_mol_identification_string(extra_stereo):
  return str(uuid.uuid1())


def _get_bad_mol_identification_string(ctab, stereo_category, extra_stereo):
  pieces = []
  ctab_str = ctab
  if ctab_str:  # make the ctab part of the key if available
    ctab_str = _fix_line_ends(ctab_str)
    ctab_str = _fix_chemdraw_header(ctab_str)
    ctab_str = '\n'.join(ctab_str.split('\n')[3:])
    pieces.append(ctab_str.replace('\n', r'\n'))  # make a handy one-line string
  else:
    pass
  if stereo_category:  # add xtra info if available
    pieces.append(f'ST={stereo_category}')
  if extra_stereo:  # add xtra info if available
    pieces.append(f'XTR={extra_stereo}')
  return '/'.join(pieces)


def _identify(err, ctab, inchi, stereo_category, extra_structure_desc=None):
  """ Compute the molecule key based on the inchi string,
    stereo category as well as extra structure
    information """
  key_string = _get_identification_string(err, ctab, inchi, stereo_category, extra_structure_desc)
  if not key_string:
    return None
  hash_key = base64.b64encode(hashlib.md5(key_string.encode('UTF-8')).digest()).decode()
  return f"{MOL_KEY_VERSION}|{hash_key}"


def _get_chiral_identification_string(n_def, n_udf):
  assert n_def >= 0
  assert n_udf >= 0

  if n_def == 0:  # no defined stereocenters
    if n_udf == 0:  # no undefined ones either
      return 'S_ACHIR'  # -> achiral
    elif n_udf == 1:  # one undefined, no defined
      return 'R_ONE'  # -> racemate by convention
    else:  # several undefined, no defined
      return 'S_UNKN'  # -> can't say anything based on the drawing
  else:  # some stereo defined
    if n_udf == 0:  # fully specified stereo
      return 'S_ABS'  # -> absolute stereo
    else:  # multiple possibilities
      return 'S_PART'  # -> assume single compound (can usually be separated)
  return 'OTHER'


def ErrorBitsToText(err):
  " returns a list of error bit descriptions for the error code provided "
  return [k for k, v in ERROR_DICT.items() if (err & v) > 0]


MolKeyResult = namedtuple(
  'MolKeyResult', ['mol_key', 'error', 'inchi', 'fixed_ctab', 'stereo_code', 'stereo_comment'])


def GetKeyForCTAB(ctab, stereo_info=None, stereo_comment=None, logger=None):
  """
    >>> from rdkit.Chem.MolKey import MolKey
    >>> from rdkit.Avalon import pyAvalonTools
    >>> res = MolKey.GetKeyForCTAB(pyAvalonTools.Generate2DCoords('c1ccccc1C(F)Cl', True))
    >>> res.mol_key
    '1|L7676nfGsSIU33wkx//NCg=='
    >>> res.stereo_code
    'R_ONE'
    >>> res = MolKey.GetKeyForCTAB(pyAvalonTools.Generate2DCoords('c1ccccc1[C@H](F)Cl', True))
    >>> res.mol_key
    '1|Aj38EIxf13RuPDQG2A0UMw=='
    >>> res.stereo_code
    'S_ABS'
    >>> res = MolKey.GetKeyForCTAB(pyAvalonTools.Generate2DCoords('c1ccccc1[C@@H](F)Cl', True))
    >>> res.mol_key
    '1|9ypfMrhxn1w0ncRooN5HXw=='
    >>> res.stereo_code
    'S_ABS'
    >>> res = MolKey.GetKeyForCTAB(pyAvalonTools.Generate2DCoords('c1cccc(C(Br)Cl)c1[C@@H](F)Cl', True))
    >>> res.mol_key
    '1|c96jMSlbn7O9GW5d5uB9Mw=='
    >>> res.stereo_code
    'S_PART'
    >>> res = MolKey.GetKeyForCTAB(pyAvalonTools.Generate2DCoords('c1cccc([C@H](Br)Cl)c1[C@@H](F)Cl', True))
    >>> res.mol_key
    '1|+B+GCEardrJteE8xzYdGLA=='
    >>> res.stereo_code
    'S_ABS'
    >>> res = MolKey.GetKeyForCTAB(pyAvalonTools.Generate2DCoords('c1cccc(C(Br)Cl)c1C(F)Cl', True))
    >>> res.mol_key
    '1|5H9R3LvclagMXHp3Clrc/g=='
    >>> res.stereo_code
    'S_UNKN'
    >>> res = MolKey.GetKeyForCTAB(pyAvalonTools.Generate2DCoords('c1cccc(C(Br)Cl)c1C(F)Cl',True), stereo_info='S_REL')
    >>> res.mol_key
    '1|cqKWVsUEY6QNpGCbDaDTYA=='
    >>> res.stereo_code
    'S_REL'
    >>> res.inchi
    'InChI=1/C8H6BrCl2F/c9-7(10)5-3-1-2-4-6(5)8(11)12/h1-4,7-8H/t7?,8?'

    """
  if logger is None:
    logger = logging
  try:
    err, inchi, fixed_mol = GetInchiForCTAB(ctab)
  except BadMoleculeException:
    logger.warning(u'Corrupt molecule substituting no-struct: --->\n{0}\n<----'.format(ctab))
    err = NULL_MOL
    key = _identify(err, '', '', None, None)
    return MolKeyResult(key, err, '', '', None, None)

  # read or estimate stereo category and/or extra structure description
  stereo_category = None
  extra_structure_desc = stereo_comment
  if stereo_info:  # check stereo_info field for coded stereo category and extra stereo info
    info_flds = stereo_info.split(' ', 1)
    code_fld = info_flds[0]
    if code_fld in stereo_code_dict:
      stereo_category = code_fld
      if not stereo_comment and len(info_flds) > 1:
        extra_structure_desc = info_flds[1].strip()
    else:
      logger.warning(f'stereo code {code_fld} not recognized. Using default value for ctab.')

  if not (err & BAD_SET):
    n_stereo, n_undef_stereo, is_meso, dummy = InchiInfo.InchiInfo(
      inchi).get_sp3_stereo()['main']['non-isotopic']
    if stereo_category is None or stereo_category == 'DEFAULT':  # compute if not set
      stereo_category = _get_chiral_identification_string(n_stereo - n_undef_stereo, n_undef_stereo)
  else:
    raise NotImplementedError(
      "currently cannot generate correct keys for molecules with struchk errors")

  key = _identify(err, fixed_mol, inchi, stereo_category, extra_structure_desc)
  return MolKeyResult(key, err, inchi, fixed_mol, stereo_category, extra_structure_desc)


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
