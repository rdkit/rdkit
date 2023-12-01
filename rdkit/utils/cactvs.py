#
#  Copyright (C) 2003  Rational Discovery LLC
#   All Rights Reserved
#
import os
import sys
import tempfile

from PIL import Image

from rdkit import RDConfig


def SmilesToGif(smiles, fileNames, size=(200, 200), cmd=None, dblSize=0, frame=0):
  if isinstance(smiles, str):
    smiles = [smiles]
    fileNames = [fileNames]
  if cmd is None:
    cmd = os.path.join('/usr', 'local', 'bin', 'csts')
  baseCmd = cmd
  assert len(smiles) == len(fileNames)
  width = size[0]
  height = size[1]
  args = ""
  nDone = 0
  while nDone < len(smiles):
    smi = smiles[nDone]
    name = fileNames[nDone]
    if not dblSize:
      args += "ens get [ens create {%(smi)s}] E_GIF {} {width %(width)d height %(height)d bgcolor white filename %(name)s format gif frame %(frame)d};" % locals(
      )
    else:
      args += "ens get [ens create {%(smi)s}] E_GIF {} {width %(width)d height %(height)d bgcolor white filename %(name)s format gif symbolfontsize 24 frame %(frame)d linewidth 2.8 linespacing 4.0};" % locals(
      )
    nDone += 1
  if args:
    with tempfile.NamedTemporaryFile('w+', suffix='.cmd', delete=False) as tmp:
      tmp.write(args + '\n')
    try:
      cmd = "%s < %s" % (baseCmd, tmp.name)
      os.system(cmd)
    except Exception:
      import traceback
      traceback.print_exc()
      sys.stderr.write('CMD: %s\n' % cmd)
      res = 0
    else:
      res = 1
      for name in fileNames:
        if not os.path.exists(name):
          res = 0
          break
    try:
      os.unlink(tmp.name)
    except Exception:
      pass
  return res


def SmilesToImage(smiles, **kwargs):
  with tempfile.NamedTemporaryFile(suffix='.gif') as tmp:
    ok = SmilesToGif(smiles, tmp.name, **kwargs)
    if ok:
      img = Image.open(tmp.name)
    else:
      img = None
  return img
