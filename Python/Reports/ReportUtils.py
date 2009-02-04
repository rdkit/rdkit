#
#  Copyright (C) 2004  Rational Discovery LLC
#   All Rights Reserved
#
import RDConfig
import tempfile,os

class TempFileHandler:
  def __init__(self):
    self.files = []

  def get(self,ext=''):
    fName = tempfile.mktemp(ext)
    self.files.append(fName)
    return fName
  def __del__(self):
    for fileN in self.files:
      try:
        os.unlink(fileN)
      except OSError:
        pass

