# $Id$
#
#  Copyright (C) 2005-2006  Greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

import sys
import traceback

from rdkit.rdBase import DisableLog, EnableLog, LogMessage

_levels = ['rdApp.debug', 'rdApp.info', 'rdApp.warning', 'rdApp.error']
DEBUG = 0
INFO = 1
WARNING = 2
ERROR = 3
CRITICAL = 4


class logger(object):

  def logIt(self, dest, msg, *args, **kwargs):
    if (args):
      msg = msg % args
    LogMessage(dest, msg + '\n')
    if kwargs.get('exc_info', False):
      exc_type, exc_val, exc_tb = sys.exc_info()
      if exc_type:
        LogMessage(dest, '\n')
        txt = ''.join(traceback.format_exception(exc_type, exc_val, exc_tb))
        LogMessage(dest, txt)

  def debug(self, msg, *args, **kwargs):
    self.logIt('rdApp.debug', 'DEBUG: ' + msg, *args, **kwargs)

  def error(self, msg, *args, **kwargs):
    self.logIt('rdApp.error', 'ERROR: ' + msg, *args, **kwargs)

  def info(self, msg, *args, **kwargs):
    self.logIt('rdApp.info', 'INFO: ' + msg, *args, **kwargs)

  def warning(self, msg, *args, **kwargs):
    self.logIt('rdApp.warning', 'WARNING: ' + msg, *args, **kwargs)

  def critical(self, msg, *args, **kwargs):
    self.logIt('rdApp.error', 'CRITICAL: ' + msg, *args, **kwargs)

  def setLevel(self, val):
    global _levels
    for i in range(val, len(_levels)):
      EnableLog(_levels[i])
    for i in range(0, val):
      DisableLog(_levels[i])
