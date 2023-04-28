#  $Id$
#
#  Copyright (C) 2000-2010 Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Configuration for the RDKit Python code

"""

import os
import sys

if 'RDBASE' in os.environ:
  RDBaseDir = os.environ['RDBASE']
  RDCodeDir = os.path.join(RDBaseDir, 'rdkit')
  RDDataDir = os.path.join(RDBaseDir, 'Data')
  RDDocsDir = os.path.join(RDBaseDir, 'Docs')
  RDDemoDir = os.path.join(RDBaseDir, 'Demo')
  RDBinDir = os.path.join(RDBaseDir, 'bin')
  RDProjDir = os.path.join(RDBaseDir, 'Projects')
  RDContribDir = os.path.join(RDBaseDir, 'Contrib')
else:
  import rdkit.RDPaths
  from rdkit.RDPaths import *
  os.environ['RDBASE'] = rdkit.RDPaths._share

rpcTestPort = 8423
pythonTestCommand = "python"

defaultDBUser = 'sysdba'
defaultDBPassword = 'masterkey'


class ObsoleteCodeError(Exception):
  pass


class UnimplementedCodeError(Exception):
  pass


# ---------------------
# the following block contains stuff used by the
# testing infrastructure
pythonExe = sys.executable

# ---------------------
# the following block contains stuff controlling database access:
usePgSQL = False
useSqlLite = False
if not os.environ.get('RD_USESQLLITE', ''):
  try:
    from pyPgSQL import PgSQL
    usePgSQL = True
  except ImportError:
    usePgSQL = False
if not usePgSQL:
  try:
    # python2.5 has this:
    import sqlite3
    useSqlLite = True
  except ImportError:
    try:
      # earlier versions of python:
      from pysqlite2 import dbapi2
      useSqlLite = True
    except ImportError:
      pass

if usePgSQL:
  RDTestDatabase = '::RDTests'
  RDDataDatabase = '::RDData'
elif useSqlLite:
  RDTestDatabase = os.path.join(RDDataDir, "RDTests.sqlt")
  RDDataDatabase = os.path.join(RDDataDir, "RDData.sqlt")
else:
  RDTestDatabase = None
  RDDataDatabase = None

# ---------------------
# the following block contains stuff controlling the program used for
#  3D molecular visualization:
molViewer = os.environ.get('RD_MOLVIEWER', 'PYMOL').upper()
