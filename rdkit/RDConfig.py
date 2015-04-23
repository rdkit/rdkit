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

import os,sys
if 'RDBASE' in os.environ:
  RDBaseDir=os.environ['RDBASE']
  RDCodeDir=os.path.join(RDBaseDir,'rdkit')
  RDDataDir=os.path.join(RDBaseDir,'Data')
  RDDocsDir=os.path.join(RDBaseDir,'Docs')
  RDDemoDir=os.path.join(RDBaseDir,'Demo')
  RDBinDir=os.path.join(RDBaseDir,'bin')
  RDProjDir=os.path.join(RDBaseDir,'Projects')
elif 'CONDA_DEFAULT_ENV' in os.environ:
  class MissingCondaDirectoryException(Exception):
    pass

  class MissingCondaDir:
    def __init__(self, directory):
      self.directory = directory
      
    def __getattr__(self, a):
      raise MissingCondaDirectoryException(
        "RDKit Anaconda builds do not ship with the RDBASE/%s directory"%self.directory)

  # we are running in a conda environ.
  RDCodeDir=os.path.dirname(__file__)
  splitdir = RDCodeDir.split(os.path.sep)
  condaDir = splitdir[:-4]
  if condaDir[0]=='':
    condaDir[0] = os.path.sep
  condaDir += ['share','RDKit']
  RDBaseDir = os.path.join(*condaDir)
  RDDataDir=os.path.join(RDBaseDir,'Data')
  RDDocsDir=os.path.join(RDBaseDir,'Docs')
  RDProjDir=os.path.join(RDBaseDir,'Projects')
  RDBinDir=MissingCondaDir('bin')
  RDDemoDir=MissingCondaDir('Demo')
  RDCodeDir=MissingCondaDir('rdkit')
else:
  from rdkit.RDPaths import *


rpcTestPort=8423
pythonTestCommand="python"

defaultDBUser='sysdba'
defaultDBPassword='masterkey'

class ObsoleteCodeError(Exception):
  pass
class UnimplementedCodeError(Exception):
  pass

# ---------------------
# the following block contains stuff used by the
# testing infrastructure
if sys.platform=='win32':
  pythonExe=sys.executable
else:
  pythonExe="python"

# ---------------------
# the following block contains stuff controlling database access:
usePgSQL=False
useSqlLite=False
if not os.environ.get('RD_USESQLLITE',''):
  try:
    from pyPgSQL import PgSQL
    usePgSQL=True
  except ImportError:
    usePgSQL=False
if not usePgSQL:
  try:
    # python2.5 has this:
    import sqlite3
    useSqlLite=True
  except ImportError:
    try:
      # earlier versions of python:
      from pysqlite2 import dbapi2
      useSqlLite=True
    except:
      pass

if usePgSQL:
  RDTestDatabase='::RDTests'
  RDDataDatabase='::RDData'
elif useSqlLite:
  RDTestDatabase=os.path.join(RDDataDir,"RDTests.sqlt")
  RDDataDatabase=os.path.join(RDDataDir,"RDData.sqlt")
else:  
  RDTestDatabase=None
  RDDataDatabase=None

# ---------------------
# the following block contains stuff controlling the program used for
#  3D molecular visualization:
molViewer=os.environ.get('RD_MOLVIEWER','PYMOL').upper()

  
