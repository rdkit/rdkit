"""
$Id$

Original author:  Thorsten Henninger <henni@brainbot.com>
Info about this file: http://mail.python.org/pipermail/python-win32/2003-December/001453.html

Use this in confuction with McMillan-Installer, if you need COM-Support in your application.
McMillan Installer as well as PythonCOM do overwrite the import Method for their needs.

_myimport advises the imports to either PythonCOM-Import or McMillan Import.

set_gen_path: use this method to specify the gencache-path, where all Python-COMObjects will be
              automatically generated in.

One needs to Patch the McMillan Installer such that the Python-nativ import has to be saved as
__oldimport__

Use this hack in your application with

import comhack
set_gen_path(mypath)

Then everything will work ...
"""

import sys
import __builtin__

#save the import-Method (either McMillan's import or Python's)
mcimport = __builtin__.__import__


def _myimport(name, globals=None, locals=None, fromlist=None):
  """
    Tell all modules to imported by McMillan's (or Python's) import.method,
    besides win32com modules automatically genrated by win32com.gencache
    """
  try:
    #fails with ImportError, if McMillan has overwritten Python's native import
    #and win32com.gen_py tries to generate a module
    return mcimport(name, globals, locals, fromlist)
  except ImportError as err:
    if name.startswith('win32com.gen_py'):
      #this is the Python-Native Import Method, if a patched McMillan-Installer exists
      return __oldimport__(name, globals, locals, fromlist)
    else:
      #win32com needs this ImportError
      raise err


def set_gen_path(path):
  """
    Set the gencache Path
    If not set, all Modules (win32com) will be generated to support/gen_py of your apllication.
    """
  import os

  import win32com
  import win32com.gen_py

  path = os.path.abspath(path)

  if not os.path.exists(path):
    os.makedirs(path)

  #check, if McMillan runs ...
  frozen = sys.__dict__.get("frozen", 0)

  #set the gencache path
  win32com.gen_py.__path__ = [path]
  win32com.__gen_path__ = path

  if not frozen:
    return

  #setup our import method, if McMillan runs.
  __builtin__.__import__ = _myimport
