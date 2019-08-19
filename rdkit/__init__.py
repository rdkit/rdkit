try:
  from .rdBase import rdkitVersion as __version__
except ImportError:
  __version__ = 'Unknown'
  raise

# if we are running in a jupyter notebook, enable the extensions
try:
  get_ipython()
  from rdkit.Chem.Draw import IPythonConsole
  from rdkit.Chem import LogWarningMsg, WrapLogs
  WrapLogs()
  LogWarningMsg("Enabling RDKit %s jupyter extensions"%__version__)
except:
  pass
