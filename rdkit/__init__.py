# Need to import rdBase to properly wrap exceptions                                                                                                         #  otherwise they will leak memory 
from . import rdBase

try:
  from .rdBase import rdkitVersion as __version__
except ImportError:
  __version__ = 'Unknown'
  raise

import logging
logger = logging.getLogger("rdkit")
# if we are running in a jupyter notebook, enable the extensions
try:
  kernel_name = get_ipython().__class__.__name__

  if kernel_name == 'ZMQInteractiveShell':
    from rdkit.Chem.Draw import IPythonConsole
    from rdkit.Chem import LogWarningMsg, WrapLogs
    WrapLogs()  
    logger.info("Enabling RDKit %s jupyter extensions"%__version__)
    
except:
  pass
