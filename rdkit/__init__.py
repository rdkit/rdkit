import logging
import sys

# Need to import rdBase to properly wrap exceptions
# otherwise they will leak memory
from . import rdBase

try:
  from .rdBase import rdkitVersion as __version__
except ImportError:
  __version__ = 'Unknown'
  raise

logger = logging.getLogger("rdkit")

# if we are running in a jupyter notebook, enable the extensions
try:
  kernel_name = get_ipython().__class__.__name__

  if kernel_name == 'ZMQInteractiveShell':
    logger.info("Enabling RDKit %s jupyter extensions" % __version__)
    from rdkit.Chem.Draw import IPythonConsole
    rdBase.LogToPythonStderr()
except Exception:
  pass

# Do logging setup at the end, so users can suppress the
# "enabling jupyter" message at the root logger.
log_handler = logging.StreamHandler(sys.stderr)
logger.addHandler(log_handler)
logger.setLevel(logging.WARN)
logger.propagate = False

# Uncomment this to use Python logging by default:
# rdBase.LogToPythonLogger()
