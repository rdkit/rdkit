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
  module_name = get_ipython().__class__.__module__

  if kernel_name == 'ZMQInteractiveShell' or module_name == 'google.colab._shell':
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

#
#  Boost vector iterators are VERY slow, but indexing is quite fast
#   this wraps the iterators into a much faster version
#

class VectIter:
    def __init__(self, vect):
        self.vect = vect
        self.l = len(vect)
        self.pos = -1

    def __iter__(self):
        return self
    
    def __next__(self):
        self.pos += 1
        if self.pos < self.l:
            return self.vect[self.pos]
        else:
            raise StopIteration()

def __vect__iter__(vect):
    return VectIter(vect)

VECT_WRAPS = {'MatchTypeVect', 'UnsignedLong_Vect', 'VectSizeT', 'VectorOfStringVectors'}

for name, object in vars(rdBase).items():
  if name.startswith("_list") or name.startswith("_vect") or name in VECT_WRAPS:
    if hasattr(object, "__iter__"):
      object.__iter__ = __vect__iter__

