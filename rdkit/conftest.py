# additional pytest configuration

collect_ignore = []
collect_ignore_glob = []
try:
  import win32com
except ImportError:
  collect_ignore.extend(['Chem/DSViewer.py'])

try:
  from rdkit.Avalon import pyAvalonTools
except ImportError:
  collect_ignore_glob.extend(['**Chem/MolKey/*.py'])

from rdkit.Chem import inchi
if not inchi.INCHI_AVAILABLE:
  collect_ignore_glob.extend(['**Chem/MolKey/*.py'])

try:
  from rdkit.Chem import rdDASHTree
except ImportError:  # not built
  collect_ignore.append('Chem/DASHTree')

try:
  import IPython
except ImportError:
  collect_ignore.extend([
    'Chem/Draw/IPythonConsole.py', 'Chem/Draw/InteractiveRenderer.py',
    'Chem/Draw/UnitTestIPython.py'
  ])
