'''

Script will test the RDkit python code for conformance with the agreed format using
yapf.

For each Python file that is found in $RDBASE (excluding the build and External
directories), yapf is used with the style configuration in $RDBASE/setup.cfg.
If a change is required, the difference is printed. At the end of the process,
all non-conformant files are listed and the required yapf command(s) printed.

If changes are found, the script will exit with error code 1, otherwise 0.

'''

import os
from yapf.yapflib.yapf_api import FormatCode
import sys

rdbase = os.environ.get('RDBASE', '')
styleConfig = os.path.join(rdbase, 'setup.cfg')

excludeDirs = [os.path.join(rdbase, 'build'),
               os.path.join(rdbase, 'External'), ]


def pythonFiles(dirname=rdbase):
  """ Find all python files below directory dirname """
  for root, _, files in os.walk(dirname):
    if any(root.startswith(d) for d in excludeDirs):
      continue
    for file in files:
      if file.endswith(".py"):
        yield os.path.join(root, file)


def yapfChanges(filename):
  """ Use yapf with the default settings to format file filename """
  try:
    with open(filename) as f:
      codeBefore = f.read()
  except UnicodeError:
    with open(filename, encoding='latin-1') as f:
      codeBefore = f.read()
  try:
    changes, changed = FormatCode(codeBefore, style_config=styleConfig, print_diff=True,
                                  filename=filename)
  except Exception:
    print(filename)
    raise
  if changed:
    print(changes)
  return changed


if __name__ == "__main__":
  changedFiles = []
  for s in pythonFiles():
    if yapfChanges(s):
      changedFiles.append(s)
  print()
  if changedFiles:
    print('yapf will make changes to the following files:')
    print('\n'.join(sorted(changedFiles)))
    print('To apply the required changes to your code use the following command(s)')
    for s in sorted(set(s.replace(rdbase, '').split(os.sep)[1] for s in changedFiles)):
      print('yapf --style $RDBASE/setup.cfg --in-place --recursive $RDBASE/{0}'.format(s))
    sys.exit(1)
  print('Code complies with the agreed formatting rules.')
  sys.exit(0)
