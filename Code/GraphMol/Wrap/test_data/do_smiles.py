
import re
splitExpr = re.compile('[\t\ ]')

from rdkit import Chem


def runit(fName):
  inLines = open(fName, 'r').readlines()
  nFailed = 0
  nPassed = 0
  nTried = 0
  for line in inLines:
    if len(line):
      smi = splitExpr.split(line)[1]
      if smi[-1] == '\n':
        smi = smi[:-1]
      if smi[-1] == '\r':
        smi = smi[:-1]
      nTried += 1
      m = Chem.MolFromSmiles(smi)
      if m:
        nPassed += 1
      else:
        print('\t%s failed' % repr(smi))
        print('\tline: %s' % (repr(line)))
        nFailed += 1
      m = None
  print('%d of %d passed' % (nPassed, nTried))


if __name__ == '__main__':
  import sys
  fName = 'ntp_smiles.txt'
  if len(sys.argv) > 1:
    fName = sys.argv[1]

  runit(fName)
