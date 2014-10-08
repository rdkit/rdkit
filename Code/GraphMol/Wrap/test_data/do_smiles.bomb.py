from __future__ import print_function

import re
splitExpr = re.compile('[\t\ ]')

from Chem import rdmol

def runit(fName):
  inLines = open(fName,'r').readlines()
  nFailed = 0
  nPassed = 0
  nTried = 0
  for line in inLines:
    if len(line):
      smi = splitExpr.split(line)[1]
      if smi[-1] == '\n': smi = smi[:-1]
      if smi[-1] == '\r': smi = smi[:-1]
      nTried += 1
      print('trying: "%s"'%smi)
      m = rdmol.MolFromSmiles(smi)
      nPassed += 1

      print('frags')
      nFrags = max(rdmol.FindMolFrags(m))+1
      print('sssr')
      rings = rdmol.FindSSSR(m)
      cyclomat = m.getNumBonds() - m.getNumAtoms() + nFrags;
      assert cyclomat == len(rings),'bad num rings for %s\n\t%s!=%s'%(smi,cyclomat,len(rings))
      print('span')
      nChords = m.getNumBonds() - len(rdmol.FindSpanningTree(m))
      assert len(rings) == nChords, 'bad num chords for %s\n\t%s!=%s'%(smi,nChords,len(rings))
      m = None        

  print('%d of %d passed'%(nPassed,nTried))

if __name__ == '__main__':
  import sys
  fName = 'ntp_smiles.txt'
  if len(sys.argv) > 1:
    fName = sys.argv[1]

  runit(fName)
