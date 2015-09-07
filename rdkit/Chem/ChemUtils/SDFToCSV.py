# $Id$
#
# Copyright (C) 2004-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from rdkit import RDConfig
from rdkit import Chem
import sys,csv

def Convert(suppl,outFile,keyCol='',stopAfter=-1,includeChirality=0,smilesFrom=''):
  w = csv.writer(outFile)
  mol = suppl[0]
  propNames = list(mol.GetPropNames())
  if keyCol and keyCol in propNames:
    propNames.remove(keyCol)

  outL = []
  if keyCol:
    outL.append(keyCol)
  outL.append('SMILES')
  outL.extend(propNames)
  w.writerow(outL)
  nDone = 0
  for mol in suppl:
    if not mol:
      continue
    if not smilesFrom or not mol.HasProp(smilesFrom):
      smi = Chem.MolToSmiles(mol,includeChirality)
    else:
      smi = mol.GetProp(smilesFrom)
      tMol = Chem.MolFromSmiles(smi)
      smi = Chem.MolToSmiles(tMol,includeChirality)
    outL = []
    if keyCol:
      outL.append(str(mol.GetProp(keyCol)))
    outL.append(smi)
    for prop in propNames:
      if mol.HasProp(prop):
        outL.append(str(mol.GetProp(prop)))
      else:
        outL.append('')
    w.writerow(outL)
    nDone += 1
    if nDone == stopAfter:
      break
  return


#-------------------
#  Testing:
import unittest
class TestCase(unittest.TestCase):
  def setUp(self):
    pass
  def tearDown(self):
    pass
  def test1(self):
    import os
    from rdkit.six.moves import cStringIO as StringIO  #@UnresolvedImport #pylint: disable=F0401
    fName = os.path.join(RDConfig.RDDataDir,'NCI','first_200.props.sdf')
    suppl = Chem.SDMolSupplier(fName)
    io = StringIO()
    try:
      Convert(suppl,io)
    except Exception:
      import traceback
      traceback.print_exc()
      self.fail('conversion failed')
    txt = io.getvalue()
    lines = txt.split('\n')
    if not lines[-1]:
      del lines[-1]
    self.assertTrue(len(lines)==201,'bad num lines: %d'%len(lines))
    line0 = lines[0].split(',')
    self.assertEqual(len(line0),20)
    self.assertTrue(line0[0]=='SMILES')
  def test2(self):
    import os
    from rdkit.six.moves import cStringIO as StringIO  #@UnresolvedImport #pylint: disable=F0401
    fName = os.path.join(RDConfig.RDDataDir,'NCI','first_200.props.sdf')
    suppl = Chem.SDMolSupplier(fName)
    io = StringIO()
    try:
      Convert(suppl,io,keyCol='AMW',stopAfter=5)
    except Exception:
      import traceback
      traceback.print_exc()
      self.fail('conversion failed')
    txt = io.getvalue()
    lines = txt.split('\n')
    if not lines[-1]:
      del lines[-1]
    self.assertTrue(len(lines)==6,'bad num lines: %d'%len(lines))
    line0 = lines[0].split(',')
    self.assertEqual(len(line0),20)
    self.assertTrue(line0[0]=='AMW')
    self.assertTrue(line0[1]=='SMILES')
    
    
    
  



#-------------------
#  CLI STuff:
def Usage():
  message = """
  Usage: SDFToCSV [-k keyCol] inFile.sdf [outFile.csv]

  """
  sys.stderr.write(message)
  sys.exit(-1)



if __name__=='__main__':
  import getopt

  try:
    args,extras = getopt.getopt(sys.argv[1:],'hk:',
                                ['test',
                                 'chiral',
                                 'smilesCol=',
                                 ])
  except Exception:
    import traceback
    traceback.print_exc()
    Usage()

  keyCol = ''
  testIt = 0
  useChirality=0
  smilesCol=''
  for arg,val in args:
    if arg=='-k':
      keyCol = val
    elif arg=='--chiral':
      useChirality=1
    elif arg=='--smilesCol':
      smilesCol=val
    elif arg=='--test':
      testIt=1
    elif arg=='-h':
      Usage()

  if not testIt and len(extras)<1:
    Usage()
      

  if not testIt:
    inFilename = extras[0]
    if len(extras)>1:
      outFilename = extras[1]
      outF = open(outFilename,'w+')
    else:
      outF = sys.stdout

    suppl = Chem.SDMolSupplier(inFilename)
    Convert(suppl,outF,keyCol=keyCol,includeChirality=useChirality,smilesFrom=smilesCol)
  else:
    sys.argv = [sys.argv[0]]
    unittest.main()

