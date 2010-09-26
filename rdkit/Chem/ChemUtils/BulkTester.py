# $Id$
#
#  Copyright (C) 2005-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from rdkit import Chem
import sys
from rdkit.Chem import Randomize

def TestMolecule(mol):
  try:
    Chem.SanitizeMol(mol)
    mol = Chem.RemoveHs(mol)
  except ValueError,msg:
    return -1
  except:
    import traceback
    traceback.print_exc()
    return -2
  if mol.GetNumAtoms():
    try:
      Randomize.CheckCanonicalization(mol,10)
    except:
      import traceback
      traceback.print_exc()
      return -3
  return 0
    
    
def TestSupplier(suppl,stopAfter=-1,reportInterval=100,reportTo=sys.stderr,
                 nameProp='_Name'):
  nDone = 0
  nFailed = 0
  while 1:
    try:
      mol = suppl.next()
    except StopIteration:
      break
    except:
      import traceback
      traceback.print_exc()
      nFailed += 1
      reportTo.flush()
      print >>reportTo,'Failure at mol %d'%nDone
    else:
      if mol:
        ok = TestMolecule(mol)
      else:
        ok = -3
      if ok<0:
        nFailed += 1
        reportTo.flush()
        if ok==-3:
          print >>reportTo,'Canonicalization',
        print >>reportTo,'Failure at mol %d'%nDone,
        if mol:
          print >>reportTo,mol.GetProp(nameProp),
        print >>reportTo,''

          
    nDone += 1
    if nDone==stopAfter:
      break
    if not nDone%reportInterval:
      print 'Done %d molecules, %d failures'%(nDone,nFailed)
  return nDone,nFailed
if __name__=='__main__':
  suppl = Chem.SDMolSupplier(sys.argv[1],False)
  if len(sys.argv)>2:
    nameProp = sys.argv[2]
  else:
    nameProp = '_Name'
    
  nDone,nFailed = TestSupplier(suppl,nameProp=nameProp)
  print '%d failures in %d mols'%(nFailed,nDone)
