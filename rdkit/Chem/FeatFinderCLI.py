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
from __future__ import print_function

_version = "$Rev$"
_splashMessage="""
-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  FeatFinderCLI version %s

  Copyright (C) 2005 Rational Discovery LLC

  This software is copyrighted.  The software may not be copied,
  reproduced, translated or reduced to any electronic medium or
  machine-readable form without the prior written consent of
  Rational Discovery LLC.
-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
"""%_version
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit import RDLogger
logger = RDLogger.logger()
import sys,os
import re
splitExpr = re.compile(r'[ \t,]')

def GetAtomFeatInfo(factory,mol):
  res = [None]*mol.GetNumAtoms()
  feats = factory.GetFeaturesForMol(mol)
  for feat in feats:
    ids = feat.GetAtomIds()
    for id in ids:
      if res[id] is None:
        res[id] = []
      res[id].append("%s-%s"%(feat.GetFamily(),feat.GetType()))
  return res

if __name__ == '__main__':
  def Usage():
    message="""
    Usage: FeatFinderCLI [-r] <fdefFilename> <smilesFilename>
    
    NOTE:
      - the smiles file should have SMILES in the first column

    """
    print(message, file=sys.stderr)


  import getopt
  args,extras = getopt.getopt(sys.argv[1:],'r')
  reverseIt=False
  for arg,val in args:
    if arg=='-r':
      reverseIt=True
      
  if len(extras)<2:
    Usage()
    sys.exit(-1)
  print(_splashMessage, file=sys.stderr)
  fdefFilename = extras[0]
  if not os.path.exists(fdefFilename):
    logger.error("Fdef file %s does not exist."%fdefFilename)
    sys.exit(-1)
  try:
    factory = ChemicalFeatures.BuildFeatureFactory(fdefFilename)
  except Exception:
    logger.error("Could not parse Fdef file %s."%fdefFilename,exc_info=True)
    sys.exit(-1)
    
  smilesFilename = extras[1]
  if not os.path.exists(smilesFilename):
    logger.error("Smiles file %s does not exist."%smilesFilename)
    sys.exit(-1)

  try:
    inF = file(smilesFilename,'r')
  except Exception:
    logger.error("Could not open smiles file %s."%smilesFilename,exc_info=True)
    sys.exit(-1)

  lineNo=0
  for line in inF.readlines():
    lineNo+=1
    line = line.strip()
    smi = splitExpr.split(line)[0].strip()
    mol = Chem.MolFromSmiles(smi)

    if mol is not None:
      print('Mol-%d\t%s'%(lineNo,smi))

      if not reverseIt:
        featInfo = GetAtomFeatInfo(factory,mol)
        for i,v in enumerate(featInfo):
          print('\t% 2s(%d)'%(mol.GetAtomWithIdx(i).GetSymbol(),i+1),end='')
          if v:
            print('\t',', '.join(v))
          else:
            print()
      else:
        feats = factory.GetFeaturesForMol(mol)
        for feat in feats:
          print('\t%s-%s: '%(feat.GetFamily(),feat.GetType()),end='')
          print(', '.join([str(x) for x in feat.GetAtomIds()]))
    else:
      logger.warning("Could not process smiles '%s' on line %d."%(smi,lineNo))


      
    
