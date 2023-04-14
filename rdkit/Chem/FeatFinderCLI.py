#
#  Copyright (C) 2005-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

import argparse
import os
import re

from rdkit import Chem, RDLogger
from rdkit.Chem import ChemicalFeatures

logger = RDLogger.logger()
splitExpr = re.compile(r'[ \t,]')


def GetAtomFeatInfo(factory, mol):
  res = [None] * mol.GetNumAtoms()
  feats = factory.GetFeaturesForMol(mol)
  for feat in feats:
    ids = feat.GetAtomIds()
    feature = "%s-%s" % (feat.GetFamily(), feat.GetType())
    for id_ in ids:
      if res[id_] is None:
        res[id_] = []
      res[id_].append(feature)
  return res


def initParser():
  """ Initialize the parser """
  parser = argparse.ArgumentParser(description='Determine pharmacophore features of molecules',
                                   epilog=_splashMessage,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument('-r', dest='reverseIt', default=False, action='store_true',
                      help='Set to get atoms lists for each feature.')
  parser.add_argument('-n', dest='maxLines', default=-1, help=argparse.SUPPRESS, type=int)
  parser.add_argument('fdefFilename', type=existingFile,
                      help='Pharmacophore feature definition file')
  parser.add_argument('smilesFilename', type=existingFile,
                      help='The smiles file should have SMILES in the first column')
  return parser


_splashMessage = """
-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  FeatFinderCLI
  Part of the RDKit (http://www.rdkit.org)
-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
"""


def existingFile(filename):
  """ 'type' for argparse - check that filename exists """
  if not os.path.exists(filename):
    raise argparse.ArgumentTypeError("{0} does not exist".format(filename))
  return filename


def processArgs(args, parser):
  try:
    factory = ChemicalFeatures.BuildFeatureFactory(args.fdefFilename)
  except Exception:
    parser.error("Could not parse Fdef file {0.fdefFilename}.".format(args))

  with open(args.smilesFilename) as inF:
    for lineNo, line in enumerate(inF, 1):
      if lineNo == args.maxLines + 1:
        break
      smi = splitExpr.split(line.strip())[0].strip()
      mol = Chem.MolFromSmiles(smi)
      if mol is None:
        logger.warning("Could not process smiles '%s' on line %d." % (smi, lineNo))
        continue

      print('Mol-%d\t%s' % (lineNo, smi))
      if args.reverseIt:
        feats = factory.GetFeaturesForMol(mol)
        for feat in feats:
          print('\t%s-%s: ' % (feat.GetFamily(), feat.GetType()), end='')
          print(', '.join([str(x) for x in feat.GetAtomIds()]))
      else:
        featInfo = GetAtomFeatInfo(factory, mol)
        for i, v in enumerate(featInfo):
          print('\t% 2s(%d)' % (mol.GetAtomWithIdx(i).GetSymbol(), i + 1), end='')
          if v:
            print('\t', ', '.join(v))
          else:
            print()


def main():
  """ Main application """
  parser = initParser()
  args = parser.parse_args()
  processArgs(args, parser)


if __name__ == '__main__':
  main()
