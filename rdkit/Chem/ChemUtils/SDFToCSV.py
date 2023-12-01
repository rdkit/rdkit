#
# Copyright (C) 2004-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import argparse
import csv
import os
import sys

from rdkit import Chem


def Convert(suppl, outFile, keyCol=None, stopAfter=-1, includeChirality=False, smilesFrom=''):
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

  for nDone, mol in enumerate(suppl, 1):
    if not mol:
      continue
    if not smilesFrom or not mol.HasProp(smilesFrom):
      smi = Chem.MolToSmiles(mol, isomericSmiles=includeChirality)
    else:
      smi = mol.GetProp(smilesFrom)
      smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi), isomericSmiles=includeChirality)

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
    if nDone == stopAfter:
      break
  return


def initParser():
  """ Initialize the parser for the CLI """
  parser = argparse.ArgumentParser(description='Convert SDF file to CSV',
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument('--key', '-k', metavar='keyCol', default=None, dest='keyCol')
  parser.add_argument('--chiral', default=False, action='store_true', dest='useChirality')
  parser.add_argument('--smilesCol', metavar='smilesCol', default='')
  parser.add_argument('inFilename', metavar='inFile.sdf', type=existingFile)
  parser.add_argument('outF', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
  return parser


def existingFile(filename):
  """ 'type' for argparse - check that filename exists """
  if not os.path.exists(filename):
    raise argparse.ArgumentTypeError("{0} does not exist".format(filename))
  return filename


def main():
  """ Main application """
  parser = initParser()
  args = parser.parse_args()
  suppl = Chem.SDMolSupplier(args.inFilename)
  Convert(suppl, args.outF, keyCol=args.keyCol, includeChirality=args.useChirality,
          smilesFrom=args.smilesCol)


if __name__ == '__main__':
  main()
