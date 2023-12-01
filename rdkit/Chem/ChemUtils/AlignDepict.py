#
#  Copyright (C) 2006 Greg Landrum
#  This file is part of RDKit and covered by $RDBASE/license.txt
#
import argparse
import sys

from rdkit import Chem, Geometry
from rdkit.Chem import rdDepictor


def AlignDepict(mol, core, corePattern=None, acceptFailure=False):
  """
  Arguments:
    - mol:          the molecule to be aligned, this will come back
                    with a single conformer.
    - core:         a molecule with the core atoms to align to;
                    this should have a depiction.
    - corePattern:  (optional) an optional molecule to be used to
                    generate the atom mapping between the molecule
                    and the core.
  """
  if core and corePattern:
    if not core.GetNumAtoms(onlyExplicit=True) == corePattern.GetNumAtoms(onlyExplicit=True):
      raise ValueError(
        'When a pattern is provided, it must have the same number of atoms as the core')
    coreMatch = core.GetSubstructMatch(corePattern)
    if not coreMatch:
      raise ValueError("Core does not map to itself")
  else:
    coreMatch = list(range(core.GetNumAtoms(onlyExplicit=True)))

  if corePattern:
    match = mol.GetSubstructMatch(corePattern)
  else:
    match = mol.GetSubstructMatch(core)

  if not match:
    if not acceptFailure:
      raise ValueError('Substructure match with core not found.')
    else:
      coordMap = {}
  else:
    conf = core.GetConformer()
    coordMap = {}
    for i, idx in enumerate(match):
      pt3 = conf.GetAtomPosition(coreMatch[i])
      coordMap[idx] = Geometry.Point2D(pt3.x, pt3.y)
  rdDepictor.Compute2DCoords(mol, clearConfs=True, coordMap=coordMap, canonOrient=False)


def initParser():
  """ Initialize the parser """
  parser = argparse.ArgumentParser(description='Create aligned depiction')
  parser.add_argument('--pattern', '-p', metavar='SMARTS', default=None, dest='patt')
  parser.add_argument('--smiles', default=False, action='store_true', dest='useSmiles',
                      help='Set if core and input are SMILES strings')
  parser.add_argument('-o', dest='outF', type=argparse.FileType('w'), default=sys.stdout,
                      metavar='OUTFILE',
                      help='Specify a file to take the output. If missing, uses stdout.')
  parser.add_argument('core', metavar='core')
  parser.add_argument('mol', metavar='molecule', help='')
  return parser


def processArgs(args):
  patt = args.patt
  if patt:
    patt = Chem.MolFromSmarts(patt)

  if args.useSmiles:
    core = Chem.MolFromSmiles(args.core)
    mol = Chem.MolFromSmiles(args.mol)
    rdDepictor.Compute2DCoords(core)
  else:
    core = Chem.MolFromMolFile(args.core)
    mol = Chem.MolFromMolFile(args.mol)

  AlignDepict(mol, core, patt)
  print(Chem.MolToMolBlock(mol), file=args.outF)


def main():
  """ Main application """
  parser = initParser()
  args = parser.parse_args()
  processArgs(args)


if __name__ == '__main__':
  main()
