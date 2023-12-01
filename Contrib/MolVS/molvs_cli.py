"""
molvs.cli
~~~~~~~~~

This module contains a command line interface for standardization.

:copyright: Copyright 2016 by Matt Swain.
:license: MIT, see LICENSE file for more details.

*Adapted for purposes of integration of MolVS into RDKit
"""

import argparse
import logging
import sys

from rdkit import Chem
from rdkit.Chem.MolStandardize import Standardizer, Validator

log = logging.getLogger(__name__)

FILETYPES = ['smi', 'mol', 'sdf']


class MolvsParser(argparse.ArgumentParser):

  def error(self, message):
    sys.stderr.write('Error: %s\n\n'.encode() % message)
    self.print_help()
    sys.exit(2)


def _read_mol(args):
  if args.smiles:
    return Chem.MolFromSmiles(args.smiles)
  elif args.intype in {'smi', 'smiles'
                       } or args.infile.name.endswith('smi') or args.infile.name.endswith('smiles'):
    return Chem.MolFromSmiles(args.infile.read())
  elif args.intype in {'mol', 'sdf'
                       } or args.infile.name.endswith('mol') or args.infile.name.endswith('sdf'):
    return Chem.MolFromMolBlock(args.infile.read())
  else:
    return Chem.MolFromSmiles(args.infile.read())


def _write_mol(mol, args):
  if args.outtype in {
      'smi', 'smiles'
  } or args.outfile.name.endswith('smi') or args.outfile.name.endswith('smiles'):
    args.outfile.write(Chem.MolToSmiles(mol))
    args.outfile.write('\n')
  elif args.outtype in {'mol', 'sdf'
                        } or args.outfile.name.endswith('mol') or args.outfile.name.endswith('sdf'):
    args.outfile.write(Chem.MolToMolBlock(mol))
  else:
    args.outfile.write(Chem.MolToSmiles(mol))
    args.outfile.write('\n')


def standardize_main(args):
  mol = _read_mol(args)
  s = Standardizer()
  mol = s.standardize(mol)
  _write_mol(mol, args)


def validate_main(args):
  mol = _read_mol(args)
  v = Validator()
  logs = v.validate(mol)
  for log in logs:
    args.outfile.write(log)
    args.outfile.write('\n')


if __name__ == '__main__':
  """Main function for molvs command line interface."""

  # Root options
  #    parser = MolvsParser(epilog='use "molvs <command> -h" to show help for a specific command')
  parser = MolvsParser(usage="usage: python cli.py [-h] {standardize,validate} ...")
  subparsers = parser.add_subparsers(title='Available commands')

  # Options common to all commands

  common_parser = MolvsParser(add_help=False)
  common_parser.add_argument('infile', nargs='?', help='input filename',
                             type=argparse.FileType('r'), default=sys.stdin)
  common_parser.add_argument('-i', '--intype', help='input filetype', choices=FILETYPES)
  common_parser.add_argument('-:', '--smiles', help='input SMILES instead of file',
                             metavar='<smiles>')
  common_parser.add_argument('-O', '--outfile', help='output filename', type=argparse.FileType('w'),
                             default=sys.stdout, metavar='<outfile>')

  # Standardize options
  standardize_parser = subparsers.add_parser('standardize', help='standardize a molecule',
                                             parents=[common_parser])
  standardize_parser.add_argument('-o', '--outtype', help='output filetype', choices=FILETYPES)
  standardize_parser.set_defaults(func=standardize_main)

  # Validate options
  validate_parser = subparsers.add_parser('validate', help='validate a molecule',
                                          parents=[common_parser])
  validate_parser.set_defaults(func=validate_main)

  args = parser.parse_args()
  try:
    args.func(args)
  except Exception as e:
    sys.stderr.write('Error: %s\n\n'.encode() % e.message)
    parser.print_help()
    sys.exit(2)
