#  Original Author: iwatobipen
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
"""
This script performs fast clustering of SMILES

Clustering method is repeated bi section, the method looks like -k-means.
To use this script, the user needs to install bayon at first.
Input format: Tab separated SMILES strings (SMILES \t molID \n ...)
Please see more details in README.

"""

import argparse
import os
import pickle
import subprocess

from rdkit import Chem
from rdkit.Chem import AllChem


def getArgParser():
  """ Create the argument parser """
  parser = argparse.ArgumentParser("Fast clustering for chemoinformatics")
  parser.add_argument("input", help="filename of input file")
  parser.add_argument("nclusters", metavar="N", help="the number of clusters")
  parser.add_argument("--output", help="filename of output, tab separated format",
                      default="clustered.tsv")
  parser.add_argument("--centroid", metavar="CENTROID",
                      help="filename of centroid information. tab separated format",
                      default="centroid.tsv")
  return parser


def smi2fp(molid, smiles):
  mol = Chem.MolFromSmiles(smiles)
  onbits = AllChem.GetMorganFingerprintAsBitVect(mol, 2).GetOnBits()
  row = molid
  for bit in onbits:
    row += "\tFP_{}\t1.0".format(bit)
  row += "\n"
  return row


if __name__ == "__main__":
  parser = getArgParser()
  args = parser.parse_args()
  with open(args.input, "r") as inputf:
    with open("fp.tsv", "w") as tempf:
      for line in inputf:
        molid, smiles = line.rstrip().split("\t")
        tempf.write(smi2fp(molid, smiles))
  res = subprocess.call(
    "time bayon -p -c {0.centroid} -n  {0.nclusters} fp.tsv > {0.output}".format(args), shell=True)

  #parse results
  parsefile = open(args.output.split(".")[0] + "_parse.tsv", "w")
  inputf = open(args.output, "r")
  for line in inputf:
    line = line.rstrip().split("\t")
    cluster_id = line[0]
    for i in range(1, len(line) - 1, 2):
      molid = line[i]
      point = line[i + 1]
      parsefile.write("{}\t{}\tCLS_ID_{}\n".format(molid, point, cluster_id))
  parsefile.close()

  if res != 0:
    parser.exit("Error running bayon")
