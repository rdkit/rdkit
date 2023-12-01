# Copyright (c) 2013, GlaxoSmithKline Research & Development Ltd.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of GlaxoSmithKline Research & Development Ltd.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Created by Jameed Hussain, July 2013

import re
import sys
from optparse import OptionParser

from rdkit import Chem
from rdkit.Chem import AllChem


def smiles_to_smarts(smi):
  mol = Chem.MolFromSmiles(smi)

  if (mol is None):
    sys.stderr.write("Can't generate mol for: %s\n" % (smi))
    return None

  #change the isotope to 42
  for atom in mol.GetAtoms():
    atom.SetIsotope(42)

  #preint out the smiles - all the atom attributes will be fully specified
  smarts = Chem.MolToSmiles(mol, isomericSmiles=True)
  #remove the 42 isotope labels
  smarts = re.sub(r'\[42', "[", smarts)
  #now have a fully specified SMARTS - simples!

  return smarts


if __name__ == '__main__':

  parser = OptionParser(
    description="Program to apply transformations to a set of input molecules",
    epilog="Example command: mol_transform.py -t TRANSFORM_FILE <SMILES_FILE\t\t"
    "Format of smiles file: SMILES ID <space or comma separated>\t\t\t"
    "Format of transform file: transform <one per line>\t\t\t"
    "Output: SMILES,ID,Transfrom,Modified_SMILES")
  parser.add_option('-f', '--file', action='store', dest='transform_file', type='string',
                    help='The file containing the transforms to apply to your input SMILES')

  (options, args) = parser.parse_args()

  #print options.transform_file
  if options.transform_file is None:
    print("Please specify the transform file.")
    sys.exit(1)

  smiles = []
  #read the STDIN
  for line in sys.stdin:
    line = line.rstrip()
    smi, id = re.split(r'\s|,', line)
    #print smiles,id
    smiles.append((smi, id))

  #read the transform file
  #all the transform must come from BioDig to guarantee they have been cansmirk'ed
  infile = open(options.transform_file, 'r')

  print("Input_SMILES,ID,RG-Transform,RG-transformedSMILES")
  for transform in infile:
    transform = transform.rstrip()

    #need to convert the smiles to smart to get rid of any potential issues
    lhs, rhs = transform.split(">>")

    if (lhs == "[*:1][H]"):
      lhs = "[*;!H0:1]"
    else:
      lhs = smiles_to_smarts(lhs)

    if (rhs == "[*:1][H]"):
      rhs = "[*:1]"
    else:
      rhs = smiles_to_smarts(rhs)

    rdkit_transform = "%s>>%s" % (lhs, rhs)
    rxn = AllChem.ReactionFromSmarts(rdkit_transform)
    #rxn = AllChem.ReactionFromSmarts(transform)

    for x in smiles:
      mol = Chem.MolFromSmiles(x[0])
      ps = rxn.RunReactants([mol])

      products = set()
      for y in range(len(ps)):
        for z in range(len(ps[y])):
          p = ps[y][z]
          Chem.SanitizeMol(p)
          products.add(Chem.MolToSmiles(p, isomericSmiles=True))

      for p in products:
        print("%s,%s,%s,%s" % (x[0], x[1], transform, p))
