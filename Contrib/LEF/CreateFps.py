#
#  Copyright (c) 2009, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
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
#     * Neither the name of Novartis Institutes for BioMedical Research Inc.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written permission.
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
# Created by Greg Landrum and Anna Vulpetti, March 2009

import pickle
import sys

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AtomPairs import Pairs, Torsions

# maxPathLength is the maximum path length in atoms
# maxPathLength=6 corresponds to F-FP-5
# maxPathLength=7 corresponds to F-FP-6
# maxPathLength=8 corresponds to F-FP-7
maxPathLength = 8

# nameField is the name of the property (from the SD file) that has molecule
# names... If the molecules have names in the first row of the file, use "_Name"
nameField = 'Compound_orig'
#nameField = '_Name'

extraQueries = (
  ('SCF3?', Chem.MolFromSmarts('SC(F)(F)F')),
  ('COCF3?', Chem.MolFromSmarts('C(=O)C(F)(F)F')),
  ('OCF3?', Chem.MolFromSmarts('OC(F)(F)F')),
  ('NCF3?', Chem.MolFromSmarts('NC(F)(F)F')),
  ('CF3?', Chem.MolFromSmarts('C(F)(F)F')),
)


def GetMolFingerprint(mol, maxPathLength):
  FQuery = Chem.MolFromSmarts('F')
  CF3Query = Chem.MolFromSmarts('[$(C(F)(F)F)]')
  CF3Rxn = AllChem.ReactionFromSmarts('[*:1]-C(F)(F)F>>[*:1]-F')
  hasCF3 = mol.HasSubstructMatch(CF3Query)
  if hasCF3:
    p = CF3Rxn.RunReactants((mol, ))[0][0]
    Chem.SanitizeMol(p)
    for nm in mol.GetPropNames():
      p.SetProp(nm, mol.GetProp(nm))
    mol = p
  match = mol.GetSubstructMatch(FQuery)
  fp = Torsions.GetHashedTopologicalTorsionFingerprint(mol, nBits=9192, targetSize=maxPathLength,
                                                       fromAtoms=match)
  for i in range(2, maxPathLength):
    nfp = Torsions.GetHashedTopologicalTorsionFingerprint(mol, nBits=9192, targetSize=i,
                                                          fromAtoms=match)
    for bit, v in nfp.GetNonzeroElements().iteritems():
      fp[bit] = fp[bit] + v
  return fp


if __name__ == '__main__':
  suppl = Chem.SDMolSupplier(sys.argv[1])
  outF = file(sys.argv[2], 'w+')
  fps = []

  for i, mol in enumerate(suppl):
    if not mol:
      continue
    smi = Chem.MolToSmiles(mol, True)
    queryMatches = [str(mol.HasSubstructMatch(y)) for x, y in extraQueries]

    fp = GetMolFingerprint(mol, maxPathLength)

    nm = mol.GetProp(nameField)
    fps.append([nm, smi, fp] + queryMatches)
  colNames = ['name', 'smiles', 'fp'] + [x for x, y in extraQueries]
  pickle.dump(colNames, outF)
  pickle.dump(fps, outF)

  print('name1 smiles1 name2 smiles2 name12 smiles12 environment_id ' +
        ' '.join([x for x, y in extraQueries]))
  if 1:
    seen = []
    smis = []
    data = []
    for row in fps:
      nm = row[0]
      smi = row[1]
      fp = row[2]
      if fp in seen and smi not in smis:
        id = seen.index(fp)
        onm, osmi = data[id]
        print(nm, smi, onm, osmi, nm + '.' + onm, smi + '.' + osmi, id + 1, ' '.join(row[3:]))
      else:
        seen.append(fp)
        smis.append(smi)
        data.append((nm, smi))
  else:
    smis = []
    for nm, smi, fp in fps:
      if smi not in smis:
        pass
