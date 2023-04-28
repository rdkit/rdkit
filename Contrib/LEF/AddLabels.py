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
import re
import sys

from rdkit import Chem
from rdkit.Chem import BRICS

inF = file(sys.argv[1], 'r')
inLs = inF.readlines()

delim = ' '

# definitions of the functional groups to look for and flag.
# format is: Name<tab>SMARTS<tab>Note
fgData = """AcidChloride	C(=O)Cl	Acid Chloride
CarboxylicAcid	C(=O)[O;H,-]	Carboxylic acid
SulfonylChloride	[$(S-!@[#6])](=O)(=O)(Cl)	Sulfonyl Chloride
Amine				[N;!H0;$(N-[#6]);!$(N-[!#6]);!$(N-C=[O,N,S])]	Amine
BoronicAcid			[$(B-!@[#6])](O)(O)		Boronic Acid
Isocyanate			[$(N-!@[#6])](=!@C=!@O)	Isocyanate
Alcohol				[O;H1;$(O-!@[#6;!$(C=!@[O,N,S])])]	Alcohol
Aldehyde			[CH;D2;!$(C-[!#6])]=O	Aldehyde
Halogen				[$([Cl,Br,I]-!@[#6]);!$([Cl,Br,I]-!@C-!@[F,Cl,Br,I]);!$([Cl,Br,I]-[C,S](=[O,S,N]))]	Halogen"""
fglines = [re.split(r'\t+', x.strip()) for x in fgData.split('\n')]
hLabels = [x[0] for x in fglines]
patts = [Chem.MolFromSmarts(x[1]) for x in fglines]

labels = inLs[0].strip().split(delim) + hLabels + ['HasBRICSBond?']
print(delim.join(labels))
for line in inLs[1:]:
  splitL = line.strip().split(delim)
  mol = Chem.MolFromSmiles(splitL[1])
  for fg in patts:
    if mol.HasSubstructMatch(fg):
      splitL.append('True')
    else:
      splitL.append('False')

  bricsRes = BRICS.BRICSDecompose(mol)
  if len(bricsRes) > 1:
    splitL.append('True')
  else:
    splitL.append('False')
  print(delim.join(splitL))
