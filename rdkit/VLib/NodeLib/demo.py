#  $Id$
#
#  Copyright (C) 2003 Rational Discovery LLC
#     All Rights Reserved
#

import os.path

from rdkit import Chem, RDConfig
from rdkit.VLib import Filter, Supply
from rdkit.VLib.NodeLib import *

# this would be a real input, from an sd file:
#fName = os.path.join(RDConfig.RDCodeDir,'VLib','NodeLib','test_data','NCI_aids.10.dupes.sdf')
#supplier = SDSupply.SDSupplyNode(fName)
# instead though, we want a simpler input:
smis = [
  'CCOC', 'CCO.Cl', 'CC(=O)[O-].[Na+]', 'CC[Cu]CC', 'OCC', 'C[N+](C)(C)C.[Cl-]', '[Na+].[Cl-]'
]
mols = [Chem.MolFromSmiles(x) for x in smis]
# name the molecules (only needed because we built them from smiles):
for i in range(len(mols)):
  mols[i].SetProp('Name', 'Mol-%d' % (i + 1))
supplier = Supply.SupplyNode(contents=mols)
# should be 7 here
print('initial:', len([x for x in supplier]))

# filter out anything with a transition metal or lanthanide:
metals = '[#21,#22,#23,#24,#25,#26,#27,#28,#29,#39,#40,#41,#42,#43,#44,#45,#46,#47,#57,#58,#59,#60,#61,#62,#63,#64,#65,#66,#67,#68,#69,#70,#71,#72,#73,#74,#75,#76,#77,#78,#79]'
smaFilter = SmartsMolFilter.SmartsFilter(patterns=[metals], counts=[1])
smaFilter.SetNegate(1)
smaFilter.AddParent(supplier)
# should be 6 here
print('post-smaFilter:', len([x for x in smaFilter]))

salts = ['[Cl;H1&X1,-]', '[Na+]', '[O;H2,H1&-,X0&-2]']
remover = SmartsRemover.SmartsRemover(patterns=salts)
remover.AddParent(smaFilter)
atsFilter = Filter.FilterNode(func=lambda x: x.GetNumAtoms() > 1)
atsFilter.AddParent(remover)
# should be 5 here
print('post-remover:', len([x for x in atsFilter]))

dupeFilter = SmilesDupeFilter.DupeFilter()
dupeFilter.AddParent(atsFilter)
# should be 4 here
print('post-dupes:', len([x for x in dupeFilter]))

from io import StringIO
# a StringIO object acts like a file:
io = StringIO()
output = SmilesOutput.OutputNode(dest=io, delim=', ', idField='Name')
output.AddParent(dupeFilter)
print('post-output:', len([x for x in output]))
print('OUTPUT:')
print(io.getvalue())
