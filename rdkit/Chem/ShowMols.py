# $Id$
#
# Copyright (C) 2004-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from rdkit import RDConfig

# change this to use another viewer:
if RDConfig.molViewer in ('WEBLAB', 'DSVIEWER'):
  from rdkit.Chem.DSViewer import *
elif RDConfig.molViewer == 'PYMOL':
  from rdkit.Chem.PyMol import *
else:
  raise ValueError('invalid RD_MOLVIEWER specified')

if __name__ == '__main__':
  import sys

  import AllChem
  if len(sys.argv) < 2:
    smi = 'c1cccc2c1cccc2CC(=O)N'
  else:
    smi = sys.argv[1]

  m = Chem.MolFromSmiles(smi)
  m = Chem.AddHs(m)
  AllChem.EmbedMolecule(m)
  v = MolViewer()
  v.ShowMol(m, 'raw')
  AllChem.UFFOptimizeMolecule(m)
  v.ShowMol(m, 'opt', showOnly=0, highlightFeatures=[(0, ), (2, ), (3, 4)])
