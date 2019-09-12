# This file is part of the RingDecomposerLib, licensed
# under BSD New license (see LICENSE in the root directory).
# Copyright (c) 2016
# University of Hamburg, ZBH - Center for Bioinformatics
# Niek Andresen, Florian Flachsenberg, Matthias Rarey
# 
# Please cite:
# 
# Kolodzik, A.; Urbaczek, S.; Rarey, M.
# Unique Ring Families: A Chemically Meaningful Description
# of Molecular Ring Topologies.
# J. Chem. Inf. Model., 2012, 52 (8), pp 2013-2021
# 
# Flachsenberg, F.; Andresen, N.; Rarey, M.
# RingDecomposerLib: An Open-Source Implementation of
# Unique Ring Families and Other Cycle Bases.
# J. Chem. Inf. Model., 2017, 57 (2), pp 122-126

from __future__ import print_function
import py_rdl
import rdkit.Chem

# p-cyclophane as URF example
molecule = rdkit.Chem.MolFromSmiles('C1c2ccc(cc2)CCc2ccc(cc2)C1')

# calculation step
data = py_rdl.Calculator.get_calculated_result(
           molecule.GetBonds(),              # pass edges
           rdkit.Chem.Bond.GetBeginAtom,     # get first atom
           rdkit.Chem.Bond.GetEndAtom,       # get second atom
           rdkit.Chem.Atom.GetIdx,           # identify atom
           rdkit.Chem.Bond.GetIdx            # identify edge
       )

# iterate over URFs
for urf in data.urfs:
    # print URF
    print(urf)
    # URF can be used as an index, get RCs
    rcs = data.get_relevant_cycles_for_urf(urf)
    # print RCs
    for rc in rcs:
        print(rc)