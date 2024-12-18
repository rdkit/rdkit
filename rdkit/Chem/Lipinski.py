#
# Copyright (C) 2001-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Calculation of Lipinski parameters for molecules

"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# -----------------------------------
# on import build the SMARTS patterns so we only have to do it once
# -----------------------------------
# The Daylight SMARTS expressions for
# recognizing H-bond donors and acceptors in the Lipinski scheme.
# HDonor     '[!#6;!H0;-0]'
# HAcceptor  '[$([!#6;+0]);!$([F,Cl,Br,I]);
#             !$([o,s,nX3]);!$([Nv5,Pv5,Sv4,Sv6])]'
# Heteroatom '[B,N,O,P,S,F,Cl,Br,I]'
# 2 definitions adapted from those in the Gobbi Paper
#  NOTE: if you want traditional Lipinski numbers, you
#  should use NOCounts (below) instead of HAcceptor
#
HDonorSmarts = Chem.MolFromSmarts('[$([N;!H0;v3]),$([N;!H0;+1;v4]),$([O,S;H1;+0]),$([n;H1;+0])]')
# changes log for HAcceptorSmarts:
#  v2, 1-Nov-2008, GL : fix amide-N exclusion; remove Fs from definition
HAcceptorSmarts = Chem.MolFromSmarts('[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),' +
                                     '$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=!@[O,N,P,S])]),' +
                                     '$([nH0,o,s;+0])]')

HeteroatomSmarts = Chem.MolFromSmarts('[!#6;!#1]')
#  NOTE: the Rotatable bond smarts here doesn't treat deuteriums (which are left in the graph
#  and therefore contribute to the degree of a carbon) the same as hydrogens (which are removed
#  from the graph). So the bond in [2H]C([2H])([2H])C([2H])([2H])[2H] *is* considered
#  rotatable.
RotatableBondSmarts = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')
NHOHSmarts = Chem.MolFromSmarts('[#8H1,#7H1,#7H2,#7H3]')
NOCountSmarts = Chem.MolFromSmarts('[#7,#8]')

# # this little trick saves duplicated code
# def _NumMatches(mol, smarts):
#   return len(mol.GetSubstructMatches(smarts, uniquify=1))

NumHDonors = lambda x: rdMolDescriptors.CalcNumHBD(x)
NumHDonors.__doc__ = "Number of Hydrogen Bond Donors"
NumHDonors.version = "1.0.0"
_HDonors = lambda x, y=HDonorSmarts: x.GetSubstructMatches(y, uniquify=1)
NumHAcceptors = lambda x: rdMolDescriptors.CalcNumHBA(x)
NumHAcceptors.__doc__ = "Number of Hydrogen Bond Acceptors"
NumHAcceptors.version = "2.0.0"
_HAcceptors = lambda x, y=HAcceptorSmarts: x.GetSubstructMatches(y, uniquify=1)
NumHeteroatoms = lambda x: rdMolDescriptors.CalcNumHeteroatoms(x)
NumHeteroatoms.__doc__ = "Number of Heteroatoms"
NumHeteroatoms.version = "1.0.0"
_Heteroatoms = lambda x, y=HeteroatomSmarts: x.GetSubstructMatches(y, uniquify=1)
NumRotatableBonds = lambda x: rdMolDescriptors.CalcNumRotatableBonds(x)
NumRotatableBonds.__doc__ = "Number of Rotatable Bonds"
NumRotatableBonds.version = "1.0.0"
_RotatableBonds = lambda x, y=RotatableBondSmarts: x.GetSubstructMatches(y, uniquify=1)
NOCount = lambda x: rdMolDescriptors.CalcNumLipinskiHBA(x)
NOCount.__doc__ = "Number of Nitrogens and Oxygens"
NOCount.version = "1.0.0"
NHOHCount = lambda x: rdMolDescriptors.CalcNumLipinskiHBD(x)
NHOHCount.__doc__ = "Number of NHs or OHs"
NHOHCount.version = "2.0.0"

RingCount = lambda x: rdMolDescriptors.CalcNumRings(x)
RingCount.version = "1.0.0"


def HeavyAtomCount(mol):
  " Number of heavy atoms a molecule."
  return mol.GetNumHeavyAtoms()


HeavyAtomCount.version = "1.0.1"

_bulkConvert = ("CalcFractionCSP3", "CalcNumAromaticRings", "CalcNumSaturatedRings",
                "CalcNumAromaticHeterocycles", "CalcNumAromaticCarbocycles",
                "CalcNumSaturatedHeterocycles", "CalcNumSaturatedCarbocycles",
                "CalcNumAliphaticRings", "CalcNumAliphaticHeterocycles",
                "CalcNumAliphaticCarbocycles", "CalcNumHeterocycles", "CalcNumBridgeheadAtoms",
                "CalcNumAmideBonds", "CalcNumAtomStereoCenters", "CalcNumHeterocycles",
                "CalcNumUnspecifiedAtomStereoCenters", "CalcNumSpiroAtoms", "CalcPhi")
for txt in _bulkConvert:
  _cfn = getattr(rdMolDescriptors, txt)
  _fn = lambda x, y=_cfn: y(x)
  try:
    _fn.version = getattr(rdMolDescriptors, "_" + txt + "_version")
  except AttributeError:
    pass
  _fn.__doc__ = _cfn.__doc__
  nm = txt.replace("Calc", "")
  locals()[nm] = _fn
