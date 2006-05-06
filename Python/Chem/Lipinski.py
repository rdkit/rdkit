# $Id: Lipinski.py 5007 2006-02-22 15:14:41Z glandrum $
#
# Copyright (C) 2001-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" Calculation of Lipinski parameters for molecules

"""
#from Chem import rdchem
import Chem

#-----------------------------------
# on import build the SMARTS patterns so we only have to do it once
#-----------------------------------

# The Daylight SMARTS expressions for
# recognizing H-bond donors and acceptors in the Lipinski scheme.
# HDonor     '[!#6;!H0;-0]'
# HAcceptor  '[$([!#6;+0]);!$([F,Cl,Br,I]);
#             !$([o,s,nX3]);!$([Nv5,Pv5,Sv4,Sv6])]'
# Heteroatom '[B,N,O,P,S,F,Cl,Br,I]'


# 2 definitions from Gobbi Paper
#  NOTE: if you want traditional Lipinski numbers, you
#  should use NOCounts (below) instead of HAcceptor
#
HDonorSmarts = Chem.MolFromSmarts('[$([N;!H0;v3]),$([N;!H0;+1;v4]),$([O,S;H1;+0]),$([n;H1;+0])]')
HAcceptorSmarts = Chem.MolFromSmarts('[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N&v3;H1,H2]-[!$(*=[O,N,P,S])]),$([N;v3;H0]),$([n,o,s;+0]),F]')
HeteroatomSmarts = Chem.MolFromSmarts('[!#6;!#1]')
RotatableBondSmarts = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')
NHOHSmarts = Chem.MolFromSmarts('[#8H1,#7H1,#7H2,#7H3]')
NOCountSmarts = Chem.MolFromSmarts('[#7,#8]')

# this little trick saves duplicated code
def _NumMatches(mol,smarts):
  return len(mol.GetSubstructMatches(smarts,uniquify=1))

NumHDonors = lambda x,y=HDonorSmarts:_NumMatches(x,y)
NumHDonors.__doc__="Number of Hydrogen Bond Donors"
NumHDonors.version="1.0.0"
_HDonors = lambda x,y=HDonorSmarts:x.GetSubstructMatches(y,uniquify=1)
NumHAcceptors = lambda x,y=HAcceptorSmarts:_NumMatches(x,y)
NumHAcceptors.__doc__="Number of Hydrogen Bond Acceptors"
NumHAcceptors.version="1.0.0"
_HAcceptors = lambda x,y=HAcceptorSmarts:x.GetSubstructMatches(y,uniquify=1)
NumHeteroatoms = lambda x,y=HeteroatomSmarts:_NumMatches(x,y)
NumHeteroatoms.__doc__="Number of Heteroatoms"
NumHeteroatoms.version="1.0.0"
_Heteroatoms = lambda x,y=HeteroatomSmarts:x.GetSubstructMatches(y,uniquify=1)
NumRotatableBonds = lambda x,y=RotatableBondSmarts:_NumMatches(x,y)
NumRotatableBonds.__doc__="Number of Rotatable Bonds"
NumRotatableBonds.version="1.0.0"
_RotatableBonds = lambda x,y=RotatableBondSmarts:x.GetSubstructMatches(y,uniquify=1)
NOCount = lambda x,y=NOCountSmarts:_NumMatches(x,y)
NOCount.__doc__="Number of Nitrogens and Oxygens"
NOCount.version="1.0.0"
NHOHCount = lambda x,y=NHOHSmarts:_NumMatches(x,y)
NHOHCount.__doc__="Number of NHs or OHs"
NHOHCount.version="1.0.0"

  
