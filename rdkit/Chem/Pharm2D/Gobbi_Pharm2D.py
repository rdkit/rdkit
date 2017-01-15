#
#  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Definitions for 2D Pharmacophores from:
  Gobbi and Poppinger, Biotech. Bioeng. _61_ 47-54 (1998)

"""
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D.SigFactory import SigFactory

fdef = """
DefineFeature Hydrophobic [$([C;H2,H1](!=*)[C;H2,H1][C;H2,H1][$([C;H1,H2,H3]);!$(C=*)]),$(C([C;H2,H3])([C;H2,H3])[C;H2,H3])]
  Family LH
  Weights 1.0
EndFeature
DefineFeature Donor [$([N;!H0;v3]),$([N;!H0;+1;v4]),$([O,S;H1;+0]),$([n;H1;+0])]
  Family HD
  Weights 1.0
EndFeature
DefineFeature Acceptor [$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N&v3;H1,H2]-[!$(*=[O,N,P,S])]),$([N;v3;H0]),$([n,o,s;+0]),F]
  Family HA
  Weights 1.0
EndFeature
DefineFeature AromaticAttachment [$([a;D3](@*)(@*)*)]
  Family AR
  Weights 1.0
EndFeature
DefineFeature AliphaticAttachment [$([A;D3](@*)(@*)*)]
  Family RR
  Weights 1.0
EndFeature
DefineFeature UnusualAtom [!#1;!#6;!#7;!#8;!#9;!#16;!#17;!#35;!#53]
  Family X
  Weights 1.0
EndFeature
DefineFeature BasicGroup [$([N;H2&+0][$([C,a]);!$([C,a](=O))]),$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))]),$([N,n;X2;+0])]
  Family BG
  Weights 1.0
EndFeature
DefineFeature AcidicGroup [$([C,S](=[O,S,P])-[O;H1])]
  Family AG
  Weights 1.0
EndFeature
"""
defaultBins = [(2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 100)]


def _init():
  global labels, patts, factory
  featFactory = ChemicalFeatures.BuildFeatureFactoryFromString(fdef)
  factory = SigFactory(featFactory, minPointCount=2, maxPointCount=3)
  factory.SetBins(defaultBins)
  factory.Init()


_init()
