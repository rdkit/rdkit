# $Id$
#
#  Copyright (C) 2003-2006 Greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" Definitions for 2D Pharmacophores from:
  Gobbi and Poppinger, Biotech. Bioeng. _61_ 47-54 (1998)

"""
from pyRDKit import Chem
from pyRDKit import RDConfig
import os,sys
from pyRDKit.Chem.Pharm2D.SigFactory import SigFactory

patts={
  'LH':('[$([C;H2,H1](!=*)[C;H2,H1][C;H2,H1][$([C;H1,H2,H3]);!$(C=*)]),$(C([C;H2,H3])([C;H2,H3])[C;H2,H3])]',
      'Hydrophobic group'),
  'HD':('[$([N;!H0;v3]),$([N;!H0;+1;v4]),$([O,S;H1;+0]),$([n;H1;+0])]',
        'Hydrogen bond donor'),
  'HA':('[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N&v3;H1,H2]-[!$(*=[O,N,P,S])]),$([N;v3;H0]),$([n,o,s;+0]),F]',
        'Hydrogen bond acceptor'),
  'AR':('[$([a;D3](@*)(@*)*)]',
        'Aromatic Attachment Point'),
  'RR':('[$([A;D3](@*)(@*)*)]',
        'Aliphatic Attachment Point'),
  'X':('[!#1;!#6;!#7;!#8;!#9;!#16;!#17;!#35;!#53]',
       'Unusual Atom'),
  'BG':('[$([N;H2&+0][$([C,a]);!$([C,a](=O))]),$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))]),$([N,n;X2;+0])]',
        'Basic Group'),
  'AG':('[C,S](=[O,S,P])-[O;H1]',
        'Acidic Group'),
  }
defaultBins = [(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,100)]

def _init():
  global labels,patts,factory

  labels = []
  smaPatts = []
  for label in patts.keys():
    sma,desc = patts[label]
    p = Chem.MolFromSmarts(sma)
    if not p:
      sys.stderr.write('feature %s failed to parse'%(label))
    else:
      labels.append(label)
      patts[label] = p,desc
      smaPatts.append(p)
  factory = SigFactory()    
  factory.SetPatterns(smaPatts)
  factory.SetLabels(labels)
  factory.SetBins(defaultBins)
  factory.SetMinCount(2)
  factory.SetMaxCount(3)
  factory.SetShortestPathsOnly(1)
  factory.SetIncludeBondOrder(0)
  
_init()
