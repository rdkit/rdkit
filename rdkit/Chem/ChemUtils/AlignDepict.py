# $Id: AlignDepict.py 736 2008-02-14 14:09:36Z landrgr1 $
#
#  Copyright (C) 2006 Greg Landrum
#  This file is part of RDKit and covered by $RDBASE/license.txt
#
from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit import Geometry

def AlignDepict(mol,core,corePattern=None,acceptFailure=False):
  """

  Arguments:
    - mol:          the molecule to be aligned, this will come back
                    with a single conformer.
    - core:         a molecule with the core atoms to align to;
                    this should have a depiction.
    - corePattern:  (optional) an optional molecule to be used to 
                    generate the atom mapping between the molecule
                    and the core.
  """
  if core and corePattern:
    if not core.GetNumAtoms(onlyExplicit=True)==corePattern.GetNumAtoms(onlyExplicit=True):
      raise ValueError('When a pattern is provided, it must have the same number of atoms as the core')
    coreMatch = core.GetSubstructMatch(corePattern)
    if not coreMatch:
      raise ValueError("Core does not map to itself")
  else:
    coreMatch = range(core.GetNumAtoms(onlyExplicit=True))
  if corePattern:
    match = mol.GetSubstructMatch(corePattern)
  else:
    match = mol.GetSubstructMatch(core)

  if not match:
    if not acceptFailure:
      raise ValueError('Substructure match with core not found.')
    else:
      coordMap={}
  else:
    conf = core.GetConformer()
    coordMap={}
    for i,idx in enumerate(match):
      pt3 = conf.GetAtomPosition(coreMatch[i])
      pt2 = Geometry.Point2D(pt3.x,pt3.y)
      coordMap[idx] = pt2
  rdDepictor.Compute2DCoords(mol,clearConfs=True,coordMap=coordMap,canonOrient=False)

if __name__=='__main__':
  import sys,getopt

  def Usage():
    pass
  
  args,extras = getopt.getopt(sys.argv[1:],'p:ho:',['smiles','pattern='])
  if len(extras)!=2:
    print('ERROR: Not enough arguments', file=sys.stderr)
    Usage()
    sys.exit(1)
  patt = None
  useSmiles = False
  outF=None
  for arg,val in args:
    if arg=='-h':
      Usage()
      sys.exit(0)
    elif arg=='-p' or arg=='--pattern':
      patt = Chem.MolFromSmarts(val)
    elif arg=='--smiles':
      useSmiles = True
    elif arg=='-o':
      outF = val
  
  if not useSmiles:
    core = Chem.MolFromMolFile(extras[0])
  else:
    core = Chem.MolFromSmiles(extras[0])
    rdDepictor.Compute2DCoords(core)

  if not useSmiles:
    mol = Chem.MolFromMolFile(extras[1])
  else:
    mol = Chem.MolFromSmiles(extras[1])

  AlignDepict(mol,core,patt)
  
  if outF:
    outF = open(outF,'w+')
  else:
    outF = sys.stdout

  print(Chem.MolToMolBlock(mol), file=outF)
  
  
  
