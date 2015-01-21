# $Id$
#
# Copyright (C) 2001-2008 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from rdmolfiles import MolFromSmarts
from rdPartialCharges import ComputeGasteigerCharges

def _get_positions(mol,confId=-1):
    conf = mol.GetConformer(confId)
    return [conf.GetAtomPosition(i) for i in range(conf.GetNumAtoms())]


def MolToMol2File(mol, filename, confId=-1):
    """Writes a Mol2 file for a molecule
      ARGUMENTS:

        - mol: the molecule
        - filename: the file to write to
        - confId: (optional) selects which conformation to output (-1 = default)
                  if set to None will return all conformers

      RETURNS:

        None
    """      
    block = MolToMol2Block(mol, confId)
    open(filename, "w").writelines(block)
    
def MolToMol2Block(mol, confId=-1):
    """Returns a Mol2 string block for a molecule
      ARGUMENTS:

        - mol: the molecule
        - confId: (optional) selects which conformation to output (-1 = default)
                  if set to None will return all conformers

      RETURNS:

        a string
    """  

    #
    # References 
    # - Format specs http://www.tripos.com/data/support/mol2.pdf
    #
    
    confIds = (confId,)
    
    if confId == None:
        confIds = Chem.Mol.GetNumConformers()
    
    blocks = []
    
    # compute charges
    ComputeGasteigerCharges(mol)
    
    for confId in confIds:
        
        molecule = """@<TRIPOS>MOLECULE
{}
{} {} 0 0 0
SMALL
GASTEIGER\n\n""".format(mol.GetProp("_Name") if mol.HasProp("_Name") else "UNK", mol.GetNumAtoms(), mol.GetNumBonds())

        # FIXME "USER_CHARGES" could become 'Gasteiger charges'
        # FIXME "SMALL" means small molecule but could become "PROTEIN"
        
        pos = _get_positions(mol, confId)
        atom_lines = ["{:>4} {:>4} {:>13.4f} {:>9.4f} {:>9.4f} {:<5} {} {} {:>7.4f}".format(a.GetIdx()+1, 
                                                                                                         a.GetSymbol(), 
                                                                                                         float(pos[a.GetIdx()][0]), float(pos[a.GetIdx()][1]), float(pos[a.GetIdx()][2]), 
                                                                                                         _sybyl_atom_type(a), 
                                                                                                         1, "UNL", 
                                                                                                         float(a.GetProp('_GasteigerCharge').replace(',','.')) if a.HasProp('_GasteigerCharge') else 0.0) for a in mol.GetAtoms()]
        atom_lines = ["@<TRIPOS>ATOM"] +  atom_lines + ["\n"]
        atom_lines = "\n".join(atom_lines)
            
        bond_lines = [ "{:>5} {:>5} {:>5} {:>2}".format(bid+1, b.GetBeginAtomIdx()+1, b.GetEndAtomIdx()+1,  "ar" if int(b.GetBondType()) > 2 else str(int(b.GetBondType())) ) for bid, (b) in enumerate(mol.GetBonds()) ]
        bond_lines = ["@<TRIPOS>BOND"] + bond_lines + [ "\n"]
        bond_lines = "\n".join(bond_lines)
        
        block = molecule + atom_lines + bond_lines
        blocks.append(block)   
    return "".join(blocks)

def _sybyl_atom_type(atom):
    """ Asign sybyl atom type
    Reference: http://www.tripos.com/mol2/atom_types.html
    """
    sybyl = None
    atom_symbol = atom.GetSymbol()
    atomic_num = atom.GetAtomicNum()
    hyb = atom.GetHybridization()-1 # -1 since 1 = sp, 2 = sp1 etc
    hyb = 3 if hyb > 3 else hyb # cap at 3
    if atomic_num == 6:
        if atom.GetIsAromatic():
            sybyl = 'C.ar'
        elif _atom_matches_smarts(atom, 'C(=N)(N)N*'): # https://github.com/rdkit/rdkit/blob/master/Data/FragmentDescriptors.csv
            sybyl = 'C.cat'
        else:
            sybyl = '%s.%i' % (atom_symbol, hyb if hyb < 3 else 3)
    elif atomic_num == 7:
        if atom.GetIsAromatic():
            sybyl = 'N.ar'
        elif _atom_matches_smarts(atom, 'C(=O)-N'): # https://github.com/rdkit/rdkit/blob/master/Data/FragmentDescriptors.csv
            sybyl = 'N.am'
        elif _atom_matches_smarts(atom, '[$([nX3](:*):*),$([nX2](:*):*),$([#7X2]=*),$([NX3](=*)=*),$([#7X3+](-*)=*),$([#7X3+H]=*)]'): # http://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html
            sybyl = 'N.pl3'
        elif hyb == 3 and atom.GetFormalCharge():
            sybyl = 'N.4'
        else:
            sybyl = '%s.%i' % (atom_symbol, hyb)
    elif atomic_num == 8:
        if _atom_matches_smarts(atom, 'C(=O)[O;H,-]'):
            sybyl = 'O.co2'
        else:
            sybyl = '%s.%i' % (atom_symbol, hyb)
    elif atomic_num == 16:
        if _atom_matches_smarts(atom, '[$([#16X3]=[OX1]),$([#16X3+][OX1-])]'): # http://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html
            sybyl = 'S.O'
        elif _atom_matches_smarts(atom, 'S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])(-[#6])-[#6]'): # https://github.com/rdkit/rdkit/blob/master/Data/FragmentDescriptors.csv
            sybyl = 'S.O2'
        else:
            sybyl = '%s.%i' % (atom_symbol, hyb)
    elif atomic_num == 15 and hyb == 3:
        sybyl = '%s.%i' % (atom_symbol, hyb)
    
    if not sybyl:
        sybyl = atom_symbol
    return sybyl

def _atom_matches_smarts(atom, smarts):
    idx = atom.GetIdx()
    patt = MolFromSmarts(smarts)
    for m in atom.GetOwningMol().GetSubstructMatches(patt):
        if idx in m:
            return True
    return False

__all__ = ["MolToMol2Block", "MolToMol2File"]    
