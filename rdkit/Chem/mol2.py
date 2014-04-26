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
    for confId in confIds:
        
        molecule = """@<TRIPOS>MOLECULE
{}
{} {} 0 0 0
SMALL
USER_CHARGES\n""".format(mol.GetProp("_Name"), mol.GetNumAtoms(), mol.GetNumBonds())

        # FIXME "USER_CHARGES" could become 'Gasteiger charges'
        # FIXME "SMALL" means small molecule but could become "PROTEIN"
        
        pos = _get_positions(mol, confId)
        atom_lines = atom_lines = ["{:>4} {:>4} {:>13.4f} {:>9.4f} {:>9.4f} {:>4} {} {} {:>7.4f}".format(a.GetIdx()+1, a.GetProp("_TriposAtomName"), float(pos[a.GetIdx()][0]), float(pos[a.GetIdx()][1]), float(pos[a.GetIdx()][2]), a.GetProp("_TriposAtomType"), 1, "UNL", float(a.GetProp("_TriposPartialCharge"))) for a in mol.GetAtoms()]
        atom_lines = ["@<TRIPOS>ATOM"] +  atom_lines + ["\n"]
        atom_lines = "\n".join(atom_lines)
            
        bond_lines = [ "{:>5} {:>5} {:>5} {:>2}".format(bid+1, b.GetBeginAtomIdx()+1, b.GetEndAtomIdx()+1,  "ar" if int(b.GetBondType()) > 2 else str(int(b.GetBondType())) ) for bid, (b) in enumerate(mol.GetBonds()) ]
        bond_lines = ["@<TRIPOS>BOND"] + bond_lines + [ "\n"]
        bond_lines = "\n".join(bond_lines)
        
        block = molecule + atom_lines + bond_lines
        blocks.append(block)   
    return "".join(blocks)
    
__all__ = ["MolToMol2Block", "MolToMol2File"]    
