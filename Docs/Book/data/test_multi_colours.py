#!/usr/bin/env python

from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from json import dumps

COLS = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0), (1.0, 0.55, 1.0)]


def get_hit_atoms_and_bonds(mol, smt):
    alist = []
    blist = []
    q = Chem.MolFromSmarts(smt)
    for match in mol.GetSubstructMatches(q):
        alist.extend(match)

    for ha1 in alist:
        for ha2 in alist:
            if ha1 > ha2:
                b = mol.GetBondBetweenAtoms(ha1, ha2)
                if b:
                    blist.append(b.GetIdx())
    
    return alist, blist


def add_colours_to_map(els, cols, col_num):
    for el in els:
        if not el in cols:
            cols[el] = []
        if COLS[col_num] not in cols[el]:
            cols[el].append(COLS[col_num])


def do_a_picture(smi, smarts, filename, label, fmt='svg'):

    rdDepictor.SetPreferCoordGen(True)
    mol = Chem.MolFromSmiles(smi)
    mol = Draw.PrepareMolForDrawing(mol)

    acols = {}
    bcols = {}
    h_rads = {}
    h_lw_mult = {}

    for i, smt in enumerate(smarts):
        alist, blist = get_hit_atoms_and_bonds(mol, smt)
        col = i % 4
        add_colours_to_map(alist, acols, col)
        add_colours_to_map(blist, bcols, col)
    
    if fmt == 'svg':
        d = rdMolDraw2D.MolDraw2DSVG(300, 300)
        mode = 'w'
    elif fmt == 'png':
        d = rdMolDraw2D.MolDraw2DCairo(300, 300)
        mode = 'wb'
    else:
        print('unknown format {}'.format(fmt))
        return
    
    d.drawOptions().fillHighlights = False
    d.DrawMoleculeWithHighlights(mol, label, acols, bcols, h_rads, h_lw_mult, -1)
    d.FinishDrawing()
        
    with open(filename, mode) as f:
        f.write(d.GetDrawingText())


smi = 'CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]'
smarts = ['CONN', 'N#CC~CO', 'C=CON', 'CONNCN']
do_a_picture(smi, smarts, 'atom_highlights_3.png', '', fmt='png')

