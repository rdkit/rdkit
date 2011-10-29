import IPython
if IPython.release.version<'0.11':
    raise ImportError,'this module requires at least v0.11 of IPython'

from rdkit.Chem import rdchem
from rdkit.Chem import Draw
from cStringIO import StringIO

molSize=(450,150)
highlightSubstructs=True
def _toPNG(mol):
    if hasattr(mol,'__sssAtoms'):
        highlightAtoms=mol.__sssAtoms
    else:
        highlightAtoms=[]
    img = Draw.MolToImage(mol,size=molSize,
                          highlightAtoms=highlightAtoms)
    sio = StringIO()
    img.save(sio,format='PNG')
    return sio.getvalue()

def _GetSubstructMatch(mol,query):
    res = mol.__GetSubstructMatch(query)
    if highlightSubstructs:
        mol.__sssAtoms=list(res)
    else:
        mol.__sssAtoms=[]
    return res

def _GetSubstructMatches(mol,query):
    res = mol.__GetSubstructMatches(query)
    mol.__sssAtoms=[]
    if highlightSubstructs:
        for entry in res:
            mol.__sssAtoms.extend(list(entry))
    return res

def InstallIPythonRenderer():
    rdchem.Mol._repr_png_=_toPNG

    if not hasattr(rdchem.Mol,'__GetSubstructMatch'):
        rdchem.Mol.__GetSubstructMatch=rdchem.Mol.GetSubstructMatch
    rdchem.Mol.GetSubstructMatch=_GetSubstructMatch
    if not hasattr(rdchem.Mol,'__GetSubstructMatches'):
        rdchem.Mol.__GetSubstructMatches=rdchem.Mol.GetSubstructMatches
    rdchem.Mol.GetSubstructMatches=_GetSubstructMatches
    
InstallIPythonRenderer()

def UninstallIPythonRenderer():
    del rdchem.Mol._repr_png_
    if hasattr(rdchem.Mol,'__GetSubstructMatch'):
        rdchem.Mol.GetSubstructMatch=rdchem.Mol.__GetSubstructMatch
        del rdchem.Mol.__GetSubstructMatch
    if hasattr(rdchem.Mol,'__GetSubstructMatches'):
        rdchem.Mol.GetSubstructMatches=rdchem.Mol.__GetSubstructMatches
        del rdchem.Mol.__GetSubstructMatches
    

