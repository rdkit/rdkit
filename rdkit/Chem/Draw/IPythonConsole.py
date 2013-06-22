import IPython
if IPython.release.version<'0.11':
    raise ImportError,'this module requires at least v0.11 of IPython'

from rdkit.Chem import rdchem
from rdkit.Chem import Draw
from cStringIO import StringIO
import copy
import numpy
try:
    import Image
except ImportError:
    from PIL import Image

molSize=(450,150)
highlightSubstructs=True
kekulizeStructures=True
def _toPNG(mol):
    if hasattr(mol,'__sssAtoms'):
        highlightAtoms=mol.__sssAtoms
    else:
        highlightAtoms=[]
    try:
        mol.GetAtomWithIdx(0).GetExplicitValence()
    except RuntimeError:
        mol.UpdatePropertyCache(False)
    
    mc = copy.deepcopy(mol)
    try:
        img = Draw.MolToImage(mc,size=molSize,kekulize=kekulizeStructures,
                            highlightAtoms=highlightAtoms)
    except ValueError:  # <- can happen on a kekulization failure
        mc = copy.deepcopy(mol)
        img = Draw.MolToImage(mc,size=molSize,kekulize=False,
                            highlightAtoms=highlightAtoms)
        
    sio = StringIO()
    img.save(sio,format='PNG')
    return sio.getvalue()

def _GetSubstructMatch(mol,query,**kwargs):
    res = mol.__GetSubstructMatch(query,**kwargs)
    if highlightSubstructs:
        mol.__sssAtoms=list(res)
    else:
        mol.__sssAtoms=[]
    return res

def _GetSubstructMatches(mol,query,**kwargs):
    res = mol.__GetSubstructMatches(query,**kwargs)
    mol.__sssAtoms=[]
    if highlightSubstructs:
        for entry in res:
            mol.__sssAtoms.extend(list(entry))
    return res

# code for displaying PIL images directly,
def display_pil_image(img):
    """displayhook function for PIL Images, rendered as PNG"""
    sio = StringIO()
    img.save(sio,format='PNG')
    return sio.getvalue()

def InstallIPythonRenderer():
    rdchem.Mol._repr_png_=_toPNG
    if not hasattr(rdchem.Mol,'__GetSubstructMatch'):
        rdchem.Mol.__GetSubstructMatch=rdchem.Mol.GetSubstructMatch
    rdchem.Mol.GetSubstructMatch=_GetSubstructMatch
    if not hasattr(rdchem.Mol,'__GetSubstructMatches'):
        rdchem.Mol.__GetSubstructMatches=rdchem.Mol.GetSubstructMatches
    rdchem.Mol.GetSubstructMatches=_GetSubstructMatches
    Image.Image._repr_png_=display_pil_image
InstallIPythonRenderer()

def UninstallIPythonRenderer():
    del rdchem.Mol._repr_png_
    if hasattr(rdchem.Mol,'__GetSubstructMatch'):
        rdchem.Mol.GetSubstructMatch=rdchem.Mol.__GetSubstructMatch
        del rdchem.Mol.__GetSubstructMatch
    if hasattr(rdchem.Mol,'__GetSubstructMatches'):
        rdchem.Mol.GetSubstructMatches=rdchem.Mol.__GetSubstructMatches
        del rdchem.Mol.__GetSubstructMatches
    del Image.Image._repr_png_

