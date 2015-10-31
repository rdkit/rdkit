import IPython

if IPython.release.version < '0.11':
    raise ImportError('this module requires at least v0.11 of IPython')
elif IPython.release.version < '2.0':
    install_nbextension=None
    _canUse3D=False
else:
    from IPython.html.nbextensions import install_nbextension
    _canUse3D=True
from rdkit import Chem
from rdkit.Chem import rdchem, rdChemReactions
from rdkit.Chem import Draw
from rdkit.six import BytesIO,StringIO
import copy
import os
import json
import uuid
import numpy
try:
    import Image
except ImportError:
    from PIL import Image

molSize = (450, 150)
highlightSubstructs = True
kekulizeStructures = True

ipython_useSVG = False
ipython_3d = False
molSize_3d = (450, 450)
drawing_type_3d = "ball and stick"
camera_type_3d = "perspective"
shader_3d = "lambert"


def _toJSON(mol):
    """For IPython notebook, renders 3D webGL objects."""

    if not ipython_3d or not mol.GetNumConformers():
        return None

    try:
        import imolecule
    except ImportError:
        raise ImportError("Cannot import 3D rendering. Please install "
                          "with `pip install imolecule`.")

    conf = mol.GetConformer()
    if not conf.Is3D():
        return None

    mol = Chem.Mol(mol)
    try:
        Chem.Kekulize(mol)
    except Exception:
        mol = Chem.Mol(mol)
    size = molSize_3d

    # center the molecule:
    atomps = numpy.array([list(conf.GetAtomPosition(x))
                          for x in range(mol.GetNumAtoms())])
    avgP = numpy.average(atomps, 0)
    atomps -= avgP

    # Convert the relevant parts of the molecule into JSON for rendering
    atoms = [{"element": atom.GetSymbol(),
              "location": list(atomps[atom.GetIdx()])}
             for atom in mol.GetAtoms()]
    bonds = [{"atoms": [bond.GetBeginAtomIdx(),
                        bond.GetEndAtomIdx()],
              "order": int(bond.GetBondTypeAsDouble())}
             for bond in mol.GetBonds()]
    mol = {"atoms": atoms, "bonds": bonds}
    return imolecule.draw({"atoms": atoms, "bonds": bonds}, format="json",
                          size=molSize_3d, drawing_type=drawing_type_3d,
                          camera_type=camera_type_3d, shader=shader_3d,
                          display_html=False)


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
    bio = BytesIO()
    img.save(bio,format='PNG')
    return bio.getvalue()

def _toSVG(mol):
    if not ipython_useSVG:
        return None
    if hasattr(mol, '__sssAtoms'):
        highlightAtoms = mol.__sssAtoms
    else:
        highlightAtoms = []
    try:
        mol.GetAtomWithIdx(0).GetExplicitValence()
    except RuntimeError:
        mol.UpdatePropertyCache(False)

    mc = copy.deepcopy(mol)
    sio = StringIO()
    try:
        Draw.MolToFile(mc, sio, size=molSize, imageType="svg",
                       kekulize=kekulizeStructures,
                       highlightAtoms=highlightAtoms)
    except ValueError:  # <- can happen on a kekulization failure
        mc = copy.deepcopy(mol)
        Draw.MolToFile(mc, sio, size=molSize, kekulize=False,
                       highlightAtoms=highlightAtoms, imageType="svg")

    # It seems that the svg renderer used doesn't quite hit the spec.
    # Here are some fixes to make it work in the notebook, although I think
    # the underlying issue needs to be resolved at the generation step
    svg_split = sio.getvalue().replace("svg:", "").splitlines()
    header = ('<svg version="1.1" xmlns="http://www.w3.org/2000/svg"'
              'width="%dpx" height="%dpx">') % molSize
    svg = "\n".join(svg_split[:1] + [header] + svg_split[5:])
    return svg


def _toReactionPNG(rxn):
    rc = copy.deepcopy(rxn)
    img = Draw.ReactionToImage(rc,subImgSize=(int(molSize[0]/3), molSize[1]))
    bio = BytesIO()
    img.save(bio,format='PNG')
    return bio.getvalue()

def _GetSubstructMatch(mol, query, **kwargs):
    res = mol.__GetSubstructMatch(query, **kwargs)
    if highlightSubstructs:
        mol.__sssAtoms = list(res)
    else:
        mol.__sssAtoms = []
    return res


def _GetSubstructMatches(mol, query, **kwargs):
    res = mol.__GetSubstructMatches(query, **kwargs)
    mol.__sssAtoms = []
    if highlightSubstructs:
        for entry in res:
            mol.__sssAtoms.extend(list(entry))
    return res


# code for displaying PIL images directly,
def display_pil_image(img):
    """displayhook function for PIL Images, rendered as PNG"""
    bio = BytesIO()
    img.save(bio,format='PNG')
    return bio.getvalue()

def InstallIPythonRenderer():
    rdchem.Mol._repr_png_ = _toPNG
    rdchem.Mol._repr_svg_ = _toSVG
    if _canUse3D:
        rdchem.Mol._repr_html_ = _toJSON
    rdChemReactions.ChemicalReaction._repr_png_ = _toReactionPNG
    if not hasattr(rdchem.Mol, '__GetSubstructMatch'):
        rdchem.Mol.__GetSubstructMatch = rdchem.Mol.GetSubstructMatch
    rdchem.Mol.GetSubstructMatch = _GetSubstructMatch
    if not hasattr(rdchem.Mol, '__GetSubstructMatches'):
        rdchem.Mol.__GetSubstructMatches = rdchem.Mol.GetSubstructMatches
    rdchem.Mol.GetSubstructMatches = _GetSubstructMatches
    Image.Image._repr_png_ = display_pil_image
InstallIPythonRenderer()


def UninstallIPythonRenderer():
    del rdchem.Mol._repr_svg_
    del rdchem.Mol._repr_png_
    if _canUse3D:
        del rdchem.Mol._repr_html_
    del rdChemReactions.ChemicalReaction._repr_png_
    if hasattr(rdchem.Mol, '__GetSubstructMatch'):
        rdchem.Mol.GetSubstructMatch = rdchem.Mol.__GetSubstructMatch
        del rdchem.Mol.__GetSubstructMatch
    if hasattr(rdchem.Mol, '__GetSubstructMatches'):
        rdchem.Mol.GetSubstructMatches = rdchem.Mol.__GetSubstructMatches
        del rdchem.Mol.__GetSubstructMatches
    del Image.Image._repr_png_
