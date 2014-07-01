import IPython
if IPython.release.version < '2.0':
    raise ImportError('this module requires at least v2.0 of IPython')
from rdkit import Chem
from rdkit.Chem import rdchem, rdChemReactions
from rdkit.Chem import Draw
from cStringIO import StringIO
import copy
import json
import numpy
try:
    import Image
except ImportError:
    from PIL import Image

molSize = (450, 150)
highlightSubstructs = True
kekulizeStructures = True

_3d_initialized = False

ipython_useSVG = False
ipython_3d = False
molSize_3d = (450, 450)
drawing_type_3d = "ball and stick"
camera_type_3d = "perspective"


def _toJSON(mol):
    """For IPython notebook, renders 3D webGL objects."""

    if not ipython_3d or not mol.GetNumConformers():
        return None

    conf = mol.GetConformer()
    if not conf.Is3D():
        return None

    mol = Chem.Mol(mol)
    try:
        Chem.Kekulize(mol)
    except:
        mol = Chem.Mol(mol)
    size = molSize_3d

    # center the molecule:
    atomps = numpy.array([list(conf.GetAtomPosition(x))
                          for x in range(mol.GetNumAtoms())])
    avgP = numpy.average(atomps, 0)
    atomps -= avgP

    lib_path = "/nbextensions/imolecule.min.js"
    # If the javascript lib has not yet been loaded, do so.
    # IPython >=2.0 does this by copying the file into ~/.ipython/nbextensions
    global _3d_initialized
    if not _3d_initialized:
        import os
        from IPython.html.nbextensions import install_nbextension
        filepath = os.path.normpath(os.path.dirname(__file__))
        install_nbextension([os.path.join(filepath, "imolecule.min.js")],
                            verbose=0)
        _3d_initialized = True

    # Convert the relevant parts of the molecule into JSON for rendering
    atoms = [{"element": atom.GetSymbol(),
              "location": list(atomps[atom.GetIdx()])}
             for atom in mol.GetAtoms()]
    bonds = [{"atoms": [bond.GetBeginAtomIdx(),
                        bond.GetEndAtomIdx()],
              "order": int(bond.GetBondTypeAsDouble())}
             for bond in mol.GetBonds()]
    mol = {"atoms": atoms, "bonds": bonds}
    json_mol = json.dumps(mol, separators=(",", ":"))

    return """require(['%s'], function () {
           var $d = $('<div/>').attr('id', 'molecule_' + utils.uuid());
           $d.width(%d); $d.height(%d);
           $d.imolecule = jQuery.extend({}, imolecule);
           $d.imolecule.create($d, {drawingType: '%s', cameraType: '%s'});
           $d.imolecule.draw(%s);
           element.append($d);
           });""" % ((lib_path,) + size + (drawing_type_3d, camera_type_3d,
                     json_mol))


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
    img = Draw.ReactionToImage(rc, subImgSize=(int(molSize[0]/3), molSize[1]))
    sio = StringIO()
    img.save(sio, format='PNG')
    return sio.getvalue()


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
    sio = StringIO()
    img.save(sio, format='PNG')
    return sio.getvalue()


def InstallIPythonRenderer():
    rdchem.Mol._repr_svg_ = _toSVG
    rdchem.Mol._repr_png_ = _toPNG
    rdchem.Mol._repr_javascript_ = _toJSON
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
    del rdchem.Mol._repr_javascript_
    del rdChemReactions.ChemicalReaction._repr_png_
    if hasattr(rdchem.Mol, '__GetSubstructMatch'):
        rdchem.Mol.GetSubstructMatch = rdchem.Mol.__GetSubstructMatch
        del rdchem.Mol.__GetSubstructMatch
    if hasattr(rdchem.Mol, '__GetSubstructMatches'):
        rdchem.Mol.GetSubstructMatches = rdchem.Mol.__GetSubstructMatches
        del rdchem.Mol.__GetSubstructMatches
    del Image.Image._repr_png_
