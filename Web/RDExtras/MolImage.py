# $Id$
#
#  Copyright (C) 2004, 2005 Rational Discovery LLC
#  All Rights Reserved
#
import Chem
from Chem.Draw.MolDrawing import MolDrawing
from Chem import TemplateAlign
from sping.SVG.pidSVG import SVGCanvas as Canvas
from mod_python import apache
from utils import cactvs
import sys, os, tempfile


def gif(req, smiles, width=100, height=100, highlight='[]', frame=0, dblSize=0, **kwargs):
  req.content_type = 'image/gif'
  width = int(width)
  height = int(height)
  frame = int(frame)
  dblSize = int(dblSize)

  # FIX: unsafe:
  highlight = eval(highlight)
  imgD = ''
  if smiles:
    fName = tempfile.NamedTemporaryFile(suffix='.gif', delete=False).name
    cactvs.SmilesToGif(smiles, fName, (width, height), dblSize=dblSize, frame=frame)
    if os.path.exists(fName):
      imgD = open(fName, 'rb').read()
      try:
        os.unlink(fName)
      except OSError:
        pass
  return imgD


def svg(req, smiles, width=100, height=100, highlight='[]', frame=0, template='', numbers=0,
        **kwargs):
  req.content_type = 'image/svg+xml'
  width = int(width)
  height = int(height)
  frame = int(frame)
  # FIX: unsafe:
  highlight = eval(highlight)
  imgD = ''
  mol = None
  if smiles:
    mol = Chem.MolFromSmiles(smiles)
  if mol:
    if kwargs.get('kekulize', True):
      Chem.Kekulize(mol)
    if template and highlight:
      try:
        patt = Chem.MolFromSmiles(template)
        Chem.Compute2DCoords(patt)
        TemplateAlign.AlignMolToTemplate2D(mol, patt, highlight)
      except Exception:
        Chem.Compute2DCoords(mol)
    else:
      Chem.Compute2DCoords(mol)
    canvas = Canvas(size=(width, height))
    drawer = MolDrawing(canvas=canvas)
    if numbers and numbers != '0':
      drawer.includeAtomNumbers = True
      drawer.atomNumberOffset = 1
    drawer.noCarbonSymbols = 1
    drawer.AddMol(mol, highlightAtoms=highlight)
    svg = canvas._txt + '</svg>'
  else:
    svg = ''
  return svg
