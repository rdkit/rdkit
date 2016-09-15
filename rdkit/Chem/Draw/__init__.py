#
# Copyright (C) 2006-2016 Greg Landrum
#  All Rights Reserved
#
from __future__ import print_function
import os,re,sys
from rdkit.six import iteritems
from rdkit.Chem import AllChem
from rdkit.Chem.Draw.MolDrawing import MolDrawing,DrawingOptions
from rdkit.Chem.Draw.rdMolDraw2D import *
import logging
def _getCanvas():
  useAGG=False
  useCairo=False
  useSping=False
  Canvas=None
  if not os.environ.get('RDKIT_CANVAS',''):
    try:
      from rdkit.Chem.Draw.cairoCanvas import Canvas
      useCairo=True
    except ImportError:
      try:
        from rdkit.Chem.Draw.aggCanvas import Canvas
        useAGG=True
      except ImportError:
        from rdkit.Chem.Draw.spingCanvas import Canvas
        useSping=True
  else:
    canv=os.environ['RDKIT_CANVAS'].lower()
    if canv =='cairo':
      from rdkit.Chem.Draw.cairoCanvas import Canvas
      useCairo=True
    elif canv =='agg':
      from rdkit.Chem.Draw.aggCanvas import Canvas
      useAGG=True
    else:
      from rdkit.Chem.Draw.spingCanvas import Canvas
      useSping=True
  if useSping:
    DrawingOptions.radicalSymbol='.' #<- the sping canvas doesn't support unicode well
  return useAGG,useCairo,Canvas

def _createCanvas(size):
  useAGG,useCairo,Canvas=_getCanvas()
  if useAGG or useCairo:
    try:
      import Image
    except ImportError:
      from PIL import Image
    img = Image.new("RGBA",size,(0,0,0,0))
    canvas = Canvas(img)
  else:
    from rdkit.Chem.Draw.spingCanvas import Canvas
    canvas = Canvas(size=size,name='MolToImageFile')
    img = canvas._image
  return img,canvas

def MolToImage(mol, size=(300,300), kekulize=True, wedgeBonds=True,
               fitImage=False, options=None, canvas=None, **kwargs):
  """Returns a PIL image containing a drawing of the molecule

      ARGUMENTS:

        - kekulize: run kekulization routine on input `mol` (default True)

        - size: final image size, in pixel (default (300,300))

        - wedgeBonds: draw wedge (stereo) bonds (default True)

        - highlightAtoms: list of atoms to highlight (default [])

        - highlightMap: dictionary of (atom, color) pairs (default None)

        - highlightBonds: list of bonds to highlight (default [])

        - highlightColor: RGB color as tuple (default [1, 0, 0])

      NOTE:

            use 'matplotlib.colors.to_rgb()' to convert string and
            HTML color codes into the RGB tuple representation, eg.

              from matplotlib.colors import ColorConverter
              img = Draw.MolToImage(m, highlightAtoms=[1,2], highlightColor=ColorConverter().to_rgb('aqua'))
              img.save("molecule.png")

      RETURNS:

        a PIL Image object
  """

  if not mol:
    raise ValueError('Null molecule provided')
  if canvas is None:
    img,canvas=_createCanvas(size)
  else:
    img=None

  if options is None:
    options = DrawingOptions()
  if fitImage:
      options.dotsPerAngstrom = int(min(size) / 10)
  options.wedgeDashedBonds = wedgeBonds
  if 'highlightColor' in kwargs:
      color = kwargs.pop('highlightColor', (1, 0, 0))
      options.selectColor = color

  drawer = MolDrawing(canvas=canvas,drawingOptions=options)

  if kekulize:
    from rdkit import Chem
    mol = Chem.Mol(mol.ToBinary())
    Chem.Kekulize(mol)

  if not mol.GetNumConformers():
    from rdkit.Chem import AllChem
    AllChem.Compute2DCoords(mol)

  if 'legend' in kwargs:
    legend = kwargs['legend']
    del kwargs['legend']
  else:
    legend=''

  drawer.AddMol(mol,**kwargs)

  if legend:
    from rdkit.Chem.Draw.MolDrawing import Font
    bbox = drawer.boundingBoxes[mol]
    pos = size[0]/2,int(.94*size[1]),0 # the 0.94 is extremely empirical
    # canvas.addCanvasPolygon(((bbox[0],bbox[1]),(bbox[2],bbox[1]),(bbox[2],bbox[3]),(bbox[0],bbox[3])),
    #                         color=(1,0,0),fill=False,stroke=True)
    # canvas.addCanvasPolygon(((0,0),(0,size[1]),(size[0],size[1]),(size[0],0)   ),
    #                         color=(0,0,1),fill=False,stroke=True)
    font=Font(face='sans',size=12)
    canvas.addCanvasText(legend,pos,font)

  if kwargs.get('returnCanvas',False):
    return img,canvas,drawer
  else:
    canvas.flush()
    return img

def MolToFile(mol,fileName,size=(300,300),kekulize=True, wedgeBonds=True,
              imageType=None, fitImage=False, options=None, **kwargs):
  """ Generates a drawing of a molecule and writes it to a file
  """
  # original contribution from Uwe Hoffmann
  if not fileName:
    raise ValueError('no fileName provided')
  if not mol:
    raise ValueError('Null molecule provided')

  if imageType is None:
    imageType=os.path.splitext(fileName)[1][1:]

  if options is None:
    options = DrawingOptions()
  useAGG,useCairo,Canvas = _getCanvas()
  if fitImage:
      options.dotsPerAngstrom = int(min(size) / 10)
  options.wedgeDashedBonds = wedgeBonds
  if useCairo or useAGG:
    canvas = Canvas(size=size,imageType=imageType,
                              fileName=fileName)
  else:
    options.radicalSymbol = '.' #<- the sping canvas doesn't support unicode well
    canvas = Canvas(size=size,name=fileName,imageType=imageType)
  drawer = MolDrawing(canvas=canvas,drawingOptions=options)
  if kekulize:
    from rdkit import Chem
    mol = Chem.Mol(mol.ToBinary())
    Chem.Kekulize(mol)

  if not mol.GetNumConformers():
    from rdkit.Chem import AllChem
    AllChem.Compute2DCoords(mol)

  drawer.AddMol(mol,**kwargs)
  if useCairo or useAGG:
    canvas.flush()
  else:
    canvas.save()

def MolToImageFile(mol,filename,size=(300,300),kekulize=True, wedgeBonds=True,
                   **kwargs):
  """  DEPRECATED:  please use MolToFile instead

  """
  img = MolToImage(mol,size=size,kekulize=kekulize,wedgeBonds=wedgeBonds,**kwargs)
  img.save(filename)

tkRoot=None
tkLabel=None
tkPI=None
def ShowMol(mol,size=(300,300),kekulize=True,wedgeBonds=True,
            title='RDKit Molecule',**kwargs):
  """ Generates a picture of a molecule and displays it in a Tkinter window
  """
  global tkRoot,tkLabel,tkPI
  try:
    import Tkinter
  except ImportError:
    import tkinter as Tkinter
  try:
    import ImageTk
  except ImportError:
    from PIL import ImageTk

  img = MolToImage(mol,size,kekulize,wedgeBonds,**kwargs)

  if not tkRoot:
    tkRoot = Tkinter.Tk()
    tkRoot.title(title)
    tkPI = ImageTk.PhotoImage(img)
    tkLabel = Tkinter.Label(tkRoot,image=tkPI)
    tkLabel.place(x=0,y=0,width=img.size[0],height=img.size[1])
  else:
    tkPI.paste(img)
  tkRoot.geometry('%dx%d'%(img.size))


def MolToMPL(mol,size=(300,300),kekulize=True, wedgeBonds=True,
             imageType=None, fitImage=False, options=None, **kwargs):
  """ Generates a drawing of a molecule on a matplotlib canvas
  """
  if not mol:
    raise ValueError('Null molecule provided')
  from rdkit.Chem.Draw.mplCanvas import Canvas
  canvas = Canvas(size)
  if options is None:
    options = DrawingOptions()
    options.bgColor=None
  if fitImage:
      drawingOptions.dotsPerAngstrom = int(min(size) / 10)
  options.wedgeDashedBonds=wedgeBonds
  drawer = MolDrawing(canvas=canvas, drawingOptions=options)
  omol=mol
  if kekulize:
    from rdkit import Chem
    mol = Chem.Mol(mol.ToBinary())
    Chem.Kekulize(mol)

  if not mol.GetNumConformers():
    from rdkit.Chem import AllChem
    AllChem.Compute2DCoords(mol)

  drawer.AddMol(mol,**kwargs)
  omol._atomPs=drawer.atomPs[mol]
  for k,v in iteritems(omol._atomPs):
    omol._atomPs[k]=canvas.rescalePt(v)
  canvas._figure.set_size_inches(float(size[0])/100,float(size[1])/100)
  return canvas._figure

def calcAtomGaussians(mol,a=0.03,step=0.02,weights=None):
  """
useful things to do with these:
fig.axes[0].imshow(z,cmap=cm.gray,interpolation='bilinear',origin='lower',extent=(0,1,0,1))
fig.axes[0].contour(x,y,z,20,colors='k')

fig=Draw.MolToMPL(m);
contribs=Crippen.rdMolDescriptors._CalcCrippenContribs(m)
logps,mrs=zip(*contribs)
x,y,z=Draw.calcAtomGaussians(m,0.03,step=0.01,weights=logps)
fig.axes[0].imshow(z,cmap=cm.jet,interpolation='bilinear',origin='lower',extent=(0,1,0,1))
fig.axes[0].contour(x,y,z,20,colors='k',alpha=0.5)
fig.savefig('coumlogps.colored.png',bbox_inches='tight')


  """
  import numpy
  from matplotlib import mlab
  x = numpy.arange(0,1,step)
  y = numpy.arange(0,1,step)
  X,Y = numpy.meshgrid(x,y)
  if weights is None:
    weights=[1.]*mol.GetNumAtoms()
  Z = mlab.bivariate_normal(X,Y,a,a,mol._atomPs[0][0], mol._atomPs[0][1])*weights[0]
  for i in range(1,mol.GetNumAtoms()):
    Zp = mlab.bivariate_normal(X,Y,a,a,mol._atomPs[i][0], mol._atomPs[i][1])
    Z += Zp*weights[i]
  return X,Y,Z


def MolsToImage(mols, subImgSize=(200,200),legends=None,**kwargs):
  """
  """
  try:
    import Image
  except ImportError:
    from PIL import Image
  if legends is None: legends = [None]*len(mols)
  res = Image.new("RGBA",(subImgSize[0]*len(mols),subImgSize[1]))
  for i,mol in enumerate(mols):
    res.paste(MolToImage(mol,subImgSize,legend=legends[i],**kwargs),(i*subImgSize[0],0))
  return res

def _moltoimg(mol,sz,highlights,legend,**kwargs):
    try:
        import Image
    except ImportError:
        from PIL import Image
    from rdkit.Chem.Draw import rdMolDraw2D
    if not hasattr(rdMolDraw2D,'MolDraw2DCairo'):
        img = MolToImage(mol,sz,legend=legend,highlightAtoms=highlights,
                             **kwargs)
    else:
        nmol = rdMolDraw2D.PrepareMolForDrawing(mol,kekulize=kwargs.get('kekulize',True))
        d2d = rdMolDraw2D.MolDraw2DCairo(sz[0],sz[1])
        d2d.DrawMolecule(nmol,legend=legend,highlightAtoms=highlights)
        from io import BytesIO
        d2d.FinishDrawing()
        sio = BytesIO(d2d.GetDrawingText())
        img = Image.open(sio)
    return img

def _MolsToGridImage(mols,molsPerRow=3,subImgSize=(200,200),legends=None,
                    highlightAtomLists=None,**kwargs):
  """ returns a PIL Image of the grid
  """
  try:
    import Image
  except ImportError:
    from PIL import Image
  if legends is None: legends = ['']*len(mols)

  nRows = len(mols)//molsPerRow
  if len(mols)%molsPerRow : nRows+=1

  res = Image.new("RGBA",(molsPerRow*subImgSize[0],nRows*subImgSize[1]),(255,255,255,0))
  for i,mol in enumerate(mols):
    row = i//molsPerRow
    col = i%molsPerRow
    highlights=None
    if highlightAtomLists and highlightAtomLists[i]:
      highlights=highlightAtomLists[i]
    if mol is not None:
      img = _moltoimg(mol,subImgSize,highlights,legends[i],**kwargs)
      res.paste(img,(col*subImgSize[0],row*subImgSize[1]))
  return res

def _MolsToGridSVG(mols,molsPerRow=3,subImgSize=(200,200),legends=None,
                    highlightAtomLists=None,stripSVGNamespace=True,**kwargs):
  """ returns an SVG of the grid
  """
  matcher = re.compile(r'^(<.*>\n)(<svg:rect .*</svg\:rect>\n)(.*)</svg\:svg>',re.DOTALL)
  if legends is None: legends = ['']*len(mols)
  hdr=''
  ftr='</svg:svg>'
  rect=''

  nRows = len(mols)//molsPerRow
  if len(mols)%molsPerRow : nRows+=1

  blocks = ['']*(nRows*molsPerRow)

  fullSize=(molsPerRow*subImgSize[0],nRows*subImgSize[1])
  for i,mol in enumerate(mols):
    highlights=None
    if highlightAtomLists and highlightAtomLists[i]:
      highlights=highlightAtomLists[i]
    if mol is not None:
        nmol = rdMolDraw2D.PrepareMolForDrawing(mol,kekulize=kwargs.get('kekulize',True))
        d2d = rdMolDraw2D.MolDraw2DSVG(subImgSize[0],subImgSize[1])
        d2d.DrawMolecule(nmol,legend=legends[i],highlightAtoms=highlights)
        d2d.FinishDrawing()
        txt = d2d.GetDrawingText()
        h,r,b = matcher.match(txt).groups()
        if not hdr:
            hdr = h.replace("width='%dpx' height='%dpx' >"%subImgSize,"width='%dpx' height='%dpx' >"%fullSize)
        if not rect:
            rect = r
        blocks[i] = b
  for i,elem in enumerate(blocks):
    row = i//molsPerRow
    col = i%molsPerRow
    elem = rect+elem
    blocks[i] = '<g transform="translate(%d,%d)" >%s</g>'%(col*subImgSize[0],row*subImgSize[1],elem)
  res = hdr + '\n'.join(blocks)+ftr
  if stripSVGNamespace:
    res = res.replace('svg:','')
  return res

def MolsToGridImage(mols,molsPerRow=3,subImgSize=(200,200),legends=None,
                    highlightAtomLists=None,useSVG=False,**kwargs):
  if useSVG:
      return _MolsToGridSVG(mols,molsPerRow=molsPerRow,subImgSize=subImgSize,
                legends=legends, highlightAtomLists=highlightAtomLists, **kwargs)
  else:
      return _MolsToGridImage(mols,molsPerRow=molsPerRow,subImgSize=subImgSize,
                legends=legends, highlightAtomLists=highlightAtomLists, **kwargs)

def ReactionToImage(rxn, subImgSize=(200,200),**kwargs):
  """
  """
  try:
    import Image
  except ImportError:
    from PIL import Image

  mols = []
  for i in range(rxn.GetNumReactantTemplates()):
    tmpl=rxn.GetReactantTemplate(i)
    tmpl.UpdatePropertyCache(False)
    mols.append(tmpl)
  mols.append(None)
  for i in range(rxn.GetNumProductTemplates()):
    tmpl = rxn.GetProductTemplate(i)
    tmpl.UpdatePropertyCache(False)
    mols.append(tmpl)

  res = Image.new("RGBA",(subImgSize[0]*len(mols),subImgSize[1]),(255,255,255,0))
  for i,mol in enumerate(mols):
    if mol is not None:
      nimg = MolToImage(mol,subImgSize,kekulize=False,**kwargs)
    else:
      nimg,canvas = _createCanvas(subImgSize)
      p0 = (10,subImgSize[1]//2)
      p1 = (subImgSize[0]-10,subImgSize[1]//2)
      p3 = (subImgSize[0]-20,subImgSize[1]//2-10)
      p4 = (subImgSize[0]-20,subImgSize[1]//2+10)
      canvas.addCanvasLine(p0,p1,lineWidth=2,color=(0,0,0))
      canvas.addCanvasLine(p3,p1,lineWidth=2,color=(0,0,0))
      canvas.addCanvasLine(p4,p1,lineWidth=2,color=(0,0,0))
      if hasattr(canvas,'flush'):
        canvas.flush()
      else:
        canvas.save()
    res.paste(nimg,(i*subImgSize[0],0))
  return res

def is2D(conf):
  """returns True if a conf has 2d coords"""
  allZeroes = True
  for atompos in range(conf.GetNumAtoms()):
    x,y,z = conf.GetAtomPosition(atompos)
    if z != 0.0:
      return False
    if x != 0.0 or y != 0.0:
      allZeroes = False
  return not allZeroes

def _prepareRxnMol(mol):
    compute2D=False
    for conf in mol.GetConformers():
      compute2D = not is2D(conf)

    try:
      AllChem.SanitizeMol(mol, 
                          sanitizeOps=(AllChem.SanitizeFlags.SANITIZE_ALL^ 
                                       AllChem.SanitizeFlags.SANITIZE_KEKULIZE))
    except:
      pass
    if compute2D:
      AllChem.Compute2DCoords(mol)
  
def _cleanupRXN(rxn):
  """Cleanup reactions so the computed coords are sane for depiction"""
  for mol in rxn.GetReactants(): _prepareRxnMol(mol)
  for mol in rxn.GetProducts(): _prepareRxnMol(mol)
  return rxn

#source https://commons.wikimedia.org/wiki/File:Tab_plus.svg
#  creative commons license (need attribution?  Is link enough? Make our own?)

svg_plus = '''<path d="M -5,95 L 25.950102,95" style="opacity:0.5;fill:#000000;fill-opacity:1 " id="p1" />
<path d="M -5.4181532,95 L 25,95" style="opacity:0.5;fill:#000000;fill-opacity:0.5" id="p2" />
<text x="6.9755859" y="14.5" style="font-size:10px;font-style:normal;font-variant:normal;font-weight:normal;font-stretch:normal;text-align:start;line-height:125%;writing-mode:lr-tb;text-anchor:start;opacity:1;fill:#000000;fill-opacity:1;stroke:#000000;stroke-width:1;stroke-linecap:round;stroke-linejoin:round;stroke-miterlimit:4;stroke-dasharray:none;stroke-dashoffset:0;stroke-opacity:1;font-family:Arial" id="t" xml:space="preserve">
<tspan x="6.9755859" y="14.5" style="font-size:14px;text-align:start;text-anchor:start" id="ts">+</tspan></text>
'''

svg_arrow = '''<path
     d="M 14.770536,27.531746 25.418045,14.261996 14.729243,0.9901797 8.8071689,0.99322045 16.990033,11.257321 
        l -16.96481647,0.0038 0.0052585,5.903094 17.07121297,0.001 -8.2708412,10.366437 5.9396892,0 z" />
'''

from xml.dom import minidom
fontSizePattern = re.compile(".*font[-]size[:](.*)px")
def getBBox(svg, width, height):
  """svg, width, height -> minx, maxx, miny, maxy"""
  doc = minidom.parseString(svg)  # parseString also exists
  # get all vector paths
  path_strings = [path.getAttribute('d') for path
                  in doc.getElementsByTagName('svg:path')]
  minx = maxx = miny = maxy = None
  for i,path in enumerate(path_strings):
    if path and path[0] == "M":
      try:
        for x,y  in [map(float,el.split(",")) for el in path[1:].split(" ") if el]:
          if minx is None or x < minx:
            minx = x
          if miny is None or y < miny:
            miny = y
          if maxx is None or x > maxx:
            maxx = x
          if maxy is None or y > maxy:
            maxy = y
      except ValueError:
        pass

  # find text areas, assume text is left justified
  text = [(text.getAttribute('x'), text.getAttribute('y'),
           text.getAttribute('style'),
           text.getElementsByTagName('svg:tspan'))
          for text in doc.getElementsByTagName("svg:text")]

  # assumes fixed width, might be bogus, should at least be greater
  for x,y,style,tspans in text:
    try:
      x = float(x)
      y = float(y)
      fontSize = max( list(map(float, fontSizePattern.findall(style))) + [0] )
      # get the text width
      count = 0
      for tspan in tspans:
        for child in tspan.childNodes:
          if child.nodeValue:
            count += len(child.nodeValue)
            
      if not count:
        count = 1
      xfontSize = fontSize*count
      if minx is None or x < minx:
        minx = x
      if miny is None or y-fontSize < miny:
        miny = y - fontSize
      if maxx is None or x+xfontSize > maxx:
        maxx = x + xfontSize
      if maxy is None or y > maxy:
        maxy = y

        
    except:
      logging.exception()
      pass
    
  if minx is None: minx = 0
  if maxx is None: maxx = width
  if miny is None: miny = 0
  if maxy is None: maxy = height
  return minx, maxx, miny, maxy

def makeRect(minx, maxx, miny, maxy):
  """Make an svg rectangle for debugging purposes"""
  rect = "<svg:rect style='opacity:0.4;fill:#FF0000;stroke:#000000' width='%s' height='%s' x='%s' y='%s'> </svg:rect>" % (
    maxx-minx, maxy-miny,
    minx, miny)
  return rect


def ReactionToSVG(rxn, subImgSize=(200,200), stripSVGNamespace=True,
                  scaleRelative=True,
                  debugRender=False,
                  **kwargs):
  """
  """
  matcher = re.compile(r'^(<.*>\n)(<svg:rect .*</svg\:rect>\n)(.*)</svg\:svg>',re.DOTALL)
  rect_matcher = re.compile("^(.*width=['])[^']*(.*)")
  num_reactants = rxn.GetNumReactantTemplates()
  num_products = rxn.GetNumProductTemplates()
  # figure out sub image sizes
  # make a copy so we don't obliterate the original
  rxn = AllChem.ChemicalReaction(rxn)
  mols = list(rxn.GetReactants()) + list(rxn.GetProducts())

  # get the relative sizes of the molecules
  relativeSizes = [1.0] * len(mols)
  
  num_mols = len(mols)
  blocks = [''] * num_mols
  hdr = ''
  ftr='</svg:svg>'
  rect = ''

  _cleanupRXN(rxn)
  xOffset = 0

  xdelta = subImgSize[0]

  if scaleRelative:
    fontSizes = []
    for col,mol in enumerate(mols):
      scale = relativeSizes[col]
      img_width = int(subImgSize[0] * scale)
      img_height = int(subImgSize[1] * scale)
      nmol = rdMolDraw2D.PrepareMolForDrawing(mol,kekulize=kwargs.get('kekulize',False))
      d2d = rdMolDraw2D.MolDraw2DSVG(img_width, img_height)
      d2d.DrawMolecule(nmol)
      d2d.FinishDrawing()
      txt = d2d.GetDrawingText()
      fontSizes.append(max( list(map(float, fontSizePattern.findall(txt))) + [0] ))

    median = list(set(fontSizes[:]))
    median.sort()
    relFont = median[int(len(median)/2)]
    relativeSizes = [(relFont/fn) for fn in fontSizes]
    fm = max(relativeSizes)
    if fm > 1.0:
      relativeSizes = [rs * 1/fm for rs in relativeSizes]
    mol_width = 0
    for rs in relativeSizes:
      mol_width += subImgSize[0] * rs
  else:
    mol_width = subImgSize[0] * len(mols)
    
  svg_size = fullSize = (mol_width + 20*num_reactants + 25 + 20*num_products,
                         subImgSize[1])
  
      
  # for each molecule make an svg and scale place it in the appropriate
  #    position
  top_layer = []
  for col,mol in enumerate(mols):
    scale = relativeSizes[col]
    img_width = int(subImgSize[0] * scale)
    img_height = int(subImgSize[1] * scale)
    nmol = rdMolDraw2D.PrepareMolForDrawing(mol,kekulize=kwargs.get('kekulize',False))
    d2d = rdMolDraw2D.MolDraw2DSVG(img_width, img_height)
    d2d.DrawMolecule(nmol)
    d2d.FinishDrawing()
    txt = d2d.GetDrawingText()
    
    # compute approx bbox
    maxFontSize = max( list(map(float, fontSizePattern.findall(txt))) + [0] )
    minx, maxx, miny, maxy = getBBox(txt, img_width, img_height)
    
    h,r,b = matcher.match(txt).groups()
    if not hdr:
      # header for the WHOLE image
      hdr = h.replace("width='%dpx' height='%dpx' >"%(img_width, img_height),
                      "width='%dpx' height='%dpx' >\n<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='%dpx' height='%dpx' x='0' y='0'> </rect>"%(
                        fullSize[0], fullSize[1],
                        fullSize[0], fullSize[1],
                      ))

    if not rect:
      rect = r

    elem = rect + b + "\n"
    if debugRender:
      elem += makeRect(minx,maxx,miny,maxy) + "\n"

    # now fit the "real" bbox into the image by subtracting the
    #  bits we don't use
    xOffset -= minx

    elem = elem.replace("opacity:1.0", "opacity:0.0")
    blocks.append('<g transform="translate(%d,%d)" >%s</g>'%(
      xOffset, subImgSize[1]/2 - miny - (maxy-miny)/2,elem))

    # add a nice offset for the next one
    xOffset += xdelta*scale

    # add the plus signs
    if col < num_reactants - 1 or col >= num_reactants and col < len(mols)-1:
      start,end = rect_matcher.match(rect.replace("opacity:1.0", "opacity:0.0")).groups()
      elem = start + repr(20) + end + svg_plus
      top_layer.append('<g transform="translate(%d,%d) scale(1.5,1.5)" >%s</g>'%(
        xOffset-10,subImgSize[1]/2.0 - 10,elem))
      xOffset += 25

    # add the arrow
    if col == num_reactants-1:
      start,end = rect_matcher.match(rect.replace("opacity:1.0", "opacity:0.0")).groups()
      elem = start + repr(20) + end + svg_arrow
      top_layer.append('<g transform="translate(%d,%d) scale(1.5,1.5)" >%s</g>'%(
        xOffset,subImgSize[1]/2.0 - 15,elem))
      xOffset += 40

  res = hdr + '\n'.join(blocks+top_layer) + ftr 
  
  if stripSVGNamespace:
    res = res.replace('svg:','')

  return res

def MolToQPixmap(mol, size=(300,300), kekulize=True,  wedgeBonds=True,
                 fitImage=False, options=None, **kwargs):
    """ Generates a drawing of a molecule on a Qt QPixmap
    """
    if not mol:
        raise ValueError('Null molecule provided')
    from rdkit.Chem.Draw.qtCanvas import Canvas
    canvas = Canvas(size)
    if options is None:
        options = DrawingOptions()
    options.bgColor = None
    if fitImage:
        options.dotsPerAngstrom = int(min(size) / 10)
    options.wedgeDashedBonds=wedgeBonds
    if kekulize:
        from rdkit import Chem
        mol = Chem.Mol(mol.ToBinary())
        Chem.Kekulize(mol)
    if not mol.GetNumConformers():
        from rdkit.Chem import AllChem
        AllChem.Compute2DCoords(mol)
    drawer = MolDrawing(canvas=canvas, drawingOptions=options)
    drawer.AddMol(mol, **kwargs)
    canvas.flush()
    return canvas.pixmap

