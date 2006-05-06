#  $Id: MolCanvas.py 5175 2006-05-06 17:58:57Z glandrum $
#
#  Copyright (C) 2005-2006 Rational Discovery LLC
#    All Rights Reserved
#
import RDConfig
from qt import *
from qtcanvas import *
from qtGui.PiddleWindowImpl import PiddleCanvasView
from Chem.Draw.MolDrawing import MolDrawing
import cPickle,os,copy,types
import Chem
from Chem import rdDepictor

from Logger import Logger
from sping.Qt.pidQt import QtCanvas as Canvas
from sping.Qt import pidQt
import traceback

class PickModes:
  pickNone=0
  pickAtoms=1
  pickBonds=2
  pickHighlights=4

class MolCanvasView(PiddleCanvasView):
  """ a canvas class for interacting with molecules

  """
  def __init__(self,*args,**kwargs):
    PiddleCanvasView.__init__(self,*args,**kwargs)
    self.drawing = MolDrawing()
    self.rescaleOnResize=False
    self.mols = []
    self.fontSize=10
    self.dotsPerAngstrom=40
    self.includeAtomNumbers=False
    self.atomNumberOffset=1
    self.highlights=[]
    if kwargs.get('interactive',False):
      self.addLeftHandler(lambda x,y:self.leftClickHandler(y))
    self.pickMode=PickModes.pickNone
    self.pickSlop=2
    if kwargs.get('enableContextMenu',True):
      self.addRightHandler(self.contextMenuRequested)
    
  def initCanvas(self,size=None):
    #print '\n\ninitCanvas',self
    #traceback.print_stack(limit=3)
    if size is None:
      minSz = self.minimumSize()
      size= max(minSz.width(),self.visibleWidth()), \
            max(minSz.height(),self.visibleHeight())
      #print '\tsz:',size
    else:
      #print '\tpolicy:',size
      self.setSizePolicy(QSizePolicy(QSizePolicy.MinimumExpanding,QSizePolicy.MinimumExpanding,
                                     0,0,False))
      self.setMinimumSize(QSize(size[0],size[1]))
    self.qCanvas = QCanvas(size[0],size[1])
    self.logCanvas = Logger.Logger(Canvas,destCanvas=self.qCanvas,
                                   size=size,
                                   loggerFlushCommand='clear',
                                   loggerIgnore=['size'])
    self.setCanvas(self.qCanvas)
    self.initDrawing()
  def initDrawing(self):
    #print 'initDrawing',self
    self.clearHighlights()
    self.drawing = MolDrawing()
    self.drawing.atomLabelFontSize=self.fontSize
    self.drawing.dotsPerAngstrom=self.dotsPerAngstrom
    self.drawing.additionalLabelPadding=(3,0)
    self.drawing.includeAtomNumbers=self.includeAtomNumbers
    self.drawing.atomNumberOffset=self.atomNumberOffset
    self.logCanvas.clear()

  def contextMenuRequested(self,view,evt):
    mol = view.getMol()
    if not mol:
      return

    pop = QPopupMenu(view)
    fn1 = lambda x,y=view:y.sendSmilesToClipboard()
    pop.insertItem('Copy &Smiles',fn1)
    fn2 = lambda x,y=view:y.sendImageToClipboard()
    pop.insertItem('Copy &Image',fn2)
    pop.exec_loop(evt.globalPos())


  def addMol(self,mol,**kwargs):
    """

    NOTE:
     if mol has no conformers, a 2D depiction will be generated

    """
    if not mol.GetNumConformers():
      rdDepictor.Compute2DCoords(mol)
    self.drawing.AddMol(mol,canvas=self.logCanvas,**kwargs)
    self.mols.append((mol,kwargs))

  def clearHighlights(self):
    if hasattr(self,'highlights'):
      for idx,ell in self.highlights:
        ell.hide()
    self.highlights = []

  def addHighlight(self,x,y,idx=-1,
                   highlightRadius=0.5,highlightColor=(1,0,1),
                   transformIt=False):
    canvas = self.canvas()
    if not canvas:
      return
    rad = self.drawing.dotsPerAngstrom*highlightRadius
    if transformIt:
      x,y = self.drawing.transformPoint((x,y))
    ellipse = QCanvasEllipse(rad*2,rad*2,canvas)
    ellipse.setX(x)
    ellipse.setY(y)
    ellipse.setZ(-1000)
    brush = self.drawing.canvas._brush
    brush.setColor(QColor(int(highlightColor[0]*255),
                        int(highlightColor[1]*255),
                        int(highlightColor[2]*255)))
    brush.setStyle(Qt.SolidPattern)
    ellipse.setBrush(brush)

    ellipse.setVisible(1)
    ellipse.highlightColor = highlightColor
    ellipse.highlightRad = highlightRadius

    self.highlights.append((idx,ellipse))
    
  def highlightAtoms(self,atomIndices,highlightRadius=.5,
                     highlightColor=(1,0,1),append=False,whichMol=-1):
    """

      highlightRadius is in Angstroms

    """  
    canvas = self.canvas()
    if not canvas:
      return
    if not append:
      self.clearHighlights()
    rad = self.drawing.dotsPerAngstrom*highlightRadius
    activeMol=self.mols[whichMol][0]
    for idx in atomIndices:
      if self.drawing.atomPs[activeMol].has_key(idx):
        x,y = self.drawing.atomPs[activeMol][idx]
        self.addHighlight(x,y,idx,highlightRadius=highlightRadius,
                          highlightColor=highlightColor)
    self.update()
    canvas.update()

  def setMol(self,mol,**kwargs):
    """

    NOTE:
     if mol has no conformers, a 2D depiction will be generated

    """
    self.mols = []
    self.clearHighlights()
    #self.initDrawing()
    self.initCanvas()
    self.addMol(mol,**kwargs)
    self.update()
    self.canvas().update()

  def getMol(self):
    if self.mols:
      return self.mols[0][0]
    else:
      return None
    
  def clearCanvas(self):
    self.initDrawing()
    self.logCanvas.clear()
  def reset(self):
    self.clearCanvas()
    self.canvas().update()
    self.mols = []

  def refresh(self):
    highlights= [(x[0],x[1].highlightColor,x[1].highlightRad,(x[1].x(),x[1].y())) for x in self.highlights]
    self.initCanvas()
    self.clearHighlights()
    for mol,kwargs in self.mols:
      apply(self.drawing.AddMol,(mol,self.logCanvas),kwargs)
    for idx,color,rad,pos in highlights:
      if idx>=0:
        self.highlightAtoms((idx,),highlightColor=color,highlightRadius=rad,
                            append=True)
      else:
        self.addHighlight(pos[0],pos[1],idx,highlightColor=color,highlightRadius=rad,
                          append=True)
        
  def getAtomAtPos(self,x,y,tol=0,whichMol=-1):
    activeMol=self.mols[whichMol][0]
    for idx,pos in self.drawing.atomPs[activeMol].iteritems():
      if abs(pos[0]-x)<=tol and abs(pos[1]-y)<=tol:
        return idx
    return -1
  def leftClickHandler(self,evt):
    if self.pickMode&PickModes.pickAtoms:
      atIdx = self.getAtomAtPos(evt.x(),evt.y(),self.pickSlop)

  def sendSmilesToClipboard(self,isomericSmiles=True,kekuleSmiles=False):
    import copy
    if not self.mols:
      return
    smis = []
    for mol,keyw in self.mols:
      #if isomericSmiles or kekuleSmiles:
      #  # when these are provided, the molecule itself
      #  # can be modified
      #  mol = copy.deepcopy(mol)
      smi = Chem.MolToSmiles(mol,isomericSmiles=isomericSmiles,
                             kekuleSmiles=kekuleSmiles)
      smis.append(smi)
    txt = '.'.join(smis)
    clipB = QApplication.clipboard()
    clipB.setText(txt)
    return

  def sendImageToClipboard(self):
    if not self.mols:
      return
    self.refresh()
    x = self.contentsX()
    y = self.contentsY()
    width = self.visibleWidth()
    height = self.visibleHeight()

    pix = QPixmap(width,height)
    painter = QPainter(pix)
    self.drawContents(painter,x,y,width,height)
    painter.end()
    
    clipB = QApplication.clipboard()
    clipB.setPixmap(pix)
    return


if __name__ == '__main__':
  import sys,getopt
  from qtGui import Gui
  
  app,widg = Gui.Launcher(MolCanvasView,None,interactive=True)

  if len(sys.argv)<2:
    #m = Chem.MolFromSmiles('c1ncc(CC(=O)O)cc1')
    #m = Chem.MolFromSmiles('Cc1ccc(O)c(Nc2nc(NCC#N)ncc2C)c1')
    m = Chem.MolFromSmiles('c1ccccc1[C@@H](Cl)[C@@H](Br)(F)')
  else:
    args,extras = getopt.getopt(sys.argv[1:],'n',['smi'])
    useSmi=False
    for arg,val in args:
      if arg=='-n':
        widg.includeAtomNumbers=True
      elif arg=='--smi':
        useSmi=True
    m = None
    if len(extras):
      if not useSmi:
        m = Chem.MolFromMolFile(extras[0])
      else:
        m = Chem.MolFromSmiles(extras[0])
    else:
      m = Chem.MolFromSmiles('c1ccccc1[C@@H](Cl)[C@@H](Br)(F)')
  assert m
  Chem.Kekulize(m)

  widg.resize(310,310)
  widg.initCanvas()
  #widg.drawing.includeAtomNumbers=True
  widg.addMol(m)
  widg.refresh()
  widg.addRightHandler(lambda x,y:x.refresh())
  widg.pickMode=PickModes.pickAtoms
  widg.pickSlop=5
  #widg.highlightAtoms((2,))

  app.exec_loop()
  widg.destroy(1)

