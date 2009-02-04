# $Id$
#
#  Copyright (C) 2008 Greg Landrum
#
#   @@ All Rights Reserved  @@
#
import Chem
import RDConfig
import math

elemDict={
  7:(0,0,1),
  8:(1,0,0),
  9:(.2,.8,.8),
  15:(1,.5,0),
  16:(.8,.8,0),
  17:(0,.8,0),
  35:(.5,.3,.1),
  0:(.5,.5,.5),
  }

class Font(object):
  face='sans'
  size='12'
  weight='normal'
  name=None
  def __init__(self,face=None,size=None,name=None,weight=None):
    if face: self.face=face
    if size: self.size=size
    if name: self.name=name
    if weight: self.weight=weight

class MolDrawing(object):
  dotsPerAngstrom = 30

  atomLabelFontFace = "sans"
  atomLabelFontSize = 12
  atomLabelMinFontSize = 10

  bondLineWidth = 1.2
  dblBondOffset = .2
  dblBondLengthFrac = .8

  defaultColor = (1,0,0)
  selectColor = (1,0,0)

  colorBonds=True
  noCarbonSymbols=True
  includeAtomNumbers=False
  atomNumberOffset=0

  dash = (1,1)
  atomPs = None
  canvas = None
  canvasSize=None

  wedgeDashedBonds=True

  def __init__(self,canvas=None):
    self.canvas = canvas
    if canvas:
      self.canvasSize=canvas.size
    self.atomPs = {}

  def transformPoint(self,pos):
    res = [0,0]
    res[0] = (pos[0] + self.molTrans[0])*self.dotsPerAngstrom + self.drawingTrans[0]
    res[1] = self.canvasSize[1]-((pos[1] + self.molTrans[1])*self.dotsPerAngstrom + \
                                 self.drawingTrans[1])
    return res
  

  def _getBondOffset(self,p1,p2):
    # get the vector between the points:
    dx = p2[0]-p1[0]
    dy = p2[1]-p1[1]

    # figure out the angle and the perpendicular:
    ang = math.atan2(dy,dx)
    perp = ang + math.pi/2.

    # here's the offset for the parallel bond:
    offsetX = math.cos(perp)*self.dblBondOffset*self.dotsPerAngstrom
    offsetY = math.sin(perp)*self.dblBondOffset*self.dotsPerAngstrom

    return perp,offsetX,offsetY

  def _getOffsetBondPts(self,p1,p2,
                        offsetX,offsetY,
                        lenFrac=None):
    if not lenFrac:
      lenFrac = self.dblBondLengthFrac
      
    dx = p2[0]-p1[0]
    dy = p2[1]-p1[1]
    # ----
    # now figure out where to start and end it:
    
    # offset the start point:
    fracP1 = p1[0]+offsetX,p1[1]+offsetY

    # now move a portion of the way along the line to the neighbor:
    frac = (1.-lenFrac)/2
    fracP1 = fracP1[0]+dx*frac,\
             fracP1[1]+dy*frac

    fracP2 = fracP1[0]+dx*lenFrac,\
             fracP1[1]+dy*lenFrac
    return fracP1,fracP2

  def _offsetDblBond(self,p1,p2,bond,a1,a2,conf,dir=1,
                     lenFrac=None):
    perp,offsetX,offsetY = self._getBondOffset(p1,p2)
    offsetX = offsetX*dir
    offsetY = offsetY*dir
    
    # if we're a ring bond, we may need to flip over to the other side:
    if bond.IsInRing():
      bondIdx = bond.GetIdx()
      a1Idx = a1.GetIdx()
      a2Idx = a2.GetIdx()
      # find a ring bond from a1 to an atom other than a2:
      for otherBond in a1.GetBonds():
        if otherBond.GetIdx()!=bondIdx and \
           otherBond.IsInRing():
          sharedRing=False
          for ring in self.bondRings:
            if bondIdx in ring and otherBond.GetIdx() in ring:
              sharedRing=True
              break
          if not sharedRing:
            continue
          a3 = otherBond.GetOtherAtom(a1)
          if a3.GetIdx() != a2Idx:
            p3 = self.transformPoint(conf.GetAtomPosition(a3.GetIdx()))
            dx2 = p3[0] - p1[0]
            dy2 = p3[1] - p1[1]
            dotP = dx2*offsetX + dy2*offsetY
            if dotP < 0:
              perp += math.pi
              offsetX = math.cos(perp)*self.dblBondOffset*self.dotsPerAngstrom
              offsetY = math.sin(perp)*self.dblBondOffset*self.dotsPerAngstrom
            
    fracP1,fracP2 = self._getOffsetBondPts(p1,p2,
                                           offsetX,offsetY,
                                           lenFrac=lenFrac)
    return fracP1,fracP2
    
  def _drawWedgedBond(self,canvas,bond,pos,nbrPos,
                      width=bondLineWidth,color=defaultColor,
                      dash=None):
    perp,offsetX,offsetY = self._getBondOffset(pos,nbrPos)
    offsetX *=.75
    offsetY *=.75
    poly = ((pos[0],pos[1]),
            (nbrPos[0]+offsetX,nbrPos[1]+offsetY),
            (nbrPos[0]-offsetX,nbrPos[1]-offsetY))
    #canvas.drawPolygon(poly,edgeColor=color,edgeWidth=1,fillColor=color,closed=1)
    if not dash:
      addCanvasPolygon(canvas,poly,color=color)
    elif self.wedgeDashedBonds and addCanvasDashedWedge:
      addCanvasDashedWedge(canvas,poly[0],poly[1],poly[2],color=color)
    else:
      addCanvasLine(canvas,pos,nbrPos,linewidth=width*2,color=color,
                    dashes=dash)
    
  def _drawBond(self,canvas,bond,atom,nbr,pos,nbrPos,conf,
                width=bondLineWidth,color=defaultColor,color2=None):
    bType=bond.GetBondType()
    if bType == Chem.BondType.SINGLE:
      bDir = bond.GetBondDir()
      if bDir in (Chem.BondDir.BEGINWEDGE,Chem.BondDir.BEGINDASH):
        # if the bond is "backwards", change the drawing direction:
        if bond.GetBeginAtom().GetChiralTag() in (Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                                                  Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
          p1,p2 = pos,nbrPos
        else:
          p2,p1 = pos,nbrPos
        if bDir==Chem.BondDir.BEGINWEDGE:
          self._drawWedgedBond(canvas,bond,p1,p2,color=(0,0,0),width=width)
        elif bDir==Chem.BondDir.BEGINDASH:
          self._drawWedgedBond(canvas,bond,p1,p2,color=(0,0,0),width=width,
                               dash=self.dash)
      else:
        addCanvasLine(canvas,pos,nbrPos,linewidth=width,color=color,color2=color2)
    elif bType == Chem.BondType.DOUBLE:
      if bond.IsInRing() or (atom.GetDegree()!=1 and bond.GetOtherAtom(atom).GetDegree()!=1):
        addCanvasLine(canvas,pos,nbrPos,linewidth=width,color=color,color2=color2)
        fp1,fp2 = self._offsetDblBond(pos,nbrPos,bond,atom,nbr,conf)
        addCanvasLine(canvas,fp1,fp2,linewidth=width,color=color,color2=color2)
      else:
        fp1,fp2 = self._offsetDblBond(pos,nbrPos,bond,atom,nbr,conf,dir=.5,
                                      lenFrac=1.0)
        addCanvasLine(canvas,fp1,fp2,linewidth=width,color=color,color2=color2)
        fp1,fp2 = self._offsetDblBond(pos,nbrPos,bond,atom,nbr,conf,dir=-.5,
                                      lenFrac=1.0)
        addCanvasLine(canvas,fp1,fp2,linewidth=width,color=color,color2=color2)
    elif bType == Chem.BondType.AROMATIC:
      addCanvasLine(canvas,pos,nbrPos,linewidth=width,color=color,color2=color2)
      fp1,fp2 = self._offsetDblBond(pos,nbrPos,bond,atom,nbr,conf)
      addCanvasLine(canvas,fp1,fp2,linewidth=width,color=color,color2=color2,
                    dash=self.dash)
    elif bType == Chem.BondType.TRIPLE:
      addCanvasLine(canvas,pos,nbrPos,linewidth=width,color=color,color2=color2)
      fp1,fp2 = self._offsetDblBond(pos,nbrPos,bond,atom,nbr,conf)
      addCanvasLine(canvas,fp1,fp2,linewidth=width,color=color,color2=color2)
      fp1,fp2 = self._offsetDblBond(pos,nbrPos,bond,atom,nbr,conf,dir=-1)
      addCanvasLine(canvas,fp1,fp2,linewidth=width,color=color,color2=color2)
      
  def _scaleAndCenter(self,mol,conf,coordCenter=False):
    canvasSize = self.canvasSize
    xAccum = 0
    yAccum = 0
    minX = 1e8
    minY = 1e8
    maxX = -1e8
    maxY = -1e8

    nAts = mol.GetNumAtoms()
    for i in range(nAts):
      pos = conf.GetAtomPosition(i)
      xAccum += pos[0]
      yAccum += pos[1]
      minX = min(minX,pos[0])
      minY = min(minY,pos[1])
      maxX = max(maxX,pos[0])
      maxY = max(maxY,pos[1])

    dx = abs(maxX-minX)
    dy = abs(maxY-minY)
    xSize = dx*self.dotsPerAngstrom
    ySize = dy*self.dotsPerAngstrom

    if coordCenter:
      molTrans = -xAccum/nAts,-yAccum/nAts
    else:
      molTrans = -(minX+(maxX-minX)/2),-(minY+(maxY-minY)/2)
    self.dotsPerAngstrom=30.0
    self.molTrans = molTrans

    if xSize>=.95*canvasSize[0]:
      scale = .9*canvasSize[0]/xSize
      xSize*=scale
      ySize*=scale
      self.dotsPerAngstrom*=scale
      self.atomLabelFontSize = max(self.atomLabelFontSize*scale,
                                   self.atomLabelMinFontSize)
    if ySize>=.95*canvasSize[1]:
      scale = .9*canvasSize[1]/ySize
      xSize*=scale
      ySize*=scale
      self.dotsPerAngstrom*=scale
      self.atomLabelFontSize = max(self.atomLabelFontSize*scale,
                                   self.atomLabelMinFontSize)
    drawingTrans = canvasSize[0]/2,canvasSize[1]/2
    self.drawingTrans = drawingTrans
    tMax = self.transformPoint((maxX,maxY))
    tMin = self.transformPoint((minX,minY))

  def _drawLabel(self,canvas,label,pos,font,color=None,
                 highlightIt=False):
    if highlightIt:
      color = self.selectColor
    elif not color:
      color = self.defaultColor
    x1 = pos[0]
    y1 = pos[1]
    labelP = x1,y1
    addCanvasText(canvas,label,(x1,y1),font,color)
    
  def AddMol(self,mol,canvas=None,centerIt=True,molTrans=(0,0),drawingTrans=(0,0),
             highlightAtoms=[],confId=-1):
    """

    Notes:
      - specifying centerIt will cause molTrans and drawingTrans to be ignored
      
    """
    try:
      dl = addCanvasLine
    except NameError:
      registerCanvas('sping')
    if canvas is None:
      canvas = self.canvas
    else:
      self.canvas = canvas
    self.canvasSize=canvas.size
      
    conf = mol.GetConformer(confId)

    if centerIt:
      self._scaleAndCenter(mol,conf)
    else:
      self.molTrans = molTrans
      self.drawingTrans = drawingTrans
    font = Font(face=self.atomLabelFontFace,size=self.atomLabelFontSize)

    if not mol.HasProp('_drawingBondsWedged'):
      Chem.WedgeMolBonds(mol,conf)
      
    self.atomPs[mol] = {}
    self.activeMol = mol
    self.bondRings = mol.GetRingInfo().BondRings()
    for atom in mol.GetAtoms():
      idx = atom.GetIdx()
      pos = self.atomPs[mol].get(idx,None)
      if pos is None:
        pos = self.transformPoint(conf.GetAtomPosition(idx))
        self.atomPs[mol][idx] = pos
      nbrSum = [0,0]
      for bond in atom.GetBonds():
        nbr = bond.GetOtherAtom(atom)
        nbrIdx = nbr.GetIdx()
        if nbrIdx > idx:
          nbrPos = self.atomPs[mol].get(nbrIdx,None)
          if nbrPos is None:
            nbrPos = self.transformPoint(conf.GetAtomPosition(nbrIdx))
            self.atomPs[mol][nbrIdx] = nbrPos
            
          if highlightAtoms and idx in highlightAtoms and nbrIdx in highlightAtoms:
            width=2.0*self.bondLineWidth
            color = self.selectColor
            color2 = self.selectColor
          else:
            width=self.bondLineWidth
            if self.colorBonds:
              color = elemDict.get(atom.GetAtomicNum(),(0,0,0))
              color2 = elemDict.get(nbr.GetAtomicNum(),(0,0,0))
            else:
              color = self.defaultColor
              color2= color
              
          # make sure we draw from the beginning to the end
          # (this was Issue400)
          if idx==bond.GetBeginAtomIdx():
            self._drawBond(canvas,bond,atom,nbr,pos,nbrPos,conf,
                           color=color,width=width,color2=color2)
          else:
            self._drawBond(canvas,bond,nbr,atom,nbrPos,pos,conf,
                           color=color2,width=width,color2=color)
        else:
          nbrPos = self.atomPs[mol][nbrIdx]
        nbrSum[0] += nbrPos[0]-pos[0]
        nbrSum[1] += nbrPos[1]-pos[1]
  

      labelIt= not self.noCarbonSymbols or \
               atom.GetAtomicNum()!=6 or \
               atom.GetFormalCharge()!=0 or \
               self.includeAtomNumbers 
      if labelIt:
        if self.includeAtomNumbers:
          symbol = str(atom.GetIdx())
        else:
          base = atom.GetSymbol()
          nHs = atom.GetTotalNumHs()
          if nHs>0:
            if nHs>1:
              hs='H%d'%nHs
            else:
              hs ='H'
          else:
            hs = ''
          chg = atom.GetFormalCharge()
          if chg!=0:
            if chg==1:
              chg = '+'
            elif chg==-1:
              chg = '-'
            elif chg>1:
              chg = '+%d'%chg
            elif chg<-1:
              chg = '-%d'%chg
          else:
            chg = ''
          if nbrSum[0]<=0:
            symbol = '%s%s%s'%(base,hs,chg)
          else:
            symbol = '%s%s%s'%(chg,hs,base)

        color = elemDict.get(atom.GetAtomicNum(),(0,0,0))
        self._drawLabel(canvas,symbol,pos,font,color=color,
                        highlightIt=(highlightAtoms and idx in highlightAtoms))

def registerCanvas(canvasNm):
  g= globals()
  if canvasNm in ('sping','SPING'):
    from spingCanvas import addCanvasLine,addCanvasText,addCanvasPolygon,addCanvasDashedWedge
  elif canvasNm in ('agg','AGG'):
    from aggCanvas import addCanvasLine,addCanvasText,addCanvasPolygon,addCanvasDashedWedge
  elif canvasNm in ('mpl','MPL'):
    from mplCanvas import addCanvasLine,addCanvasText,addCanvasPolygon
    addCanvasDashedWedge=None
  else:
    raise ValueError,'unrecognized canvas type'
  g['addCanvasLine']=addCanvasLine
  g['addCanvasText']=addCanvasText
  g['addCanvasPolygon']=addCanvasPolygon
  g['addCanvasDashedWedge']=addCanvasDashedWedge
        
if __name__=='__main__':
  import sys
  if len(sys.argv)<2:
    mol = Chem.MolFromSmiles('O=C1C([C@@H](F)C=CN[C@H](Cl)Br)C(c2c(O)c(NN)ccc2)=C1C#N')
  else:
    mol = Chem.MolFromSmiles(sys.argv[1])
    
  #mol = Chem.MolFromSmiles('OC1C(O)CC1')
  Chem.Kekulize(mol)
  from Chem import rdDepictor
  rdDepictor.Compute2DCoords(mol)

  if 1:
    from aggdraw import Draw
    registerCanvas('agg')
    from PIL import Image
    img = Image.new("RGBA",(300,300),"white")
    canvas=Draw(img)
    canvas.setantialias(True)
    drawer = MolDrawing(canvas)
    drawer.AddMol(mol)
    canvas.flush()
    img.save("foo.png")
  elif 0:
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    registerCanvas('mpl')
    fig = Figure(figsize=(3,3))
    ax = fig.add_axes((0,0,1,1),xticks=[],yticks=[],frame_on=False)
    xd = ax.get_xlim()
    xd = xd[1]-xd[0]
    yd = ax.get_ylim()
    yd = yd[1]-yd[0]
    ax.size=(xd,yd)
    drawer = MolDrawing(ax)
    drawer.AddMol(mol)
    canv = FigureCanvasAgg(fig)
    canv.print_figure("foo.png",dpi=80)
  else:
    from sping.PDF.pidPDF import PDFCanvas as Canvas
    canvas = Canvas(size=(300,300),name='test.pdf')
    registerCanvas('sping')
    drawer = MolDrawing(canvas)
    drawer.AddMol(mol)
    canvas.save()
  
  
  
    
