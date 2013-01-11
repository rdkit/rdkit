# $Id$
#
#  Copyright (C) 2008-2011 Greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from rdkit import Chem
from rdkit import RDConfig
import numpy
import math

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
  useFraction=0.85
  
  atomLabelFontFace = "sans"
  atomLabelFontSize = 12
  atomLabelMinFontSize = 7

  bondLineWidth = 1.2
  dblBondOffset = .3
  dblBondLengthFrac = .8

  defaultColor = (1,0,0)
  selectColor = (1,0,0)

  colorBonds=True
  noCarbonSymbols=True
  includeAtomNumbers=False
  atomNumberOffset=0
  radicalSymbol=u'\u2219'

  dash = (4,4)
  atomPs = None
  canvas = None
  canvasSize=None

  wedgeDashedBonds=True

  # used to adjust overall scaling for molecules that have been laid out with non-standard
  # bond lengths
  coordScale=1.0

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

  def __init__(self,canvas=None):
    self.canvas = canvas
    if canvas:
      self.canvasSize=canvas.size
    self.atomPs = {}
    self.boundingBoxes = {}

  def transformPoint(self,pos):
    res = [0,0]
    res[0] = (pos[0] + self.molTrans[0])*self.currDotsPerAngstrom*self.useFraction + self.drawingTrans[0]
    res[1] = self.canvasSize[1]-((pos[1] + self.molTrans[1])*self.currDotsPerAngstrom*self.useFraction + \
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
    offsetX = math.cos(perp)*self.dblBondOffset*self.currDotsPerAngstrom
    offsetY = math.sin(perp)*self.dblBondOffset*self.currDotsPerAngstrom

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
            p3 = self.transformPoint(conf.GetAtomPosition(a3.GetIdx())*self.coordScale)
            dx2 = p3[0] - p1[0]
            dy2 = p3[1] - p1[1]
            dotP = dx2*offsetX + dy2*offsetY
            if dotP < 0:
              perp += math.pi
              offsetX = math.cos(perp)*self.dblBondOffset*self.currDotsPerAngstrom
              offsetY = math.sin(perp)*self.dblBondOffset*self.currDotsPerAngstrom
            
    fracP1,fracP2 = self._getOffsetBondPts(p1,p2,
                                           offsetX,offsetY,
                                           lenFrac=lenFrac)
    return fracP1,fracP2
    
  def _drawWedgedBond(self,bond,pos,nbrPos,
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
      self.canvas.addCanvasPolygon(poly,color=color)
    elif self.wedgeDashedBonds and self.canvas.addCanvasDashedWedge:
      self.canvas.addCanvasDashedWedge(poly[0],poly[1],poly[2],color=color)
    else:
      self.canvas.addCanvasLine(pos,nbrPos,linewidth=width*2,color=color,
                    dashes=dash)
    
  def _drawBond(self,bond,atom,nbr,pos,nbrPos,conf,
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
          self._drawWedgedBond(bond,p1,p2,color=color,width=width)
        elif bDir==Chem.BondDir.BEGINDASH:
          self._drawWedgedBond(bond,p1,p2,color=color,width=width,
                               dash=self.dash)
      else:
        self.canvas.addCanvasLine(pos, nbrPos, linewidth=width, color=color, color2=color2)
    elif bType == Chem.BondType.DOUBLE:
      if bond.IsInRing() or (atom.GetDegree()!=1 and bond.GetOtherAtom(atom).GetDegree()!=1):
        self.canvas.addCanvasLine(pos,nbrPos,linewidth=width,color=color,color2=color2)
        fp1,fp2 = self._offsetDblBond(pos,nbrPos,bond,atom,nbr,conf)
        self.canvas.addCanvasLine(fp1,fp2,linewidth=width,color=color,color2=color2)
      else:
        fp1,fp2 = self._offsetDblBond(pos,nbrPos,bond,atom,nbr,conf,dir=.5,
                                      lenFrac=1.0)
        self.canvas.addCanvasLine(fp1,fp2,linewidth=width,color=color,color2=color2)
        fp1,fp2 = self._offsetDblBond(pos,nbrPos,bond,atom,nbr,conf,dir=-.5,
                                      lenFrac=1.0)
        self.canvas.addCanvasLine(fp1,fp2,linewidth=width,color=color,color2=color2)
    elif bType == Chem.BondType.AROMATIC:
      self.canvas.addCanvasLine(pos,nbrPos,linewidth=width,color=color,color2=color2)
      fp1,fp2 = self._offsetDblBond(pos,nbrPos,bond,atom,nbr,conf)
      self.canvas.addCanvasLine(fp1,fp2,linewidth=width,color=color,color2=color2,
                    dash=self.dash)
    elif bType == Chem.BondType.TRIPLE:
      self.canvas.addCanvasLine(pos,nbrPos,linewidth=width,color=color,color2=color2)
      fp1,fp2 = self._offsetDblBond(pos,nbrPos,bond,atom,nbr,conf)
      self.canvas.addCanvasLine(fp1,fp2,linewidth=width,color=color,color2=color2)
      fp1,fp2 = self._offsetDblBond(pos,nbrPos,bond,atom,nbr,conf,dir=-1)
      self.canvas.addCanvasLine(fp1,fp2,linewidth=width,color=color,color2=color2)
    else:
      self.canvas.addCanvasLine(pos, nbrPos, linewidth=width, color=color, color2=color2,
                                dash=(1,2))

      
      
  def scaleAndCenter(self,mol,conf,coordCenter=False,canvasSize=None,ignoreHs=False):
    self.currDotsPerAngstrom=self.dotsPerAngstrom
    self.currAtomLabelFontSize=self.atomLabelFontSize
    if canvasSize is None:
      canvasSize=self.canvasSize
    xAccum = 0
    yAccum = 0
    minX = 1e8
    minY = 1e8
    maxX = -1e8
    maxY = -1e8

    nAts = mol.GetNumAtoms()
    for i in range(nAts):
      if ignoreHs and mol.GetAtomWithIdx(i).GetAtomicNum()==1: continue
      pos = conf.GetAtomPosition(i)*self.coordScale
      xAccum += pos[0]
      yAccum += pos[1]
      minX = min(minX,pos[0])
      minY = min(minY,pos[1])
      maxX = max(maxX,pos[0])
      maxY = max(maxY,pos[1])

    dx = abs(maxX-minX)
    dy = abs(maxY-minY)
    xSize = dx*self.currDotsPerAngstrom
    ySize = dy*self.currDotsPerAngstrom

    if coordCenter:
      molTrans = -xAccum/nAts,-yAccum/nAts
    else:
      molTrans = -(minX+(maxX-minX)/2),-(minY+(maxY-minY)/2)
    self.molTrans = molTrans

    if xSize>=.95*canvasSize[0]:
      scale = .9*canvasSize[0]/xSize
      xSize*=scale
      ySize*=scale
      self.currDotsPerAngstrom*=scale
      self.currAtomLabelFontSize = max(self.currAtomLabelFontSize*scale,
                                       self.atomLabelMinFontSize)
    if ySize>=.95*canvasSize[1]:
      scale = .9*canvasSize[1]/ySize
      xSize*=scale
      ySize*=scale
      self.currDotsPerAngstrom*=scale
      self.currAtomLabelFontSize = max(self.currAtomLabelFontSize*scale,
                                       self.atomLabelMinFontSize)
    drawingTrans = canvasSize[0]/2,canvasSize[1]/2
    self.drawingTrans = drawingTrans

  def _drawLabel(self,label,pos,font,color=None,**kwargs):
    if color is None:
      color = self.defaultColor
    x1 = pos[0]
    y1 = pos[1]
    labelP = x1,y1
    self.canvas.addCanvasText(label,(x1,y1),font,color,**kwargs)
    
  def AddMol(self,mol,centerIt=True,molTrans=None,drawingTrans=None,
             highlightAtoms=[],confId=-1,flagCloseContactsDist=2,
             highlightMap=None, ignoreHs=False,highlightBonds=[],**kwargs):
    """Set the molecule to be drawn.

    Parameters:
      hightlightAtoms -- list of atoms to highlight (default [])
      highlightMap -- dictionary of (atom, color) pairs (default None)

    Notes:
      - specifying centerIt will cause molTrans and drawingTrans to be ignored
    """
    conf = mol.GetConformer(confId)
    if kwargs.has_key('coordScale'):
      self.coordScale=kwargs['coordScale']

    if centerIt:
      self.scaleAndCenter(mol,conf,ignoreHs=ignoreHs)
    else:
      if molTrans is None:
        molTrans = (0,0)
      self.molTrans = molTrans
      if drawingTrans is None:
        drawingTrans = (0,0)
      self.drawingTrans = drawingTrans

    font = Font(face=self.atomLabelFontFace,size=self.currAtomLabelFontSize)

    if not mol.HasProp('_drawingBondsWedged'):
      Chem.WedgeMolBonds(mol,conf)
      

    includeAtomNumbers = kwargs.get('includeAtomNumbers',self.includeAtomNumbers)
    self.atomPs[mol] = {}
    self.boundingBoxes[mol] = [0]*4
    self.activeMol = mol
    self.bondRings = mol.GetRingInfo().BondRings()
    for atom in mol.GetAtoms():
      if ignoreHs and atom.GetAtomicNum()==1:
        drawAtom=False
      else:
        drawAtom=True
      idx = atom.GetIdx()
      pos = self.atomPs[mol].get(idx,None)
      if pos is None:
        pos = self.transformPoint(conf.GetAtomPosition(idx)*self.coordScale)
        self.atomPs[mol][idx] = pos
        if drawAtom:
          self.boundingBoxes[mol][0]=min(self.boundingBoxes[mol][0],pos[0])
          self.boundingBoxes[mol][1]=min(self.boundingBoxes[mol][1],pos[1])
          self.boundingBoxes[mol][2]=max(self.boundingBoxes[mol][2],pos[0])
          self.boundingBoxes[mol][3]=max(self.boundingBoxes[mol][3],pos[1])
      if not drawAtom: continue
      nbrSum = [0,0]
      for bond in atom.GetBonds():
        nbr = bond.GetOtherAtom(atom)
        if ignoreHs and nbr.GetAtomicNum()==1: continue
        nbrIdx = nbr.GetIdx()
        if nbrIdx > idx:
          nbrPos = self.atomPs[mol].get(nbrIdx,None)
          if nbrPos is None:
            nbrPos = self.transformPoint(conf.GetAtomPosition(nbrIdx)*self.coordScale)
            self.atomPs[mol][nbrIdx] = nbrPos
            self.boundingBoxes[mol][0]=min(self.boundingBoxes[mol][0],nbrPos[0])
            self.boundingBoxes[mol][1]=min(self.boundingBoxes[mol][1],nbrPos[1])
            self.boundingBoxes[mol][2]=max(self.boundingBoxes[mol][2],nbrPos[0])
            self.boundingBoxes[mol][3]=max(self.boundingBoxes[mol][3],nbrPos[1])
            
          if highlightBonds and bond.GetIdx() in highlightBonds:
            width=2.0*self.bondLineWidth
            color = self.selectColor
            color2 = self.selectColor
          elif highlightAtoms and idx in highlightAtoms and nbrIdx in highlightAtoms:
            width=2.0*self.bondLineWidth
            color = self.selectColor
            color2 = self.selectColor
          elif highlightMap is not None and idx in highlightMap and nbrIdx in highlightMap:
            width=2.0*self.bondLineWidth
            color = highlightMap[idx]
            color2 = highlightMap[nbrIdx]
          else:
            width=self.bondLineWidth
            if self.colorBonds:
              color = self.elemDict.get(atom.GetAtomicNum(),(0,0,0))
              color2 = self.elemDict.get(nbr.GetAtomicNum(),(0,0,0))
            else:
              color = self.defaultColor
              color2= color
              
          # make sure we draw from the beginning to the end
          # (this was Issue400)
          if idx==bond.GetBeginAtomIdx():
            self._drawBond(bond,atom,nbr,pos,nbrPos,conf,
                           color=color,width=width,color2=color2)
          else:
            self._drawBond(bond,nbr,atom,nbrPos,pos,conf,
                           color=color2,width=width,color2=color)
        else:
          nbrPos = self.atomPs[mol][nbrIdx]
        nbrSum[0] += nbrPos[0]-pos[0]
        nbrSum[1] += nbrPos[1]-pos[1]
  
      iso = atom.GetIsotope()

      labelIt= not self.noCarbonSymbols or \
               atom.GetAtomicNum()!=6 or \
               atom.GetFormalCharge()!=0 or \
               atom.GetNumRadicalElectrons() or \
               includeAtomNumbers or \
               iso or \
               atom.HasProp('molAtomMapNumber') or \
               atom.GetDegree()==0
      orient=''
      if labelIt:
        if includeAtomNumbers:
          symbol = str(atom.GetIdx())
        else:
          base = atom.GetSymbol()
          nHs = atom.GetTotalNumHs()
          if nHs>0:
            if nHs>1:
              hs='H<sub>%d</sub>'%nHs
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
          if chg:
            chg = '<sup>%s</sup>'%chg

          if atom.GetNumRadicalElectrons():
            rad = self.radicalSymbol*atom.GetNumRadicalElectrons()
            rad = '<sup>%s</sup>'%rad
          else:
            rad = ''

          isotope=''
          if iso:
            isotope='<sup>%d</sup>'%atom.GetIsotope()
                   
          mapNum=''
          if atom.HasProp('molAtomMapNumber'):
            mapNum=':'+atom.GetProp('molAtomMapNumber')
          deg = atom.GetDegree()
          if deg>1 or nbrSum[0]<1:
            symbol = '%s%s%s%s%s%s'%(isotope,base,mapNum,hs,chg,rad)
          else:
            symbol = '%s%s%s%s%s%s'%(rad,chg,hs,isotope,base,mapNum)

          if deg==1:
            if abs(nbrSum[1])>1:
              islope=nbrSum[0]/abs(nbrSum[1])
            else:
              islope=nbrSum[0]
            if abs(islope)>.3:
              if islope>0:
                orient='W'
              else:
                orient='E'
            elif abs(nbrSum[1])>10:
              if nbrSum[1]>0:
                orient='N'
              else :
                orient='S'
          else:
            orient = 'C'
        if highlightMap and idx in highlightMap:
          color = highlightMap[idx]
        elif highlightAtoms and idx in highlightAtoms:
          color = self.selectColor
        else:
          color = self.elemDict.get(atom.GetAtomicNum(),(0,0,0))
        self._drawLabel(symbol, pos, font, color=color,orientation=orient)

    if flagCloseContactsDist>0:
      tol = flagCloseContactsDist*flagCloseContactsDist
      for i,atomi in enumerate(mol.GetAtoms()):
        pi = numpy.array(self.atomPs[mol][i])
        for j in range(i+1,mol.GetNumAtoms()):
          pj = numpy.array(self.atomPs[mol][j])
          d = pj-pi
          dist2 = d[0]*d[0]+d[1]*d[1]
          if dist2<=tol:
            self.canvas.addCanvasPolygon(((pi[0]-2*flagCloseContactsDist,
                                      pi[1]-2*flagCloseContactsDist),
                                     (pi[0]+2*flagCloseContactsDist,
                                      pi[1]-2*flagCloseContactsDist),
                                     (pi[0]+2*flagCloseContactsDist,
                                      pi[1]+2*flagCloseContactsDist),
                                     (pi[0]-2*flagCloseContactsDist,
                                      pi[1]+2*flagCloseContactsDist)),
                             color=(1.,0,0),
                             fill=False,stroke=True)
