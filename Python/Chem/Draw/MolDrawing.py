# $Id$
#
#  Copyright (C) 2004-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
import Chem
import RDConfig
from sping import pid
import math

class MolDrawing(object):
  dotsPerAngstrom = 30

  atomLabelFontFace = "helvetica"
  atomLabelFontSize = 12
  atomLabelMinFontSize = 6
  additionalLabelPadding=(.1,.1) #(as a fraction of the font size)
  minLabelPadding=(2,2) 

  bondLineWidth = 1
  dblBondOffset = .2
  dblBondLengthFrac = .8

  defaultColor = pid.Color(0,0,0)
  selectColor = pid.Color(1,0,0)

  noCarbonSymbols=True
  includeAtomNumbers=False
  atomNumberOffset=0

  dash = (1,1)
  atomPs = None
  canvas = None
  
  def __init__(self,canvas=None):
    self.canvas = canvas
    self.atomPs = {}

  def transformPoint(self,pos):
    res = [0,0]
    res[0] = (pos[0] + self.molTrans[0])*self.dotsPerAngstrom + self.drawingTrans[0]
    res[1] = self.canvas.size[1]-((pos[1] + self.molTrans[1])*self.dotsPerAngstrom + self.drawingTrans[1])
    return res
  

  def _getBondOffset(self,p1,p2):
    # get the vector between the points:
    dx = p2[0]-p1[0]
    dy = p2[1]-p1[1]
    #print '\t(%d,%d)> %.2f %.2f |'%(a1.GetIdx()+1,a2.GetIdx()+1,dx,dy),

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
          a3 = otherBond.GetOtherAtom(a1)
          if a3.GetIdx() != a2Idx:
            p3 = self.transformPoint(conf.GetAtomPosition(a3.GetIdx()))
            dx2 = p3[0] - p1[0]
            dy2 = p3[1] - p1[1]
            #print '\ntest (%d,%d,%d):'%(a1Idx+1,a2Idx+1,a3.GetIdx()+1)
            dotP = dx2*offsetX + dy2*offsetY
            if dotP < 0:
              #print '\tadd'
              perp += math.pi
              offsetX = math.cos(perp)*self.dblBondOffset*self.dotsPerAngstrom
              offsetY = math.sin(perp)*self.dblBondOffset*self.dotsPerAngstrom
            
    fracP1,fracP2 = self._getOffsetBondPts(p1,p2,
                                           offsetX,offsetY,
                                           lenFrac=lenFrac)
    #print '(%.2f %.2f)'%(fracP2[0],fracP2[1]),
    #print

    return fracP1,fracP2
    
  def _drawWedgedBond(self,canvas,bond,pos,nbrPos,
                      width=bondLineWidth,color=defaultColor):
    perp,offsetX,offsetY = self._getBondOffset(pos,nbrPos)
    poly = ((pos[0],pos[1]),
            (nbrPos[0]+offsetX,nbrPos[1]+offsetY),
            (nbrPos[0]-offsetX,nbrPos[1]-offsetY))
    canvas.drawPolygon(poly,edgeColor=color,edgeWidth=1,fillColor=color,closed=1)


  def _drawBond(self,canvas,bond,atom,nbr,pos,nbrPos,conf,
                width=bondLineWidth,color=defaultColor):
    if bond.GetBondType() == Chem.BondType.SINGLE:
      dir = bond.GetBondDir()
      if dir==Chem.BondDir.BEGINWEDGE:
        self._drawWedgedBond(canvas,bond,pos,nbrPos,color=color,width=width)
      elif dir==Chem.BondDir.BEGINDASH:
        width *= 2
        canvas.drawLine(pos[0],pos[1],nbrPos[0],nbrPos[1],
                        color=color,
                        width=width,dash=self.dash)
      else:
        canvas.drawLine(pos[0],pos[1],nbrPos[0],nbrPos[1],
                        color=color,
                        width=width)
      
    elif bond.GetBondType() == Chem.BondType.DOUBLE:
      if bond.IsInRing() or (atom.GetDegree()!=1 and bond.GetOtherAtom(atom).GetDegree()!=1):
        canvas.drawLine(pos[0],pos[1],nbrPos[0],nbrPos[1],
                        color=color,
                        width=width)
        fp1,fp2 = self._offsetDblBond(pos,nbrPos,bond,atom,nbr,conf)
        canvas.drawLine(fp1[0],fp1[1],fp2[0],fp2[1],
                        color=color,
                        width=width)
      else:
        fp1,fp2 = self._offsetDblBond(pos,nbrPos,bond,atom,nbr,conf,dir=.5,
                                      lenFrac=1.0)
        canvas.drawLine(fp1[0],fp1[1],fp2[0],fp2[1],
                        color=color,
                        width=width)
        fp1,fp2 = self._offsetDblBond(pos,nbrPos,bond,atom,nbr,conf,dir=-.5,
                                      lenFrac=1.0)
        canvas.drawLine(fp1[0],fp1[1],fp2[0],fp2[1],
                        color=color,
                        width=width)
    elif bond.GetBondType() == Chem.BondType.AROMATIC:
      canvas.drawLine(pos[0],pos[1],nbrPos[0],nbrPos[1],
                      color=color,
                      width=width)
      fp1,fp2 = self._offsetDblBond(pos,nbrPos,bond,atom,nbr,conf)
      canvas.drawLine(fp1[0],fp1[1],fp2[0],fp2[1],
                      color=color,
                      width=width,
                      dash=self.dash)
    elif bond.GetBondType() == Chem.BondType.TRIPLE:
      canvas.drawLine(pos[0],pos[1],nbrPos[0],nbrPos[1],
                      color=color,
                      width=width)
      fp1,fp2 = self._offsetDblBond(pos,nbrPos,bond,atom,nbr,conf)
      canvas.drawLine(fp1[0],fp1[1],fp2[0],fp2[1],
                      color=color,
                      width=width)
      fp1,fp2 = self._offsetDblBond(pos,nbrPos,bond,atom,nbr,conf,dir=-1)
      canvas.drawLine(fp1[0],fp1[1],fp2[0],fp2[1],
                      color=color,
                      width=width)
      
      
  def _scaleAndCenter(self,mol,conf,coordCenter=False):
    canvasSize = self.canvas.size
    xAccum = 0
    yAccum = 0
    nAts = 0
    minX = 1e8
    minY = 1e8
    maxX = -1e8
    maxY = -1e8

    # FIX: this gets horked in cases where a molecule is T or L shaped
    # e.g.: CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)CCCCCCCCCCCC
    
    for atom in mol.GetAtoms():
      pos = conf.GetAtomPosition(atom.GetIdx())
      xAccum += pos[0]
      yAccum += pos[1]
      minX = min(minX,pos[0])
      minY = min(minY,pos[1])
      maxX = max(maxX,pos[0])
      maxY = max(maxY,pos[1])
      nAts += 1

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
    self.drawingTrans=(0,0)
    tMin = self.transformPoint((minX,minY))
    tMax = self.transformPoint((maxX,maxY))

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
    if 0:
      if tMin[0] < 0:
        self.drawingTrans = self.drawingTrans[0]-tMin[0]+5,self.drawingTrans[1]
      elif tMax[0] > canvasSize[0]:
        self.drawingTrans = self.drawingTrans[0]-(canvasSize[0]-tMax[0])-5,self.drawingTrans[1]

      # this is a leetle bit confusing because of the stupid negative
      # y transformation on the canvas:
      tMin,tMax = tMax,tMin
      if tMin[1] < 0:
        self.drawingTrans = self.drawingTrans[0],self.drawingTrans[1]-tMin[1]+5
      elif tMax[1] > canvasSize[1]:
        self.drawingTrans = self.drawingTrans[0],self.drawingTrans[1]-(canvasSize[1]-tMax[1])-5

    #tMax = self.transformPoint((maxX,maxY))
    #tMin = self.transformPoint((minX,minY))
    #self.canvas.drawRect(tMax[0],tMax[1],tMin[0],tMin[1],
    #                     edgeColor=pid.Color(0,0,0),
    #                     edgeWidth=1,
    #                     fillColor=pid.transparent)

  def _drawLabel(self,canvas,label,pos,font,
                 highlightIt=False,drawBox=True):
    txtWidth=canvas.stringWidth(label,font)
    txtHeight=canvas.fontAscent(font)

    x1 = pos[0]-txtWidth/2
    y1 = pos[1]+txtHeight/2
    labelP = x1,y1
    if drawBox:
      #xPad=max(self.additionalLabelPadding[0]*font.size,self.minLabelPadding[0])
      #yPad=max(self.additionalLabelPadding[1]*font.size,self.minLabelPadding[1])
      xPad=self.minLabelPadding[0]
      yPad=self.minLabelPadding[1]

      x1 = pos[0]-txtWidth/2 - xPad
      y1 = pos[1]+txtHeight/2 + yPad
      x2 = pos[0]+txtWidth/2 + xPad
      y2 = pos[1]-txtHeight/2 - yPad
      canvas.drawRect(x1,y1,x2,y2,
                      edgeColor=pid.transparent,
                      #edgeColor=pid.Color(0,0,0),
                      edgeWidth=0,fillColor=pid.Color(1,1,1))
    if highlightIt:
      color = self.selectColor
    else:
      color = self.defaultColor
    canvas.drawString(label,labelP[0],labelP[1],font,color=color)

    
  def AddMol(self,mol,canvas=None,centerIt=True,molTrans=(0,0),drawingTrans=(0,0),
             highlightAtoms=[],confId=-1):
    """

    Notes:
      - specifying centerIt will cause molTrans and drawingTrans to be ignored
      
    """
    if canvas is None:
      canvas = self.canvas
    else:
      self.canvas = canvas
      
    conf = mol.GetConformer(confId)

    if centerIt:
      self._scaleAndCenter(mol,conf)
    else:
      self.molTrans = molTrans
      self.drawingTrans = drawingTrans
    font = pid.Font(face=self.atomLabelFontFace,size=self.atomLabelFontSize)

    Chem.WedgeMolBonds(mol,conf)
      
    self.atomPs[mol] = {}
    self.activeMol = mol
    for atom in mol.GetAtoms():
      idx = atom.GetIdx()
      #print '%s(%d):'%(atom.GetSymbol(),idx),
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
          else:
            width=self.bondLineWidth
            color = self.defaultColor
          # make sure we draw from the beginning to the end
          # (this was Issue400)
          if atom.GetIdx()==bond.GetBeginAtomIdx():
            self._drawBond(canvas,bond,atom,nbr,pos,nbrPos,conf,
                           color=color,width=width)
          else:
            self._drawBond(canvas,bond,nbr,atom,nbrPos,pos,conf,
                           color=color,width=width)
        else:
          nbrPos = self.atomPs[mol][nbrIdx]
        nbrSum[0] += nbrPos[0]-pos[0]
        nbrSum[1] += nbrPos[1]-pos[1]
  
          
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

      #symbol = str(atom.GetIdx()+1)
      labelIt= not self.noCarbonSymbols or \
               atom.GetAtomicNum()!=6 or \
               atom.GetFormalCharge()!=0
      if labelIt:
        self._drawLabel(canvas,symbol,pos,font,
                        highlightIt=(highlightAtoms and idx in highlightAtoms))

      if self.includeAtomNumbers:
        txt = str(idx+self.atomNumberOffset)
        txtWidth=canvas.stringWidth(txt,font)
        txtHeight=canvas.fontAscent(font)
        if labelIt:
          xp = pos[0]-txtWidth
          yp = pos[1]-txtHeight/4
        else:
          xp = pos[0]-txtWidth/2
          yp = pos[1]+txtHeight/2
          x1 = pos[0]-txtWidth/2 - self.additionalLabelPadding[0]
          y1 = pos[1]+txtHeight/2 + self.additionalLabelPadding[1]
          x2 = pos[0]+txtWidth/2 + self.additionalLabelPadding[0]
          y2 = pos[1]-txtHeight/2 - self.additionalLabelPadding[1]
          canvas.drawRect(x1,y1,x2,y2,
                          edgeColor=pid.transparent,
                          edgeWidth=0,fillColor=pid.Color(1,1,1))
          
        
        canvas.drawString(txt,xp,yp,font,color=self.defaultColor)
        
        
if __name__=='__main__':
  from Chem import rdDepictor
  from sping.SVG.pidSVG import SVGCanvas as Canvas
  from sping.ReportLab.pidReportLab import RLCanvas as Canvas
  canvas = Canvas(size=(300,300),name='mol.pdf')
  #canvas = Canvas(size=(300,300),name='test.svg')
  #from sping.PIL.pidPIL import PILCanvas as Canvas
  #canvas = Canvas(size=(300,300),name='test.png')
  #canvas.drawLine(0,100,200,100)
  #canvas.drawLine(100,0,100,200)
  drawer = MolDrawing(canvas=canvas)
  drawer.noCarbonSymbols = 1
  #mol = Chem.MolFromSmiles('C1C=CCC=C1')
  mol = Chem.MolFromSmiles('O=C1C(=O)C=C1')
  mol = Chem.MolFromSmiles('c1cc([C@](F)(Cl)Br)ccc1C(=O)O')
  #mol = Chem.MolFromSmiles('CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC')
  mol = Chem.MolFromSmiles('ClCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)CCCCCCCCCCCCCCCCCCCCCCCCCCF')
  #mol = Chem.MolFromSmiles('C1C=CC=CC=CC=1[C@](F)(Cl)Br')
  mol = Chem.MolFromSmiles('c1cc([C@@](F)(Cl)Br)cnc1C(=O)O')
  #mol = Chem.MolFromSmiles('c1cc([C@@](F)(Cl)Br)cnc1C(=O)OC#N')
  mol = Chem.MolFromSmiles('CN=Cc1ccccc1NCC(=O)O')
  #mol = Chem.MolFromSmiles('c1ccncc1')
  Chem.Kekulize(mol)
  rdDepictor.Compute2DCoords(mol)
  #print Chem.MolToMolBlock(mol)
  #drawer.AddMol(mol)
  #drawer.AddMol(mol,highlightAtoms=(0,1,2,7,8,9))
  drawer.AddMol(mol)
  canvas.save()
  
  
  
    
