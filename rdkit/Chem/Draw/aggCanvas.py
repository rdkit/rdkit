# $Id$
#
#  Copyright (C) 2008 Greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from aggdraw import Brush, Pen
from aggdraw import Font
import math
from rdkit import RDConfig
import os

from aggdraw import Draw
from canvasbase import CanvasBase

faceMap={'sans':os.path.join(RDConfig.RDCodeDir,'Chem','Draw','FreeSans.ttf')}

def convertColor(color):
  color = (int(color[0]*255),int(color[1]*255),int(color[2]*255))
  return color

class Canvas(CanvasBase):

  def __init__(self, img):
    self.image = Draw(img)
    self.image.setantialias(True)
    self.size = self.image.size
    
  def _doLine(self, p1, p2, pen, **kwargs):
    if kwargs.get('dashes',(0,0)) == (0,0):
      self.image.line((p1[0],p1[1],p2[0],p2[1]),pen)
    else:
      # the antialiasing makes the dashes appear too small
      dash = [x*4 for x in kwargs['dashes']]
      pts = _getLinePoints(p1,p2,dash)

      currDash = 0
      dashOn = True
      while currDash<(len(pts)-1):
        if dashOn:
          p1 = pts[currDash]
          p2 = pts[currDash+1]
          self.image.line((p1[0],p1[1],p2[0],p2[1]),pen)
        currDash+=1
        dashOn = not dashOn

  def addCanvasLine(self, p1, p2, color=(0,0,0), color2=None, **kwargs):
    if color2 and color2!=color:
      mp = (p1[0]+p2[0])/2.,(p1[1]+p2[1])/2.
      color = convertColor(color)
      self._doLine(p1,mp,Pen(color,kwargs.get('linewidth',1)),**kwargs)
      color2 = convertColor(color2)
      self._doLine(mp,p2,Pen(color2,kwargs.get('linewidth',1)),**kwargs)
    else:
      color = convertColor(color)
      self._doLine(p1,p2,Pen(color,kwargs.get('linewidth',1)),**kwargs)

  def addCanvasText(self,text,pos,font,color=(0,0,0),**kwargs):
    color = convertColor(color)
    font = Font(color,faceMap[font.face],size=font.size)
    w,h=self.image.textsize(text,font)
    bw,bh=w*1.1,h*1.1
    dPos = pos[0]-bw/2.,pos[1]-bh/2.
    bgColor=kwargs.get('bgColor',(1,1,1))
    bgColor = convertColor(bgColor)
    self.image.rectangle((dPos[0],dPos[1],dPos[0]+bw,dPos[1]+bh),
                     None,Brush(bgColor))
    dPos = pos[0]-w/2.,pos[1]-h/2.
    self.image.text(dPos,text,font)

  def addCanvasPolygon(self,ps,color=(0,0,0),fill=True,stroke=False,**kwargs):
    if not fill and not stroke: return
    dps = []
    for p in ps:
      dps.extend(p)
    color = convertColor(color)
    brush=None
    pen=None
    if fill:
      brush = Brush(color)
    if stroke:
      pen = Pen(color)
    self.image.polygon(dps,pen,brush)
 
  def addCanvasDashedWedge(self,p1,p2,p3,dash=(2,2),color=(0,0,0),
                           color2=None,**kwargs):
    pen = Pen(color,kwargs.get('linewidth',1))
    dash = (3,3)
    pts1 = self._getLinePoints(p1,p2,dash)
    pts2 = self._getLinePoints(p1,p3,dash)

    if len(pts2)<len(pts1): pts2,pts1=pts1,pts2

    for i in range(len(pts1)):
      self.image.line((pts1[i][0],pts1[i][1],pts2[i][0],pts2[i][1]),pen)

  def flush(self):
    self.image.flush()
