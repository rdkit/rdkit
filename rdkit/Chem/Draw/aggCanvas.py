# $Id$
#
#  Copyright (C) 2008 Greg Landrum
#
#   @@ All Rights Reserved  @@
#
from aggdraw import Brush, Pen
from aggdraw import Font
import math
from rdkit import RDConfig
import os

faceMap={'sans':os.path.join(RDConfig.RDCodeDir,'Chem','Draw','FreeSans.ttf')}

def convertColor(color):
  color = (int(color[0]*255),int(color[1]*255),int(color[2]*255))
  return color

def _getLinePoints(p1,p2,dash):
  x1,y1=p1
  x2,y2=p2
  dx = x2-x1
  dy = y2-y1
  lineLen = math.sqrt(dx*dx+dy*dy)
  theta = math.atan2(dy,dx)
  cosT = math.cos(theta)
  sinT = math.sin(theta)

  pos = (x1,y1)
  pts = [pos]
  dist = 0
  currDash = 0
  while dist < lineLen:
    currL = dash[currDash%len(dash)]
    if(dist+currL > lineLen): currL = lineLen-dist
    endP = (pos[0] + currL*cosT, pos[1] + currL*sinT)
    pts.append(endP)
    pos = endP
    dist += currL
    currDash += 1
  return pts

def _doLine(canvas,p1,p2,pen,**kwargs):
  if kwargs.get('dashes',(0,0)) == (0,0):
    canvas.line((p1[0],p1[1],p2[0],p2[1]),pen)
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
        canvas.line((p1[0],p1[1],p2[0],p2[1]),pen)
      currDash+=1
      dashOn = not dashOn

def addCanvasLine(canvas,p1,p2,color=(0,0,0),color2=None,**kwargs):
  if color2 and color2!=color:
    mp = (p1[0]+p2[0])/2.,(p1[1]+p2[1])/2.
    color = convertColor(color)
    _doLine(canvas,p1,mp,Pen(color,kwargs.get('linewidth',1)),**kwargs)
    color2 = convertColor(color2)
    _doLine(canvas,mp,p2,Pen(color2,kwargs.get('linewidth',1)),**kwargs)
  else:
    color = convertColor(color)
    _doLine(canvas,p1,p2,Pen(color,kwargs.get('linewidth',1)),**kwargs)

def addCanvasText(canvas,text,pos,font,color=(0,0,0),**kwargs):
  color = convertColor(color)
  font = Font(color,faceMap[font.face],size=font.size)
  w,h=canvas.textsize(text,font)
  bw,bh=w*1.1,h*1.1
  dPos = pos[0]-bw/2.,pos[1]-bh/2.
  bgColor=kwargs.get('bgColor',(1,1,1))
  bgColor = convertColor(bgColor)
  canvas.rectangle((dPos[0],dPos[1],dPos[0]+bw,dPos[1]+bh),
                   None,Brush(bgColor))
  dPos = pos[0]-w/2.,pos[1]-h/2.
  canvas.text(dPos,text,font)

def addCanvasPolygon(canvas,ps,color=(0,0,0),**kwargs):
  dps = []
  for p in ps:
    dps.extend(p)
  color = convertColor(color)
  canvas.polygon(dps,None,Brush(color));

def addCanvasDashedWedge(canvas,p1,p2,p3,dash=(2,2),color=(0,0,0),
                         color2=None,**kwargs):
  pen = Pen(color,kwargs.get('linewidth',1))
  dash = (3,3)
  pts1 = _getLinePoints(p1,p2,dash)
  pts2 = _getLinePoints(p1,p3,dash)

  if len(pts2)<len(pts1): pts2,pts1=pts1,pts2

  for i in range(len(pts1)):
    canvas.line((pts1[i][0],pts1[i][1],pts2[i][0],pts2[i][1]),pen)
