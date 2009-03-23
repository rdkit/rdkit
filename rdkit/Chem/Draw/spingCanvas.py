# $Id$
#
#  Copyright (C) 2008 Greg Landrum
#
#   @@ All Rights Reserved  @@
#
from rdkit.sping import pid
import math

faceMap={'sans':'helvetica'}

def convertColor(color):
  color = pid.Color(color[0],color[1],color[2])
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

def addCanvasLine(canvas,p1,p2,color=(0,0,0),color2=None,**kwargs):
  if color2 and color2!=color:
    mp = (p1[0]+p2[0])/2.,(p1[1]+p2[1])/2.
    color = convertColor(color)
    canvas.drawLine(p1[0],p1[1],mp[0],mp[1],
                    color=color,
                    width=int(kwargs.get('linewidth',1)),
                    dash=kwargs.get('dash',None))
    color2 = convertColor(color2)
    canvas.drawLine(mp[0],mp[1],p2[0],p2[1],
                    color=color2,
                    width=int(kwargs.get('linewidth',1)),
                    dash=kwargs.get('dash',None))
  else:
    color = convertColor(color)
    width=kwargs.get('linewidth',1)
    canvas.drawLine(p1[0],p1[1],p2[0],p2[1],
                    color=color,
                    width=int(width),
                    dash=kwargs.get('dash',None))

def addCanvasText(canvas,text,pos,font,color=(0,0,0),**kwargs):
  font = pid.Font(face=faceMap[font.face],size=font.size)
  txtWidth=canvas.stringWidth(text,font)
  txtHeight=canvas.fontAscent(font)
  labelP = pos[0]-txtWidth/2,pos[1]+txtHeight/2
  xPad = kwargs.get('xPadding',0)
  yPad = kwargs.get('yPadding',0)
  x1 = pos[0]-txtWidth/2 - xPad
  y1 = pos[1]+txtHeight/2 + yPad
  x2 = pos[0]+txtWidth/2 + xPad
  y2 = pos[1]-txtHeight/2 - yPad
  bgColor=kwargs.get('bgColor',(1,1,1))
  bgColor = convertColor(bgColor)
  canvas.drawRect(x1,y1,x2,y2,
                  edgeColor=pid.transparent,
                  edgeWidth=0,fillColor=bgColor)
  color = convertColor(color)
  canvas.drawString(text,labelP[0],labelP[1],font,color=color)

def addCanvasPolygon(canvas,ps,color=(0,0,0),**kwargs):
  edgeWidth=kwargs.get('lineWidth',0)
  edgeColor=pid.transparent
  color = convertColor(color)
  canvas.drawPolygon(ps,edgeColor=edgeColor,edgeWidth=int(edgeWidth),fillColor=color,closed=1)

def addCanvasDashedWedge(canvas,p1,p2,p3,dash=(2,2),color=(0,0,0),
                         color2=None,**kwargs):
  color = convertColor(color)
  dash = (4,4)
  pts1 = _getLinePoints(p1,p2,dash)
  pts2 = _getLinePoints(p1,p3,dash)

  if len(pts2)<len(pts1): pts2,pts1=pts1,pts2

  for i in range(len(pts1)):
    canvas.drawLine(pts1[i][0],pts1[i][1],pts2[i][0],pts2[i][1],
                    color=color,width=1)
