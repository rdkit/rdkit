# $Id$
#
#  Copyright (C) 2008 Greg Landrum
#
#   @@ All Rights Reserved  @@
#
from aggdraw import Brush, Pen
from aggdraw import Font

faceMap={'sans':'./FreeSans.ttf'}

def convertColor(color):
  color = (int(color[0]*255),int(color[1]*255),int(color[2]*255))
  return color

def addCanvasLine(canvas,p1,p2,color=(0,0,0),color2=None,**kwargs):
  if color2 and color2!=color:
    mp = (p1[0]+p2[0])/2.,(p1[1]+p2[1])/2.
    color = convertColor(color)
    canvas.line((p1[0],p1[1],mp[0],mp[1]),Pen(color,kwargs.get('linewidth',1)))
    color2 = convertColor(color2)
    canvas.line((mp[0],mp[1],p2[0],p2[1]),Pen(color2,kwargs.get('linewidth',1)))
  else:
    color = convertColor(color)
    canvas.line((p1[0],p1[1],p2[0],p2[1]),Pen(color,kwargs.get('linewidth',1)))

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
