# $Id$
#
#  Copyright (C) 2008 Greg Landrum
#
#   @@ All Rights Reserved  @@
#
from sping import pid

faceMap={'sans':'helvetica'}

def convertColor(color):
  color = pid.Color(color[0],color[1],color[2])
  return color

def addCanvasLine(canvas,p1,p2,color=(0,0,0),color2=None,**kwargs):
  if color2 and color2!=color:
    mp = (p1[0]+p2[0])/2.,(p1[1]+p2[1])/2.
    color = convertColor(color)
    canvas.drawLine(p1[0],p1[1],mp[0],mp[1],
                    color=color,
                    width=kwargs.get('linewidth',1),
                    dash=kwargs.get('dash',None))
    color2 = convertColor(color2)
    canvas.drawLine(mp[0],mp[1],p2[0],p2[1],
                    color=color2,
                    width=kwargs.get('linewidth',1),
                    dash=kwargs.get('dash',None))
  else:
    color = convertColor(color)
    width=kwargs.get('linewidth',1)
    canvas.drawLine(p1[0],p1[1],p2[0],p2[1],
                    color=color,
                    width=width,
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
  canvas.drawPolygon(ps,edgeColor=edgeColor,edgeWidth=edgeWidth,fillColor=color,closed=1)
