# $Id$
#
#  Copyright (C) 2008 Greg Landrum
#
#   @@ All Rights Reserved  @@
#
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon

def addCanvasLine(canvas,p1,p2,color=(0,0,0),color2=None,**kwargs):
  if color2 and color2!=color:
    mp = (p1[0]+p2[0])/2.,(p1[1]+p2[1])/2.
    canvas.add_line(Line2D((p1[0],mp[0]),(p1[1],mp[1]),
                           color=color,**kwargs))
    canvas.add_line(Line2D((mp[0],p2[0]),(mp[1],p2[1]),
                           color=color2,**kwargs))
  else:
    canvas.add_line(Line2D((p1[0],p2[0]),(p1[1],p2[1]),
                           color=color,**kwargs))

def addCanvasText(canvas,text,pos,font,color=(0,0,0),**kwargs):
  canvas.annotate(text,(pos[0],pos[1]),color=color,verticalalignment='center',
                  horizontalalignment='center',weight=font.weight,
                  size=font.size,family=font.face,backgroundcolor="white")

def addCanvasPolygon(canvas,ps,color=(0,0,0),**kwargs):
  canvas.add_patch(Polygon(ps,linewidth=0,facecolor=color))

