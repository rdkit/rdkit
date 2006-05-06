# $Id: PointCanvasItems.py 3196 2004-02-25 01:16:55Z glandrum $
#
#  Copyright (C) 2003  Rational Discovery LLC
#    All Rights Reserved
#
""" point classes for the point-sprinkler canvas

"""
import RDConfig
from qt import *
from qtcanvas import *
from sping import colors as SpingColors
import math

ptColors=((250,100,100),(100,100,250),(250,100,250),(50,250,50),(250,250,100))

class PointCanvasMixin:
  def __init__(self,*args,**kwargs):
    if kwargs.has_key('label'):
      self.id = kwargs['label']
      del kwargs['label']
    else:
      self.id='None'
    if kwargs.has_key('act'):
      act = kwargs['act']
      del kwargs['act']
    else:
      act = 0
    self._kwargs = kwargs  
    self._selectable=1
    self._activity = act

  def _initBrushes(self):
    self.setPen(QPen(QColor(255,255,255),0,Qt.SolidLine))
    self.setBrush(QBrush(QColor(0,0,0),Qt.SolidPattern))
    self.setActivity(self._activity)

  def setLabel(self,label):
    self.id = label
  def label(self):
    return self.id
    
  def setActivity(self,act):
    self._activity = act
    color = apply(QColor,ptColors[act%len(ptColors)])
    self.setColor(color)

  def activity(self):
    return self._activity

  def setColor(self,color):
    b = self.brush()
    b.setColor(color)
    self.setBrush(b)
    
  def spingRender(self,canvas):
    qPts = self.areaPoints()
    pts = []
    for i in range(qPts.size()):
      x,y = qPts.point(i)
      pts.append((x,y))
    color = self.brush().color()
    sColor = SpingColors.Color(color.red()/255.,color.green()/255.,color.blue()/255.)
    canvas.drawPolygon(pts,
                       edgeColor=SpingColors.transparent,
                       fillColor=sColor,
                       closed=1)

class PointCanvasCircle(PointCanvasMixin,QCanvasEllipse):
  """ origin is in center """
  def __init__(self,canvas,rad=10,**kwargs):
    PointCanvasMixin.__init__(self,**kwargs)
    args=(rad,rad,0,360*16,canvas)
    apply(QCanvasEllipse.__init__,tuple([self]+list(args)),self._kwargs)
    self._kwargs = None
    self._initBrushes()

  def spingRender(self,canvas):
    rad = self.width()/2.
    color = self.brush().color()
    sColor = SpingColors.Color(color.red()/255.,color.green()/255.,color.blue()/255.)
    canvas.drawEllipse(self.x()-rad,self.y()-rad,
                       self.x()+rad,self.y()+rad,
                       edgeColor=SpingColors.transparent,
                       fillColor=sColor)

class PointCanvasSquare(PointCanvasMixin,QCanvasRectangle):
  def __init__(self,canvas,size=10,**kwargs):
    PointCanvasMixin.__init__(self,**kwargs)
    args=(canvas,)
    apply(QCanvasRectangle.__init__,tuple([self]+list(args)),self._kwargs)
    self._kwargs = None
    self._initBrushes()
    self.setSize(size,size)
  def spingRender(self,canvas):
    side = self.width()
    color = self.brush().color()
    sColor = SpingColors.Color(color.red()/255.,color.green()/255.,color.blue()/255.)
    canvas.drawRect(self.x(),self.y(),
                       self.x()+side,self.y()+side,
                       edgeColor=SpingColors.transparent,
                       fillColor=sColor)
  def move(self,x,y):
    side = self.width()/2
    x -= side
    y -= side
    QCanvasRectangle.move(self,x,y)
    
class PointCanvasTriangle(PointCanvasMixin,QCanvasPolygon):
  """ point to the top, origin in center

  """
  def __init__(self,canvas,side=10,**kwargs):
    PointCanvasMixin.__init__(self,**kwargs)
    args=(canvas,)
    apply(QCanvasPolygon.__init__,tuple([self]+list(args)),self._kwargs)
    self._kwargs = None

    p0 = [-side/2,side/(2.*math.sqrt(3))]
    p1 = [side/2,side/(2.*math.sqrt(3))]
    p2 = [0,-side/math.sqrt(3)]
    self.setPoints(QPointArray(p0+p1+p2))
    self._initBrushes()


class PointCanvasDiamond(PointCanvasMixin,QCanvasPolygon):
  """ point to the top, origin in center

  """
  def __init__(self,canvas,side=10,**kwargs):
    PointCanvasMixin.__init__(self,**kwargs)
    args=(canvas,)
    apply(QCanvasPolygon.__init__,tuple([self]+list(args)),self._kwargs)
    self._kwargs = None

    p0 = [-side/2,0]
    p1 = [0,side/2]
    p2 = [side/2,0]
    p3 = [0,-side/2]
    self.setPoints(QPointArray(p0+p1+p2+p3))
    self._initBrushes()



shapes = [
  PointCanvasCircle,
  PointCanvasSquare,
  PointCanvasTriangle,
  PointCanvasDiamond,
          ]

