#
# Copyright (C) 2000 greg Landrum
#
""" Tools to be used with Piddle/Sping canvases

"""
import math

class ArrowLineTool:
  """ supports drawing lines with arrows on the end(s)

  """
  def draw(self,x1,y1,x2,y2,**kwargs):
    """ draws a line with arrow(s)

      **Arguments**

        - x1,y1,x2,y2: the endpoints of the line

      **Optional Keyword Arguments**

        - canvas: the canvas on which to draw

        - color: the color of the line

        - width: the width of the line

        - arrowBegin: toggles drawing an arrow at the beginning of the line

        - arrowEnd: toggles drawing an arrow at the end of the line

        - arrowAngle: the angle of the arrow's head

        - arrowLen: the length of the sides of the arrow

        - fillArrow: toggles filling the arrow head

        - arrowFillColor: the color used to fill the arrow head

    """
    optArgs = ['canvas','color','width','arrowBegin','arrowEnd','arrowAngle','arrowLen',
               'fillArrow','arrowFillColor']
    for arg in optArgs:
      try:
        exec('%s = kwargs[%s]'%(arg,arg));
      except NameError:
        exec('%s = self.%s'%(arg,arg))
    # start out with just the line
    canvas = self.canvas
    canvas.drawLine(x1,y1,x2,y2,color=color,width=width)

    arrowAngle = float(arrowAngle)/180. * math.pi

    # draw the arrow
    if arrowEnd:
      lineAngle = math.atan2(y2-y1,x2-x1)
      ang1 = -(math.pi - lineAngle - arrowAngle)
      ang2 = (math.pi + lineAngle - arrowAngle)
      p1 = (x2 + arrowLen*math.cos(ang1),y2 + arrowLen*math.sin(ang1))
      p2 = (x2,y2)
      p3 = (x2 + arrowLen*math.cos(ang2),y2 + arrowLen*math.sin(ang2))
      if fillArrow:
        canvas.drawPolygon([p1,p2,p3],edgeColor=color,
                         edgeWidth=width,fillColor=arrowFillColor,
                         closed=1)
      else:
        canvas.drawLines([p1+p2,p2+p3],color=color,width=width)

    if arrowBegin:
      lineAngle = math.atan2(y1-y2,x1-x2)
      ang1 = -(math.pi - lineAngle - arrowAngle)
      ang2 = (math.pi + lineAngle - arrowAngle)
      p1 = (x1 + arrowLen*math.cos(ang1),y1 + arrowLen*math.sin(ang1))
      p2 = (x1,y1)
      p3 = (x1 + arrowLen*math.cos(ang2),y1 + arrowLen*math.sin(ang2))
      if fillArrow:
        canvas.drawPolygon([p1,p2,p3],edgeColor=color,
                         edgeWidth=width,fillColor=arrowFillColor,
                         closed=1)
      else:
        canvas.drawLines([p1+p2,p2+p3],color=color,width=width)


    
  def __init__(self,canvas,color=None,width=None,arrowBegin=0,arrowEnd=1,
               arrowAngle=20,arrowLen=10,fillArrow=0,arrowFillColor=None):
    """ Constructor

      **Arguments**

        - canvas: the canvas on which to draw

        - color: the color of the line

        - width: the width of the line

        - arrowBegin: toggles drawing an arrow at the beginning of the line

        - arrowEnd: toggles drawing an arrow at the end of the line

        - arrowAngle: the angle of the arrow's head

        - arrowLen: the length of the sides of the arrow

        - fillArrow: toggles filling the arrow head

        - arrowFillColor: the color used to fill the arrow head


    """
    self.canvas = canvas
    self.color = color
    self.width = width
    self.arrowBegin = arrowBegin
    self.arrowEnd = arrowEnd
    self.arrowAngle = arrowAngle
    self.arrowLen = arrowLen
    self.fillArrow = fillArrow
    self.arrowFillColor = arrowFillColor
    
