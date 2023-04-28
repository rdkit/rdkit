# $Id$
#
# Copyright (C) 2008  Greg Landrum
#
#  @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""piddleQt

This module implements the PIDDLE/Sping API for a Qt4 canvas

Bits have been shamelessly cobbled from piddleSVG.py

Greg Landrum (glandrum@users.sourceforge.net)
"""
"""
  Functionality implemented:
  x drawLine
  x drawPolygon
  x drawString
  x drawImage

  Known problems:

"""

import copy
import types
from math import *

from PyQt4 import QtCore, QtGui, QtSvg

from rdkit.sping import pid


def _ColorToQt(color):
  """ convenience function for converting a sping.pid color to a Qt color

  """
  if color == pid.transparent:
    return None
  else:
    return QtGui.QColor(int(color.red * 255), int(color.green * 255), int(color.blue * 255))


#class QCanvasRotText(QCanvasText):
class QCanvasRotText:
  """ used to draw (UGLY) rotated text

  """

  def __init__(self, txt, canvas, angle=0):
    QCanvasText.__init__(self, txt, canvas)
    self._angle = angle

  def draw(self, qP):
    qP.save()
    x = self.x()
    y = self.y()
    theta = -self._angle
    qP.rotate(theta)
    qP.translate(-x, -y)
    thetaR = theta * pi / 180.
    newX = cos(-thetaR) * x - sin(-thetaR) * y
    newY = sin(-thetaR) * x + cos(-thetaR) * y
    qP.translate(newX, newY)
    QCanvasText.draw(self, qP)
    qP.restore()


class QtCanvas(pid.Canvas):

  def __init__(self, scene, size=None, name='QtCanvas'):
    if size is None:
      size = scene.width(), scene.height()
    self.size = size
    pid.Canvas.__init__(self, size, name)
    self._scene = scene
    self._brush = QtGui.QBrush()
    self._pen = QtGui.QPen()
    self._font = QtGui.QApplication.font()
    self.objs = []
    self._initOutput()
    self.nObjs = 0

  def _initOutput(self):
    for obj in self.objs:
      if type(obj) == types.TupleType:
        self._scene.removeItem(obj[0])
      else:
        self._scene.removeItem(obj)
    self.objs = []
    self.nObjs = 0

  def _adjustFont(self, font):
    if font.face:
      self._font.setFamily(font.face)
    self._font.setBold(font.bold)
    self._font.setItalic(font.italic)
    self._font.setPointSize(font.size)
    self._font.setUnderline(font.underline)

  # public functions
  def clear(self):
    self._initOutput()

  def flush(self):
    self._scene.update()

  def save(self, file=None, format=None):
    self._scene.update()

  #------------- drawing methods --------------
  def drawLine(self, x1, y1, x2, y2, color=None, width=None, dash=None, **kwargs):
    "Draw a straight line between x1,y1 and x2,y2."
    # set color...
    if color:
      if color == pid.transparent:
        return
    elif self.defaultLineColor == pid.transparent:
      return
    else:
      color = self.defaultLineColor
    qColor = _ColorToQt(color)

    if width:
      w = width
    else:
      w = self.defaultLineWidth
    self._pen.setColor(qColor)
    self._pen.setWidth(w)
    if dash is not None:
      self._pen.setStyle(QtCore.Qt.DashLine)
      #dash = [float(x)/w for x in dash]
      dash = list(dash)
      self._pen.setDashPattern(dash)
    else:
      self._pen.setStyle(QtCore.Qt.SolidLine)

    l = self._scene.addLine(x1, y1, x2, y2, self._pen)

    if dash is not None:
      self._pen.setStyle(QtCore.Qt.SolidLine)

    self.nObjs += 1
    self.objs.append(l)

  def drawPolygon(self, pointlist, edgeColor=None, edgeWidth=None, fillColor=pid.transparent,
                  closed=0, dash=None, **kwargs):
    """drawPolygon(pointlist) -- draws a polygon
    pointlist: a list of (x,y) tuples defining vertices

    """
    pts = [QtCore.QPointF(x[0], x[1]) for x in pointlist]
    poly = QtGui.QPolygonF(pts)

    # set color for fill...
    filling = 0
    if fillColor:
      if fillColor != pid.transparent:
        filling = 1
        qColor = _ColorToQt(fillColor)
        self._brush.setColor(qColor)

    if filling:
      self._brush.setStyle(QtCore.Qt.SolidPattern)
    else:
      self._brush.setStyle(QtCore.Qt.NoBrush)

    # set color for edge...
    if not edgeColor:
      edgeColor = self.defaultLineColor
    qColor = _ColorToQt(edgeColor)
    if qColor:
      self._pen.setColor(qColor)

    # set edge width...
    if edgeWidth is None:
      edgeWidth = self.defaultLineWidth
    self._pen.setWidth(edgeWidth)
    self._pen.setJoinStyle(QtCore.Qt.RoundJoin)

    if dash is not None:
      self._pen.setStyle(QtCore.Qt.DashLine)
    else:
      self._pen.setStyle(QtCore.Qt.SolidLine)
    if not qColor:
      self._pen.setStyle(QtCore.Qt.NoPen)
    poly = self._scene.addPolygon(poly, self._pen, self._brush)
    self.nObjs += 1
    poly.setZValue(self.nObjs)
    self.objs.append(poly)

    if dash is not None:
      self._pen.setStyle(QtCore.Qt.SolidLine)

    self.nObjs += 1

  def drawString(self, s, x, y, font=None, color=None, angle=0, **kwargs):
    # set color...
    if color:
      if color == pid.transparent:
        return
    elif self.defaultLineColor == pid.transparent:
      return
    else:
      color = self.defaultLineColor
    if font is None:
      font = self.defaultFont

    qColor = _ColorToQt(color)

    if font is not None:
      self._adjustFont(font)
    txt = self._scene.addText(s, self._font)
    txt.setDefaultTextColor(qColor)
    txtH = txt.boundingRect().height()
    #txt.setPos(QtCore.QPointF(x,y-txtH/2))
    txt.setPos(QtCore.QPointF(x, y - txtH))
    #txt.setPos(QtCore.QPointF(x,y))
    if angle:
      txt.rotate(angle)
    #if angle != 0:
    #  txt = QCanvasRotText(s,self._scene,angle=angle)
    #else:
    #  txt = QCanvasText(s,self._scene)

    #txt.setColor(qColor)
    #txt.setVisible(1)
    #txt.setX(x)
    #y -= font.size
    #txt.setY(y)
    txt.setZValue(self.nObjs)
    self.nObjs += 1
    self.objs.append(txt)

  def drawImage(self, image, x, y, **kwargs):
    """
    """
    from io import StringIO
    sio = StringIO()
    image.save(sio, format='png')
    base = QtGui.QPixmap()
    base.loadFromData(sio.getvalue())
    pix = self._scene.addPixmap(base)
    pix.setPos(QtCore.QPointF(x, y))
    pix.setZValue(self.nObjs)
    self.nObjs += 1
    self.objs.append(pix)

  def stringBox(self, s, font=None):
    "Return the logical width and height of the string if it were drawn \
    in the current font (defaults to self.font)."

    if not font:
      font = self.defaultFont

    if font:
      self._adjustFont(font)
    t = QtGui.QGraphicsTextItem(s)
    t.setFont(self._font)
    rect = t.boundingRect()
    return rect.width(), rect.height()

  def stringWidth(self, s, font=None):
    "Return the logical width of the string if it were drawn \
    in the current font (defaults to self.font)."

    if not font:
      font = self.defaultFont

    if font:
      self._adjustFont(font)
    t = QtGui.QGraphicsTextItem(s)
    t.setFont(self._font)
    rect = t.boundingRect()
    return rect.width()

  def fontAscent(self, font=None):
    if not font:
      font = self.defaultFont

    if font:
      self._adjustFont(font)
    t = QtGui.QGraphicsTextItem('B')
    t.setFont(self._font)
    rect = t.boundingRect()
    # FIX: this is a hack, but I can't immediately figure out how to solve the
    #  problem that the bounding rectangle includes the descent:
    return 1.0 * rect.height()

  def fontDescent(self, font=None):
    if not font:
      font = self.defaultFont
    if font:
      self._adjustFont(font)
    t = QtGui.QGraphicsTextItem('B')
    t.setFont(self._font)
    rect1 = t.boundingRect()
    t = QtGui.QGraphicsTextItem('y')
    t.setFont(self._font)
    rect2 = t.boundingRect()
    return 1. * (rect2.height() - rect1.height())


def test(canvas):
  #... for testing...
  canvas.defaultLineColor = Color(0.7, 0.7, 1.0)  # light blue
  canvas.drawLines(map(lambda i: (i * 10, 0, i * 10, 300), range(30)))
  canvas.drawLines(map(lambda i: (0, i * 10, 300, i * 10), range(30)))
  canvas.defaultLineColor = black

  canvas.drawLine(10, 200, 20, 190, color=red)

  canvas.drawEllipse(130, 30, 200, 100, fillColor=yellow, edgeWidth=4)

  canvas.drawArc(130, 30, 200, 100, 45, 50, fillColor=blue, edgeColor=navy, edgeWidth=4)

  canvas.defaultLineWidth = 4
  canvas.drawRoundRect(30, 30, 100, 100, fillColor=blue, edgeColor=maroon)
  canvas.drawCurve(20, 20, 100, 50, 50, 100, 160, 160)

  #canvas.drawString("This is a test!", 30,130, Font(face="times",size=16,bold=1),
  #                color=green, angle=-45)

  #canvas.drawString("This is a test!", 30,130, color=red, angle=-45)

  polypoints = [(160, 120), (130, 190), (210, 145), (110, 145), (190, 190)]
  canvas.drawPolygon(polypoints, fillColor=lime, edgeColor=red, edgeWidth=3, closed=1)

  canvas.drawRect(200, 200, 260, 260, edgeColor=yellow, edgeWidth=5)
  canvas.drawLine(200, 260, 260, 260, color=green, width=5)
  canvas.drawLine(260, 200, 260, 260, color=red, width=5)

  canvas.flush()


def dashtest(canvas):
  #... for testing...
  canvas.defaultLineColor = Color(0.7, 0.7, 1.0)  # light blue
  canvas.drawLines(map(lambda i: (i * 10, 0, i * 10, 300), range(30)), dash=(3, 3))
  canvas.drawLines(map(lambda i: (0, i * 10, 300, i * 10), range(30)), dash=(3, 3))
  canvas.defaultLineColor = black

  canvas.drawLine(10, 200, 20, 190, color=red, dash=(3, 3))

  canvas.drawEllipse(130, 30, 200, 100, fillColor=yellow, edgeWidth=4, dash=(3, 3))

  canvas.drawArc(130, 30, 200, 100, 45, 50, fillColor=blue, edgeColor=navy, edgeWidth=4,
                 dash=(3, 3))

  canvas.defaultLineWidth = 4
  canvas.drawRoundRect(30, 30, 100, 100, fillColor=blue, edgeColor=maroon, dash=(3, 3))
  canvas.drawCurve(20, 20, 100, 50, 50, 100, 160, 160, dash=(3, 3))

  canvas.drawString("This is a test!", 30, 130, Font(face="times", size=16, bold=1), color=green,
                    angle=-45)

  canvas.drawString("This is a test!", 30, 130, color=red, angle=-45)

  polypoints = [(160, 120), (130, 190), (210, 145), (110, 145), (190, 190)]
  canvas.drawPolygon(polypoints, fillColor=lime, edgeColor=red, edgeWidth=3, closed=1, dash=(3, 3))

  canvas.drawRect(200, 200, 260, 260, edgeColor=yellow, edgeWidth=5, dash=(3, 3))
  canvas.drawLine(200, 260, 260, 260, color=green, width=5, dash=(3, 3))
  canvas.drawLine(260, 200, 260, 260, color=red, width=5, dash=(3, 3))

  canvas.flush()


if __name__ == '__main__':
  import sys

  from PIL import Image

  from rdkit.sping.pid import *
  app = QtGui.QApplication(sys.argv)
  w = QtGui.QGraphicsView()

  scene = QtGui.QGraphicsScene(0, 0, 300, 300)
  canv = QtCanvas(scene)
  test(canv)
  w.setScene(scene)
  w.show()
  sys.exit(app.exec_())
