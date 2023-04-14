# Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""piddleQt

This module implements the PIDDLE/Sping API for a Qt canvas

Bits have been shamelessly cobbled from piddleSVG.py

Greg Landrum (Landrum@RationalDiscovery.com) 29 Octover, 2002
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
from math import *

from qt import *
from qtcanvas import *

from rdkit.sping import pid


def _ColorToQt(color):
  """ convenience function for converting a sping.pid color to a Qt color

  """
  if color == pid.transparent:
    return None
  else:
    return QColor(int(color.red * 255), int(color.green * 255), int(color.blue * 255))


class QCanvasRotText(QCanvasText):
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

  def __init__(self, destCanvas, size=(300, 300), name='QtCanvas'):
    self.size = size
    pid.Canvas.__init__(self, size, name)
    self._canvas = destCanvas
    self._brush = QBrush()
    self._pen = QPen()
    #self._font = QFont()
    self._font = QApplication.font()
    self.objs = []
    self._initOutput()
    self.nObjs = 0

  def _initOutput(self):
    for obj in self.objs:
      if type(obj) == tuple:
        obj[0].hide()
      else:
        obj.hide()
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
    self._canvas.update()

  def save(self, file=None, format=None):
    self._canvas.update()

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
    self._pen.setWidth(int(w))
    if dash is not None:
      self._pen.setStyle(Qt.DashLine)
    else:
      self._pen.setStyle(Qt.SolidLine)

    l = QCanvasLine(self._canvas)
    l.setPen(self._pen)
    l.setPoints(x1, y1, x2, y2)
    l.setVisible(1)
    l.setZ(self.nObjs)

    if dash is not None:
      self._pen.setStyle(Qt.SolidLine)

    self.nObjs += 1
    self.objs.append(l)

  def drawPolygon(self, pointlist, edgeColor=None, edgeWidth=None, fillColor=pid.transparent,
                  closed=0, dash=None, **kwargs):
    """drawPolygon(pointlist) -- draws a polygon
    pointlist: a list of (x,y) tuples defining vertices

    """
    pts = []
    for point in pointlist:
      pts += list(point)

    ptArr = QPointArray()
    ptArr.setPoints(pts)

    # set color for fill...
    filling = 0
    if fillColor:
      if fillColor != pid.transparent:
        filling = 1
        qColor = _ColorToQt(fillColor)
        self._brush.setColor(qColor)

    if filling:
      self._brush.setStyle(Qt.SolidPattern)
    else:
      self._brush.setStyle(Qt.NoBrush)

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
    self._pen.setJoinStyle(Qt.RoundJoin)

    if dash is not None:
      self._pen.setStyle(Qt.DashLine)
    else:
      self._pen.setStyle(Qt.SolidLine)

    poly = QCanvasPolygon(self._canvas)
    poly.setPen(self._pen)
    poly.setBrush(self._brush)
    poly.setPoints(ptArr)

    poly.setVisible(1)
    poly.setZ(self.nObjs)
    self.nObjs += 1
    self.objs.append(poly)

    # qt is moronic and doesn't draw the outlines of polygons
    if edgeColor != pid.transparent:
      for i in range(len(pointlist) - 1):
        l = QCanvasLine(self._canvas)
        l.setPoints(pointlist[i][0], pointlist[i][1], pointlist[i + 1][0], pointlist[i + 1][1])
        l.setPen(self._pen)
        l.setVisible(1)
        l.setZ(self.nObjs)
        self.objs.append(l)
      if closed:
        l = QCanvasLine(self._canvas)
        l.setPoints(pointlist[0][0], pointlist[0][1], pointlist[-1][0], pointlist[-1][1])
        l.setPen(self._pen)
        l.setVisible(1)
        l.setZ(self.nObjs)
        self.objs.append(l)

    if dash is not None:
      self._pen.setStyle(Qt.SolidLine)

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

    if angle != 0:
      txt = QCanvasRotText(s, self._canvas, angle=angle)
    else:
      txt = QCanvasText(s, self._canvas)

    txt.setTextFlags(Qt.AlignLeft | Qt.AlignVCenter)
    if self._font:
      txt.setFont(self._font)
    txt.setColor(qColor)
    txt.setVisible(1)
    txt.setX(x)
    y -= font.size
    txt.setY(y)
    txt.setZ(self.nObjs)
    self.nObjs += 1
    self.objs.append(txt)

  def drawImage(self, image, x1, y1, x2=None, y2=None, **kwargs):
    """
    """
    from io import StringIO
    sio = StringIO()
    image.save(sio, format='png')
    base = QPixmap()
    base.loadFromData(sio.getvalue())
    pm = QCanvasPixmap(base, QPoint(0, 0))
    pma = QCanvasPixmapArray()
    pma.setImage(0, pm)
    img = QCanvasSprite(pma, self._canvas)
    img.setVisible(1)
    img.setX(x1)
    img.setY(y1)
    self.objs.append((img, base, pm, pma))

  def stringWidth(self, s, font=None):
    "Return the logical width of the string if it were drawn \
    in the current font (defaults to self.font)."

    if not font:
      font = self.defaultFont

    if font:
      self._adjustFont(font)
    t = QCanvasText(s, self._canvas)
    t.setFont(self._font)
    rect = t.boundingRect()
    return rect.width()

  def fontAscent(self, font=None):
    if not font:
      font = self.defaultFont

    if font:
      self._adjustFont(font)
    t = QCanvasText('B', self._canvas)
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
    t = QCanvasText('B', self._canvas)
    t.setFont(self._font)
    rect1 = t.boundingRect()
    t = QCanvasText('y', self._canvas)
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

  from rdkit.sping.pid import *
  app = QApplication(sys.argv)
  w = QCanvasView()
  qCanv = QCanvas(300, 300)
  w.setCanvas(qCanv)
  canv = QtCanvas(qCanv)
  dashtest(canv)
  w.show()
  w.adjustSize()
  app.setMainWidget(w)
  app.exec_loop()
