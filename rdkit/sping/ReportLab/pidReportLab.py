# $Id$
# Copyright (C) 2005  Greg Landrum and Rational Discovery LLC
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
"""pidReportLab

Bits have been shamelessly cobbled from piddlePDF.py and/or
piddlePS.py

Greg Landrum (greg.landrum@gmail.com) 3/28/2005
"""
"""
  Functionality implemented:
  -drawLine
  -drawPolygon
  -drawEllipse
  -drawArc
  -drawCurve
  -drawString (rotated text is, mostly, fine... see below)
  -drawImage

"""

import os
import types
from math import *

from reportlab.graphics import shapes
from reportlab.lib import colors

from rdkit.sping.PDF import pdfmetrics, pidPDF
from rdkit.sping.pid import *


def colorToRL(color):
  if color != transparent:
    return colors.Color(color.red, color.green, color.blue)
  else:
    return None


class RLCanvas(Canvas):

  def __init__(self, size=(300, 300), name='RLCanvas'):
    self.size = size
    self._initOutput()
    Canvas.__init__(self, size, name)
    self.drawing = shapes.Drawing(size[0], size[1])

  def _initOutput(self):
    pass

  # public functions
  def clear(self):
    self._initOutput()

  def flush(self):
    # self.save('svg')
    pass  # to fit new definition of flush() -cwl

  def save(self, file=None, format=None):
    """Hand this either a file= <filename> or
    file = <an open file object>.  
    """
    if not file:
      file = self.name
    from reportlab.graphics import renderPDF
    renderPDF.drawToFile(self.drawing, file, self.name)

  def fixY(self, y):
    return self.size[1] - y

  # taken from pidPDF.py
  def _findPostScriptFontName(self, font):
    """Attempts to return proper font name."""
    #maps a piddle font to a postscript one.
    #step 1 - no face ends up serif, others are lowercased
    if not font.face:
      face = 'serif'
    else:
      face = font.face.lower()
    while face in pidPDF.font_face_map:
      face = pidPDF.font_face_map[face]
    #step 2, - resolve bold/italic to get the right PS font name
    psname = pidPDF.ps_font_map[(face, font.bold, font.italic)]
    return psname

  #------------- drawing methods --------------
  def drawLine(self, x1, y1, x2, y2, color=None, width=None, dash=None, **kwargs):
    "Draw a straight line between x1,y1 and x2,y2."
    # set color...
    if color:
      if color == transparent:
        return
    elif self.defaultLineColor == transparent:
      return
    else:
      color = self.defaultLineColor
    color = colorToRL(color)
    if width:
      w = width
    else:
      w = self.defaultLineWidth

    self.drawing.add(
      shapes.Line(x1, self.fixY(y1), x2, self.fixY(y2), strokeColor=color, strokeWidth=w,
                  strokeDashArray=dash))
    return

  def drawCurve(self, x1, y1, x2, y2, x3, y3, x4, y4, closed=0, **kwargs):
    "Draw a Bezier curve with control points x1,y1 to x4,y4."
    pts = self.curvePoints(x1, y1, x2, y2, x3, y3, x4, y4)
    if not closed:
      pointlist = [(pts[x][0], pts[x][1], pts[x + 1][0], pts[x + 1][1])
                   for x in range(len(pts) - 1)]
      self.drawLines(pointlist, **kwargs)
    else:
      self.drawPolygon(pointlist, closed=1, **kwargs)

  def drawArc(self, x1, y1, x2, y2, startAng=0, extent=360, edgeColor=None, edgeWidth=None,
              fillColor=None, dash=None, **kwargs):
    """Draw a partial ellipse inscribed within the rectangle x1,y1,x2,y2, 
    starting at startAng degrees and covering extent degrees.   Angles 
    start with 0 to the right (+x) and increase counter-clockwise. 
    These should have x1<x2 and y1<y2."""

    center = (x1 + x2) / 2, (y1 + y2) / 2
    pointlist = self.arcPoints(x1, y1, x2, y2, startAng, extent)

    # Fill...
    self.drawPolygon(pointlist + [center], edgeColor=transparent, edgeWidth=0, fillColor=fillColor)

    # Outline...
    pts = pointlist
    pointlist = [(pts[x][0], pts[x][1], pts[x + 1][0], pts[x + 1][1]) for x in range(len(pts) - 1)]
    self.drawLines(pointlist, edgeColor, edgeWidth, dash=dash, **kwargs)

  def drawPolygon(self, pointlist, edgeColor=None, edgeWidth=None, fillColor=transparent, closed=0,
                  dash=None, **kwargs):
    """drawPolygon(pointlist) -- draws a polygon
    pointlist: a list of (x,y) tuples defining vertices
    """
    if not edgeColor:
      edgeColor = self.defaultLineColor

    edgeColor = colorToRL(edgeColor)
    if not fillColor or fillColor == transparent:
      fillColor = None
    else:
      fillColor = colorToRL(fillColor)

    if edgeWidth:
      w = edgeWidth
    else:
      w = self.defaultLineWidth

    points = []
    for x, y in pointlist:
      points.append(x)
      points.append(self.fixY(y))
    self.drawing.add(
      shapes.Polygon(points, strokeColor=edgeColor, strokeWidth=w, strokeDashArray=dash,
                     fillColor=fillColor))

  def drawString(self, s, x, y, font=None, color=None, angle=0, **kwargs):
    # set color...
    if color:
      if color == transparent:
        return
    elif self.defaultLineColor == transparent:
      return
    else:
      color = self.defaultLineColor
    color = colorToRL(color)
    if font is None:
      font = self.defaultFont

    txt = shapes.String(0, 0, s, fillColor=color)
    txt.fontName = self._findPostScriptFontName(font)
    txt.fontSize = font.size

    g = shapes.Group(txt)
    g.translate(x, self.fixY(y))
    g.rotate(angle)
    self.drawing.add(g)
    return

  def drawImage(self, image, x1, y1, x2=None, y2=None, **kwargs):
    """
      to the best of my knowledge, the only real way to get an image
    """
    return

  def stringWidth(self, s, font=None):
    "Return the logical width of the string if it were drawn \
    in the current font (defaults to self.font)."

    if not font:
      font = self.defaultFont
    fontName = self._findPostScriptFontName(font)
    return pdfmetrics.stringwidth(s, fontName) * font.size * 0.001

  def fontAscent(self, font=None):
    if not font:
      font = self.defaultFont
    #return -font.size
    fontName = self._findPostScriptFontName(font)
    return pdfmetrics.ascent_descent[fontName][0] * 0.001 * font.size

  def fontDescent(self, font=None):
    if not font:
      font = self.defaultFont
    fontName = self._findPostScriptFontName(font)
    return -pdfmetrics.ascent_descent[fontName][1] * 0.001 * font.size


def test():
  #... for testing...
  canvas = RLCanvas(name="test")

  canvas.defaultLineColor = Color(0.7, 0.7, 1.0)  # light blue
  canvas.drawLines(map(lambda i: (i * 10, 0, i * 10, 300), range(30)))
  canvas.drawLines(map(lambda i: (0, i * 10, 300, i * 10), range(30)))
  canvas.defaultLineColor = black

  canvas.drawLine(10, 200, 20, 190, color=red)

  canvas.drawEllipse(130, 30, 200, 100, fillColor=yellow, edgeWidth=4)

  canvas.drawArc(130, 30, 200, 100, 45, 50, fillColor=blue, edgeColor=navy, edgeWidth=4)

  canvas.defaultLineWidth = 4
  canvas.drawRoundRect(30, 30, 100, 100, fillColor=blue, edgeColor=maroon, dash=(3, 3))
  canvas.drawCurve(20, 20, 100, 50, 50, 100, 160, 160)

  canvas.drawString("This is a test!", 30, 130, Font(face="times", size=16, bold=1), color=green,
                    angle=-45)

  canvas.drawString("This is a test!", 30, 130, color=red)

  polypoints = [(160, 120), (130, 190), (210, 145), (110, 145), (190, 190)]
  canvas.drawPolygon(polypoints, fillColor=lime, edgeColor=red, edgeWidth=3, closed=1)

  canvas.drawRect(200, 200, 260, 260, edgeColor=yellow, edgeWidth=5)
  canvas.drawLine(200, 260, 260, 260, color=green, width=5)
  canvas.drawLine(260, 200, 260, 260, color=red, width=5)

  canvas.save('test.pdf')


def dashtest():
  #... for testing...
  canvas = RLCanvas(name="test.pdf")

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

  canvas.save()


def testit(canvas, s, x, y, font=None):
  canvas.defaultLineColor = black
  canvas.drawString(s, x, y, font=font)
  canvas.defaultLineColor = blue
  w = canvas.stringWidth(s, font=font)
  canvas.drawLine(x, y, x + w, y)
  canvas.drawLine(x, y - canvas.fontAscent(font=font), x + w, y - canvas.fontAscent(font=font))
  canvas.drawLine(x, y + canvas.fontDescent(font=font), x + w, y + canvas.fontDescent(font=font))


def test2():

  canvas = RLCanvas(name="Foogar.pdf")
  testit(canvas, "Foogar", 20, 30)

  testit(canvas, "Foogar", 20, 90, font=Font(size=24))

  testit(canvas, "Foogar", 20, 150, font=Font(face='courier', size=24))

  testit(canvas, "Foogar", 20, 240, font=Font(face='courier'))
  canvas.flush()
  canvas.save()


if __name__ == '__main__':
  test()
  #dashtest()
  #test2()
