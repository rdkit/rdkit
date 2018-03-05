#
#  Copyright (C) 2008 Greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import re

from rdkit.Chem.Draw.canvasbase import CanvasBase
from rdkit.sping import pid

faceMap = {'sans': 'helvetica', 'serif': 'times'}


def convertColor(color):
  color = pid.Color(color[0], color[1], color[2])
  return color


class Canvas(CanvasBase):

  def __init__(self, size, name, imageType='png'):
    if imageType == "pdf":
      from rdkit.sping.PDF.pidPDF import PDFCanvas as _Canvas
    elif imageType == "ps":
      from rdkit.sping.PS.pidPS import PSCanvas as _Canvas
    elif imageType == "svg":
      from rdkit.sping.SVG.pidSVG import SVGCanvas as _Canvas
    elif imageType == "png":
      from rdkit.sping.PIL.pidPIL import PILCanvas as _Canvas
    else:
      raise ValueError('unrecognized format: %s' % imageType)
    self.canvas = _Canvas(size=size, name=name)
    if hasattr(self.canvas, '_image'):
      self._image = self.canvas._image
    else:
      self._image = None
    self.size = size

  def addCanvasLine(self, p1, p2, color=(0, 0, 0), color2=None, **kwargs):
    if color2 and color2 != color:
      mp = (p1[0] + p2[0]) / 2., (p1[1] + p2[1]) / 2.
      color = convertColor(color)
      self.canvas.drawLine(p1[0], p1[1], mp[0], mp[1], color=color,
                           width=int(kwargs.get('linewidth', 1)), dash=kwargs.get('dash', None))
      color2 = convertColor(color2)
      self.canvas.drawLine(mp[0], mp[1], p2[0], p2[1], color=color2,
                           width=int(kwargs.get('linewidth', 1)), dash=kwargs.get('dash', None))
    else:
      color = convertColor(color)
      width = kwargs.get('linewidth', 1)
      self.canvas.drawLine(p1[0], p1[1], p2[0], p2[1], color=color, width=int(width),
                           dash=kwargs.get('dash', None))

  def addCanvasText(self, text, pos, font, color=(0, 0, 0), **kwargs):
    text = re.sub(r'\<.+?\>', '', text)
    font = pid.Font(face=faceMap[font.face], size=font.size)
    txtWidth, txtHeight = self.canvas.stringBox(text, font)
    bw, bh = txtWidth + txtHeight * 0.4, txtHeight * 1.4
    offset = txtWidth * pos[2]
    labelP = pos[0] - txtWidth / 2 + offset, pos[1] + txtHeight / 2
    color = convertColor(color)
    self.canvas.drawString(text, labelP[0], labelP[1], font, color=color)
    return (bw, bh, offset)

  def addCanvasPolygon(self, ps, color=(0, 0, 0), fill=True, stroke=False, **kwargs):
    if not fill and not stroke:
      return
    edgeWidth = kwargs.get('lineWidth', 0)
    edgeColor = pid.transparent
    color = convertColor(color)
    if not stroke:
      edgeWidth = 0
      edgeColor = pid.transparent
    else:
      edgeWidth = kwargs.get('lineWidth', 1)
      edgeColor = color
    if not fill:
      fillColor = pid.transparent
    else:
      fillColor = color
    self.canvas.drawPolygon(ps, edgeColor=edgeColor, edgeWidth=int(edgeWidth), fillColor=fillColor,
                            closed=1)

  def addCanvasDashedWedge(self, p1, p2, p3, dash=(2, 2), color=(0, 0, 0), color2=None, **kwargs):
    color = convertColor(color)
    dash = (4, 4)
    pts1 = self._getLinePoints(p1, p2, dash)
    pts2 = self._getLinePoints(p1, p3, dash)

    if len(pts2) < len(pts1):
      pts2, pts1 = pts1, pts2

    for i in range(len(pts1)):
      self.canvas.drawLine(pts1[i][0], pts1[i][1], pts2[i][0], pts2[i][1], color=color, width=1)

  def flush(self):
    self.canvas.flush()

  def save(self):
    self.canvas.save()
