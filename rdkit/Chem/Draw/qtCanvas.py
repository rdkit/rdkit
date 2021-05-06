#
#  Copyright (C) 2014 Seiji Matsuoka
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

from rdkit.Chem.Draw.canvasbase import CanvasBase
try:
  from PySide import QtGui, QtCore
except ImportError:
  try:
    from PyQt5 import QtGui, QtCore
  except ImportError:
    from PySide2 import QtGui, QtCore


class Canvas(CanvasBase):

  def __init__(self, size):
    self.size = size
    self.qsize = QtCore.QSize(*size)
    self.pixmap = QtGui.QPixmap(self.qsize)
    self.painter = QtGui.QPainter(self.pixmap)
    self.painter.setRenderHint(QtGui.QPainter.Antialiasing, True)
    self.painter.setRenderHint(QtGui.QPainter.SmoothPixmapTransform, True)
    self.painter.fillRect(0, 0, size[0], size[1], QtCore.Qt.white)

  def addCanvasLine(self, p1, p2, color=(0, 0, 0), color2=None, **kwargs):
    if 'dash' in kwargs:
      line_type = QtCore.Qt.DashLine
    else:
      line_type = QtCore.Qt.SolidLine
    qp1 = QtCore.QPointF(*p1)
    qp2 = QtCore.QPointF(*p2)
    qpm = QtCore.QPointF((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)
    if color2 and color2 != color:
      rgb = [int(c * 255) for c in color]
      pen = QtGui.QPen(QtGui.QColor(*rgb), 1, line_type)
      self.painter.setPen(pen)
      self.painter.drawLine(qp1, qpm)
      rgb2 = [int(c * 255) for c in color2]
      pen.setColor(QtGui.QColor(*rgb2))
      self.painter.setPen(pen)
      self.painter.drawLine(qpm, qp2)
    else:
      rgb = [int(c * 255) for c in color]
      pen = QtGui.QPen(QtGui.QColor(*rgb), 1, line_type)
      self.painter.setPen(pen)
      self.painter.drawLine(qp1, qp2)

  def addCanvasText(self, text, pos, font, color=(0, 0, 0), **kwargs):
    orientation = kwargs.get('orientation', 'E')
    qfont = QtGui.QFont("Helvetica", font.size * 1.5)
    qtext = QtGui.QTextDocument()
    qtext.setDefaultFont(qfont)
    colored = [int(c * 255) for c in color]
    colored.append(text)
    html_format = "<span style='color:rgb({},{},{})'>{}</span>"
    formatted = html_format.format(*colored)
    qtext.setHtml(formatted)
    if orientation == 'N':
      qpos = QtCore.QPointF(pos[0] - qtext.idealWidth() / 2, pos[1] - font.size)
    elif orientation == 'W':
      qpos = QtCore.QPointF(pos[0] - qtext.idealWidth() + font.size, pos[1] - font.size)
    else:
      qpos = QtCore.QPointF(pos[0] - font.size, pos[1] - font.size)
    self.painter.save()
    self.painter.translate(qpos)
    qtext.drawContents(self.painter)
    self.painter.restore()
    return font.size * 1.8, font.size * 1.8, 0

  def addCanvasPolygon(self, ps, color=(0, 0, 0), fill=True, stroke=False, **kwargs):
    polygon = QtGui.QPolygonF()
    for ver in ps:
      polygon.append(QtCore.QPointF(*ver))
    pen = QtGui.QPen(QtGui.QColor(*color), 1, QtCore.Qt.SolidLine)
    self.painter.setPen(pen)
    self.painter.setBrush(QtGui.QColor(0, 0, 0))
    self.painter.drawPolygon(polygon)

  def addCanvasDashedWedge(self, p1, p2, p3, dash=(2, 2), color=(0, 0, 0), color2=None, **kwargs):
    rgb = [int(c * 255) for c in color]
    pen = QtGui.QPen(QtGui.QColor(*rgb), 1, QtCore.Qt.SolidLine)
    self.painter.setPen(pen)
    dash = (4, 4)
    pts1 = self._getLinePoints(p1, p2, dash)
    pts2 = self._getLinePoints(p1, p3, dash)
    if len(pts2) < len(pts1):
      pts2, pts1 = pts1, pts2
    for i in range(len(pts1)):
      qp1 = QtCore.QPointF(pts1[i][0], pts1[i][1])
      qp2 = QtCore.QPointF(pts2[i][0], pts2[i][1])
      self.painter.drawLine(qp1, qp2)

  def flush(self):
    self.painter.end()
