# $Id$
#
# Copyright (C) 2001-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""Cluster tree visualization using Sping

"""

try:
  from rdkit.sping import pid
  piddle = pid
except ImportError:
  from rdkit.piddle import piddle
import numpy

from . import ClusterUtils


class VisOpts(object):
  """ stores visualization options for cluster viewing

    **Instance variables**

      - x/yOffset: amount by which the drawing is offset from the edges of the canvas

      - lineColor: default color for drawing the cluster tree

      - lineWidth: the width of the lines used to draw the tree

  """
  xOffset = 20
  yOffset = 20
  lineColor = piddle.Color(0, 0, 0)
  hideColor = piddle.Color(.8, .8, .8)
  terminalColors = [
    piddle.Color(1, 0, 0),
    piddle.Color(0, 0, 1),
    piddle.Color(1, 1, 0),
    piddle.Color(0, .5, .5),
    piddle.Color(0, .8, 0),
    piddle.Color(.5, .5, .5),
    piddle.Color(.8, .3, .3),
    piddle.Color(.3, .3, .8),
    piddle.Color(.8, .8, .3),
    piddle.Color(.3, .8, .8)
  ]
  lineWidth = 2
  hideWidth = 1.1
  nodeRad = 15
  nodeColor = piddle.Color(1., .4, .4)
  highlightColor = piddle.Color(1., 1., .4)
  highlightRad = 10


def _scaleMetric(val, power=2, min=1e-4):
  val = float(val)
  nval = pow(val, power)
  if nval < min:
    return 0.0
  else:
    return numpy.log(nval / min)


class ClusterRenderer(object):

  def __init__(self, canvas, size, ptColors=[], lineWidth=None, showIndices=0, showNodes=1,
               stopAtCentroids=0, logScale=0, tooClose=-1):
    self.canvas = canvas
    self.size = size
    self.ptColors = ptColors
    self.lineWidth = lineWidth
    self.showIndices = showIndices
    self.showNodes = showNodes
    self.stopAtCentroids = stopAtCentroids
    self.logScale = logScale
    self.tooClose = tooClose

  def _AssignPointLocations(self, cluster, terminalOffset=4):
    self.pts = cluster.GetPoints()
    self.nPts = len(self.pts)
    self.xSpace = float(self.size[0] - 2 * VisOpts.xOffset) / float(self.nPts - 1)
    ySize = self.size[1]
    for i in range(self.nPts):
      pt = self.pts[i]
      if self.logScale > 0:
        v = _scaleMetric(pt.GetMetric(), self.logScale)
      else:
        v = float(pt.GetMetric())
      pt._drawPos = (VisOpts.xOffset + i * self.xSpace,
                     ySize - (v * self.ySpace + VisOpts.yOffset) + terminalOffset)

  def _AssignClusterLocations(self, cluster):
    # first get the search order (top down)
    toDo = [cluster]
    examine = cluster.GetChildren()[:]
    while len(examine):
      node = examine.pop(0)
      children = node.GetChildren()
      if len(children):
        toDo.append(node)
        for child in children:
          if not child.IsTerminal():
            examine.append(child)
    # and reverse it (to run from bottom up)
    toDo.reverse()
    for node in toDo:
      if self.logScale > 0:
        v = _scaleMetric(node.GetMetric(), self.logScale)
      else:
        v = float(node.GetMetric())
      # average our children's x positions
      childLocs = [x._drawPos[0] for x in node.GetChildren()]
      if len(childLocs):
        xp = sum(childLocs) / float(len(childLocs))
        yp = self.size[1] - (v * self.ySpace + VisOpts.yOffset)
        node._drawPos = (xp, yp)

  def _DrawToLimit(self, cluster):
    """
      we assume that _drawPos settings have been done already
    """
    if self.lineWidth is None:
      lineWidth = VisOpts.lineWidth
    else:
      lineWidth = self.lineWidth

    examine = [cluster]
    while len(examine):
      node = examine.pop(0)
      xp, yp = node._drawPos
      children = node.GetChildren()
      if abs(children[1]._drawPos[0] - children[0]._drawPos[0]) > self.tooClose:
        # draw the horizontal line connecting things
        drawColor = VisOpts.lineColor
        self.canvas.drawLine(children[0]._drawPos[0], yp, children[-1]._drawPos[0], yp, drawColor,
                             lineWidth)
        # and draw the lines down to the children
        for child in children:
          if self.ptColors and child.GetData() is not None:
            drawColor = self.ptColors[child.GetData()]
          else:
            drawColor = VisOpts.lineColor
          cxp, cyp = child._drawPos
          self.canvas.drawLine(cxp, yp, cxp, cyp, drawColor, lineWidth)
          if not child.IsTerminal():
            examine.append(child)
          else:
            if self.showIndices and not self.stopAtCentroids:
              try:
                txt = str(child.GetName())
              except Exception:
                txt = str(child.GetIndex())
              self.canvas.drawString(txt, cxp - self.canvas.stringWidth(txt) / 2, cyp)

      else:
        # draw a "hidden" line to the bottom
        self.canvas.drawLine(xp, yp, xp, self.size[1] - VisOpts.yOffset, VisOpts.hideColor,
                             lineWidth)

  def DrawTree(self, cluster, minHeight=2.0):
    if self.logScale > 0:
      v = _scaleMetric(cluster.GetMetric(), self.logScale)
    else:
      v = float(cluster.GetMetric())
    if v <= 0:
      v = minHeight
    self.ySpace = float(self.size[1] - 2 * VisOpts.yOffset) / v

    self._AssignPointLocations(cluster)
    self._AssignClusterLocations(cluster)
    if not self.stopAtCentroids:
      self._DrawToLimit(cluster)
    else:
      raise NotImplementedError('stopAtCentroids drawing not yet implemented')


def DrawClusterTree(cluster, canvas, size, ptColors=[], lineWidth=None, showIndices=0, showNodes=1,
                    stopAtCentroids=0, logScale=0, tooClose=-1):
  """ handles the work of drawing a cluster tree on a Sping canvas

    **Arguments**

      - cluster: the cluster tree to be drawn

      - canvas:  the Sping canvas on which to draw

      - size: the size of _canvas_

      - ptColors: if this is specified, the _colors_ will be used to color
        the terminal nodes of the cluster tree.  (color == _pid.Color_)

      - lineWidth: if specified, it will be used for the widths of the lines
        used to draw the tree

   **Notes**

     - _Canvas_ is neither _save_d nor _flush_ed at the end of this

     - if _ptColors_ is the wrong length for the number of possible terminal
       node types, this will throw an IndexError

     - terminal node types are determined using their _GetData()_ methods

  """
  renderer = ClusterRenderer(canvas, size, ptColors, lineWidth, showIndices, showNodes,
                             stopAtCentroids, logScale, tooClose)
  renderer.DrawTree(cluster)


def _DrawClusterTree(cluster, canvas, size, ptColors=[], lineWidth=None, showIndices=0, showNodes=1,
                     stopAtCentroids=0, logScale=0, tooClose=-1):
  """ handles the work of drawing a cluster tree on a Sping canvas

    **Arguments**

      - cluster: the cluster tree to be drawn

      - canvas:  the Sping canvas on which to draw

      - size: the size of _canvas_

      - ptColors: if this is specified, the _colors_ will be used to color
        the terminal nodes of the cluster tree.  (color == _pid.Color_)

      - lineWidth: if specified, it will be used for the widths of the lines
        used to draw the tree

   **Notes**

     - _Canvas_ is neither _save_d nor _flush_ed at the end of this

     - if _ptColors_ is the wrong length for the number of possible terminal
       node types, this will throw an IndexError

     - terminal node types are determined using their _GetData()_ methods

  """
  if lineWidth is None:
    lineWidth = VisOpts.lineWidth
  pts = cluster.GetPoints()
  nPts = len(pts)
  if nPts <= 1:
    return
  xSpace = float(size[0] - 2 * VisOpts.xOffset) / float(nPts - 1)
  if logScale > 0:
    v = _scaleMetric(cluster.GetMetric(), logScale)
  else:
    v = float(cluster.GetMetric())
  ySpace = float(size[1] - 2 * VisOpts.yOffset) / v

  for i in range(nPts):
    pt = pts[i]
    if logScale > 0:
      v = _scaleMetric(pt.GetMetric(), logScale)
    else:
      v = float(pt.GetMetric())
    pt._drawPos = (VisOpts.xOffset + i * xSpace, size[1] - (v * ySpace + VisOpts.yOffset))


#     if not stopAtCentroids or not hasattr(pt, '_isCentroid'):
#       allNodes.remove(pt)  # allNodes not defined

  if not stopAtCentroids:
    allNodes = ClusterUtils.GetNodeList(cluster)
  else:
    allNodes = ClusterUtils.GetNodesDownToCentroids(cluster)

  while len(allNodes):
    node = allNodes.pop(0)
    children = node.GetChildren()
    if len(children):
      if logScale > 0:
        v = _scaleMetric(node.GetMetric(), logScale)
      else:
        v = float(node.GetMetric())
      yp = size[1] - (v * ySpace + VisOpts.yOffset)
      childLocs = [x._drawPos[0] for x in children]
      xp = sum(childLocs) / float(len(childLocs))
      node._drawPos = (xp, yp)
      if not stopAtCentroids or node._aboveCentroid > 0:
        for child in children:
          if ptColors != [] and child.GetData() is not None:
            drawColor = ptColors[child.GetData()]
          else:
            drawColor = VisOpts.lineColor
          if showNodes and hasattr(child, '_isCentroid'):
            canvas.drawLine(child._drawPos[0], child._drawPos[1] - VisOpts.nodeRad / 2,
                            child._drawPos[0], node._drawPos[1], drawColor, lineWidth)
          else:
            canvas.drawLine(child._drawPos[0], child._drawPos[1], child._drawPos[0],
                            node._drawPos[1], drawColor, lineWidth)
        canvas.drawLine(children[0]._drawPos[0], node._drawPos[1], children[-1]._drawPos[0],
                        node._drawPos[1], VisOpts.lineColor, lineWidth)
      else:
        for child in children:
          drawColor = VisOpts.hideColor
          canvas.drawLine(child._drawPos[0], child._drawPos[1], child._drawPos[0], node._drawPos[1],
                          drawColor, VisOpts.hideWidth)
        canvas.drawLine(children[0]._drawPos[0], node._drawPos[1], children[-1]._drawPos[0],
                        node._drawPos[1], VisOpts.hideColor, VisOpts.hideWidth)

    if showIndices and (not stopAtCentroids or node._aboveCentroid >= 0):
      txt = str(node.GetIndex())
      if hasattr(node, '_isCentroid'):
        txtColor = piddle.Color(1, .2, .2)
      else:
        txtColor = piddle.Color(0, 0, 0)

      canvas.drawString(txt, node._drawPos[0] - canvas.stringWidth(txt) / 2,
                        node._drawPos[1] + canvas.fontHeight() / 4, color=txtColor)

    if showNodes and hasattr(node, '_isCentroid'):
      rad = VisOpts.nodeRad
      canvas.drawEllipse(node._drawPos[0] - rad / 2, node._drawPos[1] - rad / 2,
                         node._drawPos[0] + rad / 2, node._drawPos[1] + rad / 2, piddle.transparent,
                         fillColor=VisOpts.nodeColor)
      txt = str(node._clustID)
      canvas.drawString(txt, node._drawPos[0] - canvas.stringWidth(txt) / 2,
                        node._drawPos[1] + canvas.fontHeight() / 4, color=piddle.Color(0, 0, 0))

  if showIndices and not stopAtCentroids:
    for pt in pts:
      txt = str(pt.GetIndex())
      canvas.drawString(str(pt.GetIndex()), pt._drawPos[0] - canvas.stringWidth(txt) / 2,
                        pt._drawPos[1])


def ClusterToPDF(cluster, fileName, size=(300, 300), ptColors=[], lineWidth=None, showIndices=0,
                 stopAtCentroids=0, logScale=0):
  """ handles the work of drawing a cluster tree to an PDF file

    **Arguments**

      - cluster: the cluster tree to be drawn

      - fileName: the name of the file to be created

      - size: the size of output canvas

      - ptColors: if this is specified, the _colors_ will be used to color
        the terminal nodes of the cluster tree.  (color == _pid.Color_)

      - lineWidth: if specified, it will be used for the widths of the lines
        used to draw the tree

   **Notes**

     - if _ptColors_ is the wrong length for the number of possible terminal
       node types, this will throw an IndexError

     - terminal node types are determined using their _GetData()_ methods

  """
  try:
    from rdkit.sping.PDF import pidPDF
  except ImportError:
    from rdkit.piddle import piddlePDF
    pidPDF = piddlePDF

  canvas = pidPDF.PDFCanvas(size, fileName)
  if lineWidth is None:
    lineWidth = VisOpts.lineWidth
  DrawClusterTree(cluster, canvas, size, ptColors=ptColors, lineWidth=lineWidth,
                  showIndices=showIndices, stopAtCentroids=stopAtCentroids, logScale=logScale)
  if fileName:
    canvas.save()
  return canvas


def ClusterToSVG(cluster, fileName, size=(300, 300), ptColors=[], lineWidth=None, showIndices=0,
                 stopAtCentroids=0, logScale=0):
  """ handles the work of drawing a cluster tree to an SVG file

    **Arguments**

      - cluster: the cluster tree to be drawn

      - fileName: the name of the file to be created

      - size: the size of output canvas

      - ptColors: if this is specified, the _colors_ will be used to color
        the terminal nodes of the cluster tree.  (color == _pid.Color_)

      - lineWidth: if specified, it will be used for the widths of the lines
        used to draw the tree

   **Notes**

     - if _ptColors_ is the wrong length for the number of possible terminal
       node types, this will throw an IndexError

     - terminal node types are determined using their _GetData()_ methods

  """
  try:
    from rdkit.sping.SVG import pidSVG
  except ImportError:
    from rdkit.piddle.piddleSVG import piddleSVG
    pidSVG = piddleSVG

  canvas = pidSVG.SVGCanvas(size, fileName)

  if lineWidth is None:
    lineWidth = VisOpts.lineWidth
  DrawClusterTree(cluster, canvas, size, ptColors=ptColors, lineWidth=lineWidth,
                  showIndices=showIndices, stopAtCentroids=stopAtCentroids, logScale=logScale)
  if fileName:
    canvas.save()
  return canvas


def ClusterToImg(cluster, fileName, size=(300, 300), ptColors=[], lineWidth=None, showIndices=0,
                 stopAtCentroids=0, logScale=0):
  """ handles the work of drawing a cluster tree to an image file

    **Arguments**

      - cluster: the cluster tree to be drawn

      - fileName: the name of the file to be created

      - size: the size of output canvas

      - ptColors: if this is specified, the _colors_ will be used to color
        the terminal nodes of the cluster tree.  (color == _pid.Color_)

      - lineWidth: if specified, it will be used for the widths of the lines
        used to draw the tree

   **Notes**

     - The extension on  _fileName_ determines the type of image file created.
       All formats supported by PIL can be used.

     - if _ptColors_ is the wrong length for the number of possible terminal
       node types, this will throw an IndexError

     - terminal node types are determined using their _GetData()_ methods

  """
  try:
    from rdkit.sping.PIL import pidPIL
  except ImportError:
    from rdkit.piddle import piddlePIL
    pidPIL = piddlePIL
  canvas = pidPIL.PILCanvas(size, fileName)
  if lineWidth is None:
    lineWidth = VisOpts.lineWidth
  DrawClusterTree(cluster, canvas, size, ptColors=ptColors, lineWidth=lineWidth,
                  showIndices=showIndices, stopAtCentroids=stopAtCentroids, logScale=logScale)
  if fileName:
    canvas.save()
  return canvas
