#
#  Copyright (C) 2010 Gianluca Sforna
#
#  All Rights Reserved
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.

import math


class CanvasBase:
  """Base class for specialized canvas backends"""

  def addCanvasLine(self, p1, p2, color=(0, 0, 0), color2=None, **kwargs):
    """Draw a single line on the canvas

            This function will draw a line between `p1` and `p2` with the
            given `color`.
            If `color2` is specified, it will be used to draw the second half
            of the segment
        """
    raise NotImplementedError('This should be implemented')

  def addCanvasText(self, text, pos, font, color=(0, 0, 0), **kwargs):
    """Draw some text

           The provided `text` is drawn at position `pos` using the given
           `font` and the chosen `color`.
        """
    raise NotImplementedError('This should be implemented')

  def addCanvasPolygon(self, ps, color=(0, 0, 0), **kwargs):
    """Draw a polygon

           Draw a polygon identified by vertexes given in `ps` using
           the given `color`
        """
    raise NotImplementedError('This should be implemented')

  def addCanvasDashedWedge(self, p1, p2, p3, dash=(2, 2), color=(0, 0, 0), color2=None, **kwargs):
    """Draw a dashed wedge

           The wedge is identified by the three points `p1`, `p2`, and `p3`.
           It will be drawn using the given `color`; if `color2` is specified
           it will be used for the second half of the wedge

           TODO: fix comment, I'm not sure what `dash` does

        """
    raise NotImplementedError('This should be implemented')

  def flush(self):
    """Complete any remaining draw operation

           This is supposed to be the last operation on the canvas before
           saving it
        """
    raise NotImplementedError('This should be implemented')

  def _getLinePoints(self, p1, p2, dash):
    x1, y1 = p1
    x2, y2 = p2
    dx = x2 - x1
    dy = y2 - y1
    lineLen = math.sqrt(dx * dx + dy * dy)
    theta = math.atan2(dy, dx)
    cosT = math.cos(theta)
    sinT = math.sin(theta)

    pos = (x1, y1)
    pts = [pos]
    dist = 0
    currDash = 0
    while dist < lineLen:
      currL = dash[currDash % len(dash)]
      if dist + currL > lineLen:
        currL = lineLen - dist
      endP = (pos[0] + currL * cosT, pos[1] + currL * sinT)
      pts.append(endP)
      pos = endP
      dist += currL
      currDash += 1
    return pts
