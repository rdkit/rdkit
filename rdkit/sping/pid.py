# $Id$
# based upon
# piddle.py -- Plug In Drawing, Does Little Else
# Copyright (C) 1999  Joseph J. Strout
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

#  Progress Reports...
#  JJS, 2/10/99: as discussed, I've removed the Shape classes and moved
#       the drawing methods into the Canvas class.  Numerous other changes
#       as discussed by email as well.

#  JJS, 2/11/99: removed Canvas default access functions; added fontHeight
#       etc. functions; fixed numerous typos; added drawRect and drawRoundRect
#       (how could I forget those?).  Added StateSaver utility class.

#       2/11/99 (later): minor fixes.

#  JJS, 2/12/99: removed scaling/sizing references.  Changed event handler
#       mechanism per Magnus's idea.  Changed drawCurve into a fillable
#       drawing function (needs default implementation).  Removed edgeList
#       from drawPolygon.  Added drawFigure.  Changed drawLines to draw
#       a set of disconnected lines (of uniform color and width).

#       2/12/99 (later): added HexColor function and WWW color constants.
#       Fixed bug in StateSaver.  Changed params to drawArc.

#  JJS, 2/17/99: added operator methods to Color; added default implementation
#       of drawRoundRect in terms of Line, Rect, and Arc.

#  JJS, 2/18/99: added isInteractive and canUpdate methods to Canvas.

#  JJS, 2/19/99: added drawImage method; added angle parameter to drawString.

#  JJS, 3/01/99: nailed down drawFigure interface (and added needed constants).

#  JJS, 3/08/99: added arcPoints and curvePoints methods; added default
#       implementations for drawRect, drawRoundRect, drawArc, drawCurve,
#       drawEllipse, and drawFigure (!), mostly thanks to Magnus.

#  JJS, 3/09/99: added 'closed' parameter to drawPolygon, drawCurve, and
#       drawFigure.  Made use of this in several default implementations.

#  JJS, 3/11/99: added 'onKey' callback and associated constants; also added
#       Canvas.setInfoLine(s) method.

#  JJS, 3/12/99: typo in drawFigure.__doc__ corrected (thanks to Magnus).

#  JJS, 3/19/99: fixed bug in drawArc (also thanks to Magnus).

#  JJS, 5/30/99: fixed bug in arcPoints.

#  JJS, 6/10/99: added __repr__ method to Font.

#  JJS, 6/22/99: added additional WWW colors thanks to Rafal Smotrzyk

#  JJS, 6/29/99: added inch and cm units

#  JJS, 6/30/99: added size & name parameters to Canvas.__init__

#  JJS, 9/21/99: fixed bug in arcPoints

#  JJS, 9/29/99: added drawMultiLineStrings, updated fontHeight with new definition

#  JJS, 10/21/99: made Color immutable; fixed bugs in default fontHeight,
#               drawMultiLineString
"""
PIDDLE (Plug-In Drawing, Does Little Else)
2D Plug-In Drawing System

Magnus Lie Hetland
Andy Robinson
Joseph J. Strout
and others

February-March 1999

On coordinates: units are Big Points, approximately 1/72 inch.
The origin is at the top-left, and coordinates increase down (y)
and to the right (x).

"""

__version_maj_number__ = "1.0"  # if release should match "1.0"
__version_min_number__ = "0"  # should match "12"
__version__ = __version_maj_number__ + "." + __version_min_number__  # c.f. "1.0.12"

from rdkit.sping.colors import *

inch = 72  # 1 PIDDLE drawing unit == 1/72 imperial inch
cm = inch / 2.54  # more sensible measurement unit


# -------------------------------------------------------------------------
# StateSaver
# -------------------------------------------------------------------------
class StateSaver:
  """This is a little utility class for saving and restoring the
          default drawing parameters of a canvas.  To use it, add a line
          like this before changing any of the parameters:

                  saver = StateSaver(myCanvas)

          then, when "saver" goes out of scope, it will automagically
          restore the drawing parameters of myCanvas."""

  def __init__(self, canvas):
    self.canvas = canvas
    self.defaultLineColor = canvas.defaultLineColor
    self.defaultFillColor = canvas.defaultFillColor
    self.defaultLineWidth = canvas.defaultLineWidth
    self.defaultFont = canvas.defaultFont

  def __del__(self):
    self.canvas.defaultLineColor = self.defaultLineColor
    self.canvas.defaultFillColor = self.defaultFillColor
    self.canvas.defaultLineWidth = self.defaultLineWidth
    self.canvas.defaultFont = self.defaultFont


# -------------------------------------------------------------------------
# Font
# -------------------------------------------------------------------------
class Font:
  "This class represents font typeface, size, and style."

  def __init__(self, size=12, bold=0, italic=0, underline=0, face=None):
    # public mode variables
    d = self.__dict__
    d["bold"] = bold
    d["italic"] = italic
    d["underline"] = underline

    # public font size (points)
    d["size"] = size

    # typeface -- a name or set of names, interpreted by the Canvas,
    # or "None" to indicate the Canvas-specific default typeface
    d["face"] = face

  def __cmp__(self, other):
    """Compare two fonts to see if they're the same."""
    if self.face == other.face and self.size == other.size and \
       self.bold == other.bold and self.italic == other.italic \
       and self.underline == other.underline:
      return 0
    else:
      return 1

  def __repr__(self):
    return "Font(%d,%d,%d,%d,%s)" % (self.size, self.bold, self.italic, self.underline,
                                     repr(self.face))

  def __setattr__(self, name, value):
    raise TypeError("piddle.Font has read-only attributes")


# -------------------------------------------------------------------------
# constants needed for Canvas.drawFigure
# -------------------------------------------------------------------------
figureLine = 1
figureArc = 2
figureCurve = 3

# -------------------------------------------------------------------------
# key constants used for special keys in the onKey callback
# -------------------------------------------------------------------------
keyBksp = '\010'  # (erases char to left of cursor)
keyDel = '\177'  # (erases char to right of cursor)
keyLeft = '\034'
keyRight = '\035'
keyUp = '\036'
keyDown = '\037'
keyPgUp = '\013'
keyPgDn = '\014'
keyHome = '\001'
keyEnd = '\004'
keyClear = '\033'
keyTab = '\t'

modShift = 1  # shift key was also held
modControl = 2  # control key was also held


# -------------------------------------------------------------------------
# Canvas
# -------------------------------------------------------------------------
class Canvas:
  """This is the base class for a drawing canvas.  The 'plug-in renderers'
          we speak of are really just classes derived from this one, which implement
          the various drawing methods."""

  def __init__(self, size=(300, 300), name="PIDDLE"):
    """Initialize the canvas, and set default drawing parameters.
                    Derived classes should be sure to call this method."""
    # defaults used when drawing
    self.defaultLineColor = black
    self.defaultFillColor = transparent
    self.defaultLineWidth = 1
    self.defaultFont = Font()

    # set up null event handlers

    # onClick: x,y is Canvas coordinates of mouseclick
    def ignoreClick(canvas, x, y):
      pass

    self.onClick = ignoreClick

    # onOver: x,y is Canvas location of mouse
    def ignoreOver(canvas, x, y):
      pass

    self.onOver = ignoreOver

    # onKey: key is printable character or one of the constants above;
    #       modifiers is a tuple containing any of (modShift, modControl)
    def ignoreKey(canvas, key, modifiers):
      pass

    self.onKey = ignoreKey

    # size and name, for user's reference
    self.size, self.name = size, name

  def getSize(self):
    # gL
    return self.size

  # ------------ canvas capabilities -------------
  def isInteractive(self):
    "Returns 1 if onClick, onOver, and onKey events are possible, 0 otherwise."
    return 0

  def canUpdate(self):
    "Returns 1 if the drawing can be meaningfully updated over time \
                    (e.g., screen graphics), 0 otherwise (e.g., drawing to a file)."

    return 0

  # ------------ general management -------------
  def clear(self):
    "Call this to clear and reset the graphics context."
    pass

  def flush(self):
    "Call this to indicate that any comamnds that have been issued \
                    but which might be buffered should be flushed to the screen"

    pass

  def save(self, file=None, format=None):
    """For backends that can be save to a file or sent to a
                    stream, create a valid file out of what's currently been
                    drawn on the canvas.  Trigger any finalization here.
                    Though some backends may allow further drawing after this call,
                    presume that this is not possible for maximum portability

                    file may be either a string or a file object with a write method
                         if left as the default, the canvas's current name will be used

                    format may be used to specify the type of file format to use as
                         well as any corresponding extension to use for the filename
                         This is an optional argument and backends may ignore it if
                         they only produce one file format."""
    pass

  def setInfoLine(self, s):
    "For interactive Canvases, displays the given string in the \
                    'info line' somewhere where the user can probably see it."

    pass

  # ------------ string/font info ------------

  def stringBox(self, s, font=None):
    return self.stringWidth(s, font), self.fontHeight(font)

  def stringWidth(self, s, font=None):
    "Return the logical width of the string if it were drawn \
                    in the current font (defaults to self.font)."

    raise NotImplementedError('stringWidth')

  def fontHeight(self, font=None):
    "Find the height of one line of text (baseline to baseline) of the given font."
    # the following approximation is correct for PostScript fonts,
    # and should be close for most others:
    if not font:
      font = self.defaultFont
    return 1.2 * font.size

  def fontAscent(self, font=None):
    "Find the ascent (height above base) of the given font."
    raise NotImplementedError('fontAscent')

  def fontDescent(self, font=None):
    "Find the descent (extent below base) of the given font."
    raise NotImplementedError('fontDescent')

  # ------------- drawing helpers --------------

  def arcPoints(self, x1, y1, x2, y2, startAng=0, extent=360):
    "Return a list of points approximating the given arc."
    # Note: this implementation is simple and not particularly efficient.
    xScale = abs((x2 - x1) / 2.0)
    yScale = abs((y2 - y1) / 2.0)

    x = min(x1, x2) + xScale
    y = min(y1, y2) + yScale

    # "Guesstimate" a proper number of points for the arc:
    steps = min(max(xScale, yScale) * (extent / 10.0) / 10, 200)
    if steps < 5:
      steps = 5

    from math import cos, pi, sin

    pointlist = []
    step = float(extent) / steps
    angle = startAng
    for i in range(int(steps + 1)):
      point = (x + xScale * cos((angle / 180.0) * pi), y - yScale * sin((angle / 180.0) * pi))
      pointlist.append(point)
      angle = angle + step

    return pointlist

  def curvePoints(self, x1, y1, x2, y2, x3, y3, x4, y4):
    "Return a list of points approximating the given Bezier curve."

    # Adapted from BEZGEN3.HTML, one of the many
    # Bezier utilities found on Don Lancaster's Guru's Lair at
    # <URL: http://www.tinaja.com/cubic01.html>
    bezierSteps = min(
      max(max(x1, x2, x3, x4) - min(x1, x2, x3, x3),
          max(y1, y2, y3, y4) - min(y1, y2, y3, y4)), 200)

    dt1 = 1. / bezierSteps
    dt2 = dt1 * dt1
    dt3 = dt2 * dt1

    xx = x1
    yy = y1
    ux = uy = vx = vy = 0

    ax = x4 - 3 * x3 + 3 * x2 - x1
    ay = y4 - 3 * y3 + 3 * y2 - y1
    bx = 3 * x3 - 6 * x2 + 3 * x1
    by = 3 * y3 - 6 * y2 + 3 * y1
    cx = 3 * x2 - 3 * x1
    cy = 3 * y2 - 3 * y1

    mx1 = ax * dt3
    my1 = ay * dt3

    lx1 = bx * dt2
    ly1 = by * dt2

    kx = mx1 + lx1 + cx * dt1
    ky = my1 + ly1 + cy * dt1

    mx = 6 * mx1
    my = 6 * my1

    lx = mx + 2 * lx1
    ly = my + 2 * ly1

    pointList = [(xx, yy)]

    for i in range(bezierSteps):
      xx = xx + ux + kx
      yy = yy + uy + ky
      ux = ux + vx + lx
      uy = uy + vy + ly
      vx = vx + mx
      vy = vy + my
      pointList.append((xx, yy))

    return pointList

  def drawMultiLineString(self, s, x, y, font=None, color=None, angle=0, **kwargs):
    "Breaks string into lines (on \n, \r, \n\r, or \r\n), and calls drawString on each."
    import math
    h = self.fontHeight(font)
    dy = h * math.cos(angle * math.pi / 180.0)
    dx = h * math.sin(angle * math.pi / 180.0)
    s = s.replace('\r\n', '\n')
    s = s.replace('\n\r', '\n')
    s = s.replace('\r', '\n')
    lines = s.split('\n')
    for line in lines:
      self.drawString(line, x, y, font, color, angle)
      x = x + dx
      y = y + dy

  # ------------- drawing methods --------------

  # Note default parameters "=None" means use the defaults
  # set in the Canvas method: defaultLineColor, etc.

  def drawLine(self, x1, y1, x2, y2, color=None, width=None, dash=None, **kwargs):
    "Draw a straight line between x1,y1 and x2,y2."
    raise NotImplementedError('drawLine')

  def drawLines(self, lineList, color=None, width=None, dash=None, **kwargs):
    "Draw a set of lines of uniform color and width.  \
                    lineList: a list of (x1,y1,x2,y2) line coordinates."

    # default implementation:
    for x1, y1, x2, y2 in lineList:
      self.drawLine(x1, y1, x2, y2, color, width, dash=dash, **kwargs)

    #       For text, color defaults to self.lineColor.

  def drawString(self, s, x, y, font=None, color=None, angle=0, **kwargs):
    "Draw a string starting at location x,y."
    # NOTE: the baseline goes on y; drawing covers (y-ascent,y+descent)
    raise NotImplementedError('drawString')

  #       For fillable shapes, edgeColor defaults to self.defaultLineColor,
  #       edgeWidth defaults to self.defaultLineWidth, and
  #       fillColor defaults to self.defaultFillColor.
  #       Specify "don't fill" by passing fillColor=transparent.

  def drawCurve(self, x1, y1, x2, y2, x3, y3, x4, y4, edgeColor=None, edgeWidth=None,
                fillColor=None, closed=0, dash=None, **kwargs):
    "Draw a Bezier curve with control points x1,y1 to x4,y4."
    pointlist = self.curvePoints(x1, y1, x2, y2, x3, y3, x4, y4)
    self.drawPolygon(pointlist, edgeColor=edgeColor, edgeWidth=edgeWidth, fillColor=fillColor,
                     closed=closed, dash=dash, **kwargs)

  def drawRect(self, x1, y1, x2, y2, edgeColor=None, edgeWidth=None, fillColor=None, dash=None,
               **kwargs):
    "Draw the rectangle between x1,y1, and x2,y2. \
                    These should have x1<x2 and y1<y2."

    pointList = [(x1, y1), (x2, y1), (x2, y2), (x1, y2)]
    self.drawPolygon(pointList, edgeColor, edgeWidth, fillColor, closed=1, dash=dash, **kwargs)

  def drawRoundRect(self, x1, y1, x2, y2, rx=8, ry=8, edgeColor=None, edgeWidth=None,
                    fillColor=None, dash=None, **kwargs):
    "Draw a rounded rectangle between x1,y1, and x2,y2, \
                    with corners inset as ellipses with x radius rx and y radius ry. \
                    These should have x1<x2, y1<y2, rx>0, and ry>0."

    x1, x2 = min(x1, x2), max(x1, x2)
    y1, y2 = min(y1, y2), max(y1, y2)

    dx = rx * 2
    dy = ry * 2

    partList = [(figureArc, x1, y1, x1 + dx, y1 + dy, 180, -90),
                (figureLine, x1 + rx, y1, x2 - rx, y1),
                (figureArc, x2 - dx, y1, x2, y1 + dy, 90, -90),
                (figureLine, x2, y1 + ry, x2, y2 - ry),
                (figureArc, x2 - dx, y2, x2, y2 - dy, 0, -90),
                (figureLine, x2 - rx, y2, x1 + rx, y2),
                (figureArc, x1, y2, x1 + dx, y2 - dy, -90, -90),
                (figureLine, x1, y2 - ry, x1, y1 + rx)]

    self.drawFigure(partList, edgeColor, edgeWidth, fillColor, closed=1, dash=dash, **kwargs)

  def drawEllipse(self, x1, y1, x2, y2, edgeColor=None, edgeWidth=None, fillColor=None, dash=None,
                  **kwargs):
    "Draw an orthogonal ellipse inscribed within the rectangle x1,y1,x2,y2. \
                    These should have x1<x2 and y1<y2."

    pointlist = self.arcPoints(x1, y1, x2, y2, 0, 360)
    self.drawPolygon(pointlist, edgeColor, edgeWidth, fillColor, closed=1, dash=dash, **kwargs)

  def drawArc(self, x1, y1, x2, y2, startAng=0, extent=360, edgeColor=None, edgeWidth=None,
              fillColor=None, dash=None, **kwargs):
    "Draw a partial ellipse inscribed within the rectangle x1,y1,x2,y2, \
                    starting at startAng degrees and covering extent degrees.   Angles \
                    start with 0 to the right (+x) and increase counter-clockwise. \
                    These should have x1<x2 and y1<y2."

    center = (x1 + x2) / 2, (y1 + y2) / 2
    pointlist = self.arcPoints(x1, y1, x2, y2, startAng, extent)

    # Fill...
    self.drawPolygon(pointlist + [center] + [pointlist[0]], transparent, 0, fillColor)

    # Outline...
    self.drawPolygon(pointlist, edgeColor, edgeWidth, transparent, dash=dash, **kwargs)

  def drawPolygon(self, pointlist, edgeColor=None, edgeWidth=None, fillColor=None, closed=0,
                  dash=None, **kwargs):
    """drawPolygon(pointlist) -- draws a polygon
                    pointlist: a list of (x,y) tuples defining vertices
                    closed: if 1, adds an extra segment connecting the last point to the first
                    """
    raise NotImplementedError('drawPolygon')

  def drawFigure(self, partList, edgeColor=None, edgeWidth=None, fillColor=None, closed=0,
                 dash=None, **kwargs):
    """drawFigure(partList) -- draws a complex figure
                    partlist: a set of lines, curves, and arcs defined by a tuple whose
                                      first element is one of figureLine, figureArc, figureCurve
                                      and whose remaining 4, 6, or 8 elements are parameters."""

    pointList = []

    for tuple in partList:
      op = tuple[0]
      args = list(tuple[1:])

      if op == figureLine:
        pointList.extend([args[:2], args[2:]])
      elif op == figureArc:
        pointList.extend(self.arcPoints(*args))
      elif op == figureCurve:
        pointList.extend(self.curvePoints(*args))
      else:
        raise TypeError("unknown figure operator: " + op)

      self.drawPolygon(pointList, edgeColor, edgeWidth, fillColor, closed=closed, dash=dash,
                       **kwargs)

  # no colors apply to drawImage; the image is drawn as-is

  def drawImage(self, image, x1, y1, x2=None, y2=None, **kwargs):
    """Draw a PIL Image into the specified rectangle.  If x2 and y2 are
                    omitted, they are calculated from the image size."""
    raise NotImplementedError('drawImage')


# -------------------------------------------------------------------------

# utility functions #


def getFileObject(file, openFlags="wb"):
  """Common code for every Canvas.save() operation takes a string
          or a potential file object and assures that a valid fileobj is returned"""

  if file:
    if isinstance(file, str):
      fileobj = open(file, openFlags)
    else:
      if hasattr(file, "write"):
        fileobj = file
      else:
        raise ValueError('Invalid file argument to save')
  else:
    raise ValueError('Invalid file argument to save')

  return fileobj


class AffineMatrix:
  # A = [ a c e]
  #     [ b d f]
  #     [ 0 0 1]
  # self.A = [a b c d e f] = " [ A[0] A[1] A[2] A[3] A[4] A[5] ]"
  def __init__(self, init=None):
    if init:
      if len(init) == 6:
        self.A = init
      if type(init) == type(self):  # erpht!!! this seems so wrong
        self.A = init.A
    else:
      self.A = [1.0, 0, 0, 1.0, 0.0, 0.0]  # set to identity

  def scale(self, sx, sy):
    self.A = [sx * self.A[0], sx * self.A[1], sy * self.A[2], sy * self.A[3], self.A[4], self.A[5]]

  def rotate(self, theta):
    "counter clockwise rotation in standard SVG/libart coordinate system"
    # clockwise in postscript "y-upward" coordinate system
    # R = [ c  -s  0 ]
    #     [ s   c  0 ]
    #     [ 0   0  1 ]

    co = math.cos(PI * theta / 180.0)
    si = math.sin(PI * theta / 180.0)
    self.A = [
      self.A[0] * co + self.A[2] * si, self.A[1] * co + self.A[3] * si,
      -self.A[0] * si + self.A[2] * co, -self.A[1] * si + self.A[3] * co, self.A[4], self.A[5]
    ]

  def translate(self, tx, ty):
    self.A = [
      self.A[0], self.A[1], self.A[2], self.A[3], self.A[0] * tx + self.A[2] * ty + self.A[4],
      self.A[1] * tx + self.A[3] * ty + self.A[5]
    ]
