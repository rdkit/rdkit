# piddlePIL.py -- a Python Imaging Library backend for PIDDLE
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
"""pidPIL

This module implements a Python Imaging Library PIDDLE canvas.
In other words, this is a PIDDLE backend that renders into a
PIL Image object.  From there, you can save as GIF, plot into
another PIDDLE canvas, etc.

Joe Strout (joe@strout.net), 10/26/99
modified for use with sping.
This requires Imaging to be installed as a package PIL
"""

# 6/22/99: updated drawString to handle non-integer x and y

from rdkit.sping.pid import *

try:
  import Image
  import ImageDraw
  import ImageFont
except ImportError:
  from PIL import Image, ImageFont, ImageDraw

import math
import os

if __name__ == '__main__':
  _fontprefix = os.path.join(os.curdir, 'pilfonts')
else:
  _fontprefix = os.path.join(os.path.split(__file__)[0], 'pilfonts')

# load font metrics
try:
  f = open(os.path.join(_fontprefix, 'metrics.dat'), 'rb')
  import pickle
  _widthmaps = pickle.load(f)
  _ascents = pickle.load(f)
  _descents = pickle.load(f)
  f.close()
except Exception:
  print("Warning: unable to load font metrics from dir {0}".format(_fontprefix))
  _widthmaps = {}
  _ascents = {}
  _descents = {}
# finally:
#	pass	# (just here so we can comment out the except clause for debugging)


def _closestSize(size):
  supported = [8, 10, 12, 14, 18, 24]  # list of supported sizes
  if size in supported:
    return size
  best = supported[0]
  bestdist = abs(size - best)
  for trial in supported[1:]:
    dist = abs(size - trial)
    if dist < bestdist:
      best = trial
      bestdist = dist
  return best


def _pilFontPath(face, size, bold=0):
  if face == 'monospaced':
    face = 'courier'
  elif face == 'serif':
    face = 'times'
  elif face == 'sansserif' or face == 'system':
    face = 'helvetica'

  if bold and face != 'symbol':
    fname = "%s-bold-%d.pil" % (face, size)
  else:
    fname = "%s-%d.pil" % (face, size)
  path = os.path.join(_fontprefix, fname)
  return path


def _matchingFontPath(font):
  # returns a font path which matches info in our font metrics
  if font.face:
    face = font.face
  else:
    face = 'times'

  size = _closestSize(font.size)
  if isinstance(face, str):
    path = _pilFontPath(face, size, font.bold)
    path = path.split(os.sep)[-1]
    if path in _widthmaps.keys():
      return path
  else:
    for item in font.face:
      path = _pilFontPath(item, size, font.bold)
      path = path.split(os.sep)[-1]
      if path in _widthmaps.keys():
        return path
  # not found?  Try it with courier, which should always be there
  path = _pilFontPath('courier', size, font.bold)
  return path.split(os.sep)[-1]


def _pilFont(font):
  if font.face:
    face = font.face
  else:
    face = 'times'

  size = _closestSize(font.size)
  if isinstance(face, str):
    try:
      pilfont = ImageFont.load_path(_pilFontPath(face, size, font.bold))
    except Exception:
      return 0  # font not found!
  else:
    for item in font.face:
      pilfont = None
      try:
        pilfont = ImageFont.load_path(_pilFontPath(item, size, font.bold))
        break
      except Exception:
        pass
    if pilfont is None:
      return 0  # font not found!
  return pilfont


class PILCanvas(Canvas):

  def __init__(self, size=(300, 300), name='piddlePIL'):
    self._image = Image.new('RGB', (int(size[0]), int(size[1])), (255, 255, 255))
    self._pen = ImageDraw.ImageDraw(self._image)
    self._setFont(Font())
    Canvas.__init__(self, size, name)

  def __setattr__(self, attribute, value):
    self.__dict__[attribute] = value
    if attribute == "defaultLineColor":
      self._setColor(self.defaultLineColor)

  # utility functions
  def _setColor(self, c):
    "Set the pen color from a piddle color."
    self._color = (int(c.red * 255), int(c.green * 255), int(c.blue * 255))

  def _setFont(self, font):
    self._font = _pilFont(font)

  # public functions

  def getImage(self):
    return self._image

  def save(self, file=None, format=None):
    """format may be a string specifying a file extension corresponding to
        an image file format. Ex: 'png', 'jpeg', 'gif', 'tif' etc.
        These are defined by PIL, not by us so you need to check the docs.
        In general, I just specify an extension and let format default to None"""
    file = file or self.name
    if hasattr(file, 'write'):
      raise ValueError('fileobj not implemented for piddlePIL')
    # below here, file is guaranteed to be a string
    if format is None:
      if '.' not in file:
        filename = file + '.png'  # default to producing jpg
      else:
        filename = file
        # format = file[-3:] # os.path.splitext(..)
    else:
      filename = file + '.' + format

    self._image.save(filename, format=format)

  def clear(self):
    # why is edgeColor yellow ???
    self.drawRect(0, 0, self.size[0], self.size[1], edgeColor=yellow, fillColor=white)
    # FIXME: need to reset canvas as well to defaults ???

  # ------------ string/font info ------------

  def stringWidth(self, s, font=None):
    "Return the logical width of the string if it were drawn \
        in the current font (defaults to self.defaultFont)."

    if not font:
      font = self.defaultFont
    if not _widthmaps:
      return font.size * len(s)

    path = _matchingFontPath(font)
    map = _widthmaps[path]
    out = 0
    for c in s:
      out += map.get(c, font.size)
    return out

  def fontAscent(self, font=None):
    "Find the ascent (height above base) of the given font."

    if not font:
      font = self.defaultFont
    if not _ascents:
      return font.size

    path = _matchingFontPath(font)
    return _ascents[path]

  def fontDescent(self, font=None):
    "Find the descent (extent below base) of the given font."

    if not font:
      font = self.defaultFont
    if not _descents:
      return font.size / 2

    path = _matchingFontPath(font)
    return _descents[path]

  # ------------- drawing methods --------------
  def drawLine(self, x1, y1, x2, y2, color=None, width=None, dash=None, **kwargs):
    "Draw a straight line between x1,y1 and x2,y2."
    # set color...
    if width is None:
      w = self.defaultLineWidth
    elif width:
      w = width
    else:
      return

    if color:
      if color == transparent:
        return
      self._setColor(color)
    elif self.defaultLineColor == transparent:
      return

    if not dash:
      self._pen.line((x1, y1, x2, y2), fill=self._color, width=w)
    else:
      dx = x2 - x1
      dy = y2 - y1
      lineLen = math.sqrt(dx * dx + dy * dy)
      theta = math.atan2(dy, dx)
      cosT = math.cos(theta)
      sinT = math.sin(theta)

      pos = (x1, y1)
      dist = 0
      currDash = 0
      dashOn = 1
      while dist < lineLen:
        currL = dash[currDash % len(dash)]
        if (dist + currL > lineLen):
          currL = lineLen - dist
        endP = (pos[0] + currL * cosT, pos[1] + currL * sinT)
        if dashOn:
          self.drawLine(pos[0], pos[1], endP[0], endP[1], color=color, width=width, dash=None,
                        **kwargs)
        pos = endP
        dist += currL
        currDash += 1
        dashOn = not dashOn

  def drawPolygon(self, pointlist, edgeColor=None, edgeWidth=None, fillColor=None, closed=0,
                  dash=None, **kwargs):
    """drawPolygon(pointlist) -- draws a polygon
        pointlist: a list of (x,y) tuples defining vertices
        """
    # PIL's routine requires a sequence of tuples...
    # the input is not so restricted, so fix it
    pts = list(pointlist)
    for i in range(len(pts)):
      pts[i] = tuple(pts[i])

    # set color for fill...
    filling = 0
    if fillColor:
      if fillColor != transparent:
        self._setColor(fillColor)
        filling = 1
    elif self.defaultFillColor != transparent:
      self._setColor(self.defaultFillColor)
      filling = 1

    # do the fill
    if filling:
      pts = [(int(x[0]), int(x[1])) for x in pts]
      self._pen.polygon(pts, fill=self._color)

    # set edge width...
    if edgeWidth is None:
      edgeWidth = self.defaultLineWidth
    elif not edgeWidth:
      return

    # set color for edge...
    if edgeColor:
      self._setColor(edgeColor)
    else:
      self._setColor(self.defaultLineColor)

    # draw the outline
    if (closed or (pts[0][0] == pts[-1][0] and pts[0][1] == pts[-1][1])) \
            and edgeWidth <= 1:
      self._pen.polygon(pts, outline=self._color)
    else:
      # ...since PIL's polygon routine insists on closing,
      # and does not support thick edges, we'll use our drawLine instead
      # OFI: use default color/width to speed this up!
      oldp = pts[0]
      if closed:
        pts.append(oldp)
      for p in pts[1:]:
        self.drawLine(oldp[0], oldp[1], p[0], p[1], edgeColor, edgeWidth, dash=dash, **kwargs)
        oldp = p

  def drawString(self, s, x, y, font=None, color=None, angle=0, **kwargs):
    "Draw a string starting at location x,y."
    x = int(x)
    y = int(y)
    if '\n' in s or '\r' in s:
      self.drawMultiLineString(s, x, y, font, color, angle, **kwargs)
      return
    if not font:
      font = self.defaultFont

    if not color:
      color = self.defaultLineColor
    if color == transparent:
      return

    # draw into an offscreen Image
    # tmpsize was originally 1.2* stringWidth, added code to give enough room for single character strings (piddle bug#121995)
    sHeight = (self.fontAscent(font) + self.fontDescent(font))
    sWidth = self.stringWidth(s, font)
    tempsize = max(sWidth * 1.2, sHeight * 2.0)
    tempimg = Image.new('RGB', (int(tempsize), int(tempsize)), (0, 0, 0))

    temppen = ImageDraw.ImageDraw(tempimg)

    pilfont = _pilFont(font)
    if not pilfont:
      raise ValueError("bad font: %s" % font)
    pos = [4, int(tempsize / 2 - self.fontAscent(font)) - self.fontDescent(font)]
    temppen.text(pos, s, font=pilfont, fill=(255, 255, 255))
    pos[1] = int(tempsize / 2)

    if font.underline:
      ydown = (0.5 * self.fontDescent(font))
      # thickness = 0.08 * font.size # may need to ceil this
      temppen.line((pos[0], pos[1] + ydown, pos[0] + sWidth, pos[1] + ydown))

    # rotate
    if angle:
      from math import cos, pi, sin
      tempimg = tempimg.rotate(angle, Image.BILINEAR)
      temppen = ImageDraw.ImageDraw(tempimg)
      radians = -angle * pi / 180.0
      r = tempsize / 2 - pos[0]
      pos[0] = int(tempsize / 2 - r * cos(radians))
      pos[1] = int(pos[1] - r * sin(radians))

    # temppen.rectangle( (pos[0],pos[1],pos[0]+2,pos[1]+2) ) # PATCH for debugging
    # colorize, and copy it in
    mask = tempimg.convert('L').point(lambda c: c)
    clr = (int(color.red * 255), int(color.green * 255), int(color.blue * 255))
    temppen.rectangle((0, 0, tempsize, tempsize), fill=clr)
    self._image.paste(tempimg, (int(x) - pos[0], int(y) - pos[1]), mask)

  def drawImage(self, image, x1, y1, x2=None, y2=None, **kwargs):
    """Draw a PIL Image into the specified rectangle.  If x2 and y2 are
        omitted, they are calculated from the image size."""

    if x2 and y2:
      bbox = image.getbbox()
      if x2 - x1 != bbox[2] - bbox[0] or y2 - y1 != bbox[3] - bbox[1]:
        image = image.resize((x2 - x1, y2 - y1))
    self._image.paste(image, (x1, y1))


def test():
  # ... for testing...
  canvas = PILCanvas()

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

  canvas.drawString("This is a test!", 30, 130, Font(face="times", size=16, bold=1), color=green,
                    angle=-45)

  canvas.drawString("This is a test!", 30, 130, color=red, angle=45)

  polypoints = [(160, 120), (130, 190), (210, 145), (110, 145), (190, 190)]
  canvas.drawPolygon(polypoints, fillColor=lime, edgeColor=red, edgeWidth=3, closed=1)

  canvas.drawRect(200, 200, 260, 260, edgeColor=yellow, edgeWidth=5)
  canvas.drawLine(200, 260, 260, 260, color=green, width=5)
  canvas.drawLine(260, 200, 260, 260, color=red, width=5)

  # now, for testing, save the image as a PNG file
  canvas.flush()
  canvas.getImage().save("test.png")

  return canvas


def testit(canvas, s, x, y, font=None):
  canvas.defaultLineColor = black
  canvas.drawString(s, x, y, font=font)
  canvas.defaultLineColor = blue
  w = canvas.stringWidth(s, font=font)
  canvas.drawLine(x, y, x + w, y)
  canvas.drawLine(x, y - canvas.fontAscent(font=font), x + w, y - canvas.fontAscent(font=font))
  canvas.drawLine(x, y + canvas.fontDescent(font=font), x + w, y + canvas.fontDescent(font=font))


def test2():

  canvas = PILCanvas()
  testit(canvas, "Foogar", 20, 30)

  testit(canvas, "Foogar", 20, 90, font=Font(size=24))
  global dammit
  dammit = _pilFont(Font(size=24))

  testit(canvas, "Foogar", 20, 150, font=Font(face='courier', size=24))

  testit(canvas, "Foogar", 20, 240, font=Font(face='courier'))

  import piddleQD
  global qdcanvas
  try:
    qdcanvas.close()
  except Exception:
    pass
  qdcanvas = piddleQD.QDCanvas()
  qdcanvas.drawImage(canvas.getImage(), 0, 0)


if __name__ == '__main__':
  test()
