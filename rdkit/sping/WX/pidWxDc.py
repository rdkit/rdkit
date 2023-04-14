'''\
a piddle wrapper for wxPython DeviceContexts
piddleWxDc.py

By Paul and Kevin Jacobs

History -

   1.0  Many fixes and the rest of required piddle functionality added.
   0.5  Much work done by Jeffrey Kunce on image support and code factoring.

PiddleWxDc adds Piddle-compatible methods to any wxPython DeviceContext.
It can be used any where a wxDC is used (onPaint, onDraw, etc).

Code factoring and pil image support by Jeffrey Kunce

see also piddleWxDcDemo.py
'''

from wxPython.wx import *

from rdkit.sping import pid as sping_pid


class WxCanvasError(RuntimeError):
  pass


class PiddleWxDc(sping_pid.Canvas):

  def __init__(self, aWxDc, size=(300, 300), name="piddleWX"):
    sping_pid.Canvas.__init__(self, size, name)
    self.dc = aWxDc
    self.dc.BeginDrawing()

  def __del__(self):
    self.dc.EndDrawing()

  def _getWXcolor(self, color, default=None):
    '''Converts PIDDLE colors to wx colors'''
    if color is not None:
      if color == sping_pid.transparent:
        return None
      elif color.red >= 0 and color.green >= 0 and color.blue >= 0:
        return wxColour(color.red * 255, color.green * 255, color.blue * 255)

    if default is not None:
      return self._getWXcolor(default)
    else:
      return None  # End of the line

  def _getWXbrush(self, color, default_color=None):
    '''Converts PIDDLE colors to a wx brush'''

    if color == sping_pid.transparent:
      return wxTRANSPARENT_BRUSH

    wxcolor = self._getWXcolor(color)

    if wxcolor is None:
      if default_color is not None:
        return self._getWXbrush(default_color)
      else:
        raise WxCanvasError("Cannot create brush.")

    return wxBrush(wxcolor)

  def _getWXpen(self, width, color, default_color=None):
    '''Converts PIDDLE colors to a wx pen'''

    if width is None or width < 0:
      width = self.defaultLineWidth

    if color == sping_pid.transparent:
      return wxTRANSPARENT_PEN

    wxcolor = self._getWXcolor(color)

    if wxcolor is None:
      if default_color is not None:
        return self._getWXpen(width, default_color)
      else:
        raise WxCanvasError("Cannot create pen.")

    return wxPen(wxcolor, width)

  def _getWXfont(self, font):
    '''Returns a wxFont roughly equivalent to the requested PIDDLE font'''
    if font is None:
      font = self.defaultFont
    #  PIDDLE fonts are matched to wxFont families.  While it is possible to
    #  match them to individual fonts, this is difficult to do in a platform
    #  independent way
    if font.face is None or font.face == 'times':
      family = wxDEFAULT
    elif font.face == 'courier' or font.face == 'monospaced':
      family = wxMODERN
    elif font.face == 'helvetica' or font.face == 'sansserif':
      family = wxSWISS
    elif font.face == 'serif' or font.face == 'symbol':
      family = wxDEFAULT
    else:
      family = wxDEFAULT
    weight = wxNORMAL
    style = wxNORMAL
    underline = 0
    if font.bold == 1:
      weight = wxBOLD
    if font.underline == 1:
      underline = 1
    if font.italic == 1:
      style = wxITALIC
    return wxFont(font.size, family, style, weight, underline)

  def _setWXfont(self, font=None):
    '''set/return the current font for the dc
        jjk  10/28/99'''
    wx_font = self._getWXfont(font)
    self.dc.SetFont(wx_font)
    return (wx_font)

  def isInteractive(self):
    return (0)

  def canUpdate(self):
    return 1

  def clear(self):
    self.dc.Clear()

  #------------ string/font info ------------

  def stringWidth(self, s, font=None):
    '''Return the logical width of the string if it were drawn \
        in the current font (defaults to self.font).'''
    wx_font = self._setWXfont(font)
    return self.dc.GetTextExtent(s)[0]

  def fontHeight(self, font=None):
    '''Find the total height (ascent + descent) of the given font.'''
    return self.fontAscent(font) + self.fontDescent(font)

  def fontAscent(self, font=None):
    '''Find the ascent (height above base) of the given font.'''
    wx_font = self._setWXfont(font)
    return self.dc.GetCharHeight() - self.fontDescent(font)

  def fontDescent(self, font=None):
    '''Find the descent (extent below base) of the given font.'''
    wx_font = self._setWXfont(font)
    extents = self.dc.GetFullTextExtent(' ', wx_font)
    return extents[2]


#------------- drawing methods --------------
# Note default parameters "=None" means use the defaults set in the
# Canvas method: defaultLineColor, etc.

  def drawLine(self, x1, y1, x2, y2, color=None, width=None):
    '''Draw a straight line between x1,y1 and x2,y2.'''

    if width is None or width < 0:
      width = self.defaultLineWidth
    self.dc.SetPen(self._getWXpen(width, color, self.defaultLineColor))
    self.dc.DrawLine(x1, y1, x2, y2)

  def drawString(self, s, x, y, font=None, color=None, angle=None):
    '''Draw a string starting at location x,y.
        NOTE: the baseline goes on y; drawing covers (y-ascent,y+descent)
        Text rotation (angle%360 != 0) is not supported.'''

    self._setWXfont(font)

    if color == sping_pid.transparent:
      return

    # No defaultFontColor?
    wx_color = self._getWXcolor(color, self.defaultLineColor)

    if wx_color is None:
      wx_color = wxBLACK

    self.dc.SetTextForeground(wx_color)

    if '\n' in s or '\r' in s:
      #normalize line ends
      s = s.replace('\r\n', '\n')
      s = s.replace('\n\r', '\n')
      lines = s.split('\n')
    else:
      lines = [s]

    if angle is not None:
      self._drawRotatedString(lines, x, y, font, wx_color, angle)
    else:
      line_height = self.fontHeight(font)
      for l in range(0, len(lines)):
        self.dc.DrawText(lines[l], x, y - self.fontAscent(font) + l * line_height)

  def _drawRotatedString(self, lines, x, y, font=None, color=None, angle=0):

    import math

    # [kbj] Hack since the default system font may not be able to rotate.
    if font is None:
      font = sping_pid.Font(face='helvetica')
      self._setWXfont(font)

    ascent = self.fontAscent(font)
    height = self.fontHeight(font)

    rad = angle * math.pi / 180.
    s = math.sin(rad)
    c = math.cos(rad)
    dx = s * height
    dy = c * height
    lx = x - dx
    ly = y - c * ascent

    for i in range(0, len(lines)):
      self.dc.DrawRotatedText(lines[i], lx + i * dx, ly + i * dy, angle)

    # drawPolygon: For fillable shapes, edgeColor defaults to
    # self.defaultLineColor, edgeWidth defaults to self.defaultLineWidth, and
    # fillColor defaults to self.defaultFillColor.  Specify "don't fill" by
    # passing fillColor=transparent.

  def drawPolygon(self, pointlist, edgeColor=None, edgeWidth=None, fillColor=None, closed=0):
    """drawPolygon(pointlist) -- draws a polygon 
        pointlist: a list of (x,y) tuples defining vertices
        closed:    if 1, adds an extra segment connecting the last point 
            to the first
        """

    # Because wxPython automatically closes polygons, the polygon fill and the border
    # are drawn separately, so open polygons will display correctly
    self.dc.SetPen(wxTRANSPARENT_PEN)
    self.dc.SetBrush(self._getWXbrush(fillColor, self.defaultFillColor))

    #  Workaround : PIDDLE will pass mixed lists of lists and 2-tuples
    #  instead of just 2-tuples.  Therefore, pointlist must be re-created as
    #  only 2-tuples

    pointlist = map(lambda i: tuple(i), pointlist)
    if closed == 1:
      pointlist.append(pointlist[0])

    self.dc.DrawPolygon(pointlist)

    # Create a list of lines (4-tuples) to pass to drawLines
    linelist = []
    if len(pointlist) > 1:
      for i in range(1, len(pointlist)):
        linelist.append(
          (pointlist[i - 1][0], pointlist[i - 1][1], pointlist[i][0], pointlist[i][1]))
      else:
        linelist.append((pointlist[0][0], pointlist[0][1], pointlist[0][0], pointlist[0][1]))
    self.drawLines(linelist, edgeColor, edgeWidth)

  # no colors apply to drawImage; the image is drawn as-is
  def drawImage(self, image, x1, y1, x2=None, y2=None):
    """Draw a PIL Image into the specified rectangle.  If x2 and y2 are
        omitted, they are calculated from the image size.
        jjk  11/03/99"""

    try:
      from PIL import Image
    except ImportError:
      print('PIL not installed as package')
      try:
        import Image
      except ImportError:
        raise RuntimeError("PIL not available!")

    if (x2 and y2 and x2 > x1 and y2 > y1):
      imgPil = image.resize((x2 - x1, y2 - y1))
    else:
      imgPil = image
    if (imgPil.mode != 'RGB'):
      imgPil = imgPil.convert('RGB')
    imgData = imgPil.tobytes()
    imgWx = wxEmptyImage(imgPil.size[0], imgPil.size[1])
    imgWx.SetData(imgData)
    self.dc.DrawBitmap(imgWx.ConvertToBitmap(), x1, y1)
