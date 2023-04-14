"""
A Tkinter based backend for piddle.

Perry A. Stoll

Created: February 15, 1999

Requires PIL for rotated string support.

Known Problems:
     - Doesn't handle the interactive commands yet.
     - PIL based canvas inherits lack of underlining strings from piddlePIL

You can find the latest version of this file:
    via http://piddle.sourceforge.net
"""
# we depend on PIL for rotated strings so watch for changes in PIL

import tkFont
import Tkinter

tk = Tkinter
import rdkit.sping.pid

__version__ = "0.3"
__date__ = "April 8, 1999"
__author__ = "Perry Stoll,  perry.stoll@mail.com "

# fixups by chris lee, cwlee@artsci.wustl.edu
# $Id$
# - added drawImage scaling support
# - shifted baseline y parameter in drawString to work around font metric
#         shift due to Tkinter's Canvas text_item object
# - fixed argument names so that argument keywords agreed with piddle.py (passes discipline.py)
#
#
# ToDo: for TKCanvas
#         make sure that fontHeight() is returnng appropriate measure.  Where is this info?
#
# $Log: pidTK.py,v $
# Revision 1.1  2002/07/12 18:34:47  glandrum
# added
#
# Revision 1.6  2000/11/03 00:56:57  clee
# fixed sizing error in TKCanvas
#
# Revision 1.5  2000/11/03 00:25:37  clee
# removed reference to "BaseTKCanvas" (should just use TKCanvas as default)
#
# Revision 1.4  2000/10/29 19:35:31  clee
# eliminated BaseTKCanvas in favor of straightforward "TKCanvas" name
#
# Revision 1.3  2000/10/29 01:57:41  clee
# - added scrollbar support to both TKCanvas and TKCanvasPIL
# - added getTKCanvas() access method to TKCanvasPIL
#
# Revision 1.2  2000/10/15 00:47:17  clee
# commit before continuing after getting pil to work as package
#
# Revision 1.1.1.1  2000/09/27 03:53:15  clee
# Simple Platform Independent Graphics
#
# Revision 1.6  2000/04/06 01:55:34  pmagwene
# - TKCanvas now uses multiple inheritance from Tkinter.Canvas and piddle.Canvas
# 	* for the most part works much like a normal Tkinter.Canvas object
# - TKCanvas draws rotated strings using PIL image, other objects using normal Tk calls
# - Minor fixes to FontManager and TKCanvas so can specify root window other than Tk()
# - Removed Quit/Clear buttons from default canvas
#
# Revision 1.5  2000/03/12 07:07:42  clee
# sync with 1_x
#
# Revision 1.4  2000/02/26 23:12:42  clee
# turn off compression by default on piddlePDF
# add doc string to new pil-based piddleTK
#
# Revision 1.3  2000/02/26 21:23:19  clee
# update that makes PIL based TKCanvas the default Canvas for TK.
# Updated piddletest.py.  Also, added clear() methdo to piddlePIL's
# canvas it clears to "white" is this correct behavior?  Not well
# specified in current documents.
#


class FontManager:
  __alt_faces = {"serif": "Times", "sansserif": "Helvetica", "monospaced": "Courier"}

  def __init__(self, master):
    self.master = master
    self.font_cache = {}

  # the main interface
  def stringWidth(self, s, font):
    tkfont = self.piddleToTkFont(font)
    return tkfont.measure(s)

  def fontHeight(self, font):
    tkfont = self.piddleToTkFont(font)
    return self._tkfontHeight(tkfont)

  def fontAscent(self, font):
    tkfont = self.piddleToTkFont(font)
    return self._tkfontAscent(tkfont)

  def fontDescent(self, font):
    tkfont = self.piddleToTkFont(font)
    return self._tkfontDescent(tkfont)

  def getTkFontString(self, font):
    """Return a string suitable to pass as the -font option to
        to a Tk widget based on the piddle-style FONT"""
    tkfont = self.piddleToTkFont(font)
    # XXX: should just return the internal tk font name?
    # return str(tkfont)
    return ('-family %(family)s -size %(size)s '
            '-weight %(weight)s -slant %(slant)s '
            '-underline %(underline)s' % tkfont.config())

  def getTkFontName(self, font):
    """Return a name associated with the piddle-style FONT"""
    tkfont = self.piddleToTkFont(font)
    return str(tkfont)

  def piddleToTkFont(self, font):
    """Return a tkFont instance based on the pid-style FONT"""
    if font is None:
      return ''

    #default 12 pt, "Times", non-bold, non-italic
    size = 12
    family = "Times"
    weight = "normal"
    slant = "roman"
    underline = "false"

    if font.face:
      # check if the user specified a generic face type
      # like serif or monospaced. check is case-insensitive.
      f = font.face.lower()
      if f in self.__alt_faces:
        family = self.__alt_faces[f]
      else:
        family = font.face

    size = font.size or 12

    if font.bold:
      weight = "bold"
    if font.italic:
      slant = "italic"
    if font.underline:
      underline = 'true'

    # ugh... is there a better way to do this?
    key = (family, size, weight, slant, underline)

    # check if we've already seen this font.
    if key in self.font_cache:
      # yep, don't bother creating a new one. just fetch it.
      font = self.font_cache[key]
    else:
      # nope, let's create a new tk font.
      # this way we will return info about the actual font
      # selected by Tk, which may be different than what we ask
      # for if it's not available.
      font = tkFont.Font(self.master, family=family, size=size, weight=weight, slant=slant,
                         underline=underline)
      self.font_cache[(family, size, weight, slant, underline)] = font

    return font

  def _tkfontAscent(self, tkfont):
    return tkfont.metrics("ascent")

  def _tkfontDescent(self, tkfont):
    return tkfont.metrics("descent")


class TKCanvas(tk.Canvas, rdkit.sping.pid.Canvas):

  __TRANSPARENT = ''  # transparent for Tk color

  def __init__(
      self,
      size=(300, 300),
      name="sping.TK",
      master=None,
      scrollingViewPortSize=None,  # a 2-tuple to define the size of the viewport
      **kw):
    """This canvas allows you to add a tk.Canvas with a sping API for drawing.
        To add scrollbars, the simpliest method is to set the 'scrollingViewPortSize'
        equal to a tuple that describes the width and height of the visible porition
        of the canvas on screen.  This sets scrollregion=(0,0, size[0], size[1]).
        Then you can add scrollbars as you would any tk.Canvas.

        Note, because this is a subclass of tk.Canvas, you can use the normal keywords
        to specify a tk.Canvas with scrollbars, however, you should then be careful to
        set the "scrollregion" option to the same size as the 'size' passed to __init__.
        Tkinter's scrollregion option essentially makes 'size' ignored.  """

    rdkit.sping.pid.Canvas.__init__(self, size=size, name=size)
    if scrollingViewPortSize:  # turn on ability to scroll
      kw["scrollregion"] = (0, 0, size[0], size[1])
      kw["height"] = scrollingViewPortSize[0]
      kw["width"] = scrollingViewPortSize[1]
    else:
      kw["width"] = size[0]
      kw["height"] = size[1]

    tk.Canvas.__init__(self, master, **kw)
    self.config(background="white")
    self.width, self.height = size

    self._font_manager = FontManager(self)
    self._configure()
    self._item_ids = []
    self._images = []

  def _configure(self):
    pass

  def _display(self):
    self.flush()
    self.mainloop()

  def _quit(self):
    self.quit()

  # Hmmm...the postscript generated by this causes my Ghostscript to barf...
  def _to_ps_file(self, filename):
    self.postscript(file=filename)

  def isInteractive(self):
    return 0

  def onOver(self, event):
    pass

  def onClick(self, event):
    pass

  def onKey(self, event):
    pass

  def flush(self):
    tk.Canvas.update(self)

  def clear(self):
    map(self.delete, self._item_ids)
    self._item_ids = []

  def _colorToTkColor(self, c):
    return "#%02X%02X%02X" % (int(c.red * 255), int(c.green * 255), int(c.blue * 255))

  def _getTkColor(self, color, defaultColor):
    if color is None:
      color = defaultColor
    if color is rdkit.sping.pid.transparent:
      color = self.__TRANSPARENT
    else:
      color = self._colorToTkColor(color)
    return color

  def drawLine(self, x1, y1, x2, y2, color=None, width=None):
    color = self._getTkColor(color, self.defaultLineColor)
    if width is None:
      width = self.defaultLineWidth
    new_item = self.create_line(x1, y1, x2, y2, fill=color, width=width)
    self._item_ids.append(new_item)

# NYI: curve with fill
#def drawCurve(self, x1, y1, x2, y2, x3, y3, x4, y4,
#              edgeColor=None, edgeWidth=None, fillColor=None, closed=0):
#

  def stringWidth(self, s, font=None):
    return self._font_manager.stringWidth(s, font or self.defaultFont)

  def fontAscent(self, font=None):
    return self._font_manager.fontAscent(font or self.defaultFont)

  def fontDescent(self, font=None):
    return self._font_manager.fontDescent(font or self.defaultFont)

  def drawString(self, s, x, y, font=None, color=None, angle=None):
    if angle:
      try:
        self._drawRotatedString(s, x, y, font, color, angle)
        return
      except ImportError:
        print("PIL not available. Using unrotated strings.")
    # fudge factor for TK on linux (at least)
    # strings are being drawn using create_text in canvas
    y = y - self.fontHeight(font) * .28  # empirical
    #y = y - self.fontDescent(font)

    color = self._getTkColor(color, self.defaultLineColor)
    font = self._font_manager.getTkFontString(font or self.defaultFont)
    new_item = self.create_text(x, y, text=s, font=font, fill=color, anchor=Tkinter.W)
    self._item_ids.append(new_item)

  def _drawRotatedString(self, s, x, y, font=None, color=None, angle=0):
    # we depend on PIL for rotated strings so watch for changes in PIL
    try:
      from PIL import Image, ImageTk

      import rdkit.sping.PIL.pidPIL
      pp = rdkit.sping.PIL.pidPIL
    except ImportError:
      raise ImportError("Rotated strings only possible with PIL support")

    pilCan = pp.PILCanvas(size=(self.width, self.height))
    pilCan.defaultFont = self.defaultFont
    pilCan.defaultLineColor = self.defaultLineColor

    if '\n' in s or '\r' in s:
      self.drawMultiLineString(s, x, y, font, color, angle)
      return
    if not font:
      font = pilCan.defaultFont

    if not color:
      color = self.defaultLineColor
    if color == rdkit.sping.pid.transparent:
      return

    # draw into an offscreen Image
    tempsize = pilCan.stringWidth(s, font) * 1.2
    tempimg = Image.new('RGB', (tempsize, tempsize), (0, 0, 0))
    txtimg = Image.new('RGB', (tempsize, tempsize), (255, 255, 255))

    from PIL import ImageDraw
    temppen = ImageDraw.ImageDraw(tempimg)
    temppen.setink((255, 255, 255))
    pilfont = pp._pilFont(font)

    if not pilfont:
      raise ValueError("Bad font: %s" % font)

    temppen.setfont(pilfont)
    pos = [4, int(tempsize / 2 - pilCan.fontAscent(font)) - pilCan.fontDescent(font)]
    temppen.text(pos, s)
    pos[1] = int(tempsize / 2)

    # rotate
    if angle:
      from math import cos, pi, sin
      tempimg = tempimg.rotate(angle, Image.BILINEAR)
      temppen = ImageDraw.ImageDraw(tempimg)
      radians = -angle * pi / 180.0
      r = tempsize / 2 - pos[0]
      pos[0] = int(tempsize / 2 - r * cos(radians))
      pos[1] = int(pos[1] - r * sin(radians))

    ###temppen.rectangle( (pos[0],pos[1],pos[0]+2,pos[1]+2) ) # PATCH for debugging
    # colorize, and copy it in
    mask = tempimg.convert('L').point(lambda c: c)
    temppen.setink((color.red * 255, color.green * 255, color.blue * 255))
    temppen.setfill(1)
    temppen.rectangle((0, 0, tempsize, tempsize))

    txtimg.paste(tempimg, (0, 0), mask)

    ##Based on code posted by John Michelson in the PIL SIG
    transp = txtimg.convert("RGBA")
    source = transp.split()
    R, G, B, A = 0, 1, 2, 3
    mask = transp.point(lambda i: i < 255 and 255)  # use white as transparent
    source[A].paste(mask)
    transp = Image.merge(transp.mode, source)  # build a new multiband image

    self.drawImage(transp, x - pos[0], y - pos[1])

  def drawRect(self, x1, y1, x2, y2, edgeColor=None, edgeWidth=None, fillColor=None):
    fillColor = self._getTkColor(fillColor, self.defaultFillColor)
    edgeColor = self._getTkColor(edgeColor, self.defaultLineColor)
    if edgeWidth is None:
      edgeWidth = self.defaultLineWidth
    new_item = self.create_rectangle(x1, y1, x2, y2, fill=fillColor, width=edgeWidth,
                                     outline=edgeColor)
    self._item_ids.append(new_item)

  # NYI:
  #def drawRoundRect(self, x1,y1, x2,y2, rx=5, ry=5,
  #                  edgeColor=None, edgeWidth=None, fillColor=None):

  def drawEllipse(self, x1, y1, x2, y2, edgeColor=None, edgeWidth=None, fillColor=None):
    fillColor = self._getTkColor(fillColor, self.defaultFillColor)
    edgeColor = self._getTkColor(edgeColor, self.defaultLineColor)
    if edgeWidth is None:
      edgeWidth = self.defaultLineWidth
    new_item = self.create_oval(x1, y1, x2, y2, fill=fillColor, outline=edgeColor, width=edgeWidth)

    self._item_ids.append(new_item)

  def drawArc(self, x1, y1, x2, y2, startAng=0, extent=360, edgeColor=None, edgeWidth=None,
              fillColor=None):
    fillColor = self._getTkColor(fillColor, self.defaultFillColor)
    edgeColor = self._getTkColor(edgeColor, self.defaultLineColor)
    if edgeWidth is None:
      edgeWidth = self.defaultLineWidth
    new_item = self.create_arc(x1, y1, x2, y2, start=startAng, extent=extent, fill=fillColor,
                               width=edgeWidth, outline=edgeColor)
    self._item_ids.append(new_item)

  def drawPolygon(self, pointlist, edgeColor=None, edgeWidth=None, fillColor=None, closed=0):
    fillColor = self._getTkColor(fillColor, self.defaultFillColor)
    edgeColor = self._getTkColor(edgeColor, self.defaultLineColor)
    if edgeWidth is None:
      edgeWidth = self.defaultLineWidth
    if closed:
      # draw a closed shape
      new_item = self.create_polygon(pointlist, fill=fillColor, width=edgeWidth, outline=edgeColor)
    else:
      if fillColor == self.__TRANSPARENT:
        # draw open-ended set of lines
        d = {'fill': edgeColor, 'width': edgeWidth}
        new_item = self.create_line(pointlist, **d)
      else:
        # open filled shape.
        # draw it twice:
        #    once as a polygon with no edge outline with the fill color
        #    and once as an open set of lines of the appropriate color
        new_item = self.create_polygon(pointlist, fill=fillColor, outline=self.__TRANSPARENT)

        self._item_ids.append(new_item)

        d = {'fill': edgeColor, 'width': edgeWidth}
        new_item = self.create_line(pointlist, **d)

    self._item_ids.append(new_item)


#def drawFigure(self, partList,
#               edgeColor=None, edgeWidth=None, fillColor=None):
# use default implementation

  def drawImage(self, image, x1, y1, x2=None, y2=None):

    try:
      from PIL import ImageTk
    except ImportError:
      raise NotImplementedError('drawImage - require the ImageTk module')

    w, h = image.size
    if not x2:
      x2 = w + x1
    if not y2:
      y2 = h + y1

    if (w != x2 - x1) or (h != y2 - y1):  # need to scale image
      myimage = image.resize((x2 - x1, y2 - y1))
    else:
      myimage = image

    # unless I keep a copy of this PhotoImage, it seems to be garbage collected
    # and the image is removed from the display after this function. weird
    itk = ImageTk.PhotoImage(myimage, master=self)
    new_item = self.create_image(x1, y1, image=itk, anchor=Tkinter.NW)
    self._item_ids.append(new_item)
    self._images.append(itk)

try:
  import rdkit.sping.PIL

  class TKCanvasPIL(rdkit.sping.PIL.PILCanvas):
    """This canvas maintains a PILCanvas as its backbuffer.  Drawing calls
        are made to the backbuffer and flush() sends the image to the screen
        using TKCanvas.
           You can also save what is displayed to a file in any of the formats
        supported by PIL"""

    def __init__(self, size=(300, 300), name='TKCanvas', master=None, **kw):
      rdkit.sping.PIL.PILCanvas.__init__(self, size=size, name=name)
      self._tkcanvas = TKCanvas(size, name, master, **kw)

    def flush(self):
      rdkit.sping.PIL.PILCanvas.flush(self)  # call inherited one first
      self._tkcanvas.drawImage(self._image, 0, 0)  # self._image should be a PIL image
      self._tkcanvas.flush()

    def getTKCanvas(self):
      return self._tkcanvas

except ImportError:
  raise ImportError("TKCanvasPIL requires sping  PIL Canvas, PIL may not be installed")
