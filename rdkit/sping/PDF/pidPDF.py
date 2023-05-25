# piddlePDF.py - a PDF backend for PIDDLE
"""This is the PIDDLE back end for PDF.  It acts as a wrapper
over the pdfgen.Canvas class, and translates between the
PIDDLE graphics state and the PDF/PostScript one. It only
exposes PIDDLE methods; however, it has an attribute
self.pdf which offers numerous lower-level drawing routines.
"""

# Copyright Robinson Analytics 1998-9
# 6/10/99 - near total rewrite.  This is just a wrapper
# on top of a lower-level graphics canvas which speaks the
# Postscript metaphor.  Offers native methods for everything
# except drawFigure, which doesn't behave like PostScript
# paths so I left it unchanged.

import os
from math import cos, pi, sin

from rdkit.sping import pagesizes
from rdkit.sping.pid import *

from . import pdfgen, pdfgeom, pdfmetrics

# edit this is the setting offends you, or set it in the constructor
DEFAULT_PAGE_SIZE = pagesizes.A4
#DEFAULT_PAGE_SIZE = pagesizes.letter

##########################################################################
#
#      Font Mappings
# The brute force approach to finding the correct postscript font name;
# much safer than the rule-based ones we tried.
# preprocessor to reduce font face names to the shortest list
# possible.  Add any aliases you wish; it keeps looking up
# until it finds no more translations to do.  Any input
# will be lowercased before checking.
font_face_map = {
  'serif': 'times',
  'sansserif': 'helvetica',
  'monospaced': 'courier',
  'arial': 'helvetica'
}

# maps a piddle font to a postscript one.
ps_font_map = {
  # face, bold, italic -> ps name
  ('times', 0, 0): 'Times-Roman',
  ('times', 1, 0): 'Times-Bold',
  ('times', 0, 1): 'Times-Italic',
  ('times', 1, 1): 'Times-BoldItalic',
  ('courier', 0, 0): 'Courier',
  ('courier', 1, 0): 'Courier-Bold',
  ('courier', 0, 1): 'Courier-Oblique',
  ('courier', 1, 1): 'Courier-BoldOblique',
  ('helvetica', 0, 0): 'Helvetica',
  ('helvetica', 1, 0): 'Helvetica-Bold',
  ('helvetica', 0, 1): 'Helvetica-Oblique',
  ('helvetica', 1, 1): 'Helvetica-BoldOblique',

  # there is only one Symbol font
  ('symbol', 0, 0): 'Symbol',
  ('symbol', 1, 0): 'Symbol',
  ('symbol', 0, 1): 'Symbol',
  ('symbol', 1, 1): 'Symbol',

  # ditto for dingbats
  ('zapfdingbats', 0, 0): 'ZapfDingbats',
  ('zapfdingbats', 1, 0): 'ZapfDingbats',
  ('zapfdingbats', 0, 1): 'ZapfDingbats',
  ('zapfdingbats', 1, 1): 'ZapfDingbats',
}

######################################################################
#
#     Canvas class
#
######################################################################


class PDFCanvas(Canvas):
  """This works by accumulating a list of strings containing
      PDF page marking operators, as you call its methods.  We could
      use a big string but this is more efficient - only concatenate
      it once, with control over line ends.  When
      done, it hands off the stream to a PDFPage object."""

  def __init__(self, size=None, name="pidPDF.pdf", pagesize=DEFAULT_PAGE_SIZE):
    # if no extension, add .PDF
    root, ext = os.path.splitext(name)
    if ext == '':
      name = root + '.pdf'

    # create the underlying pdfgen canvas and set some attributes
    self.pdf = pdfgen.Canvas(name, pagesize=pagesize, bottomup=0)  # add pagesize to constructor
    # by default do not use comrpression (mod by cwl, may not be necessary w/ newer pdfgen)
    self.pdf.setPageCompression(0)

    self.pdf.setLineCap(2)
    # now call super init, which will trigger
    # calls into self.pdf

    Canvas.__init__(self, size=size, name=name)

    # memorize stuff
    self.pagesize = pagesize
    self.filename = name

    # self.pdf.setPageSize(pagesize) # This doesn't seem to work correctly -cwl
    if size is None:
      # take the page size, which might not be default
      self.drawingsize = self.pagesize
    else:
      # convenience for other platformslike GUI views
      # we let them centre a smaller drawing in a page
      self.drawingsize = size

    self.pageTransitionString = ''
    self.pageNumber = 1  # keep a count

    # if they specified a size smaller than page,
    # be helpful and centre their diagram
    if self.pagesize != self.drawingsize:
      dx = 0.5 * (self.pagesize[0] - self.drawingsize[0])
      dy = 0.5 * (self.pagesize[1] - self.drawingsize[1])
      self.pdf.translate(dx, dy)

  def _resetDefaults(self):
    """Only used in setup - persist from page to page"""
    self.defaultLineColor = black
    self.defaultFillColor = transparent
    self.defaultLineWidth = 1
    self.defaultFont = Font()
    self.pdf.setLineCap(2)

  def showPage(self):
    """ensure basic settings are the same after a page break"""
    self.pdf.showPage()
    self.defaultFont = self.defaultFont
    self.defaultLineColor = self.defaultLineColor
    self.defaultFillColor = self.defaultFillColor
    self.defaultLineWidth = self.defaultLineWidth
    self.pdf.setLineCap(2)

  # ------------ canvas capabilities -------------

  def isInteractive(self):
    return 0

  def canUpdate(self):
    return 0

  # ------------ general management -------------
  def clear(self):
    "Not wll defined for file formats, use same as ShowPage"
    self.showPage()

  def flush(self):
    pass

  def save(self, file=None, format=None):
    """Saves the file.  If holding data, do
            a showPage() to save them having to."""

    if self.pdf.pageHasData():
      self.pdf.showPage()

    if hasattr(file, 'write'):
      self.pdf.save(fileobj=file)
    elif isinstance(file, str):
      self.pdf.save(filename=file)
    else:
      self.pdf.save()

  def setInfoLine(self, s):
    self.pdf.setTitle(s)

  # -------------handle assignment to defaultXXX-------

  def __setattr__(self, key, value):
    # we let it happen...
    self.__dict__[key] = value
    # ...but take action if needed
    if key == "defaultLineColor":
      self._updateLineColor(value)
    elif key == "defaultLineWidth":
      self._updateLineWidth(value)
    elif key == "defaultFillColor":
      self._updateFillColor(value)
    elif key == "defaultFont":
      self._updateFont(value)

  def _updateLineColor(self, color):
    """Triggered when someone assigns to defaultLineColor"""
    self.pdf.setStrokeColorRGB(color.red, color.green, color.blue)

  def _updateFillColor(self, color):
    """Triggered when someone assigns to defaultFillColor"""
    self.pdf.setFillColorRGB(color.red, color.green, color.blue)

  def _updateLineWidth(self, width):
    """Triggered when someone assigns to defaultLineWidth"""
    self.pdf.setLineWidth(width)

  def _updateFont(self, font):
    """Triggered when someone assigns to defaultFont"""
    psfont = self._findPostScriptFontName(font)
    self.pdf.setFont(psfont, font.size)

  def _findPostScriptFontName(self, font):
    """Attempts to return proper font name."""

    # step 1 - no face ends up serif, others are lowercased
    if not font.face:
      face = 'serif'
    else:
      face = font.face.lower()
    while face in font_face_map:
      face = font_face_map[face]
    # step 2, - resolve bold/italic to get the right PS font name
    psname = ps_font_map[(face, font.bold, font.italic)]
    return psname

  def _escape(self, s):
    """PDF escapes are like Python ones, but brackets need slashes before them too.
            Use Python's repr function and chop off the quotes first"""
    s = repr(s)[1:-1]
    s = s.replace('(', r'\(')
    s = s.replace(')', r'\)')
    return s

  def resetDefaults(self):
    """If you drop down to a lower level, PIDDLE can lose
            track of the current graphics state.  Calling this after
            wards ensures that the canvas is updated to the same
            defaults as PIDDLE thinks they should be."""
    self.defaultFont = self.defaultFont
    self.defaultLineColor = self.defaultLineColor
    self.defaultFillColor = self.defaultFillColor
    self.defaultLineWidth = self.defaultLineWidth

  # ------------ string/font info ------------

  def stringWidth(self, s, font=None):
    "Return the logical width of the string if it were drawn \
            in the current font (defaults to self.font)."

    if not font:
      font = self.defaultFont
    fontname = self._findPostScriptFontName(font)
    return pdfmetrics.stringwidth(s, fontname) * font.size * 0.001

  def fontHeight(self, font=None):
    if not font:
      font = self.defaultFont
    return font.size

  def fontAscent(self, font=None):
    if not font:
      font = self.defaultFont
    fontname = self._findPostScriptFontName(font)
    return pdfmetrics.ascent_descent[fontname][0] * 0.001 * font.size

  def fontDescent(self, font=None):
    if not font:
      font = self.defaultFont
    fontname = self._findPostScriptFontName(font)
    return -pdfmetrics.ascent_descent[fontname][1] * 0.001 * font.size

  # ------------- drawing helpers --------------
  def _endPath(self, path, edgeColor, fillColor):
    """in PIDDLE, the edge and fil colors might be transparent,
            and might also be None, in which case they should be taken
            from the defaults.  This leads to a standard 10 lines of code
            when closing each shape, which are wrapped up here.  Use
            these if you implement new PIDDLE shapes."""
    # allow for transparent fills and lines
    fill = fillColor or self.defaultFillColor
    edge = edgeColor or self.defaultLineColor
    if (fill == transparent and edge == transparent):
      pass
    else:
      self.pdf.drawPath(
        path,
        (edge != transparent),  # whether to stroke
        (fill != transparent)  # whether to fill
      )

    # ------------- drawing methods --------------

  def drawLine(self, x1, y1, x2, y2, color=None, width=None, dash=None, **kwargs):
    """Calls the underlying methods in pdfgen.canvas.  For the
            highest performance, use canvas.setDefaultFont and
            canvas.setLineWidth, and draw batches of similar
            lines together."""
    # set the state if needed
    if color:
      self._updateLineColor(color)
    if width:
      self._updateLineWidth(width)

    # now do the work
    self.pdf.line(x1, y1, x2, y2)

    # now reset state if needed
    if color:
      self._updateLineColor(self.defaultLineColor)
    if width:
      self._updateLineWidth(self.defaultLineWidth)

  def drawLines(self, lineList, color=None, width=None, dash=None, **kwargs):
    """Draws several distinct lines, all with same color
            and width, efficiently"""
    if color:
      self._updateLineColor(color)
    if width:
      self._updateLineWidth(width)

    self.pdf.lines(lineList)

    if color:
      self._updateLineColor(self.defaultLineColor)
    if width:
      self._updateLineWidth(self.defaultLineWidth)

  def drawString(self, s, x, y, font=None, color=None, angle=0, **kwargs):
    """As it says, but many options to process.  It translates
            user space rather than text space, in case underlining is
            needed on rotated text.  It cheats and does literals
            for efficiency, avoiding changing the python graphics state."""
    self.pdf.addLiteral('%begin drawString')
    col = color or self.defaultLineColor
    if col != transparent:
      if '\n' in s or '\r' in s:
        # normalize line ends
        s = s.replace('\r\n', '\n')
        s = s.replace('\n\r', '\n')
        lines = s.split('\n')
      else:
        lines = [s]
      fnt = font or self.defaultFont
      self._updateFont(fnt)
      text = self.pdf._escape(s)

      # start of Chris's hacking
      # inserting basic commands here  to see if can get working
      textobj = self.pdf.beginText()

      if col != self.defaultFillColor:
        textobj.setFillColorRGB(col.red, col.green, col.blue)

      if angle != 0:
        co = cos(angle * pi / 180.0)
        si = sin(angle * pi / 180.0)
        textobj.setTextTransform(co, -si, si, co, x, y)  # top down coords so reverse angle
      else:
        textobj.setTextOrigin(x, y)

      for line in lines:
        # keep underlining separate - it is slow and unusual anyway
        if fnt.underline:
          # breaks on angled text - FIXME
          ycursor = textobj.getY()  # returns offset from last set origin
          dy = 0.5 * self.fontDescent(fnt)
          width = self.stringWidth(line, fnt)
          linewidth = fnt.size * 0.1

          self.pdf.saveState()
          self.pdf.setLineWidth(linewidth)
          self.pdf.translate(x, y)  # need to translate first before rotate
          if angle != 0:
            self.pdf.rotate(-angle)
          self.pdf.translate(0, ycursor - y)  # move down to start of current text line
          self.pdf.line(0, dy, width, dy)
          self.pdf.restoreState()
          lasty = ycursor
        textobj.textLine(line)  # adds text to textobj, advances getY's cursor
      # finally actually send text object to the page
      self.pdf.drawText(textobj)  # draw all the text afterwards? Doesn't seem right
    self.pdf.addLiteral('%end drawString')
    # done wth drawString()

  def drawCurve(self, x1, y1, x2, y2, x3, y3, x4, y4, edgeColor=None, edgeWidth=None,
                fillColor=None, closed=0, dash=None, **kwargs):
    """This could do two totally different things.  If not closed,
            just does a bezier curve so fill is irrelevant.  If closed,
            it is actually a filled shape."""
    if closed:
      if edgeColor:
        self._updateLineColor(edgeColor)
      if edgeWidth:
        self._updateLineWidth(edgeWidth)
      if fillColor:
        self._updateFillColor(fillColor)

      p = self.pdf.beginPath()
      p.moveTo(x1, y1)
      p.curveTo(x2, y2, x3, y3, x4, y4)
      p.close()
      self._endPath(p, edgeColor, fillColor)  # handles case of transparency

      if edgeColor:
        self._updateLineColor(self.defaultLineColor)
      if edgeWidth:
        self._updateLineWidth(self.defaultLineWidth)
      if fillColor:
        self._updateFillColor(self.defaultFillColor)
    else:
      # just a plain old line segment
      if edgeColor:
        self._updateLineColor(edgeColor)
      if edgeWidth:
        self._updateLineWidth(edgeWidth)

      self.pdf.bezier(x1, y1, x2, y2, x3, y3, x4, y4)

      if edgeColor:
        self._updateLineColor(self.defaultLineColor)
      if edgeWidth:
        self._updateLineWidth(self.defaultLineWidth)

  def drawRect(self, x1, y1, x2, y2, edgeColor=None, edgeWidth=None, fillColor=None, dash=None,
               **kwargs):
    if edgeColor:
      self._updateLineColor(edgeColor)
    if edgeWidth:
      self._updateLineWidth(edgeWidth)
    if fillColor:
      self._updateFillColor(fillColor)

    p = self.pdf.beginPath()
    p.rect(x1, y1, x2 - x1, y2 - y1)
    self._endPath(p, edgeColor, fillColor)  # handles case of transparency

    if edgeColor:
      self._updateLineColor(self.defaultLineColor)
    if edgeWidth:
      self._updateLineWidth(self.defaultLineWidth)
    if fillColor:
      self._updateFillColor(self.defaultFillColor)

  # drawRoundRect is inherited - cannot really improve on that one,
  # and figures are quite efficient now.
  def drawEllipse(self, x1, y1, x2, y2, edgeColor=None, edgeWidth=None, fillColor=None, dash=None,
                  **kwargs):
    if edgeColor:
      self._updateLineColor(edgeColor)
    if edgeWidth:
      self._updateLineWidth(edgeWidth)
    if fillColor:
      self._updateFillColor(fillColor)

    p = self.pdf.beginPath()
    p.ellipse(x1, y1, x2 - x1, y2 - y1)
    self._endPath(p, edgeColor, fillColor)  # handles case of transparency

    if edgeColor:
      self._updateLineColor(self.defaultLineColor)
    if edgeWidth:
      self._updateLineWidth(self.defaultLineWidth)
    if fillColor:
      self._updateFillColor(self.defaultFillColor)

  def drawArc(self, x1, y1, x2, y2, startAng=0, extent=90, edgeColor=None, edgeWidth=None,
              fillColor=None, dash=None, **kwargs):
    """This draws a PacMan-type shape connected to the centre.  One
            idiosyncrasy - if you specify an edge color, it apples to the
            outer curved rim but not the radial edges."""
    if edgeColor:
      self._updateLineColor(edgeColor)
    if edgeWidth:
      self._updateLineWidth(edgeWidth)
    if fillColor:
      self._updateFillColor(fillColor)
    # I need to do some more work on flipping the coordinate system -
    # in pdfgen - note the angle reversal needed when drawing top-down.
    pointList = pdfgeom.bezierArc(x1, y1, x2, y2, -startAng, -extent)
    start = pointList[0]
    end = pointList[-1]
    x_cen = 0.5 * (x1 + x2)
    y_cen = 0.5 * (y1 + y2)

    # first the fill
    p = self.pdf.beginPath()
    p.moveTo(x_cen, y_cen)
    p.lineTo(start[0], start[1])
    for curve in pointList:
      p.curveTo(curve[2], curve[3], curve[4], curve[5], curve[6], curve[7])
    p.close()  # back to centre
    self._endPath(p, transparent, fillColor)  # handles case of transparency
    # now the outer rim
    p2 = self.pdf.beginPath()
    p2.moveTo(start[0], start[1])
    for curve in pointList:
      p2.curveTo(curve[2], curve[3], curve[4], curve[5], curve[6], curve[7])
    self._endPath(p2, edgeColor, transparent)  # handles case of transparency

    if edgeColor:
      self._updateLineColor(self.defaultLineColor)
    if edgeWidth:
      self._updateLineWidth(self.defaultLineWidth)
    if fillColor:
      self._updateFillColor(self.defaultFillColor)

  def drawPolygon(self, pointlist, edgeColor=None, edgeWidth=None, fillColor=None, closed=0,
                  dash=None, **kwargs):
    """As it says.  Easy with paths!"""
    if edgeColor:
      self._updateLineColor(edgeColor)
    if edgeWidth:
      self._updateLineWidth(edgeWidth)
    if fillColor:
      self._updateFillColor(fillColor)

    p = self.pdf.beginPath()
    p.moveTo(pointlist[0][0], pointlist[0][1])
    for point in pointlist[1:]:
      p.lineTo(point[0], point[1])
    if closed:
      p.close()

    self._endPath(p, edgeColor, fillColor)  # handles case of transparency

    if edgeColor:
      self._updateLineColor(self.defaultLineColor)
    if edgeWidth:
      self._updateLineWidth(self.defaultLineWidth)
    if fillColor:
      self._updateFillColor(self.defaultFillColor)


# def drawFigure(self, partList,
# edgeColor=None, edgeWidth=None, fillColor=None, closed=0):
# """This is PIDDLE's attempt at Postscript paths.  Due to
# necessary limitations in the algorithm, if the start and end
# points are not connected but closed=1, then you get the full shape
# filled but the final line segment does not join up.  I have to
# do extra work to simulate this."""
##
# if edgeColor:
# self._updateLineColor(edgeColor)
# if edgeWidth:
# self._updateLineWidth(edgeWidth)
# if fillColor:
# self._updateFillColor(fillColor)
##
# p1 = self.pdf.beginPath()  #use for the fill (i.e. closed)
# p2 = self.pdf.beginPath()  #use for the edge (may not be closed)
##
# move to first point
##        start = (partList[0][1:3])
##        end = None
##        p1.moveTo(start[0], start[1])
##        p2.moveTo(start[0], start[1])
##
# for tuple in partList:
##            op = tuple[0]
##            args = list(tuple[1:])
##            start = args[0:2]
# lineTo the start if not coincident with end of last segment
# if start != end:
##                p1.lineTo(start[0], start[1])
##                p2.lineTo(start[0], start[1])
##
# now draw appropriate segment
# if op == figureLine:
##                p1.lineTo(args[2], args[3])
##                p2.lineTo(args[2], args[3])
##                end = args[2:4]
# elif op == figureArc:
# p1.arcTo(args[0], args[1], args[2], args[3], args[4], args[5])
# p2.arcTo(args[0], args[1], args[2], args[3], args[4], args[5])
##                p1.arc(args[0], args[1], args[2], args[3], args[4], args[5])
##                p2.arc(args[0], args[1], args[2], args[3], args[4], args[5])
##                end = args[2:4]
# elif op == figureCurve:
##                p1.curveTo(args[2], args[3], args[4], args[5], args[6], args[7])
##                p2.curveTo(args[2], args[3], args[4], args[5], args[6], args[7])
##                end = args[6:8]
# else:
# raise TypeError, "unknown figure operator: " + op
##
# now for the weirdness
# p1.close()
# if closed:
# p2.close()
# print 'closed edge path'
# print 'inner path p1:' + p1.getCode()
# print 'outer path p2:' + p2.getCode()
##
##        self._endPath(p1, transparent, fillColor)
##        self._endPath(p2, edgeColor, transparent)
##
# if edgeColor:
# self._updateLineColor(self.defaultLineColor)
# if edgeWidth:
# self._updateLineWidth(self.defaultLineWidth)
# if fillColor:
# self._updateFillColor(self.defaultFillColor)

  def drawImage(self, image, x1, y1, x2=None, y2=None, **kwargs):
    """Draw a PIL Image or image filename into the specified rectangle.
            If x2 and y2 are omitted, they are calculated from the image size.
            """
    # chris starts meddling here -cwl
    # piddle only takes PIL images
    im_width, im_height = image.size
    if not x2:
      x2 = x1 + im_width
    if not y2:
      y2 = y1 + im_height

    self.pdf.saveState()  # I'm changing coordinates to isolate the problem -cwl

    self.pdf.translate(x1, y1)
    self.pdf.drawInlineImage(image, 0, 0, abs(x1 - x2), abs(y1 - y2))

    self.pdf.restoreState()

  # original code below -cwl
  # the underlying canvas uses a bott-up coord system, so flips things
  # if x2:
  #width = abs(x2 - x1)
  #x = min(x1, x2)
  # if y2:
  #height = abs(y2 - y1)
  #y = min(y1, y2)
  #self.pdf.drawInlineImage(image, x, y, width, height)

  ##########################################################
  #
  #   non-standard extensions outside Piddle API
  #
  ##########################################################

  def drawLiteral(self, literal):
    # adds a chunk of raw stuff to the PDF contents stream
    self.code.append(literal)


def test():
  # ... for testing...
  canvas = PDFCanvas(name="test")

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


if __name__ == '__main__':
  test()
  # dashtest()
  # test2()
