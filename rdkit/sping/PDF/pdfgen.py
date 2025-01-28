## Automatically adapted for numpy.oldnumeric Jun 27, 2008 by -c

#pdfgen.py
"""
PDFgen is a library to generate PDF files containing text and graphics.  It is the
foundation for a complete reporting solution in Python.  It is also the
foundation for piddlePDF, the PDF back end for PIDDLE.

Documentation is a little slim right now; run then look at testpdfgen.py
to get a clue.

---------- Licence Terms (same as the Python license) -----------------
(C) Copyright Robinson Analytics 1998-1999.

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted, provided
that the above copyright notice appear in all copies and that both that
copyright notice and this permission notice appear in supporting
documentation, and that the name of Robinson Analytics not be used
in advertising or publicity pertaining to distribution of the software
without specific, written prior permission.

ROBINSON ANALYTICS LTD. DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS,
IN NO EVENT SHALL ROBINSON ANALYTICS BE LIABLE FOR ANY SPECIAL, INDIRECT
OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE.

Progress Reports:
0.82, 1999-10-27, AR:
        Fixed some bugs on printing to Postscript.  Added 'Text Object'
        analogous to Path Object to control entry and exit from text mode.
        Much simpler clipping API.  All verified to export Postscript and
        redistill.
        One limitation still - clipping to text paths is fine in Acrobat
        but not in Postscript (any level)

0.81,1999-10-13, AR:
        Adding RoundRect; changed all format strings to use %0.2f instead of %s,
        so we don't get exponentials in the output.
0.8,1999-10-07, AR:  all changed!

2000-02-07: changed all %0.2f's to %0.4f in order to allow precise ploting of graphs that
        range between 0 and 1 -cwl

"""
##  0.81    1999-10-13:
##
##
##

import os
import sys
import tempfile
import time
from io import StringIO
from math import ceil, cos, pi, sin, tan
from types import *

from . import pdfdoc, pdfgeom, pdfmetrics, pdfutils


class PDFError(ValueError):
  pass


# Robert Kern
# Constants for closing paths.
# May be useful if one changes 'arc' and 'rect' to take a
# default argument that tells how to close the path.
# That way we can draw filled shapes.

FILL_EVEN_ODD = 0
FILL_NON_ZERO = 1
#this is used by path-closing routines.
#map stroke, fill, fillmode -> operator
# fillmode: 1 = non-Zero (obviously), 0 = evenOdd
PATH_OPS = {
  (0, 0, FILL_EVEN_ODD): 'n',  #no op
  (0, 0, FILL_NON_ZERO): 'n',  #no op
  (1, 0, FILL_EVEN_ODD): 'S',  #stroke only
  (1, 0, FILL_NON_ZERO): 'S',  #stroke only
  (0, 1, FILL_EVEN_ODD): 'f*',  #Fill only
  (0, 1, FILL_NON_ZERO): 'f',  #Fill only
  (1, 1, FILL_EVEN_ODD): 'B*',  #Stroke and Fill
  (1, 1, FILL_NON_ZERO): 'B',  #Stroke and Fill
}

close = 'h'
newpath = 'n'
stroke = 'S'
closeStroke = 's'
nzFill = 'f'
eoFill = 'f*'
fillStroke = 'B'
closeFillStroke = 'b'
eoFillStroke = 'B*'
closeEoFillStroke = 'b*'


class Canvas:
  """This is a low-level interface to the PDF file format.  The plan is to
    expose the whole pdfgen API through this.  Its drawing functions should have a
    one-to-one correspondence with PDF functionality.  Unlike PIDDLE, it thinks
    in terms of RGB values, Postscript font names, paths, and a 'current graphics
    state'.  Just started development at 5/9/99, not in use yet.

    """

  def __init__(self, filename, pagesize=(595.27, 841.89), bottomup=1):
    """Most of the attributes are private - we will use set/get methods
        as the preferred interface.  Default page size is A4."""
    self._filename = filename
    self._doc = pdfdoc.PDFDocument()
    self._pagesize = pagesize
    self._currentPageHasImages = 1
    self._pageTransitionString = ''

    self._pageCompression = 1  #on by default - turn off when debugging!
    self._pageNumber = 1  # keep a count
    self._code = []  #where the current page's marking operators accumulate

    #PostScript has the origin at bottom left. It is easy to achieve a top-
    #down coord system by translating to the top of the page and setting y
    #scale to -1, but then text is inverted.  So self.bottomup is used
    #to also set the text matrix accordingly.  You can now choose your
    #drawing coordinates.
    self.bottomup = bottomup
    if self.bottomup:
      #set initial font
      #self._preamble = 'BT /F9 12 Tf 14.4 TL ET'
      self._preamble = '1 0 0 1 0 0 cm BT /F9 12 Tf 14.4 TL ET'
    else:
      #switch coordinates, flip text and set font
      #self._preamble = '1 0 0 -1 0 %0.4f cm BT /F9 12 Tf 14.4 TL ET' % self._pagesize[1]
      self._preamble = '1 0 0 -1 0 %0.4f cm BT /F9 12 Tf 14.4 TL ET' % self._pagesize[1]

    #initial graphics state
    self._x = 0
    self._y = 0
    self._fontname = 'Times-Roman'
    self._fontsize = 12
    self._textMode = 0  #track if between BT/ET
    self._leading = 14.4
    self._currentMatrix = (1., 0., 0., 1., 0., 0.)
    self._fillMode = 0  #even-odd

    #text state
    self._charSpace = 0
    self._wordSpace = 0
    self._horizScale = 100
    self._textRenderMode = 0
    self._rise = 0
    self._textLineMatrix = (1., 0., 0., 1., 0., 0.)
    self._textMatrix = (1., 0., 0., 1., 0., 0.)

    # line drawing
    self._lineCap = 0
    self._lineJoin = 0
    self._lineDash = None  #not done
    self._lineWidth = 0
    self._mitreLimit = 0

    self._fillColorRGB = (0, 0, 0)
    self._strokeColorRGB = (0, 0, 0)

  def _escape(self, s):
    """PDF escapes are like Python ones, but brackets need slashes before them too.
        Use Python's repr function and chop off the quotes first"""
    s = repr(s)[1:-1]
    s = s.replace('(', r'\(')
    s = s.replace(')', r'\)')
    return s

  #info functions - non-standard
  def setAuthor(self, author):
    self._doc.setAuthor(author)

  def setTitle(self, title):
    self._doc.setTitle(title)

  def setSubject(self, subject):
    self._doc.setSubject(subject)

  def pageHasData(self):
    "Info function - app can call it after showPage to see if it needs a save"
    return len(self._code) == 0

  def showPage(self):
    """This is where the fun happens"""
    page = pdfdoc.PDFPage()
    page.pagewidth = self._pagesize[0]
    page.pageheight = self._pagesize[1]
    page.hasImages = self._currentPageHasImages
    page.pageTransitionString = self._pageTransitionString
    page.setCompression(self._pageCompression)
    #print stream
    page.setStream([self._preamble] + self._code)
    self._doc.addPage(page)

    #now get ready for the next one
    self._pageNumber = self._pageNumber + 1
    self._code = []  # ready for more...
    self._currentPageHasImages = 0

  def getPageNumber(self):
    return self._pageNumber

  def save(self, filename=None, fileobj=None):
    """Saves the pdf document to fileobj or to file with name filename.
        If holding data, do a showPage() to save them having to."""

    if len(self._code):
      self.showPage()  # what's the effect of multiple 'showPage's
    if fileobj:
      self._doc.SaveToFileObject(fileobj)
    elif filename:
      self._doc.SaveToFile(filename)
    else:
      self._doc.SaveToFile(self._filename)

  def setPageSize(self, size):
    """accepts a 2-tuple in points for paper size for this
        and subsequent pages"""
    self._pagesize = size

  def addLiteral(self, s, escaped=1):
    if escaped == 0:
      s = self._escape(s)
    self._code.append(s)

    ######################################################################
    #
    #      coordinate transformations
    #
    ######################################################################

  def transform(self, a, b, c, d, e, f):
    """How can Python track this?"""
    a0, b0, c0, d0, e0, f0 = self._currentMatrix
    self._currentMatrix = (a0 * a + c0 * b, b0 * a + d0 * b, a0 * c + c0 * d, b0 * c + d0 * d,
                           a0 * e + c0 * f + e0, b0 * e + d0 * f + f0)
    self._code.append('%0.4f %0.4f %0.4f %0.4f %0.4f %0.4f cm' % (a, b, c, d, e, f))

  def translate(self, dx, dy):
    self.transform(1, 0, 0, 1, dx, dy)

  def scale(self, x, y):
    self.transform(x, 0, 0, y, 0, 0)

  def rotate(self, theta):
    """Canvas.rotate(theta)

        theta is in degrees."""
    c = cos(theta * pi / 180)
    s = sin(theta * pi / 180)
    self.transform(c, s, -s, c, 0, 0)

  def skew(self, alpha, beta):
    tanAlpha = tan(alpha * pi / 180)
    tanBeta = tan(beta * pi / 180)
    self.transform(1, tanAlpha, tanBeta, 1, 0, 0)

    ######################################################################
    #
    #      graphics state management
    #
    ######################################################################

  def saveState(self):
    """These need expanding to save/restore Python's state tracking too"""
    self._code.append('q')

  def restoreState(self):
    """These need expanding to save/restore Python's state tracking too"""
    self._code.append('Q')

    ###############################################################
    #
    #   Drawing methods.  These draw things directly without
    #   fiddling around with Path objects.  We can add any geometry
    #   methods we wish as long as their meaning is precise and
    #   they are of general use.
    #
    #   In general there are two patterns.  Closed shapes
    #   have the pattern shape(self, args, stroke=1, fill=0);
    #   by default they draw an outline only. Line segments come
    #   in three flavours: line, bezier, arc (which is a segment
    #   of an elliptical arc, approximated by up to four bezier
    #   curves, one for each quadrant.
    #
    #   In the case of lines, we provide a 'plural' to unroll
    #   the inner loop; it is useful for drawing big grids
    ################################################################

    #--------first the line drawing methods-----------------------
  def line(self, x1, y1, x2, y2):
    "As it says"
    self._code.append('n %0.4f %0.4f m %0.4f %0.4f l S' % (x1, y1, x2, y2))

  def lines(self, linelist):
    """As line(), but slightly more efficient for lots of them -
        one stroke operation and one less function call"""
    self._code.append('n')
    for (x1, y1, x2, y2) in linelist:
      self._code.append('%0.4f %0.4f m %0.4f %0.4f l' % (x1, y1, x2, y2))
    self._code.append('S')

  def grid(self, xlist, ylist):
    """Lays out a grid in current line style.  Suuply list of
        x an y positions."""
    assert len(xlist) > 1, "x coordinate list must have 2+ items"
    assert len(ylist) > 1, "y coordinate list must have 2+ items"
    lines = []
    y0, y1 = ylist[0], ylist[-1]
    x0, x1 = xlist[0], xlist[-1]
    for x in xlist:
      lines.append(x, y0, x, y1)
    for y in ylist:
      lines.append(x0, y, x1, y)
    self.lines(lines)

  def bezier(self, x1, y1, x2, y2, x3, y3, x4, y4):
    "Bezier curve with the four given control points"
    self._code.append('n %0.4f %0.4f m %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f c S' %
                      (x1, y1, x2, y2, x3, y3, x4, y4))

  def arc(self, x1, y1, x2, y2, startAng=0, extent=90):
    """Contributed to piddlePDF by Robert Kern, 28/7/99.
        Trimmed down by AR to remove color stuff for pdfgen.canvas and
        revert to positive coordinates.

        Draw a partial ellipse inscribed within the rectangle x1,y1,x2,y2,
        starting at startAng degrees and covering extent degrees.   Angles
        start with 0 to the right (+x) and increase counter-clockwise.
        These should have x1<x2 and y1<y2.

        The algorithm is an elliptical generalization of the formulae in
        Jim Fitzsimmon's TeX tutorial <URL: http://www.tinaja.com/bezarc1.pdf>."""

    pointList = pdfgeom.bezierArc(x1, y1, x2, y2, startAng, extent)
    #move to first point
    self._code.append('n %0.4f %0.4f m' % pointList[0][:2])
    for curve in pointList:
      self._code.append('%0.4f %0.4f %0.4f %0.4f %0.4f %0.4f c' % curve[2:])
    # stroke
    self._code.append('S')

    #--------now the shape drawing methods-----------------------
  def rect(self, x, y, width, height, stroke=1, fill=0):
    "draws a rectangle"
    self._code.append('n %0.4f %0.4f %0.4f %0.4f re ' % (x, y, width, height) +
                      PATH_OPS[stroke, fill, self._fillMode])

  def ellipse(self, x1, y1, x2, y2, stroke=1, fill=0):
    """Uses bezierArc, which conveniently handles 360 degrees -
        nice touch Robert"""
    pointList = pdfgeom.bezierArc(x1, y1, x2, y2, 0, 360)
    #move to first point
    self._code.append('n %0.4f %0.4f m' % pointList[0][:2])
    for curve in pointList:
      self._code.append('%0.4f %0.4f %0.4f %0.4f %0.4f %0.4f c' % curve[2:])
    #finish
    self._code.append(PATH_OPS[stroke, fill, self._fillMode])

  def wedge(self, x1, y1, x2, y2, startAng, extent, stroke=1, fill=0):
    """Like arc, but connects to the centre of the ellipse.
        Most useful for pie charts and PacMan!"""

    x_cen = (x1 + x2) / 2.
    y_cen = (y1 + y2) / 2.
    pointList = pdfgeom.bezierArc(x1, y1, x2, y2, startAng, extent)

    self._code.append('n %0.4f %0.4f m' % (x_cen, y_cen))
    # Move the pen to the center of the rectangle
    self._code.append('%0.4f %0.4f l' % pointList[0][:2])
    for curve in pointList:
      self._code.append('%0.4f %0.4f %0.4f %0.4f %0.4f %0.4f c' % curve[2:])
    # finish the wedge
    self._code.append('%0.4f %0.4f l ' % (x_cen, y_cen))
    # final operator
    self._code.append(PATH_OPS[stroke, fill, self._fillMode])

  def circle(self, x_cen, y_cen, r, stroke=1, fill=0):
    """special case of ellipse"""

    x1 = x_cen - r
    x2 = x_cen + r
    y1 = y_cen - r
    y2 = y_cen + r
    self.ellipse(x1, y1, x2, y2, stroke, fill)

  def roundRect(self, x, y, width, height, radius, stroke=1, fill=0):
    """Draws a rectangle with rounded corners.  The corners are
        approximately quadrants of a circle, with the given radius."""
    #use a precomputed set of factors for the bezier approximation
    #to a circle. There are six relevant points on the x axis and y axis.
    #sketch them and it should all make sense!
    t = 0.4472 * radius

    x0 = x
    x1 = x0 + t
    x2 = x0 + radius
    x3 = x0 + width - radius
    x4 = x0 + width - t
    x5 = x0 + width

    y0 = y
    y1 = y0 + t
    y2 = y0 + radius
    y3 = y0 + height - radius
    y4 = y0 + height - t
    y5 = y0 + height

    self._code.append('n %0.4f %0.4f m' % (x2, y0))
    self._code.append('%0.4f %0.4f l' % (x3, y0))  # bottom row
    self._code.append('%0.4f %0.4f %0.4f %0.4f %0.4f %0.4f c' %
                      (x4, y0, x5, y1, x5, y2))  # bottom right

    self._code.append('%0.4f %0.4f l' % (x5, y3))  # right edge
    self._code.append('%0.4f %0.4f %0.4f %0.4f %0.4f %0.4f c' %
                      (x5, y4, x4, y5, x3, y5))  # top right

    self._code.append('%0.4f %0.4f l' % (x2, y5))  # top row
    self._code.append('%0.4f %0.4f %0.4f %0.4f %0.4f %0.4f c' %
                      (x1, y5, x0, y4, x0, y3))  # top left

    self._code.append('%0.4f %0.4f l' % (x0, y2))  # left edge
    self._code.append('%0.4f %0.4f %0.4f %0.4f %0.4f %0.4f c' %
                      (x0, y1, x1, y0, x2, y0))  # bottom left

    self._code.append('h')  #close off, although it should be where it started anyway

    self._code.append(PATH_OPS[stroke, fill, self._fillMode])
    ##################################################
    #
    #  Text methods
    #
    # As with graphics, a separate object ensures that
    # everything is bracketed between  text operators.
    # The methods below are a high-level convenience.
    # use PDFTextObject for multi-line text.
    ##################################################

  def drawString(self, x, y, text):
    """Draws a string in the current text styles."""
    #we could inline this for speed if needed
    t = self.beginText(x, y)
    t.textLine(text)
    self.drawText(t)

  def drawRightString(self, x, y, text):
    """Draws a string right-aligned with the y coordinate"""
    width = self.stringWidth(text, self._fontname, self._fontsize)
    t = self.beginText(x - width, y)
    t.textLine(text)
    self.drawText(t)

  def drawCentredString(self, x, y, text):
    """Draws a string right-aligned with the y coordinate.  I
        am British so the spelling is correct, OK?"""
    width = self.stringWidth(text, self._fontname, self._fontsize)
    t = self.beginText(x - 0.5 * width, y)
    t.textLine(text)
    self.drawText(t)

  def getAvailableFonts(self):
    """Returns the list of PostScript font names available.
        Standard set now, but may grow in future with font embedding."""
    fontnames = self._doc.getAvailableFonts()
    fontnames.sort()
    return fontnames

  def setFont(self, psfontname, size, leading=None):
    """Sets the font.  If leading not specified, defaults to 1.2 x
        font size. Raises a readable exception if an illegal font
        is supplied.  Font names are case-sensitive! Keeps track
        of font anme and size for metrics."""
    self._fontname = psfontname
    self._fontsize = size
    pdffontname = self._doc.getInternalFontName(psfontname)
    if leading is None:
      leading = size * 1.2
    self._leading = leading
    self._code.append('BT %s %0.1f Tf %0.1f TL ET' % (pdffontname, size, leading))

  def stringWidth(self, text, fontname, fontsize):
    "gets width of a string in the given font and size"
    return pdfmetrics.stringwidth(text, fontname) * 0.001 * fontsize

  # basic graphics modes
  def setLineWidth(self, width):
    self._lineWidth = width
    self._code.append('%0.4f w' % width)

  def setLineCap(self, mode):
    """0=butt,1=round,2=square"""
    assert mode in (0, 1, 2), "Line caps allowed: 0=butt,1=round,2=square"
    self._lineCap = mode
    self._code.append('%d J' % mode)

  def setLineJoin(self, mode):
    """0=mitre, 1=round, 2=bevel"""
    assert mode in (0, 1, 2), "Line Joins allowed: 0=mitre, 1=round, 2=bevel"
    self._lineJoin = mode
    self._code.append('%d j' % mode)

  def setMiterLimit(self, limit):
    self._miterLimit = limit
    self._code.append('%0.4f M' % limit)

  def setDash(self, array=[], phase=0):
    """Two notations.  pass two numbers, or an array and phase"""
    if isinstance(array, (int, float)):
      self._code.append('[%s %s] 0 d' % (array, phase))
    elif isinstance(array, (list, tuple)):
      assert phase <= len(array), "setDash phase must be l.t.e. length of array"
      textarray = ' '.join(map(str, array))
      self._code.append('[%s] %s d' % (textarray, phase))

  def setFillColorRGB(self, r, g, b):
    self._fillColorRGB = (r, g, b)
    self._code.append('%0.4f %0.4f %0.4f rg' % (r, g, b))

  def setStrokeColorRGB(self, r, g, b):
    self._strokeColorRGB = (r, g, b)
    self._code.append('%0.4f %0.4f %0.4f RG' % (r, g, b))

  # path stuff - the separate path object builds it
  def beginPath(self):
    """Returns a fresh path object"""
    return PDFPathObject()

  def drawPath(self, aPath, stroke=1, fill=0):
    "Draw in the mode indicated"
    op = PATH_OPS[stroke, fill, self._fillMode]
    self._code.append(aPath.getCode() + ' ' + op)

  def clipPath(self, aPath, stroke=1, fill=0):
    "clip as well as drawing"
    op = PATH_OPS[stroke, fill, self._fillMode]
    self._code.append(aPath.getCode() + ' W ' + op)

  def beginText(self, x=0, y=0):
    """Returns a fresh text object"""
    return PDFTextObject(self, x, y)

  def drawText(self, aTextObject):
    """Draws a text object"""
    self._code.append(aTextObject.getCode())

    ######################################################
    #
    #   Image routines
    #
    ######################################################
  def drawInlineImage(self, image, x, y, width=None, height=None):
    """Draw a PIL Image into the specified rectangle.  If width and
        height are omitted, they are calculated from the image size.
        Also allow file names as well as images.  This allows a
        caching mechanism"""
    # print "drawInlineImage: x=%s, y=%s, width = %s, height=%s " % (x,y, width, height)
    try:
      from PIL import Image
    except ImportError:
      print('Python Imaging Library not available')
      return
    try:
      import zlib
    except ImportError:
      print('zlib not available')
      return

    self._currentPageHasImages = 1
    if isinstance(image, str):
      if os.path.splitext(image)[1] in ['.jpg', '.JPG']:
        #directly process JPEG files
        #open file, needs some error handling!!
        imageFile = open(image, 'rb')
        info = self.readJPEGInfo(imageFile)
        imgwidth, imgheight = info[0], info[1]
        if info[2] == 1:
          colorSpace = 'DeviceGray'
        elif info[2] == 3:
          colorSpace = 'DeviceRGB'
        else:  #maybe should generate an error, is this right for CMYK?
          colorSpace = 'DeviceCMYK'
        imageFile.seek(0)  #reset file pointer
        imagedata = []
        imagedata.append('BI')  # begin image
        # this describes what is in the image itself
        imagedata.append('/Width %0.4f /Height %0.4f' % (info[0], info[1]))
        imagedata.append('/BitsPerComponent 8')
        imagedata.append('/ColorSpace /%s' % colorSpace)
        imagedata.append('/Filter [ /ASCII85Decode /DCTDecode]')
        imagedata.append('ID')
        #write in blocks of (??) 60 characters per line to a list
        compressed = imageFile.read()
        encoded = pdfutils._AsciiBase85Encode(compressed)
        outstream = StringIO(encoded)
        dataline = outstream.read(60)
        while dataline != "":
          imagedata.append(dataline)
          dataline = outstream.read(60)
        imagedata.append('EI')
      else:
        if not pdfutils.cachedImageExists(image):
          pdfutils.cacheImageFile(image)
        #now we have one cached, slurp it in
        cachedname = os.path.splitext(image)[0] + '.a85'
        imagedata = open(cachedname, 'rb').readlines()
        #trim off newlines...
        imagedata = [s.strip() for s in imagedata]

        #parse line two for width, height
        words = imagedata[1].split()
        imgwidth = int(words[1])
        imgheight = int(words[3])
    else:
      #PIL Image
      #work out all dimensions
      myimage = image.convert('RGB')
      imgwidth, imgheight = myimage.size
      imagedata = []
      imagedata.append('BI')  # begin image

      # this describes what is in the image itself
      imagedata.append('/W %0.4f /H %0.4f /BPC 8 /CS /RGB /F [/A85 /Fl]' % (imgwidth, imgheight))
      imagedata.append('ID')

      #use a flate filter and Ascii Base 85 to compress
      raw = myimage.tobytes()
      assert len(raw) == imgwidth * imgheight, "Wrong amount of data for image"
      compressed = zlib.compress(raw)  #this bit is very fast...
      encoded = pdfutils._AsciiBase85Encode(compressed)  #...sadly this isn't

      #write in blocks of (??) 60 characters per line to a list
      outstream = StringIO(encoded)
      dataline = outstream.read(60)
      while dataline != "":
        imagedata.append(dataline)
        dataline = outstream.read(60)
      imagedata.append('EI')

    #now build the PDF for the image.
    if not width:
      width = imgwidth
    if not height:
      height = imgheight

    # this says where and how big to draw it
    #self._code.append('ET')
    #self._code.append('q %0.4f 0 0 %0.4f %0.4f %0.4f cm' % (width, height, x, y+height))
    if self.bottomup:
      self._code.append('q %0.4f 0 0 %0.4f %0.4f %0.4f cm' % (width, height, x, y))
    else:
      # multiply  height by (-1) to overcome flip in coordinate system -cwl
      self._code.append('q %0.4f 0 0 %0.4f %0.4f %0.4f cm' % (width, -height, x, y + height))
    self._code.extend(imagedata)
    self._code.append('Q')
    #self._code.append('BT')

  #########################################################################
  #
  #  JPEG processing code - contributed by Eric Johnson
  #
  #########################################################################

  # Read data from the JPEG file. We should probably be using PIL to
  # get this information for us -- but this way is more fun!
  # Returns (width, height, color components) as a triple
  # This is based on Thomas Merz's code from GhostScript (viewjpeg.ps)

  def readJPEGInfo(self, image):
    "Read width, height and number of components from JPEG file"
    import struct

    #Acceptable JPEG Markers:
    #  SROF0=baseline, SOF1=extended sequential or SOF2=progressive
    validMarkers = [0xC0, 0xC1, 0xC2]

    #JPEG markers without additional parameters
    noParamMarkers = \
        [ 0xD0, 0xD1, 0xD2, 0xD3, 0xD4, 0xD5, 0xD6, 0xD7, 0xD8, 0x01 ]

    #Unsupported JPEG Markers
    unsupportedMarkers = \
        [ 0xC3, 0xC5, 0xC6, 0xC7, 0xC8, 0xC9, 0xCA, 0xCB, 0xCD, 0xCE, 0xCF ]

    #read JPEG marker segments until we find SOFn marker or EOF
    done = 0
    while not done:
      x = struct.unpack('B', image.read(1))
      if x[0] == 0xFF:  #found marker
        x = struct.unpack('B', image.read(1))
        #print "Marker: ", '%0.2x' % x[0]
        #check marker type is acceptable and process it
        if x[0] in validMarkers:
          image.seek(2, 1)  #skip segment length
          x = struct.unpack('B', image.read(1))  #data precision
          if x[0] != 8:
            raise PDFError(' JPEG must have 8 bits per component')
          y = struct.unpack('BB', image.read(2))
          height = (y[0] << 8) + y[1]
          y = struct.unpack('BB', image.read(2))
          width = (y[0] << 8) + y[1]
          y = struct.unpack('B', image.read(1))
          color = y[0]
          return width, height, color
          done = 1
        elif x[0] in unsupportedMarkers:
          raise PDFError(' Unsupported JPEG marker: {%0.2x}'.format(x[0]))
        elif x[0] not in noParamMarkers:
          #skip segments with parameters
          #read length and skip the data
          x = struct.unpack('BB', image.read(2))
          image.seek((x[0] << 8) + x[1] - 2, 1)

  def setPageCompression(self, onoff=1):
    """Possible values 1 or 0 (1 for 'on' is the default).
        If on, the page data will be compressed, leading to much
        smaller files, but takes a little longer to create the files.
        This applies to all subsequent pages, or until setPageCompression()
        is next called."""
    self._pageCompression = onoff

  def setPageTransition(self, effectname=None, duration=1, direction=0, dimension='H', motion='I'):
    """PDF allows page transition effects for use when giving
        presentations.  There are six possible effects.  You can
        just guive the effect name, or supply more advanced options
        to refine the way it works.  There are three types of extra
        argument permitted, and here are the allowed values:
            direction_arg = [0,90,180,270]
            dimension_arg = ['H', 'V']
            motion_arg = ['I','O'] (start at inside or outside)

        This table says which ones take which arguments:

        PageTransitionEffects = {
            'Split': [direction_arg, motion_arg],
            'Blinds': [dimension_arg],
            'Box': [motion_arg],
            'Wipe' : [direction_arg],
            'Dissolve' : [],
            'Glitter':[direction_arg]
            }
        Have fun!
"""
    if not effectname:
      self._pageTransitionString = ''
      return

    #first check each optional argument has an allowed value
    if direction in [0, 90, 180, 270]:
      direction_arg = '/Di /%d' % direction
    else:
      raise PDFError(' directions allowed are 0,90,180,270')

    if dimension in ['H', 'V']:
      dimension_arg = '/Dm /%s' % dimension
    else:
      raise PDFError('dimension values allowed are H and V')

    if motion in ['I', 'O']:
      motion_arg = '/M /%s' % motion
    else:
      raise PDFError('motion values allowed are I and O')

    # this says which effects require which argument types from above
    PageTransitionEffects = {
      'Split': [direction_arg, motion_arg],
      'Blinds': [dimension_arg],
      'Box': [motion_arg],
      'Wipe': [direction_arg],
      'Dissolve': [],
      'Glitter': [direction_arg]
    }

    try:
      args = PageTransitionEffects[effectname]
    except KeyError:
      raise PDFError('Unknown Effect Name "{%s}"'.format(effectname))
      self._pageTransitionString = ''
      return

    self._pageTransitionString = (('/Trans <</D %d /S /%s ' % (duration, effectname)) +
                                  ' '.join(args) + ' >>')


class PDFPathObject:
  """Represents a graphic path.  There are certain 'modes' to PDF
    drawing, and making a separate object to expose Path operations
    ensures they are completed with no run-time overhead.  Ask
    the Canvas for a PDFPath with getNewPathObject(); moveto/lineto/
    curveto wherever you want; add whole shapes; and then add it back
    into the canvas with one of the relevant operators.

    Path objects are probably not long, so we pack onto one line"""

  def __init__(self):
    self._code = []
    self._code.append('n')  #newpath

  def getCode(self):
    "pack onto one line; used internally"
    return ' '.join(self._code)

  def moveTo(self, x, y):
    self._code.append('%0.4f %0.4f m' % (x, y))

  def lineTo(self, x, y):
    self._code.append('%0.4f %0.4f l' % (x, y))

  def curveTo(self, x1, y1, x2, y2, x3, y3):
    self._code.append('%0.4f %0.4f %0.4f %0.4f %0.4f %0.4f c' % (x1, y1, x2, y2, x3, y3))

  def arc(self, x1, y1, x2, y2, startAng=0, extent=90):
    """Contributed to piddlePDF by Robert Kern, 28/7/99.
        Draw a partial ellipse inscribed within the rectangle x1,y1,x2,y2,
        starting at startAng degrees and covering extent degrees.   Angles
        start with 0 to the right (+x) and increase counter-clockwise.
        These should have x1<x2 and y1<y2.

        The algorithm is an elliptical generalization of the formulae in
        Jim Fitzsimmon's TeX tutorial <URL: http://www.tinaja.com/bezarc1.pdf>."""

    pointList = pdfgeom.bezierArc(x1, y1, x2, y2, startAng, extent)
    #move to first point
    self._code.append('%0.4f %0.4f m' % pointList[0][:2])
    for curve in pointList:
      self._code.append('%0.4f %0.4f %0.4f %0.4f %0.4f %0.4f c' % curve[2:])

  def arcTo(self, x1, y1, x2, y2, startAng=0, extent=90):
    """Like arc, but draws a line from the current point to
        the start if the start is not the current point."""
    pointList = pdfgeom.bezierArc(x1, y1, x2, y2, startAng, extent)
    self._code.append('%0.4f %0.4f l' % pointList[0][:2])
    for curve in pointList:
      self._code.append('%0.4f %0.4f %0.4f %0.4f %0.4f %0.4f c' % curve[2:])

  def rect(self, x, y, width, height):
    """Adds a rectangle to the path"""
    self._code.append('%0.4f %0.4f %0.4f %0.4f re' % (x, y, width, height))

  def ellipse(self, x, y, width, height):
    """adds an ellipse to the path"""
    pointList = pdfgeom.bezierArc(x, y, x + width, y + height, 0, 360)
    self._code.append('%0.4f %0.4f m' % pointList[0][:2])
    for curve in pointList:
      self._code.append('%0.4f %0.4f %0.4f %0.4f %0.4f %0.4f c' % curve[2:])

  def circle(self, x_cen, y_cen, r):
    """adds a circle to the path"""
    x1 = x_cen - r
    x2 = x_cen + r
    y1 = y_cen - r
    y2 = y_cen + r
    self.ellipse(x_cen - r, y_cen - r, x_cen + r, y_cen + r)

  def close(self):
    "draws a line back to where it started"
    self._code.append('h')


class PDFTextObject:
  """PDF logically separates text and graphics drawing; you can
    change the coordinate systems for text and graphics independently.
    If you do drawings while in text mode, they appear in the right places
    on the page in Acrobat Reader, bur when you export Postscript to
    a printer the graphics appear relative to the text coordinate
    system.  I regard this as a bug in how Acrobat exports to PostScript,
    but this is the workaround.  It forces the user to separate text
    and graphics.  To output text, ask te canvas for a text object
    with beginText(x, y).  Do not construct one directly. It keeps
    track of x and y coordinates relative to its origin."""

  def __init__(self, canvas, x=0, y=0):
    self._code = []
    self._code.append('BT')
    self._canvas = canvas  #canvas sets this so it has access to size info
    self._fontname = self._canvas._fontname
    self._fontsize = self._canvas._fontsize
    self._leading = self._canvas._leading

    self.setTextOrigin(x, y)

  def getCode(self):
    "pack onto one line; used internally"
    self._code.append('ET')
    return ' '.join(self._code)

  def setTextOrigin(self, x, y):
    if self._canvas.bottomup:
      self._code.append('1 0 0 1 %0.4f %0.4f Tm' % (x, y))  #bottom up
    else:
      self._code.append('1 0 0 -1 %0.4f %0.4f Tm' % (x, y))  #top down
    self._x = x
    self._y = y
    self._x0 = x  #the margin

  def setTextTransform(self, a, b, c, d, e, f):
    "Like setTextOrigin, but does rotation, scaling etc."
    # flip "y" coordinate for top down coordinate system -cwl
    # (1  0)   (a  b)      ( a   b)
    # (0 -1)   (c  d)   =  (-c  -d)
    self._code.append('%0.4f %0.4f %0.4f %0.4f %0.4f %0.4f Tm' % (a, b, -c, -d, e, f))  #top down

    #self._code.append('%0.4f %0.4f %0.4f %0.4f %0.4f %0.4f Tm' % (a, b, c, d, e, f)) #bottom up
    #we only measure coords relative to present text matrix
    self._x = e
    self._y = f

  def moveCursor(self, dx, dy):
    """Moves to a point dx, dy away from the start of the
        current line - NOT from the current point! So if
        you call it in mid-sentence, watch out."""
    self._code.append('%s %s Td' % (dx, -dy))

  def getCursor(self):
    """Returns current text position relative to the last origin."""
    return (self._x, self._y)

  def getX(self):
    """Returns current x position relative to the last origin."""
    return self._x

  def getY(self):
    """Returns current y position relative to the last origin."""
    return self._y

  def setFont(self, psfontname, size, leading=None):
    """Sets the font.  If leading not specified, defaults to 1.2 x
        font size. Raises a readable exception if an illegal font
        is supplied.  Font names are case-sensitive! Keeps track
        of font anme and size for metrics."""
    self._fontname = psfontname
    self._fontsize = size
    pdffontname = self._canvas._doc.getInternalFontName(psfontname)
    if leading is None:
      leading = size * 1.2
    self._leading = leading
    self._code.append('%s %0.1f Tf %0.1f TL' % (pdffontname, size, leading))

  def setCharSpace(self, charSpace):
    """Adjusts inter-character spacing"""
    self._charSpace = charSpace
    self._code.append('%0.4f Tc' % charSpace)

  def setWordSpace(self, wordSpace):
    """Adjust inter-word spacing.  This can be used
        to flush-justify text - you get the width of the
        words, and add some space between them."""
    self._wordSpace = wordSpace
    self._code.append('%0.4f Tw' % wordSpace)

  def setHorizScale(self, horizScale):
    "Stretches text out horizontally"
    self._horizScale = 100 + horizScale
    self._code.append('%0.4f Tz' % horizScale)

  def setLeading(self, leading):
    "How far to move down at the end of a line."
    self._leading = leading
    self._code.append('%0.4f TL' % leading)

  def setTextRenderMode(self, mode):
    """Set the text rendering mode.

        0 = Fill text
        1 = Stroke text
        2 = Fill then stroke
        3 = Invisible
        4 = Fill text and add to clipping path
        5 = Stroke text and add to clipping path
        6 = Fill then stroke and add to clipping path
        7 = Add to clipping path"""

    assert mode in (0, 1, 2, 3, 4, 5, 6, 7), "mode must be in (0,1,2,3,4,5,6,7)"
    self._textRenderMode = mode
    self._code.append('%d Tr' % mode)

  def setRise(self, rise):
    "Move text baseline up or down to allow superscrip/subscripts"
    self._rise = rise
    self._y = self._y - rise  # + ?  _textLineMatrix?
    self._code.append('%0.4f Ts' % rise)

  def setStrokeColorRGB(self, r, g, b):
    self._strokeColorRGB = (r, g, b)
    self._code.append('%0.4f %0.4f %0.4f RG' % (r, g, b))

  def setFillColorRGB(self, r, g, b):
    self._fillColorRGB = (r, g, b)
    self._code.append('%0.4f %0.4f %0.4f rg' % (r, g, b))

  def textOut(self, text):
    "prints string at current point, text cursor moves across"
    text = self._canvas._escape(text)
    self._x = self._x + self._canvas.stringWidth(text, self._fontname, self._fontsize)
    self._code.append('(%s) Tj' % text)

  def textLine(self, text=''):
    """prints string at current point, text cursor moves down.
        Can work with no argument to simply move the cursor down."""
    text = self._canvas._escape(text)
    self._x = self._x0
    if self._canvas.bottomup:
      self._y = self._y - self._leading
    else:
      self._y = self._y + self._leading
    self._code.append('(%s) Tj T*' % text)

  def textLines(self, stuff, trim=1):
    """prints multi-line or newlined strings, moving down.  One
        common use is to quote a multi-line block in your Python code;
        since this may be indented, by default it trims whitespace
        off each line and from the beginning; set trim=0 to preserve
        whitespace."""
    if isinstance(stuff, str):
      lines = stuff.strip().split('\n')
      if trim == 1:
        lines = [s.strip() for s in lines]
    elif isinstance(stuff, (tuple, list)):
      lines = stuff
    else:
      raise ValueError("argument to textlines must be string, list or tuple")

    for line in lines:
      escaped_text = self._canvas._escape(line)
      self._code.append('(%s) Tj T*' % escaped_text)
      if self._canvas.bottomup:
        self._y = self._y - self._leading
      else:
        self._y = self._y + self._leading
    self._x = self._x0


if __name__ == '__main__':
  print('For test scripts, run testpdfgen.py')
