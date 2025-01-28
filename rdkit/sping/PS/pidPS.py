"""
piddlePS - a PostScript backend for the PIDDLE drawing module

   Magnus Lie Hetland

   1999
"""
# some fixups by Chris Lee (cwlee@artsci.wustl.edu)
# help from Paul S. and Bernhard Herzog and others

# $Id$
#
# For each page, the coordinate system is initialized with "0 canvasHeight translate" so
# that coordinate (0,0) is at the top left of the page
# Therefore all y coordinates must be entered in opposite sign to go down
# Also, remember that angles are reversed relative to postscript standard -cwl

# TODO for version 1.0:
#     X add greyscale image support for smaller images (added to level 2)
#     X add Level2 postscript options
#     X add underline implementation from piddlePDF
#     X ** fix \n breaking in drawString
#     ? Bezier curve thing (Robert Kern's fix is in piddlePDF)
#          Hmm..drawCurve(..) already uses postscript's bezier curveto-- is there anything else
#          to be done?
#     X drawArc stuff from Eric

# In the Future:
#    _ Base85 encoding just use hex encoding involves 1:2 expansion of image data vs 4:5
#    _ Don't see a flate/deflate filter for Postscript, jpeg DCTEncode could be added.
#         PIL may have a LZW encoder
#    _ check Adobe Document structuring conventions (half done)
#    _ look at symbol's font metrics--they appear to be a little off
#    X postscript native implementation of drawRoundRect
#    _ improve underlining placement...doesn't look good for courier and symbol

#  DSC: plan uses flags for keeping track of BeginX/EndX pairs.
#            convention: use flag _inXFlag

import math
from io import StringIO

from rdkit.sping.pid import *

from . import psmetrics  # for font info


class PostScriptLevelException(ValueError):
  pass


linesep = '\n'
### constants for fonts ###

# This is actually a mapping between legal font names and PSFontMapXXX keys
# note: all piddle font names are lower-cased before running against this mapping)
PiddleLegalFonts = {
  "helvetica": "helvetica",  # note: keys are lowercased
  "times": "times",
  "courier": "courier",
  "serif": "times",
  "sansserif": "helvetica",
  "monospaced": "courier",
  "symbol": "symbol"
}  # Could add more...

Roman = "Roman"
Bold = "Bold"
Italic = "Italic"

# This is starting to look like a class
PSFontMapStdEnc = {
  ("helvetica", Roman): "Helvetica-Roman",
  ("helvetica", Bold): "Helvetica-Bold",
  ("helvetica", Italic): "Helvetica-Oblique",
  ("times", Roman): "Times-Roman",
  ("times", Bold): "Times-Bold",
  ("times", Italic): "Times-Italic",
  ("courier", Roman): "Courier-Roman",
  ("courier", Bold): "Courier-Bold",
  ("courier", Italic): "Courier-Oblique",
  ("symbol", Roman): "Symbol",
  ("symbol", Bold): "Symbol",
  ("symbol", Italic): "Symbol",
  "EncodingName": 'StandardEncoding'
}

PSFontMapLatin1Enc = {
  ("helvetica", Roman): "Helvetica-Roman-ISOLatin1",
  ("helvetica", Bold): "Helvetica-Bold-ISOLatin1",
  ("helvetica", Italic): "Helvetica-Oblique-ISOLatin1",
  ("times", Roman): "Times-Roman-ISOLatin1",
  ("times", Bold): "Times-Bold-ISOLatin1",
  ("times", Italic): "Times-Italic-ISOLatin1",
  ("courier", Roman): "Courier-Roman-ISOLatin1",
  ("courier", Bold): "Courier-Bold-ISOLatin1",
  ("courier", Italic): "Courier-Oblique-ISOLatin1",
  ("symbol", Roman): "Symbol",
  ("symbol", Bold): "Symbol",
  ("symbol", Italic): "Symbol",
  "EncodingName": 'Latin1Encoding'
}

###############################################################


def latin1FontEncoding(fontname):
  """use this to generating PS code for re-encoding a font as ISOLatin1
    from font with name 'fontname' defines reencoded font, 'fontname-ISOLatin1'"""

  latin1FontTemplate = """/%s findfont
dup length dict begin
  {1 index /FID ne
        {def}
        {pop pop}
      ifelse
   } forall
   /Encoding ISOLatin1Encoding  def
   currentdict
end
/%s-ISOLatin1 exch definefont pop
"""
  #
  return latin1FontTemplate % (fontname, fontname)


def dashLineDefinition():
  res = r"""
%% This is hacked straight out of the Blue book (from Adobe)
/centerdash
  { /pattern exch def
    /pathlen pathlength def
    /patternlength 0 def
    pattern
      { patternlength add /patternlength exch def
      } forall
    pattern length 2 mod 0 ne
      { /patternlength patternlength 2 mul def } if
    /first pattern 0 get def
    /last patternlength first sub def
    /n pathlen last sub cvi patternlength idiv def
    /endpart pathlen patternlength n mul sub
       last sub 2 div def
    /offset first endpart sub def
    pattern offset setdash
  } def

/pathlength
    { flattenpath
      /dist 0 def

      { /yfirst exch def /xfirst exch def
        /ymoveto yfirst def /xmoveto xfirst def }
      { /ynext exch def /xnext exch def
        /dist dist ynext yfirst sub dup mul
          xnext xfirst sub dup mul add sqrt add def
        /yfirst ynext def /xfirst xnext def }
      {}

      { /ynext ymoveto def /xnext xmoveto def
        /dist dist ynext yfirst sub dup mul
          xnext xfirst sub dup mul add sqrt add def
        /yfirst ynext def /xfirst xnext def }
      pathforall
      dist
    } def
"""
  return res


class PsDSC:

  # remember %% will be reduced to % when using string substitution
  # returned strings do not end with \n

  def __init__(self):
    pass

  ## Genral DSC conventions
  def documentHeader(self):
    return "%!PS-Adobe-3.0"

  def boundingBoxStr(self, x0, y0, x1, y1):
    "coordinates of bbox in default PS coordinates"
    return "%%BoundingBox: " + "%s %s %s %s" % (x0, y0, x1, y1)

  def BeginPageStr(self, pageSetupStr, pageName=None):
    """Use this at the beginning of each page, feed it your setup code
        in the form of a string of postscript.  pageName is the "number" of the
        page.  By default it will be 0."""

    self.inPageFlag = 1  # keep track

    if not pageName:
      pageDeclaration = r"%%Page: %d %d" % (1, 1)  # default 1
    else:
      pageDeclaration = "%%Page: " + pageName

    ret = pageDeclaration + "\n" + \
"""%%BeginPageSetup
/pgsave save def
"""
    # print pageSetupStr ???
    return ret + pageSetupStr + "\n%%EndPageSetup"

  def EndPageStr(self):
    self.inPageFlag = 0
    return ""


class EpsDSC(PsDSC):

  def __init__(self):
    PsDSC.__init__(self)

  ## Genral DSC conventions
  def documentHeader(self):
    return "%!PS-Adobe-3.0 EPSF-3.0"


##########################################################################


class PSCanvas(Canvas):
  """This canvas is meant for generating encapsulated PostScript files
    (EPS) used for inclusion in other documents; thus really only
    single-page documents are supported.  For historical reasons and
    because they can be printed (a showpage is included), the files are
    given a .ps extension by default, and a primitive sort of multipage
    document can be generated using nextPage() or clear().  Use at your own
    risk!  Future versions of piddlePS will include an EPSCanvas and a
    PSCanvas which will clearly delineate between single and multipage
    documents.

    Note: All font encodings must be taken care in __init__, you can't add
          more after this"""

  def __init__(self, size=(300, 300), name='piddlePS', PostScriptLevel=2,
               fontMapEncoding=PSFontMapLatin1Enc):

    Canvas.__init__(self, size, name)
    width, height = self.size = size
    self.filename = name
    if len(name) < 3 or name[-3:].lower() != '.ps':
      self.filename = name + ".ps"

    # select between postscript level 1 or level 2
    if PostScriptLevel == 1:
      self.drawImage = self._drawImageLevel1
    elif PostScriptLevel == 2:
      self.drawImage = self._drawImageLevel2
    else:
      raise PostScriptLevelException

    self.code = []
    self.dsc = PsDSC()  # handle document structing conventions

    c = self._currentColor = self.defaultLineColor
    r, g, b = c.red, c.green, c.blue

    w = self._currentWidth = self.defaultLineWidth

    self.defaultFont = Font(face='serif')
    self.fontMapEncoding = fontMapEncoding

    self._currentFont = self.defaultFont
    f = self._findFont(self._currentFont)
    s = self._currentFont.size

    # Page Structure State
    #----------------------
    self._inDocumentFlag = 0  # this is set in psBeginDocument
    self._inPageFlag = 0  # we haven't started a page

    self.pageNum = 1  # User is free to reset this or even make this a string

    self.psBeginDocument()

    ## finally begin the first page ##
    self.psBeginPage()  # each psBeginPage() needs to be closed w/ a psEndPage()

  def psBeginDocument(self):
    # General DSC Prolog := <header> [<defaults>] <procedures>
    self.code.append(self.dsc.documentHeader())
    self.code.append(self.dsc.boundingBoxStr(0, 0, self.size[0], self.size[1]))
    self.code.append("%%Pages: (atend)")
    self._inDocumentFlag = 1  # we need a Trailer to fix this up

    ######  Defaults Procedures Prolog Setup  ??? ######
    # are these procedures??? check this chris

    # Now create Latin1ISO font encodings for non-symbol fonts (need to add pdf fonts too)
    shapes = {
      "Helvetica": ["Roman", "Bold", "Oblique"],
      "Times": ["Roman", "Bold", "Italic"],
      "Courier": ["Roman", "Bold", "Oblique"]
    }
    fntnames = []

    for basename in ['Helvetica', 'Times', 'Courier']:
      for mys in shapes[basename]:
        fntnames.append(basename + '-' + mys)

    # assign the default FontMapping (which also determines the encoding)

    for fontname in fntnames:
      self.code.append(latin1FontEncoding(fontname))

    self.code.append(dashLineDefinition())

  def psEndDocument(self):

    if self._inDocumentFlag:
      # Take care of Trailer
      self.code.append("%%Trailer")
      self.code.append("%%%%Pages: %d" % self.pageNum)

      # signal end of file
      # check on need for %%EOF
      self.code.append("%%EOF")  # remove \n at end of EOF

  def psBeginPage(self, pageName=None):
    # call this function when beginning a new page but before any piddle drawing commands
    # pagesetup contains code to insert into DSC page 'header'
    if not pageName:
      pageName = "%s %s" % (self.pageNum, self.pageNum)
    pagesetup = self._psPageSetupStr(self.size[1], self.defaultLineColor,
                                     self._findFont(self.defaultFont), self.defaultFont.size,
                                     self.defaultLineWidth)
    self.code.append(self.dsc.BeginPageStr(pageSetupStr=pagesetup, pageName=pageName))
    self._inPageFlag = 1

  def _psPageSetupStr(self, pageheight, initialColor, font_family, font_size, line_width):
    "ps code for settin up coordinate system for page in accords w/ piddle standards"
    r, g, b = initialColor.red, initialColor.green, initialColor.blue
    return '''
%% initialize

2 setlinecap

0 %d
translate

%s %s %s setrgbcolor
(%s) findfont %s scalefont setfont
%s setlinewidth''' % (pageheight, r, g, b, font_family, font_size, line_width)

  def psEndPage(self):
    self.code.append("pgsave restore")
    self.code.append("showpage")
    self._inPageFlag = 0

  def _findFont(self, font):

    requested = font.face or "Serif"  # Serif is the default
    if isinstance(requested, str):
      requested = [requested]

    # once again, fall back to default, redundant, no?
    face = PiddleLegalFonts["serif"].lower()
    for reqFace in requested:
      if reqFace.lower() in PiddleLegalFonts:
        face = PiddleLegalFonts[reqFace.lower()].lower()
        break

    if font.bold:
      shape = Bold
    elif font.italic:
      shape = Italic
    else:
      shape = Roman

    return self.fontMapEncoding[(face, shape)]

  def _findExternalFontName(self, font):  #copied from piddlePDF by cwl- hack away!
    """Attempts to return proper font name.
        PDF uses a standard 14 fonts referred to
        by name. Default to self.defaultFont('Helvetica').
        The dictionary allows a layer of indirection to
        support a standard set of PIDDLE font names."""

    piddle_font_map = {
      'Times': 'Times',
      'times': 'Times',
      'Courier': 'Courier',
      'courier': 'Courier',
      'helvetica': 'Helvetica',
      'Helvetica': 'Helvetica',
      'symbol': 'Symbol',
      'Symbol': 'Symbol',
      'monospaced': 'Courier',
      'serif': 'Times',
      'sansserif': 'Helvetica',
      'ZapfDingbats': 'ZapfDingbats',
      'zapfdingbats': 'ZapfDingbats',
      'arial': 'Helvetica'
    }

    try:
      face = piddle_font_map[font.facereqFace.lower()]
    except Exception:
      return 'Helvetica'

    name = face + '-'
    if font.bold and face in ['Courier', 'Helvetica', 'Times']:
      name = name + 'Bold'
    if font.italic and face in ['Courier', 'Helvetica']:
      name = name + 'Oblique'
    elif font.italic and face == 'Times':
      name = name + 'Italic'

    if name == 'Times-':
      name = name + 'Roman'
    # symbol and ZapfDingbats cannot be modified!

    #trim and return
    if name[-1] == '-':
      name = name[0:-1]
    return name

  def _psNextPage(self):
    "advance to next page of document.  "
    self.psEndPage()
    self.pageNum = self.pageNum + 1
    self.psBeginPage()

  def nextPage(self):
    self.clear()

  def clear(self):
    """clear resets the canvas to it's default state.  Though this
        canvas is really only meant to be an EPS canvas, i.e., single page,
        for historical reasons we will allow multipage documents.  Thus
        clear will end the page, clear the canvas state back to default,
        and start a new page.  In the future, this PSCanvas will become
        EPSCanvas and will not support multipage documents.  In that case,
        the canvas will be reset to its default state and the file will be
        emptied of all previous drawing commands"""

    self.resetToDefaults()
    self._psNextPage()

  def resetToDefaults(self):
    self._currentColor = self.defaultLineColor
    self._currentWidth = self.defaultLineWidth
    self._currentFont = self.defaultFont

  def flush(self):
    pass

    # Comment the following out to make flush() a null function
    # file = open(self.filename,'w')
    # from string import join
    # file.write(join(self.code,linesep))
    # file.write('\nshowpage\n')
    # file.close()

  def save(self, file=None, format=None):
    """Write the current document to a file or stream and close the file
        Computes any final trailers, etc. that need to be done in order to
        produce a well formed postscript file.  At least for now though,
        it still allows you to add to the file after a save by not actually
        inserting the finalization code into self.code

        the format argument is not used"""

    # save() will now become part of the spec.
    file = file or self.filename
    fileobj = getFileObject(file, 'wb')
    fileobj.write(linesep.join(self.code).encode('UTF-8'))
    # here's a hack. we might want to be able to add more after saving so
    # preserve the current code ???
    preserveCode = self.code
    self.code = finalizationCode = [""]

    # might be able to move these things into DSC class & use a save state call
    preserve_inPageFlag = self._inPageFlag
    preserve_inDocumentFlag = self._inDocumentFlag

    # now do finalization code in via self.code
    # first question: are we in the middle of a page?

    if self._inPageFlag:
      self.psEndPage()

    self.psEndDocument()  # depends on _inDocumentFlag :(

    fileobj.write(linesep.join(finalizationCode).encode('UTF-8'))
    fileobj.close()
    #  fileobj.close()  ### avoid this for now
    ## clean up my mess: This is not a good way to do things FIXME!!! ???
    self.code = preserveCode
    self._inPageFlag = preserve_inPageFlag
    self._inDocumentFlag = preserve_inDocumentFlag

  def stringWidth(self, s, font=None):
    "Return the logical width of the string if it were drawn \
        in the current font (defaults to self.font)."

    if not font:
      font = self.defaultFont
    fontname = self._findExternalFontName(font)
    return psmetrics.stringwidth(s, fontname,
                                 self.fontMapEncoding["EncodingName"]) * font.size * 0.001

  ### def fontHeight(self, font=None) ### use default piddle method
  # distance between baseline of two lines of text 1.2 * font.size for piddle.py
  # Thus is 1.2 times larger than postscript minimum baseline to baseline separation
  # for single-spaced text.  fontHeight() used internally for drawString() multi-line text

  def fontAscent(self, font=None):
    if not font:
      font = self.defaultFont
    fontname = self._findExternalFontName(font)
    return psmetrics.ascent_descent[fontname][0] * 0.001 * font.size

  def fontDescent(self, font=None):
    if not font:
      font = self.defaultFont
    fontname = self._findExternalFontName(font)
    return -psmetrics.ascent_descent[fontname][1] * 0.001 * font.size

  def _updateLineColor(self, color):
    color = color or self.defaultLineColor
    if color != self._currentColor:
      self._currentColor = color
      r, g, b = color.red, color.green, color.blue
      self.code.append('%s %s %s setrgbcolor' % (r, g, b))

  def _updateFillColor(self, color):
    color = color or self.defaultFillColor
    if color != self._currentColor:
      self._currentColor = color
      r, g, b = color.red, color.green, color.blue
      self.code.append('%s %s %s setrgbcolor' % (r, g, b))

  def _updateLineWidth(self, width):
    if width is None:
      width = self.defaultLineWidth
    if width != self._currentWidth:
      self._currentWidth = width
      self.code.append('%s setlinewidth' % width)

  def _updateFont(self, font):
    font = font or self.defaultFont
    if font != self._currentFont:
      self._currentFont = font
      f = self._findFont(font)
      s = font.size
      self.code.append('(%s) findfont %s scalefont setfont' % (f, s))

  def drawLine(self, x1, y1, x2, y2, color=None, width=None, dash=None, **kwargs):
    self._updateLineColor(color)
    self._updateLineWidth(width)
    if dash is not None:
      dashTxt = '[' + '%d ' * len(dash) + ']'
      dashTxt = dashTxt % dash
      self.code.append('%s centerdash' % dashTxt)
    if self._currentColor != transparent:
      self.code.append('%s %s neg moveto %s %s neg lineto stroke' % (x1, y1, x2, y2))
    if dash is not None:
      self.code.append('[1 0] centerdash')

  def drawLines(self, lineList, color=None, width=None, dash=None, **kwargs):
    self._updateLineColor(color)
    self._updateLineWidth(width)
    codeline = '%s %s neg moveto %s %s neg lineto stroke'
    if self._currentColor != transparent:
      for line in lineList:
        self.code.append(codeline % line)

  def _escape(self, s):
    # return a copy of string s with special characters in postscript strings
    # escaped" with backslashes."""
    # Have not handled characters that are converted normally in python strings
    # i.e. \n -> newline
    str = s.replace(chr(0x5C), r'\\')
    str = str.replace('(', r'\(')
    str = str.replace(')', r'\)')
    return str

  # ??? check to see if \n response is handled correctly (should move cursor down)

  def _drawStringOneLineNoRot(self, s, x, y, font=None, **kwargs):
    # PRE: x and y and position at which to place text
    # PRE: helper function, only called from drawString(..)
    text = self._escape(s)
    self.code.append('%s %s neg moveto (%s) show' % (x, y, text))
    if self._currentFont.underline:
      swidth = self.stringWidth(s, self._currentFont)
      ypos = (0.5 * self.fontDescent(self._currentFont))
      thickness = 0.08 * self._currentFont.size  # relate to font.descent?
      self.code.extend([
        '%s setlinewidth' % thickness,
        '0 %s neg rmoveto' % (ypos),
        '%s 0 rlineto stroke' % -swidth
      ])

  def _drawStringOneLine(self, s, x, y, font=None, color=None, angle=0, **kwargs):
    # PRE: assumes that coordinate system has already been set for rotated strings
    # PRE: only meant to be called by drawString(..)
    text = self._escape(s)
    self.code.extend(['%f %f neg moveto (%s) show' % (x, y, text)])

    if self._currentFont.underline:
      swidth = self.stringWidth(s, self._currentFont)
      dy = (0.5 * self.fontDescent(self._currentFont))
      thickness = 0.08 * self._currentFont.size  # relate to font.descent?
      self.code.extend([
        '%s setlinewidth' % thickness,
        '%f %f neg moveto' % (x, dy + y),
        '%f 0 rlineto stroke' % swidth
      ])

  def drawString(self, s, x, y, font=None, color=None, angle=0, **kwargs):
    """drawString(self, s, x, y, font=None, color=None, angle=0)
        draw a string s at position x,y"""
    self._updateLineColor(color)
    self._updateFont(font)
    if self._currentColor != transparent:

      lines = s.split('\n')
      lineHeight = self.fontHeight(font)

      if angle == 0:  # do special case of angle = 0 first. Avoids a bunch of gsave/grestore ops
        for line in lines:
          self._drawStringOneLineNoRot(line, x, y, font=font, **kwargs)
      else:  # general case, rotated text
        self.code.extend(['gsave', '%s %s neg translate' % (x, y), repr(angle) + ' rotate'])
        down = 0
        for line in lines:
          self._drawStringOneLine(line, 0, 0 + down, font, color, angle, **kwargs)
          down = down + lineHeight
        self.code.extend(['grestore'])

  def drawCurve(self, x1, y1, x2, y2, x3, y3, x4, y4, edgeColor=None, edgeWidth=None,
                fillColor=None, closed=0, dash=None, **kwargs):
    codeline = '%s %s neg moveto %s %s neg %s %s neg %s %s neg curveto'
    data = (x1, y1, x2, y2, x3, y3, x4, y4)
    self._updateFillColor(fillColor)
    if self._currentColor != transparent:
      self.code.append((codeline % data) + ' eofill')
    self._updateLineWidth(edgeWidth)
    self._updateLineColor(edgeColor)
    if self._currentColor != transparent:
      self.code.append((codeline % data) + ((closed and ' closepath') or '') + ' stroke')

  ########################################################################################

  def drawRoundRect(self, x1, y1, x2, y2, rx=8, ry=8, edgeColor=None, edgeWidth=None,
                    fillColor=None, dash=None, **kwargs):
    "Draw a rounded rectangle between x1,y1, and x2,y2, \
        with corners inset as ellipses with x radius rx and y radius ry. \
        These should have x1<x2, y1<y2, rx>0, and ry>0."

    # Path is drawn in counter-clockwise direction"

    x1, x2 = min(x1, x2), max(x1, x2)  # from piddle.py
    y1, y2 = min(y1, y2), max(y1, y2)

    # Note: arcto command draws a line from current point to beginning of arc
    # save current matrix, translate to center of ellipse, scale by rx ry, and draw
    # a circle of unit radius in counterclockwise dir, return to original matrix
    # arguments are (cx, cy, rx, ry, startAngle, endAngle)
    ellipsePath = 'matrix currentmatrix %s %s neg translate %s %s scale 0 0 1 %s %s arc setmatrix'

    # choice between newpath and moveto beginning of arc
    # go with newpath for precision, does this violate any assumptions in code???
    # rrcode = ['%s %s neg moveto' % (x1+rx, y1)]  # this also works
    rrcode = ['newpath']  # Round Rect code path
    # upper left corner ellipse is first
    rrcode.append(ellipsePath % (x1 + rx, y1 + ry, rx, ry, 90, 180))
    rrcode.append(ellipsePath % (x1 + rx, y2 - ry, rx, ry, 180, 270))
    rrcode.append(ellipsePath % (x2 - rx, y2 - ry, rx, ry, 270, 360))
    rrcode.append(ellipsePath % (x2 - rx, y1 + ry, rx, ry, 0, 90))
    rrcode.append('closepath')

    # This is what you are required to do to take care of all color cases
    # should fix this so it doesn't define path twice, just use gsave if need
    # to fill and stroke path-need to figure out this system
    self._updateFillColor(fillColor)
    if self._currentColor != transparent:
      self.code.extend(rrcode)
      self.code.append("eofill")
    self._updateLineWidth(edgeWidth)
    self._updateLineColor(edgeColor)
    if self._currentColor != transparent:
      self.code.extend(rrcode)
      self.code.append("stroke")

  def drawEllipse(self, x1, y1, x2, y2, edgeColor=None, edgeWidth=None, fillColor=None, dash=None,
                  **kwargs):
    "Draw an orthogonal ellipse inscribed within the rectangle x1,y1,x2,y2. \
        These should have x1<x2 and y1<y2."

    #Just invoke drawArc to actually draw the ellipse
    self.drawArc(x1, y1, x2, y2, edgeColor=edgeColor, edgeWidth=edgeWidth, fillColor=fillColor)

  def drawArc(self, x1, y1, x2, y2, startAng=0, extent=360, edgeColor=None, edgeWidth=None,
              fillColor=None, dash=None, **kwargs):
    "Draw a partial ellipse inscribed within the rectangle x1,y1,x2,y2, \
        starting at startAng degrees and covering extent degrees.   Angles \
        start with 0 to the right (+x) and increase counter-clockwise. \
        These should have x1<x2 and y1<y2."

    #calculate centre of ellipse
    cx, cy = (x1 + x2) / 2.0, (y1 + y2) / 2.0
    rx, ry = (x2 - x1) / 2.0, (y2 - y1) / 2.0

    codeline = self._genArcCode(x1, y1, x2, y2, startAng, extent)

    # fill portion
    self._updateFillColor(fillColor)

    if self._currentColor != transparent:
      self.code.append('%s %s neg moveto' % (cx, cy))  # move to center of circle
      self.code.append(codeline + ' eofill')

    # stroke portion
    self._updateLineWidth(edgeWidth)
    self._updateLineColor(edgeColor)

    if self._currentColor != transparent:
      # move current point to start of arc, note negative angle because y increases down
      self.code.append('%s %s neg moveto' %
                       (cx + rx * math.cos(-startAng), cy + ry * math.sin(-startAng)))
      self.code.append(codeline + ' stroke')

  def _genArcCode(self, x1, y1, x2, y2, startAng, extent):
    "Calculate the path for an arc inscribed in rectangle defined by (x1,y1),(x2,y2)"
    #calculate semi-minor and semi-major axes of ellipse
    xScale = abs((x2 - x1) / 2.0)
    yScale = abs((y2 - y1) / 2.0)
    #calculate centre of ellipse
    x, y = (x1 + x2) / 2.0, (y1 + y2) / 2.0

    codeline = 'matrix currentmatrix '+\
            '%s %s neg translate %s %s scale 0 0 1 %s %s %s '+\
            'setmatrix'

    if extent >= 0:
      arc = 'arc'
    else:
      arc = 'arcn'
    data = (x, y, xScale, yScale, startAng, startAng + extent, arc)

    return codeline % data

  def drawPolygon(self, pointlist, edgeColor=None, edgeWidth=None, fillColor=None, closed=0,
                  dash=None, **kwargs):

    start = pointlist[0]
    pointlist = pointlist[1:]

    polyCode = []
    polyCode.append("%s %s neg moveto" % start)
    for point in pointlist:
      polyCode.append("%s %s neg lineto" % point)
    if closed:
      polyCode.append("closepath")

    self._updateFillColor(fillColor)
    if self._currentColor != transparent:
      self.code.extend(polyCode)
      self.code.append("eofill")
    self._updateLineWidth(edgeWidth)
    self._updateLineColor(edgeColor)
    if self._currentColor != transparent:
      self.code.extend(polyCode)
      self.code.append("stroke")

  def drawFigure(self, partList, edgeColor=None, edgeWidth=None, fillColor=None, closed=0,
                 dash=None, **kwargs):

    figureCode = []
    first = 1

    for part in partList:
      op = part[0]
      args = list(part[1:])

      if op == figureLine:
        if first:
          first = 0
          figureCode.append("%s %s neg moveto" % tuple(args[:2]))
        else:
          figureCode.append("%s %s neg lineto" % tuple(args[:2]))
        figureCode.append("%s %s neg lineto" % tuple(args[2:]))

      elif op == figureArc:
        first = 0
        x1, y1, x2, y2, startAngle, extent = args[:6]
        figureCode.append(self._genArcCode(x1, y1, x2, y2, startAngle, extent))

      elif op == figureCurve:
        if first:
          first = 0
          figureCode.append("%s %s neg moveto" % tuple(args[:2]))
        else:
          figureCode.append("%s %s neg lineto" % tuple(args[:2]))
        figureCode.append("%s %s neg %s %s neg %s %s neg curveto" % tuple(args[2:]))
      else:
        raise TypeError("unknown figure operator: " + op)

    if closed:
      figureCode.append("closepath")

    self._updateFillColor(fillColor)
    if self._currentColor != transparent:
      self.code.extend(figureCode)
      self.code.append("eofill")
    self._updateLineWidth(edgeWidth)
    self._updateLineColor(edgeColor)
    if self._currentColor != transparent:
      self.code.extend(figureCode)
      self.code.append("stroke")

  ############################################################################################
  # drawImage(self. image, x1, y1, x2=None, y2=None) is now defined by either _drawImageLevel1
  #    ._drawImageLevel2, the choice is made in .__init__ depending on option

  def _drawImageLevel1(self, image, x1, y1, x2=None, y2=None, **kwargs):
    # Postscript Level1 version available for fallback mode when Level2 doesn't work
    """drawImage(self,image,x1,y1,x2=None,y2=None) : If x2 and y2 are omitted, they are
       calculated from image size.  (x1,y1) is upper left of image, (x2,y2) is lower right of
       image in piddle coordinates."""
    try:
      from PIL import Image
    except ImportError:
      print('Python Imaging Library not available')
      return
    # For now let's start with 24 bit RGB images (following piddlePDF again)
    print("Trying to drawImage in piddlePS")
    component_depth = 8
    myimage = image.convert('RGB')
    imgwidth, imgheight = myimage.size
    if not x2:
      x2 = imgwidth + x1
    if not y2:
      y2 = y1 + imgheight
    drawwidth = x2 - x1
    drawheight = y2 - y1
    print('Image size (%d, %d); Draw size (%d, %d)' % (imgwidth, imgheight, drawwidth, drawheight))
    # now I need to tell postscript how big image is

    # "image operators assume that they receive sample data from
    # their data source in x-axis major index order.  The coordinate
    # of the lower-left corner of the first sample is (0,0), of the
    # second (1,0) and so on" -PS2 ref manual p. 215
    #
    # The ImageMatrix maps unit square of user space to boundary of the source image
    #

    # The CurrentTransformationMatrix (CTM) maps the unit square of
    # user space to the rect...on the page that is to receive the
    # image. A common ImageMatrix is [width 0 0 -height 0 height]
    # (for a left to right, top to bottom image )

    # first let's map the user coordinates start at offset x1,y1 on page

    self.code.extend([
      'gsave',
      '%s %s translate' % (x1, -y1 - drawheight),  # need to start are lower left of image
      '%s %s scale' % (drawwidth, drawheight),
      '/scanline %d 3 mul string def' % imgwidth  # scanline by multiples of image width
    ])

    # now push the dimensions and depth info onto the stack
    # and push the ImageMatrix to map the source to the target rectangle (see above)
    # finally specify source (PS2 pp. 225 ) and by exmample
    self.code.extend([
      '%s %s %s' % (imgwidth, imgheight, component_depth),
      '[%s %s %s %s %s %s]' % (imgwidth, 0, 0, -imgheight, 0, imgheight),
      '{ currentfile scanline readhexstring pop } false 3', 'colorimage '
    ])

    # data source output--now we just need to deliver a hex encode
    # series of lines of the right overall size can follow
    # piddlePDF again

    rawimage = myimage.tobytes()
    assert len(rawimage) == imgwidth * imgheight, 'Wrong amount of data for image'
    #compressed = zlib.compress(rawimage) # no zlib at moment
    hex_encoded = self._AsciiHexEncode(rawimage)

    # write in blocks of 78 chars per line
    outstream = StringIO(hex_encoded)

    dataline = outstream.read(78)
    while dataline != "":
      self.code.append(dataline)
      dataline = outstream.read(78)
    self.code.append('% end of image data')  # for clarity
    self.code.append('grestore')  # return coordinates to normal

  # end of drawImage

  def _AsciiHexEncode(self, input):  # also based on piddlePDF
    "Helper function used by images"
    output = StringIO()
    for char in input:
      output.write('%02x' % ord(char))
    output.reset()
    return output.read()

  def _drawImageLevel2(self, image, x1, y1, x2=None, y2=None):  # Postscript Level2 version
    try:
      from PIL import Image
    except ImportError:
      print('Python Imaging Library not available')
      return
      # I don't have zlib -cwl
    #         try:
    #             import zlib
    #         except ImportError:
    #             print 'zlib not available'
    #             return

    ### what sort of image are we to draw
    if image.mode == 'L':
      print('found image.mode= L')
      imBitsPerComponent = 8
      imNumComponents = 1
      myimage = image
    elif image.mode == '1':
      print('found image.mode= 1')
      myimage = image.convert('L')
      imNumComponents = 1
      myimage = image
    else:
      myimage = image.convert('RGB')
      imNumComponents = 3
      imBitsPerComponent = 8

    imwidth, imheight = myimage.size
    # print 'imwidth = %s, imheight = %s' % myimage.size
    if not x2:
      x2 = imwidth + x1
    if not y2:
      y2 = y1 + imheight
    drawwidth = x2 - x1
    drawheight = y2 - y1
    self.code.extend([
      'gsave',
      '%s %s translate' % (x1, -y1 - drawheight),  # need to start are lower left of image
      '%s %s scale' % (drawwidth, drawheight)
    ])

    if imNumComponents == 3:
      self.code.append('/DeviceRGB setcolorspace')
    elif imNumComponents == 1:
      self.code.append('/DeviceGray setcolorspace')
      print('setting colorspace gray')
    # create the image dictionary
    self.code.append("""
<<
/ImageType 1
/Width %d /Height %d  %% dimensions of source image
/BitsPerComponent %d""" % (imwidth, imheight, imBitsPerComponent))

    if imNumComponents == 1:
      self.code.append('/Decode [0 1]')
    if imNumComponents == 3:
      self.code.append('/Decode [0 1 0 1 0 1]  %% decode color values normally')

    self.code.extend([
      '/ImageMatrix [%s 0 0 %s 0 %s]' % (imwidth, -imheight, imheight),
      '/DataSource currentfile /ASCIIHexDecode filter', '>> % End image dictionary', 'image'
    ])
    # after image operator just need to dump image dat to file as hexstring
    rawimage = myimage.tobytes()
    assert len(rawimage) == imwidth * imheight, 'Wrong amount of data for image'
    #compressed = zlib.compress(rawimage) # no zlib at moment
    hex_encoded = self._AsciiHexEncode(rawimage)

    # write in blocks of 78 chars per line
    outstream = StringIO(hex_encoded)

    dataline = outstream.read(78)
    while dataline != "":
      self.code.append(dataline)
      dataline = outstream.read(78)
    self.code.append('> % end of image data')  # > is EOD for hex encoded filterfor clarity
    self.code.append('grestore')  # return coordinates to normal
