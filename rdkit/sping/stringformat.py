"""
Module StringFormat
The StringFormat module allows for character-by-character formatting of
strings. It imitates the SPING string drawing and string metrics
interface. The string formatting is done with specialized XML syntax
within the string. Therefore, the interface for the StringFormat module
consists of wrapper functions for the SPING string interface and
various XML tags and characters.

StringFormat functions

       drawString(canvas, s, x, y, [font], [color], [angle])
       stringWidth(canvas, s, [font])
       fontHeight(canvas, [font])
       fontAscent(canvas, [font])
       fontDescent(canvas, [font])
StringFormat XML tags

       <b> </b> - bold
       <i> </i> - italics
       <u> </u> - underline
       <super> </super> - superscript
       <sub> </sub> - subscript

StringFormat XML characters

       Greek Letter Symbols as specified in MathML
"""

#       How it works: Each tag grouping <b></b> sets a flag upon entry and
#       clears the flag upon exit.  Each call to handle_data creates a
#       StringSegment which takes on all of the characteristics specified
#       by flags currently set.  The greek letters can be specified as either
#       &alpha; or <alpha/>.  The are essentially transformed into <alpha/>
#       no matter what and then there is a handler for each greek letter.
#       To add or change greek letter to symbol font mappings only
#       the greekchars map needs to change.

import math

import xmllib

from rdkit.sping.pid import Font

#------------------------------------------------------------------------
# constants
sizedelta = 2  # amount to reduce font size by for super and sub script
subFraction = 0.5  # fraction of font size that a sub script should be lowered
superFraction = 0.5  # fraction of font size that a super script should be raised

#------------------------------------------------------------------------
# greek mapping dictionary

# characters not supported: epsi, Gammad, gammad, kappav, rhov
#    Upsi, upsi

greekchars = {
  'alpha': 'a',
  'beta': 'b',
  'chi': 'c',
  'Delta': 'D',
  'delta': 'd',
  'epsiv': 'e',
  'eta': 'h',
  'Gamma': 'G',
  'gamma': 'g',
  'iota': 'i',
  'kappa': 'k',
  'Lambda': 'L',
  'lambda': 'l',
  'mu': 'm',
  'nu': 'n',
  'Omega': 'W',
  'omega': 'w',
  'omicron': 'x',
  'Phi': 'F',
  'phi': 'f',
  'phiv': 'j',
  'Pi': 'P',
  'pi': 'p',
  'piv': 'v',
  'Psi': 'Y',
  'psi': 'y',
  'rho': 'r',
  'Sigma': 'S',
  'sigma': 's',
  'sigmav': 'V',
  'tau': 't',
  'Theta': 'Q',
  'theta': 'q',
  'thetav': 'j',
  'Xi': 'X',
  'xi': 'x',
  'zeta': 'z'
}


#------------------------------------------------------------------------
class StringSegment:
  """class StringSegment contains the intermediate representation of string
        segments as they are being parsed by the XMLParser.
        """

  def __init__(self):
    self.super = 0
    self.sub = 0
    self.bold = 0
    self.italic = 0
    self.underline = 0
    self.s = ""
    self.width = 0
    self.greek = 0

  def calcNewFont(self, font):
    "Given a font (does not accept font==None), creates a \
                new font based on the format of this text segment."

    # if we are a greek character we need to pick a different fontface
    if self.greek:
      face = "symbol"
    else:
      face = font.face

    # want to make sure that we don't lose any of the base
    # font formatting
    return Font(face=face, size=font.size - (self.super * sizedelta) - (self.sub * sizedelta),
                underline=self.underline or font.underline, bold=self.bold or font.bold,
                italic=self.italic or font.italic)

  def calcNewY(self, font, y):
    "Returns a new y coordinate depending on its \
                whether the string is a sub and super script."

    # should this take into account angle, I think probably not
    if self.sub == 1:
      return y + (font.size * subFraction)
    elif self.super == 1:
      return y - (font.size * superFraction)
    else:
      return y

  def dump(self):
    print("StringSegment: ]%s[" % self.s)
    print("\tsuper = ", self.super)
    print("\tsub = ", self.sub)
    print("\tbold = ", self.bold)
    print("\titalic = ", self.italic)
    print("\tunderline = ", self.underline)
    print("\twidth = ", self.width)
    print("\tgreek = ", self.greek)


#------------------------------------------------------------------
# The StringFormatter will be able to format the following xml
# tags:
#      < b > < /b > - bold
#      < i > < /i > - italics
#      < u > < /u > - underline
#      < super > < /super > - superscript
#      < sub > < /sub > - subscript
#
# It will also be able to handle any MathML specified Greek characters.
#
# Possible future additions: changing color and font
#       character-by-character
#------------------------------------------------------------------
class StringFormatter(xmllib.XMLParser):

  #----------------------------------------------------------
  # First we will define all of the xml tag handler functions.
  #
  # start_<tag>(attributes)
  # end_<tag>()
  #
  # While parsing the xml StringFormatter will call these
  # functions to handle the string formatting tags.
  # At the start of each tag the corresponding field will
  # be set to 1 and at the end tag the corresponding field will
  # be set to 0.  Then when handle_data is called the options
  # for that data will be apparent by the current settings.
  #----------------------------------------------------------

  #### bold
  def start_b(self, attributes):
    self.bold = 1

  def end_b(self):
    self.bold = 0

  #### italics
  def start_i(self, attributes):
    self.italic = 1

  def end_i(self):
    self.italic = 0

  #### underline
  def start_u(self, attributes):
    self.underline = 1

  def end_u(self):
    self.underline = 0

  #### super script
  def start_super(self, attributes):
    self.super = 1

  def end_super(self):
    self.super = 0

  #### sub script
  def start_sub(self, attributes):
    self.sub = 1

  def end_sub(self):
    self.sub = 0

  #### greek script
  def start_greek(self, attributes, letter):
    #               print("creating a greek letter... ", letter)
    self.greek = 1
    self.handle_data(letter)

  def end_greek(self):
    self.greek = 0

  #----------------------------------------------------------------

  def __init__(self):
    xmllib.XMLParser.__init__(self)

    # initialize list of string segments to empty
    self.segmentlist = []

    # initialize tag values
    self.sub = 0
    self.super = 0
    self.bold = 0
    self.italic = 0
    self.underline = 0

    # set up handlers for various tags
    self.elements = {
      'b': (self.start_b, self.end_b),
      'u': (self.start_u, self.end_u),
      'i': (self.start_i, self.end_i),
      'super': (self.start_super, self.end_super),
      'sub': (self.start_sub, self.end_sub)
    }

    # automatically add handlers for all of the greek characters
    for item in greekchars.keys():
      self.elements[item] = (
        lambda attr, self=self, letter=greekchars[item]: self.start_greek(attr, letter),
        self.end_greek)

    # flag for greek characters
    self.greek = 0
    # set up dictionary for greek characters, this is a class variable
    # should I copy it and then update it?
    for item in greekchars.keys():
      self.entitydefs[item] = '<%s/>' % item

  #----------------------------------------------------------------
  #       def syntax_error(self,message):
  #               print(message)

  #----------------------------------------------------------------

  def handle_data(self, data):
    "Creates an intermediate representation of string segments."

    # segment first has data
    segment = StringSegment()
    segment.s = data

    # if sub and super are both one they will cancel each other out
    if self.sub == 1 and self.super == 1:
      segment.sub = 0
      segment.super = 0
    else:
      segment.sub = self.sub
      segment.super = self.super

    # bold, italic, and underline
    segment.bold = self.bold
    segment.italic = self.italic
    segment.underline = self.underline

    # greek character
    segment.greek = self.greek

    self.segmentlist.append(segment)

  #----------------------------------------------------------------
  def parseSegments(self, s):
    "Given a formatted string will return a list of \
                StringSegment objects with their calculated widths."

    # the xmlparser requires that all text be surrounded by xml
    # tags, therefore we must throw some unused flags around the
    # given string
    self.feed("<formattedstring>" + s + "</formattedstring>")
    self.close()  # force parsing to complete
    self.reset()  # get rid of any previous data
    segmentlist = self.segmentlist
    self.segmentlist = []
    return segmentlist

  #------------------------------------------------------------------------
  # These functions just implement an interface layer to SPING


def fontHeight(canvas, font=None):
  "Find the total height (ascent + descent) of the given font."
  return canvas.fontHeight(font)


def fontAscent(canvas, font=None):
  "Find the ascent (height above base) of the given font."
  return canvas.fontAscent(font)


def fontDescent(canvas, font=None):
  "Find the descent (extent below base) of the given font."
  return canvas.fontDescent(font)


#------------------------------------------------------------------------
# create an instantiation of the StringFormatter
#sformatter = StringFormatter()

#------------------------------------------------------------------------
# stringWidth and drawString both have to parse the formatted strings


def stringWidth(canvas, s, font=None):
  "Return the logical width of the string if it were drawn \
        in the current font (defaults to canvas.font)."

  sformatter = StringFormatter()

  segmentlist = sformatter.parseSegments(s)

  # to calculate a new font the segments must be given an actual font
  if not font:
    font = canvas.defaultFont

  # sum up the string widths of each formatted segment
  sum = 0
  for seg in segmentlist:
    sum = sum + canvas.stringWidth(seg.s, seg.calcNewFont(font))
  return sum


def rotateXY(x, y, theta):
  "Rotate (x,y) by theta degrees.  Got transformation \
        from page 299 in linear algebra book."

  radians = theta * math.pi / 180.0
  # had to change the signs to deal with the fact that the y coordinate
  # is positive going down the screen
  return (math.cos(radians) * x + math.sin(radians) * y,
          -(math.sin(radians) * x - math.cos(radians) * y))


def drawString(canvas, s, x, y, font=None, color=None, angle=0):
  "Draw a formatted string starting at location x,y in canvas."
  sformatter = StringFormatter()

  segmentlist = sformatter.parseSegments(s)

  # to calculate a new font the segments must be given an actual font
  if not font:
    font = canvas.defaultFont

  # have each formatted string segment specify its own font
  startpos = x
  for seg in segmentlist:
    # calculate x and y for this segment based on the angle
    # if the string wasn't at an angle then
    # (draw_x,draw_y) = (startpos, seg.calcNewY(font, y)) want to
    # rotate around original x and y
    (delta_x, delta_y) = rotateXY(startpos - x, seg.calcNewY(font, y) - y, angle)

    canvas.drawString(seg.s, x + delta_x, y + delta_y, seg.calcNewFont(font), color, angle)

    # new x start position, startpos is calculated assuming no angle
    startpos = startpos + canvas.stringWidth(seg.s, seg.calcNewFont(font))

  #------------------------------------------------------------------------
  # Testing
  #------------------------------------------------------------------------


from sping.PDF import PDFCanvas


def test1():
  canvas = PDFCanvas('test1.pdf')
  drawString(canvas, "<u><b>hello there</b></u><super>hi</super>", 10, 20)
  drawString(canvas, "hello!", 10, 40)

  print("'hello!' width = ", stringWidth(canvas, "hello!"))
  print("'hello!' SPING width = ", canvas.stringWidth("hello!"))

  drawString(canvas, "<b>hello!</b> goodbye", 10, 60)
  print("'<b>hello!</b> goodbye' width = ", stringWidth(canvas, "<b>hello!</b> goodbye"))
  drawString(canvas, "hello!", 10, 80, Font(bold=1))
  print("'hello!' Font(bold=1) SPING width = ", canvas.stringWidth("hello!", Font(bold=1)))
  drawString(canvas, " goodbye", 10, 100)
  print("' goodbye' SPING width = ", canvas.stringWidth(" goodbye"))
  canvas.flush()


def test2():
  canvas = PDFCanvas('test2.pdf')

  drawString(canvas, "<alpha/>", 10, 10)

  #       drawString(canvas, "&amp;", 10, 10)
  drawString(canvas, "&alpha;", 10, 30)
  #       drawString(canvas, "a", 10, 50, Font(face="symbol"))
  #       drawString(canvas, "hello there!", 30, 90, angle= -90)
  #       drawString(canvas, "<b>goodbye!</b> <u>yall</u>", 100, 90, angle= 45)
  #       drawString(canvas, "there is a <u>time</u> and a <b>place</b><super>2</super>",
  #               100, 90, angle= -75)
  canvas.flush()


def allTagCombos(canvas, x, y, font=None, color=None, angle=0):
  """Try out all tags and various combinations of them.
        Starts at given x,y and returns possible next (x,y)."""

  oldDefault = canvas.defaultFont
  if font:
    canvas.defaultFont = font

  oldx = x
  dx = stringWidth(canvas, " ")
  dy = canvas.defaultFont.size * 1.5

  drawString(canvas, "<b>bold</b>", x, y, color=color, angle=angle)
  x = x + stringWidth(canvas, "<b>bold</b>") + dx

  drawString(canvas, "<i>italic</i>", x, y, color=color, angle=angle)
  x = x + stringWidth(canvas, "<i>italic</i>") + dx

  drawString(canvas, "<u>underline</u>", x, y, color=color, angle=angle)
  x = x + stringWidth(canvas, "<u>underline</u>") + dx

  drawString(canvas, "<super>super</super>", x, y, color=color, angle=angle)
  x = x + stringWidth(canvas, "<super>super</super>") + dx

  drawString(canvas, "<sub>sub</sub>", x, y, color=color, angle=angle)
  y = y + dy

  drawString(canvas, "<b><u>bold+underline</u></b>", oldx, y, color=color, angle=angle)
  x = oldx + stringWidth(canvas, "<b><u>bold+underline</u></b>") + dx

  drawString(canvas, "<super><i>super+italic</i></super>", x, y, color=color, angle=angle)
  x = x + stringWidth(canvas, "<super><i>super+italic</i></super>") + dx

  drawString(canvas, "<b><sub>bold+sub</sub></b>", x, y, color=color, angle=angle)
  #       x = x + stringWidth(canvas,"<b><sub>bold+sub</sub></b>") + dx
  y = y + dy

  canvas.defaultFont = oldDefault
  return (oldx, y)


def stringformatTest():

  # change the following line only to try a different SPING backend
  canvas = PDFCanvas('bigtest1.pdf')

  ################################################### testing drawString tags
  #      < b > < /b > - bold
  #      < i > < /i > - italics
  #      < u > < /u > - underline
  #      < super > < /super > - superscript
  #      < sub > < /sub > - subscript

  x = 10
  y = canvas.defaultFont.size * 1.5

  ##### try out each possible tags and all combos
  (x, y) = allTagCombos(canvas, x, y)

  ##### now try various fonts
  (x, y) = allTagCombos(canvas, x, y + 30, Font(face="serif"))
  (x, y) = allTagCombos(canvas, x, y + 30, Font(face="monospaced"))

  # what about rotated
  (x, y) = allTagCombos(canvas, x, y + 30, Font(face="serif"), angle=-30)

  ##### now try a couple of different font sizes
  (x, y) = allTagCombos(canvas, x, y + 30, Font(size=16))
  (x, y) = allTagCombos(canvas, x, y + 30, Font(size=9))

  ##### now try a different default style setting
  (x, y) = allTagCombos(canvas, x, y + 30, Font(underline=1))

  ##### now try a combo of the above 4 and a different color
  (x, y) = allTagCombos(canvas, x, y + 30, color=red)

  ################################################### testing stringWidth tags
  sfwidth = stringWidth(canvas,
                        "<b><sub>bold+sub</sub></b> hello <u><super>underline+super</super></u>")

  # break down the various string widths
  print('sw("<b><sub>bold+sub</sub></b>") = ', stringWidth(canvas, "<b><sub>bold+sub</sub></b>"))
  print('sw(" hello ") = ', stringWidth(canvas, " hello "))
  print('sw("<u><super>underline+super</super></u>") = ',
        stringWidth(canvas, "<u><super>underline+super</super></u>"))

  pwidth1 = canvas.stringWidth("bold+sub", Font(size=canvas.defaultFont.size - sizedelta, bold=1))
  print("pwidth1 = ", pwidth1)
  pwidth2 = canvas.stringWidth(" hello ")
  print("pwidth2 = ", pwidth2)
  pwidth3 = canvas.stringWidth("underline+super",
                               Font(size=canvas.defaultFont.size - sizedelta, underline=1))
  print("pwidth3 = ", pwidth3)

  # these should be the same
  print("sfwidth = ", sfwidth, " pwidth = ", pwidth1 + pwidth2 + pwidth3)

  ################################################### testing greek characters
  # looks better in a larger font
  canvas = PDFCanvas('bigtest2.pdf')
  x = 10
  y = canvas.defaultFont.size * 1.5
  drawString(canvas, "&alpha; &beta; <chi/> &Delta; <delta/>", x, y, Font(size=16), color=blue)
  print("line starting with alpha should be font size 16")
  y = y + 30
  drawString(canvas, "&epsiv; &eta; &Gamma; <gamma/>", x, y, color=green)
  y = y + 30
  drawString(canvas, "&iota; &kappa; &Lambda; <lambda/>", x, y, color=blue)
  y = y + 30
  drawString(canvas, "<u>&mu;</u> &nu; <b>&Omega;</b> <omega/>", x, y, color=green)
  print("mu should be underlined, Omega should be big and bold")
  y = y + 30
  drawString(canvas, "&omicron; &Phi; &phi; <phiv/>", x, y, color=blue)
  y = y + 30
  drawString(canvas, "&Pi; &pi; &piv; <Psi/> &psi; &rho;", x, y, color=green)
  y = y + 30
  drawString(canvas, "<u>&Sigma; &sigma; &sigmav; <tau/></u>", x, y, color=blue)
  print("line starting with sigma should be completely underlined")
  y = y + 30
  drawString(canvas, "&Theta; &theta; &thetav; <Xi/> &xi; &zeta;", x, y, color=green)
  y = y + 30
  drawString(canvas, "That's &alpha;ll <u>folks</u><super>&omega;</super>", x, y)

  canvas.flush()


#test1()
#test2()
#stringformatTest()
