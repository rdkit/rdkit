# -*- coding: latin-1 -*-
#pdfdoc.py
""" 
PDFgen is a library to generate PDF files containing text and graphics.  It is the 
foundation for a complete reporting solution in Python.  

The module pdfdoc.py handles the 'outer structure' of PDF documents, ensuring that
all objects are properly cross-referenced and indexed to the nearest byte.  The 
'inner structure' - the page descriptions - are presumed to be generated before 
each page is saved.
pdfgen.py calls this and provides a 'canvas' object to handle page marking operators.
piddlePDF calls pdfgen and offers a high-level interface.

(C) Copyright Andy Robinson 1998-1999
"""

import os
import sys
import tempfile
import time
from math import ceil, cos, pi, sin
from types import *

try:
  import zlib
except ImportError:
  print("zlib not available, page compression not available")

from . import pdfmetrics, pdfutils
from .pdfgeom import bezierArc
from .pdfutils import LINEEND  # this constant needed in both

##############################################################
#
#            Constants and declarations
#
##############################################################

StandardEnglishFonts = [
  'Courier', 'Courier-Bold', 'Courier-Oblique', 'Courier-BoldOblique', 'Helvetica',
  'Helvetica-Bold', 'Helvetica-Oblique', 'Helvetica-BoldOblique', 'Times-Roman', 'Times-Bold',
  'Times-Italic', 'Times-BoldItalic', 'Symbol', 'ZapfDingbats'
]

PDFError = 'PDFError'
AFMDIR = '.'

A4 = (595.27, 841.89)  #default page size


class PDFDocument:
  """Responsible for linking and writing out the whole document.
    Builds up a list of objects using add(key, object).  Each of these
    must inherit from PDFObject and be able to write itself into the file.
    For cross-linking, it provides getPosition(key) which tells you where
    another object is, or raises a KeyError if not found.  The rule is that
    objects should only refer ones previously written to file.
    """

  def __init__(self):
    self.objects = []
    self.objectPositions = {}

    self.fonts = MakeType1Fonts()

    #mapping of Postscriptfont names to internal ones;
    #needs to be dynamically built once we start adding
    #fonts in.
    self.fontMapping = {}
    for i in range(len(StandardEnglishFonts)):
      psname = StandardEnglishFonts[i]
      pdfname = '/F%d' % (i + 1)
      self.fontMapping[psname] = pdfname

    self.pages = []
    self.pagepositions = []

    # position 1
    cat = PDFCatalog()
    cat.RefPages = 3
    cat.RefOutlines = 2
    self.add('Catalog', cat)

    # position 2 - outlines
    outl = PDFOutline()
    self.add('Outline', outl)

    # position 3 - pages collection
    self.PageCol = PDFPageCollection()
    self.add('PagesTreeRoot', self.PageCol)

    # positions 4-16 - fonts
    fontstartpos = len(self.objects) + 1
    for font in self.fonts:
      self.add('Font.' + font.keyname, font)
    self.fontdict = MakeFontDictionary(fontstartpos, len(self.fonts))

    # position 17 - Info
    self.info = PDFInfo()  #hang onto it!
    self.add('Info', self.info)
    self.infopos = len(self.objects)  #1-based, this gives its position

  def add(self, key, obj):
    self.objectPositions[key] = len(self.objects)  # its position
    self.objects.append(obj)
    obj.doc = self
    return len(self.objects) - 1  # give its position

  def getPosition(self, key):
    """Tell you where the given object is in the file - used for
        cross-linking; an object can call self.doc.getPosition("Page001")
        to find out where the object keyed under "Page001" is stored."""
    return self.objectPositions[key]

  def setTitle(self, title):
    "embeds in PDF file"
    self.info.title = title

  def setAuthor(self, author):
    "embedded in PDF file"
    self.info.author = author

  def setSubject(self, subject):
    "embeds in PDF file"
    self.info.subject = subject

  def printXref(self):
    self.startxref = sys.stdout.tell()
    print('xref')
    print(0, len(self.objects) + 1)
    print('0000000000 65535 f')
    for pos in self.xref:
      print('%0.10d 00000 n' % pos)

  def writeXref(self, f):
    self.startxref = f.tell()
    f.write('xref' + LINEEND)
    f.write('0 %d' % (len(self.objects) + 1) + LINEEND)
    f.write('0000000000 65535 f' + LINEEND)
    for pos in self.xref:
      f.write('%0.10d 00000 n' % pos + LINEEND)

  def printTrailer(self):
    print('trailer')
    print('<< /Size %d /Root %d 0 R /Info %d 0 R>>' % (len(self.objects) + 1, 1, self.infopos))
    print('startxref')
    print(self.startxref)

  def writeTrailer(self, f):
    f.write('trailer' + LINEEND)
    f.write('<< /Size %d /Root %d 0 R /Info %d 0 R>>' % (len(self.objects) + 1, 1, self.infopos) +
            LINEEND)
    f.write('startxref' + LINEEND)
    f.write(str(self.startxref) + LINEEND)

  def SaveToFile(self, filename):
    fileobj = open(filename, 'wb')
    self.SaveToFileObject(fileobj)
    fileobj.close()

  def SaveToFileObject(self, fileobj):
    """Open a file, and ask each object in turn to write itself to
        the file.  Keep track of the file position at each point for
        use in the index at the end"""
    f = fileobj
    i = 1
    self.xref = []
    f.write("%PDF-1.2" + LINEEND)  # for CID support
    f.write("%�춾" + LINEEND)
    for obj in self.objects:
      pos = f.tell()
      self.xref.append(pos)
      f.write(str(i) + ' 0 obj' + LINEEND)
      obj.save(f)
      f.write('endobj' + LINEEND)
      i = i + 1
    self.writeXref(f)
    self.writeTrailer(f)
    f.write('%%EOF')  # no lineend needed on this one!

    # with the Mac, we need to tag the file in a special
    #way so the system knows it is a PDF file.
    #This supplied by Joe Strout
    if os.name == 'mac':
      import macfs
      try:
        macfs.FSSpec(filename).SetCreatorType('CARO', 'PDF ')
      except Exception:
        pass

  def printPDF(self):
    "prints it to standard output.  Logs positions for doing trailer"
    print("%PDF-1.0")
    print("%�춾")
    i = 1
    self.xref = []
    for obj in self.objects:
      pos = sys.stdout.tell()
      self.xref.append(pos)
      print(i, '0 obj')
      obj.printPDF()
      print('endobj')
      i = i + 1
    self.printXref()
    self.printTrailer()
    print("%%EOF", end='')

  def addPage(self, page):
    """adds page and stream at end.  Maintains pages list"""
    #page.buildstream()
    pos = len(self.objects)  # work out where added

    page.ParentPos = 3  #pages collection
    page.info = {
      'parentpos': 3,
      'fontdict': self.fontdict,
      'contentspos': pos + 2,
    }

    self.PageCol.PageList.append(pos + 1)
    self.add('Page%06d' % len(self.PageCol.PageList), page)
    #self.objects.append(page)
    self.add('PageStream%06d' % len(self.PageCol.PageList), page.stream)
    #self.objects.append(page.stream)

  def hasFont(self, psfontname):
    return psfontname in self.fontMapping

  def getInternalFontName(self, psfontname):
    try:
      return self.fontMapping[psfontname]
    except KeyError:
      raise PDFError("Font {0} not available in document".format(psfontname))

  def getAvailableFonts(self):
    fontnames = self.fontMapping.keys()
    fontnames.sort()
    return fontnames

  ##############################################################
  #
  #            Utilities
  #
  ##############################################################


class OutputGrabber:
  """At times we need to put something in the place of standard
    output.  This grabs stdout, keeps the data, and releases stdout
    when done.
    
    NOT working well enough!"""

  def __init__(self):
    self.oldoutput = sys.stdout
    sys.stdout = self
    self.closed = 0
    self.data = []

  def write(self, x):
    if not self.closed:
      self.data.append(x)

  def getData(self):
    return ' '.join(self.data)

  def close(self):
    sys.stdout = self.oldoutput
    self.closed = 1

  def __del__(self):
    if not self.closed:
      self.close()


def testOutputGrabber():
  gr = OutputGrabber()
  for i in range(10):
    print('line', i)
  data = gr.getData()
  gr.close()
  print('Data...', data)


##############################################################
#
#            PDF Object Hierarchy
#
##############################################################


class PDFObject:
  """Base class for all PDF objects.  In PDF, precise measurement
    of file offsets is essential, so the usual trick of just printing
    and redirecting output has proved to give different behaviour on
    Mac and Windows.  While it might be soluble, I'm taking charge
    of line ends at the binary level and explicitly writing to a file.
    The LINEEND constant lets me try CR, LF and CRLF easily to help
    pin down the problem."""

  def save(self, file):
    "Save its content to an open file"
    file.write('% base PDF object' + LINEEND)

  def printPDF(self):
    self.save(sys.stdout)


class PDFLiteral(PDFObject):
  " a ready-made one you wish to quote"

  def __init__(self, text):
    self.text = text

  def save(self, file):
    file.write(self.text + LINEEND)


class PDFCatalog(PDFObject):
  "requires RefPages and RefOutlines set"

  def __init__(self):
    self.template = LINEEND.join(
      ['<<', '/Type /Catalog', '/Pages %d 0 R', '/Outlines %d 0 R', '>>'])

  def save(self, file):
    file.write(self.template % (self.RefPages, self.RefOutlines) + LINEEND)


class PDFInfo(PDFObject):
  """PDF documents can have basic information embedded, viewable from
    File | Document Info in Acrobat Reader.  If this is wrong, you get
    Postscript errors while printing, even though it does not print."""

  def __init__(self):
    self.title = "untitled"
    self.author = "anonymous"
    self.subject = "unspecified"

    now = time.localtime(time.time())
    self.datestr = '%04d%02d%02d%02d%02d%02d' % tuple(now[0:6])

  def save(self, file):
    file.write(
      LINEEND.join([
        "<</Title (%s)", "/Author (%s)", "/CreationDate (D:%s)", "/Producer (PDFgen)",
        "/Subject (%s)", ">>"
      ]) % (pdfutils._escape(self.title), pdfutils._escape(self.author), self.datestr,
            pdfutils._escape(self.subject)) + LINEEND)


class PDFOutline(PDFObject):
  "null outline, does nothing yet"

  def __init__(self):
    self.template = LINEEND.join(['<<', '/Type /Outlines', '/Count 0', '>>'])

  def save(self, file):
    file.write(self.template + LINEEND)


class PDFPageCollection(PDFObject):
  "presumes PageList attribute set (list of integers)"

  def __init__(self):
    self.PageList = []

  def save(self, file):
    lines = ['<<', '/Type /Pages', '/Count %d' % len(self.PageList), '/Kids [']
    for page in self.PageList:
      lines.append(str(page) + ' 0 R ')
    lines.append(']')
    lines.append('>>')
    text = LINEEND.join(lines)
    file.write(text + LINEEND)


class PDFPage(PDFObject):
  """The Bastard.  Needs list of Resources etc. Use a standard one for now.
    It manages a PDFStream object which must be added to the document's list
    of objects as well."""

  def __init__(self):
    self.drawables = []
    self.pagewidth = 595  #these are overridden by piddlePDF
    self.pageheight = 842
    self.stream = PDFStream()
    self.hasImages = 0
    self.pageTransitionString = ''  # presentation effects
    # editors on different systems may put different things in the line end
    # without me noticing.  No triple-quoted strings allowed!
    self.template = LINEEND.join([
      '<<',
      '/Type /Page',
      '/Parent %(parentpos)d 0 R',
      '/Resources',
      '   <<',
      '   /Font %(fontdict)s',
      '   /ProcSet %(procsettext)s',
      '   >>',
      '/MediaBox [0 0 %(pagewidth)d %(pageheight)d]',  #A4 by default
      '/Contents %(contentspos)d 0 R',
      '%(transitionString)s',
      '>>'
    ])

  def setCompression(self, onoff=0):
    "Turns page compression on or off"
    assert onoff in [0, 1], "Page compression options are 1=on, 2=off"
    self.stream.compression = onoff

  def save(self, file):
    self.info['pagewidth'] = self.pagewidth
    self.info['pageheight'] = self.pageheight
    # check for image support
    if self.hasImages:
      self.info['procsettext'] = '[/PDF /Text /ImageC]'
    else:
      self.info['procsettext'] = '[/PDF /Text]'
    self.info['transitionString'] = self.pageTransitionString

    file.write(self.template % self.info + LINEEND)

  def clear(self):
    self.drawables = []

  def setStream(self, data):
    if isinstance(data, (list, tuple)):
      data = LINEEND.join(data)
    self.stream.setStream(data)


TestStream = "BT /F6 24 Tf 80 672 Td 24 TL (   ) Tj T* ET"


class PDFStream(PDFObject):
  "Used for the contents of a page"

  def __init__(self):
    self.data = None
    self.compression = 0

  def setStream(self, data):
    self.data = data

  def save(self, file):
    #avoid crashes if they wrote nothing in the page
    if self.data is None:
      self.data = TestStream

    if self.compression == 1:
      comp = zlib.compress(self.data)  #this bit is very fast...
      base85 = pdfutils._AsciiBase85Encode(comp)  #...sadly this isn't
      wrapped = pdfutils._wrap(base85)
      data_to_write = wrapped
    else:
      data_to_write = self.data
    # the PDF length key should contain the length including
    # any extra LF pairs added by Print on DOS.

    #lines = len(string.split(self.data,'\n'))
    #length = len(self.data) + lines   # one extra LF each
    length = len(data_to_write) + len(LINEEND)  #AR 19980202
    if self.compression:
      file.write('<< /Length %d /Filter [/ASCII85Decode /FlateDecode]>>' % length + LINEEND)
    else:
      file.write('<< /Length %d >>' % length + LINEEND)
    file.write('stream' + LINEEND)
    file.write(data_to_write + LINEEND)
    file.write('endstream' + LINEEND)


class PDFImage(PDFObject):
  # sample one while developing.  Currently, images go in a literals
  def save(self, file):
    file.write(
      LINEEND.join([
        '<<', '/Type /XObject', '/Subtype /Image', '/Name /Im0', '/Width 24', '/Height 23',
        '/BitsPerComponent 1', '/ColorSpace /DeviceGray', '/Filter /ASCIIHexDecode', '/Length 174',
        '>>', 'stream', '003B00 002700 002480 0E4940 114920 14B220 3CB650',
        '75FE88 17FF8C 175F14 1C07E2 3803C4 703182 F8EDFC',
        'B2BBC2 BB6F84 31BFC2 18EA3C 0E3E00 07FC00 03F800', '1E1800 1FF800>', 'endstream', 'endobj'
      ]) + LINEEND)


class PDFType1Font(PDFObject):

  def __init__(self, key, font):
    self.fontname = font
    self.keyname = key
    self.template = LINEEND.join([
      '<<', '/Type /Font', '/Subtype /Type1', '/Name /%s', '/BaseFont /%s',
      '/Encoding /MacRomanEncoding', '>>'
    ])

  def save(self, file):
    file.write(self.template % (self.keyname, self.fontname) + LINEEND)

  ##############################################################
  #
  #            some helpers
  #
  ##############################################################


def MakeType1Fonts():
  "returns a list of all the standard font objects"
  fonts = []
  pos = 1
  for fontname in StandardEnglishFonts:
    font = PDFType1Font('F' + str(pos), fontname)
    fonts.append(font)
    pos = pos + 1
  return fonts


def MakeFontDictionary(startpos, count):
  "returns a font dictionary assuming they are all in the file from startpos"
  dict = "  <<" + LINEEND
  pos = startpos
  for i in range(count):
    dict = dict + '\t\t/F%d %d 0 R ' % (i + 1, startpos + i) + LINEEND
  dict = dict + "\t\t>>" + LINEEND
  return dict


if __name__ == '__main__':
  print('For test scripts, run test1.py to test6.py')
