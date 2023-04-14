# lots to do:
#   __ native drawLines
#   __ add native drawCurve method
#   __ native rectangle/round rect method
#   __ native drawEllipse
#   __ native drawArc
#   __ drawImage support (work on Pyart side of things)

import Fontmapping  # helps by mapping pid font classes to Pyart font names
import pyart

from rdkit.sping.PDF import pdfmetrics
from rdkit.sping.pid import *

# note for now I'm just going to do the standard PDF fonts & forget the rest


class PyartCanvas(Canvas):
  "note the default face is 'times' and is set in Fontmapping.py"

  def __init__(self, size=(300, 300), name='PyartCanvas.png'):

    self._pycan = pyart.Canvas(size[0], size[1], dpi=72)
    self.filename = name

    Canvas.__init__(self, size, name)

  # self.defaultFillColor = transparent

  # now we need to setup our tracking of the defaults vs the current state

  # see if the __setattr__ approach is any better than the _updateXX strategy

  def __setattr__(self, name, value):
    if name == 'defaultLineColor':
      if value:
        # print('setting defaultLineColor to %s, 0x%x' % (value, value.toHexRGB()))
        if value != transparent:
          self._pycan.gstate.stroke = value.toHexRGB()
        self.__dict__[name] = value
    elif name == 'defaultFillColor':
      if value:
        if value != transparent:
          self._pycan.gstate.fill = value.toHexRGB()
        self.__dict__[name] = value
    elif name == 'defaultLineWidth':
      if value:
        self._pycan.gstate.stroke_width = value
        self.__dict__[name] = value
    elif name == 'defaultFont':
      if value:
        self.__dict__[name] = value
        self._setPyartFont(value)
      else:  # received None so set to default font face & size=12
        self.__dict__[name] = Font(face='times')
        self._setPyartFont(self.__dict__[name])

    else:
      self.__dict__[name] = value

    ## Private methods ##

  def _protectArtState(self, bool):
    if bool:
      self._pycan.gsave()
    return bool

  def _restoreArtState(self, bool):
    if bool:
      self._pycan.grestore()

  def _setPyartFont(self, fontInstance):
    # accounts for "None" option
    # does not act on self.defaultFont at all

    fontsize = fontInstance.size
    self._pycan.gstate.font_size = fontsize
    # map pid name for font to Pyart name
    pyartname = Fontmapping.getPyartName(fontInstance)
    self._pycan.gstate.setfont(pyartname)

  # # # # #

  ### public PID Canvas methods ##

  def clear(self):
    pass

  def flush(self):
    pass

  def save(self, file=None, format=None):
    # fileobj = getFileObject(file)

    if not file:
      file = self.filename

    if isinstance(file, StringType):
      self._pycan.save(file)
    else:
      raise NotImplementedError

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
      face = piddle_font_map[font.face.lower()]
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

  def stringWidth(self, s, font=None):
    if not font:
      font = self.defaultFont
    fontname = Fontmapping.getPdfName(font)
    return pdfmetrics.stringwidth(s, fontname) * font.size * 0.001

  def fontAscent(self, font=None):
    if not font:
      font = self.defaultFont
    fontname = Fontmapping.getPdfName(font)
    return pdfmetrics.ascent_descent[fontname][0] * 0.001 * font.size

  def fontDescent(self, font=None):
    if not font:
      font = self.defaultFont
    fontname = Fontmapping.getPdfName(font)
    return -pdfmetrics.ascent_descent[fontname][1] * 0.001 * font.size

  def drawLine(self, x1, y1, x2, y2, color=None, width=None):
    ## standard code ##
    color = color or self.defaultLineColor
    width = width or self.defaultLineWidth
    if color != transparent:
      changed = self._protectArtState((color != self.defaultLineColor)
                                      or (width != self.defaultLineWidth))
      if color != self.defaultLineColor:
        self._pycan.gstate.stroke = color.toHexRGB()
        # print("color is %s <-> %s" % (color, color.toHexStr()))
      if width != self.defaultLineWidth:
        self._pycan.gstate.stroke_width = width
      ###################

      # actual drawing
      p = pyart.VectorPath(3)
      p.moveto_open(x1, y1)
      p.lineto(x2, y2)
      self._pycan.stroke(p)

      ## standard code ##
      if changed:
        self._pycan.grestore()
      ###################

      #  def drawLines(self, lineList, color=None, width=None):
      #    pass

  def drawString(self, s, x, y, font=None, color=None, angle=0):
    # start w/ the basics
    self._pycan.drawString(x, y, s)

  def drawPolygon(self, pointlist, edgeColor=None, edgeWidth=None, fillColor=None, closed=0):

    eColor = edgeColor or self.defaultLineColor
    fColor = fillColor or self.defaultFillColor
    eWidth = edgeWidth or self.defaultLineWidth

    changed = self._protectArtState((eColor != self.defaultLineColor)
                                    or (eWidth != self.defaultLineWidth)
                                    or (fColor != self.defaultFillColor))

    if eColor != self.defaultLineColor:
      self._pycan.gstate.stroke = eColor.toHexRGB()

    if fColor != self.defaultFillColor:
      self._pycan.gstate.fill = fColor.toHexRGB()

    if eWidth != self.defaultLineWidth:
      self._pycan.gstate.stroke_width = eWidth

    path = pyart.VectorPath(len(pointlist) + 1)
    if closed:
      path.moveto_closed(pointlist[0][0], pointlist[0][1])
    else:
      path.moveto_open(pointlist[0][0], pointlist[0][1])

    for pt in pointlist[1:]:
      path.lineto(pt[0], pt[1])

    if closed:
      path.close()

    if fColor != transparent and closed:
      self._pycan.fill(path)

    if eColor != transparent:
      self._pycan.stroke(path)

    self._restoreArtState(changed)

  #def drawCurve(self, x1, y1, x2, y2, x3, y3, x4, y4,
  #                  edgeColor=None, edgeWidth=None, fillColor=None, closed=0):
  #   pass

  #  def drawRoundRect(self, x1,y1, x2,y2, rx=8, ry=8,
  #                        edgeColor=None, edgeWidth=None, fillColor=None):
  #          pass

  #      def drawEllipse(self, x1,y1, x2,y2, edgeColor=None, edgeWidth=None,
  #                      fillColor=None):
  #          pass

  #      def drawArc(self, x1,y1, x2,y2, startAng=0, extent=360, edgeColor=None,
  #                  edgeWidth=None, fillColor=None):

  #          pass

  #      def drawFigure(self, partList,
  #                     edgeColor=None, edgeWidth=None, fillColor=None, closed=0):
  #          pass

  #      def drawImage(self, image, x1, y1, x2=None,y2=None):
  #          pass

  ## basic tests ##


if __name__ == '__main__':
  import rdkit.sping.tests.pidtest
  can = PyartCanvas(size=(300, 300), name='basictest.png')

  #can.defaultLineColor = Color(0.7, 0.7, 1.0)
  #can.drawLine(10,10, 290,290)
  #can.drawLine(10,10, 50, 10, color=green, width = 4.5)
  rdkit.sping.tests.pidtest.drawBasics(can)
  can.save(file='basicTest.png')
  print('saving basicTest.png')

  can = PyartCanvas(size=(400, 400), name='test-strings.png')
  rdkit.sping.tests.pidtest.drawStrings(can)
  can.save()
