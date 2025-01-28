## font mapping stuff ##
# by default, I'm using the standard acrobat fonts so that pdfmetrics should work well
#
DefaultFace = 'times'

# these are the required fonts for piddle

PidLegalFonts = {
  "courier": "courier",  # note: keys are lowercased 
  "helvetica": "helvetica",
  "monospaced": "courier",
  "sansserif": "helvetica",
  "serif": "times",
  "times": "times",
  "symbol": "symbol",
  "zapfdingbats": "zapfdingbats"
}  # Could add more...

#PidLegalShapes:
Roman = "Roman"
Bold = "Bold"
Italic = "Italic"
BoldItalic = "BoldOblique"

MapPid2PyartFontName = {
  ("helvetica", Roman): "Helvetica",
  ("helvetica", Bold): "Helvetica-Bold",
  ("helvetica", Italic): "Helvetica-Oblique",
  ("helvetica", BoldItalic): "Helvetica-BoldOblique",
  ("times", Roman): "Times-Roman",
  ("times", Bold): "Times-Bold",
  ("times", Italic): "Times-Italic",
  ("times", BoldItalic): "Times-BoldItalic",
  ("courier", Roman): "Courier",
  ("courier", Bold): "Courier-Bold",
  ("courier", Italic): "Courier-Oblique",
  ("courier", BoldItalic): "Courier-BoldOblique",
  ("symbol", Roman): "Symbol",
  ("symbol", Bold): "Symbol",
  ("symbol", Italic): "Symbol",
  ("symbol", BoldItalic): "Symbol",
  ("zapfdingbats", Roman): "ZapfDingbats",
  ("zapfdingbats", Bold): "ZapfDingbats",
  ("zapfdingbats", Italic): "ZapfDingbats",
  ("zapfdingbats", BoldItalic): "ZapfDingbats"
}

#  PDFFontMapping = { ("helvetica", Roman): "Helvetica",
#                     ("helvetica", Bold): "Helvetica-Bold",
#                     ("helvetica", Italic): "Helvetica-Oblique",
#                     ("times", Roman) : "Times-Roman",
#                     ("times", Bold) : "Times-Bold",
#                     ("times", Italic) : "Times-Italic",
#                     ("courier", Roman) : "Courier",
#                     ("courier", Bold) : "Courier-Bold",
#                     ("courier", Italic) : "Courier-Oblique",
#                     ("symbol", Roman) : "Symbol",
#                     ("symbol", Bold) :  "Symbol",
#                     ("symbol", Italic) : "Symbol",
#                     ("zapfdingbats", Roman ) : "ZapfDingbats",
#                     ("zapfdingbats", Bold ) : "ZapfDingbats",
#                     ("zapfdingbats", Italic ) : "ZapfDingbats" }


def getPyartName(pidfont):
  if not pidfont.bold:
    if pidfont.italic:
      shape = Italic
    else:
      shape = Roman
  else:
    if pidfont.italic:
      shape = BoldItalic
    else:
      shape = Bold

  face = pidfont.face or DefaultFace
  # print "pidfont.face = %s" % pidfont.face

  face = face.lower()
  if face in PidLegalFonts:
    return MapPid2PyartFontName[(PidLegalFonts[face], shape)]
  else:
    raise ValueError("Illegal Font")


getPdfName = getPyartName
#  def getPdfName(pidfont):
#      if pidfont.bold:
#          shape = Bold
#      elif pidfont.italic:
#          shape = Italic
#      else:
#          shape = Roman

#      face = pidfont.face or DefaultFace
#      # print "pidfont.face = %s" % pidfont.face

#      face = face.lower()
#      if face in PidLegalFonts:
#          return PDFFontMapping[ ( PidLegalFonts[face], shape) ]
#      else:
#          raise "Illegal Font"
