#
#  Copyright (C) 2003,2004  Rational Discovery LLC
#   All Rights Reserved
#
import RDConfig
from reportlab import platypus
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab import rl_config
import PIL
import os
from reportlab.platypus import *
import ReportUtils

styles = getSampleStyleSheet()


def FormatTable(data,rowHeadings=0,colHeadings=0,border=0,
                fontName='',fontSize=0,centerRows=0):
  """ data should be a sequence of sequences
  """
  data = list(data)
  for i,row in enumerate(data):
    row = list(row)
    for j,cell in enumerate(row):
      if isinstance(cell,basestring):
        row[j] = cell.replace('<br>','\n')
        #pass
    data[i] = row  
  tbl = platypus.Table(data)
  style = []
  if border:
    style.append(('GRID',(0,0),(-1,-1),1,colors.black))
  if not fontName:
    fontName = styles['BodyText'].fontName
  if not fontSize:
    fontSize = styles['BodyText'].fontSize
  style.append(('FONT',(0,0),(-1,-1),fontName,fontSize))

  # FIX: font specification should be more general:
  if(colHeadings):
    style.append(('FONT',(0,0),(-1,0),'Times-Bold',fontSize))
  if(rowHeadings):
    style.append(('FONT',(0,0),(0,-1),'Times-Bold',fontSize))
    style.append(('ALIGN',(0,0),(0,-1),'RIGHT'))
  if centerRows:
    style.append(('VALIGN',(0,0),(-1,-1),'MIDDLE'))

  tbl.setStyle(platypus.TableStyle(style))
  return tbl

    
               

  

#--------------------
#
#  These are included to support headings and the like that get bookmarked
#  automagically
#
#--------------------
bookmarkID=1
class HeadingParagraph(Paragraph):
  def __init__(self,text,style,level=0):
    self._label = text
    self._level = level
    Paragraph.__init__(self,text,style)
  def drawPara(self,debug=0):
    global bookmarkID
    key='bk-%d'%bookmarkID
    bookmarkID+=1

    canvas = self.canv
    canvas.bookmarkPage(key,fitType="XYZ")
    lev = self._level
    while lev>=0:
      try:
        canvas.addOutlineEntry(self._label,key,lev,0)
      except:
        lev -=1
      else:
        break
    Paragraph.drawPara(self,debug=debug)
    
class Heading1(HeadingParagraph):
  def __init__(self,text):
    HeadingParagraph.__init__(self,text,styles['Heading1'],level=0)
    
class Heading2(HeadingParagraph):
  def __init__(self,text):
    HeadingParagraph.__init__(self,text,styles['Heading2'],level=1)



class PDFReport(object):
  pageHeader = ""
  pageSize = rl_config.defaultPageSize
  includeLogo = 1
  headerLogoFilename = os.path.join(RDConfig.RDCodeDir,'Reports',
                                    'RD-Logo.jpg')
  footerLogoFilename = os.path.join(RDConfig.RDCodeDir,'Reports',
                                    'RDOnline-Logo.jpg')
  def __init__(self,verbose=0):
    self.verbose = verbose
  def onPage(self,canvas,doc):
    if self.verbose:
      print '----- page ----'
    canvas.saveState()
    if self.includeLogo:
      if self.headerLogoFilename:
        img = platypus.Image(self.headerLogoFilename)
        img.drawHeight=1*inch
        img.drawWidth=1*inch
        img.drawOn(canvas,0.5*inch,self.pageSize[1]-1.0*inch)

      if self.footerLogoFilename:
        pImg = PIL.Image.open(self.footerLogoFilename)
        aspect = float(pImg.size[1]) / pImg.size[0]
        pImg = None
        img = platypus.Image(self.footerLogoFilename)
        img.drawHeight=.33*inch
        img.drawWidth = img.drawHeight / aspect
        img.drawOn(canvas,self.pageSize[0]-2.0*inch,0.33*inch)
    canvas.setStrokeColor(colors.black)
    canvas.line(.75*inch,self.pageSize[1]-1.0*inch,
                self.pageSize[0]-.75*inch,self.pageSize[1]-1.0*inch)
    canvas.line(.75*inch,.75*inch,
                self.pageSize[0]-.75*inch,.75*inch)
    canvas.setFont('Times-Roman',9)
    canvas.drawString(inch, 0.5 * inch, "Page %d %s" % (doc.page,
                                                         self.pageHeader))
    canvas.restoreState()
    
if __name__=='__main__':
  styles = getSampleStyleSheet()
  title = 'Test'
  reportTemplate = PDFReport()
  reportTemplate.pageHeader = title
  doc = platypus.SimpleDocTemplate('test.pdf')
  elements = [platypus.Paragraph("This is a test.",
                                 styles['Normal'])]
  elements.append(platypus.Paragraph("lt: 3 &lt; 4",
                                 styles['Normal']))
  elements.append(platypus.Paragraph("lt: 3&lt;4",
                                 styles['Normal']))
  elements.append(platypus.Paragraph("gt: 4&gt;3",
                                 styles['Normal']))
  elements.append(platypus.Paragraph("both: 3&lt;4 4&gt;3",
                                 styles['Normal']))

  elements.append(platypus.Paragraph("""| GHS_LD50_Class in mg/kg |
  | 0 | &lt;=5 |
  | 1 | 5&lt;x&lt;=50 |
  | 2 | 50&lt;x&lt;=300 |
  | 3 | 300&lt;x&lt;=2000 |
  | 4 | 2000&lt;x&lt;=5000 |
  | 5 | &gt;5000 | 
  """,styles['Normal']))
  doc.build(elements,onFirstPage=reportTemplate.onPage,
            onLaterPages=reportTemplate.onPage)

