# $Id$
#
#  Copyright (C) 2003-2005  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" defines classes required for using 2D pharmacophore signatures in QTables

"""    
from rdkit import RDConfig
from qt import *
from qttable import *
from rdkit.qtGui import GuiTable,qtUtils
from rdkit.qtGui.GuiLib import MolTable
from rdkit.Chem.Pharm2D import Matcher
from rdkit.qtGui.GuiTextViewer import GuiTextViewer
from rdkit.Chem import CDXMLWriter
import cStringIO

try:
  from rdkit.utils import chemdraw
except ImportError:
  hasCDX=0
else:
  try:
    tmp = chemdraw.CDXDisplay
  except AttributeError:
    hasCDX=0
  else:
    hasCDX=1
from rdkit.utils import PilTools
import StringIO
  

class SignatureTableItem(MolTable.MolTableItem):
  """ #DOC

  """
  def __init__(self,*args,**kwargs):
    MolTable.MolTableItem.__init__(self,*args,**kwargs)
    self._target = None
    self._mol = None
    self._img = None
    self._sig = None
    self._bit = None
    self._matches = None
    self._cdxml = None
    self._active = 0

  def active(self):
    return self._active
  def setActive(self,which):
    self._active = which
    
  def setSig(self,sig):
    if not hasattr(sig,'GetShortestPathsOnly'):
      raise ValueError,'bad signature'
    self._sig = sig
  def sig(self):
    return self._sig

  def setBit(self,bit):
    self._bit = bit
  def bit(self):
    return self._bit

  def matches(self):
    return self._matches
  def setMatches(self,matches):
    self._matches = matches

  def cdxml(self):
    # FIX: make this use more than just the first match
    if self._cdxml is None:
      matches = self.matches()
      if matches is None:
        matches = Matcher.GetAtomsMatchingBit(self.sig(),self.bit(),self.mol())
        self.setMatches(matches)
      if len(matches):
        match = []
        for entry in matches[0]:
          match += list(entry)
        io = cStringIO.StringIO()
        CDXMLWriter.MolToCDXML(self.mol(),io,highlightAtoms=match)
        self._cdxml = io.getvalue()
      else:
        self._cdxml = ''
    return self._cdxml
      
    
    
  def draw(self,where=None):
    if not self._bit:
      MolTable.MolTableItem.draw(self,where=where)
      return
    
    if where is None:
      if not hasCDX: return
      if MolTable.drawInChemdraw:
        chemdraw.CDXDisplay(str(self.cdxml()),inFormat='chemical/cdx',
                            clear=1)
      else:
        if where is None:
          where = self.drawTarget()
        if where:
          img = self.image()
          if img is None:
            cdxml = self.cdxml()
            if cdxml:
              img = chemdraw.MolToPilImage(cdxml,inFormat='text/xml',
                                           outFormat='image/gif')
              tgtSize = where.canvas().size
              img = PilTools.FitImage(img,tgtSize)
              self.setImage(img)
              canv = where.canvas()  
          if img is not None:
            canv.clear()
            canv.drawImage(img,0,0)
            canv.flush()
          where.show()

  def details(self):
    if self.bit() is not None:
      return self.sig().GetBitDescription(self.bit(),includeBins=1)
    else:
      return ''
      
def insertSignatureTable(where,*args,**kwargs):
  return GuiTable.insertTable(where,SignatureTable,*args,**kwargs)

  
class SignatureTable(MolTable.MolTable):
  """ #DOC

  """
  def __init__(self,*args,**kwargs):
    MolTable.MolTable.__init__(self,*args,**kwargs)
    self._hypoDetailsWin = GuiTextViewer()
    self._hypoDetailsWin.hide()
    self.connect(self,SIGNAL("clicked(int,int,int,const QPoint&)"),self.clicked)
  def contentsToDelimText(self,delim=',',skipCols=None):
    """ converts the contents of the table to a delimited text string
       suitable for saving or putting on the clipboard
       
    **Arguments**

      - delim: (optional) the delimiter string to use

      - skipCols: (optional) a sequence containing columns which
        should not be included in the output.  If this is not
        provided, a list will be automatically generated to ensure
        that molecule columns do not get dumped.

    **Returns**

       a string

    """
    # FIX: update this
    if skipCols is None:
      skipCols = [0]
    return GuiTable.GuiTable.contentsToDelimText(self,delim=delim,skipCols=skipCols)

  def clicked(self,row,col,button,pos):
    if button == 2:
      itm = self.item(row,col)
      try:
        html = itm.details()
      except:
        qtUtils.logger.info('no details available',exc_info=True)
      else:
        if html:
          self._hypoDetailsWin.setText(html)
          self._hypoDetailsWin.show()
  #
  #  signals/slots
  #
  def currentChanged(self,row,col):
    itm = self.item(row,col)
    if itm:
      if isinstance(itm,SignatureTableItem):
        itm.draw()
      else:
        qtUtils.logger.info('not drawable')


    

    
