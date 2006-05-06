# $Id: MolTable.py 5092 2006-03-15 19:36:19Z NightlyBuild $
#
#  Copyright (C) 2002-2005  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" defines classes required for using molecules in QTables

"""    
import RDConfig
import Chem
from Chem import rdDepictor
from qt import *
from qttable import *
from utils.PilTools import PilImgToQPixmap
from qtGui import GuiTable,qtUtils
tryChemdraw=False
hasCDX=0

if tryChemdraw:
  try:
    from utils import chemdraw
  except ImportError:
    hasCDX=0
  else:
    try:
      tmp = chemdraw.CDXDisplay
    except AttributeError:
      hasCDX=0
    else:
      hasCDX=1
  try:
    from utils import chemdraw_qax
  except ImportError:
    hasCDX_ax=0
  else:
    try:
      tmp = chemdraw_qax.ChemdrawPanel
    except AttributeError:
      hasCDX_ax=0
    else:
      hasCDX_ax=1
from utils import PilTools
import StringIO,types

drawInChemdraw=1
if hasCDX and not hasCDX_ax:
  drawInChemdraw=1


class MolTableItem(QTableItem):
  """ container for holding molecules in QTables

    Important Attributes:

     - target:  The destination where the molecule will draw itself
       (not yet implemented)

     - mol: the molecule associated with the item  

  """
  def __init__(self,*args,**kwargs):
    QTableItem.__init__(self,*args)
    self._target = None
    self._mol = None
    self._img = None
    self._flagged = 0
    self._flagColor = Qt.yellow

  def smiles(self):
    """  Returns the molecule's SMILES

    """
    return self.text()
  def setDrawTarget(self,tgt):
    """ setter for the target attribute

    """
    self._target = tgt
  def drawTarget(self):
    """ getter for the target attribute

    """
    return self._target
  def setMol(self,mol):
    """ setter for the molecule attribute

    """
    self._mol = mol
    #img = self.toImage(size=(200,200))
    #pm = PilTools.PilImgToQPixmap(img)
    #self.setPixmap(pm)
    
  def mol(self):
    """ getter for the molecule attribute

    """
    return self._mol

  def setImage(self,img):
    """ img is a PIL image

    """
    self._img = img
  def image(self):
    return self._img

  def flagged(self):
    return self._flagged
  def setFlagged(self,val,flagColor=None):
    self._flagged=val
    if flagColor is not None:
      self.setFlagColor(flagColor)
      
  def flagColor(self):
    return self._flagColor
  def setFlagColor(self,flagColor):
    self._flagColor=flagColor

  
  def paint(self,painter,colorGroup,rect,selected):
    if self._flagged:
      colorGroup.setColor(QColorGroup.Base,self._flagColor)
    QTableItem.paint(self,painter,colorGroup,rect,selected)
    
  def drawCDX(self,where=None):
    if not hasCDX: return
    if drawInChemdraw:
      chemdraw.CDXDisplay(str(self.smiles()),inFormat='chemical/daylight-smiles',
                          clear=1)
    else:
      if where is None:
        where = self.drawTarget()
      if where:
        img = self.image()
        if img is None:
          smi = str(self.smiles())
          img = chemdraw.SmilesToPilImage(smi)
          tgtSize = where.canvas().size
          img = PilTools.FitImage(img,tgtSize)
          self.setImage(img)
        canv = where.canvas()  
        canv.clear()
        canv.drawImage(img,0,0)
        canv.flush()
        where.show()

  def draw(self,where=None):
    if hasCDX:
      self.drawCDX(where=where)
    else:
      if where is None:
        where = self.drawTarget()
      if where and self.mol():
        mol = self.mol()
        where.setMol(mol)
        if mol.HasProp('_Name'):
          where.setCaption(mol.GetProp('_Name'))
        if not where.isVisible(): where.show()

  def toImage(self,fName=None,size=(100,100),fontSize=12,lineWidth=.5):
    from Chem.Draw.MolDrawing import MolDrawing
    from sping.ReportLab.pidReportLab import RLCanvas as Canvas
    canv = Canvas(size)
    d = MolDrawing(canvas=canv)
    d.atomLabelFontSize=fontSize
    d.atomLabelMinFontSize=4
    d.bondLineWidth=lineWidth
    mol = self.mol()
    if not mol.GetNumConformers():
      rdDepictor.Compute2DCoords(mol)
    d.AddMol(mol)
    return canv.drawing

    
  def __cmp__(self,other):
    """ by default we'll compare using smiles (assume it's canonical) """
    if not isinstance(other,MolTableItem):
      return -1
    return cmp(str(self.smiles()),str(other.smiles()))
    
def insertMolTable(where,*args,**kwargs):
  return GuiTable.insertTable(where,MolTable,*args,**kwargs)
  
class MolTable(GuiTable.GuiTable):
  """ a table class which can include a row of molecules

  """
  def __init__(self,*args,**kwargs):
    if kwargs.has_key('molCanvas'):
      self.molCanvas = kwargs['molCanvas']
      del kwargs['molCanvas']
    else:
      self.molCanvas=None
    GuiTable.GuiTable.__init__(self,*args,**kwargs)

    self.setSelectionMode(QTable.Single)
    self._includeCheckboxes=0
    self._checkboxCol=0
    self._molCol=0
    self._allowInserts=0
    self._allowSaves=0
    self._dir = "."
    self._displayableTypes = (types.IntType,types.FloatType,
                              types.StringType,types.UnicodeType,
                              types.LongType,types.ComplexType)
    self._changedCallbacks = []
    
  def allowInserts(self):
    return self._allowInserts
  def setAllowInserts(self,val):
    self._allowInserts = val
  def allowSaves(self):
    return self._allowSaves
  def setAllowSaves(self,val):
    self._allowSaves = val

  def addChangedCallback(self,cb):
    """ adds a callback for when the active cell is changed

    callback functions take 3 args:
       1) this table
       2) the row
       3) the column

    """
    self._changedCallbacks.append(cb)
  def clearChangedCallbacks(self):
    self._changedCallbacks = []
    
  def contentsToDelimText(self,delim=',',skipCols=None,replaceMolName=0):
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
    if skipCols is None:
      skipCols = []
      for col in range(self.numCols()):
        if isinstance(self.item(0,col),MolTableItem):
          skipCols.append(col)
    if replaceMolName:
      hdr = self.horizontalHeader()
      origLabel = str(hdr.label(self._molCol))
      hdr.setLabel(self._molCol,'SMILES')
    try:
      res= GuiTable.GuiTable.contentsToDelimText(self,delim=delim,skipCols=skipCols)
    except:
      qtUtils.logger.info('problems encountered converting table to text',exc_info=True)
      res = ""
    if replaceMolName:
      hdr = self.horizontalHeader()
      hdr.setLabel(self._molCol,origLabel)
    return res


  def contentsToPDF(self,fName=None,skipCols=None,replaceMolName=0,sketchWidth=1.5):
    from Reports import PDFReport, ReportUtils
    from cStringIO import StringIO
    import os
    
    if not fName:
      fName = str(QFileDialog.getSaveFileName(self._dir,'Acrobat files (*.pdf);;All files (*.*)'))
      if not fName:
        return None
      splitName = os.path.split(fName)
      self._dir = os.sep.join(splitName[0:-1])


    try:
      outF =open(fName,'wb+')
    except IOError:
      qtUtils.error('Could not open output file %s for writing.'%fName)
      return

    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    tblD= []
    if skipCols is None:
      skipCols = []
    cols = self._limitCols(skipCols)
    hdr = self.horizontalHeader()
    nCols = len(cols)
    tmp = ['']*nCols
    for i in range(nCols):
      label = str(hdr.label(cols[i]))
      if len(label)>=10:
        label =label[:10]
      tmp[i] = label
    tblD.append(tmp)

    tempHandler = ReportUtils.TempFileHandler()

    nRows = self.numRows()
    dlg = qtUtils.ProgressDialogHolder('Generating PDF',nRows)
    nCols = len(cols)
    for row in xrange(nRows):
    #for row in xrange(6):
      tmp = ['']*nCols
      for i,col in enumerate(cols):
        itm = self.item(row,col)
        if not itm:
          tmp[i] = ''
        else:
          if hasattr(itm,'toImage'):
            fN = tempHandler.get('.jpg')
            #itm.toImage(fN)
            #img = PDFReport.platypus.Image(fN,sketchWidth*PDFReport.inch,sketchWidth*PDFReport.inch)
            img = itm.toImage()
            #v = 'mol: '+str(itm.text())
            tmp[i] = img
          else:
            try:
              v = '%8.6g'%float(str(itm.text()))
            except ValueError:
              v = str(itm.text())
              
            tmp[i] = v
      dlg(row)
      tblD.append(tmp)
          
    elements = []
    tbl = PDFReport.LongTable(tblD,repeatRows=1)
    tbl.setStyle(PDFReport.TableStyle([('GRID',(0,0),(-1,-1),1,
                                        PDFReport.colors.black),
                                       ('FONT',(0,0),(-1,0),
                                        'Times-Bold',10),
                                       ('FONT',(0,1),(-1,-1),
                                        'Times-Roman',10),
                                       ]))
    tbl._argW[self._molCol] = sketchWidth*PDFReport.inch*1.1
    elements.append(tbl)
    template = PDFReport.PDFReport()
    template.pageHeader = fName
    doc = PDFReport.SimpleDocTemplate(outF)
    doc.build(elements,onFirstPage=template.onPage,
              onLaterPages=template.onPage)
    outF.close()
    QApplication.restoreOverrideCursor()
    
  def loadFromMolSupplier(self,mols,drawTarget=None,
                          includeCheckboxes=0,adjustCols=1,
                          processFunc=None,kekulize=False):
    # clear the current contents:
    self.setNumRows(0)
    self.setNumCols(0)

    if drawTarget is None:
      drawTarget=self.molCanvas

    self._includeCheckboxes = includeCheckboxes
    colNames=mols.GetColumnNames()
    nCols = len(colNames)
    if not self._includeCheckboxes:
      self._molCol=0
      nExtraCols=1
    else:
      nExtraCols=2
      self._molCol=1
    self.setNumCols(nCols+nExtraCols)
    self.setUpdatesEnabled(0)
    # set column headers
    hdr = self.horizontalHeader()
    hdr.setLabel(self._molCol,'Molecule')
    if self._includeCheckboxes:
      hdr.setLabel(self._checkboxCol,'  ')
      self.adjustColumn(self._checkboxCol)
    for i in range(nCols):
      hdr.setLabel(i+nExtraCols,colNames[i])

    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    if hasattr(mols,'__len__'):
      dlg = qtUtils.ProgressDialogHolder('Loading Molecules',len(mols))
    else:
      dlg = None

    nFailed = 0
    row = 0
    for mol in mols:
      if mol:
        if processFunc:
          processFunc(mol)
        
        self.insertRows(row,1)
        if self._includeCheckboxes:
          itm = QCheckTableItem(self,'')
          itm.setChecked(1)
          self.setItem(row,self._checkboxCol,itm)
        d = mol._fieldsFromDb
        if not mol.HasProp('_Name') and type(d[0]) in self._displayableTypes:
          mol.SetProp('_Name',str(d[0]))

        if kekulize:
          Chem.Kekulize(mol)
        smi = Chem.MolToSmiles(mol).strip()
        itm = MolTableItem(self,QTableItem.Never,smi)
        itm.setMol(mol)
        itm.setDrawTarget(drawTarget)
        self.setItem(row,self._molCol,itm)
        for col in range(nCols):
          if type(d[col]) in self._displayableTypes:
            self.setText(row,col+nExtraCols,str(d[col]))
          else:
            self.setText(row,col+nExtraCols,"<Undisplayable>")
        row += 1
        if dlg:
          try:
            dlg(row+nFailed)
          except qtUtils.ProgressStop:
            break
      else:
        nFailed += 1
    if nFailed>0:
      nRows = row+nFailed
      rowsToRemove = range(row,nRows)
      for entry in rowsToRemove:
        self.removeRow(self.numRows()-1)

    QApplication.restoreOverrideCursor()
    if adjustCols:
      for i in range(nCols):
        if i!=self._molCol:
          self.adjustColumn(i)
    self.setUpdatesEnabled(1)
    self.redraw()

  def loadFromVLib(self,mols,drawTarget=None,
                   includeCheckboxes=0,adjustCols=1,
                   processFunc=None):
    # clear the current contents:
    self.setNumRows(0)
    self.setNumCols(0)

    mols.reset()
    self._includeCheckboxes = includeCheckboxes
    colNames=mols.next().GetPropNames()
    nCols = len(colNames)
    if not self._includeCheckboxes:
      self._molCol=0
      nExtraCols=1
    else:
      nExtraCols=2
      self._molCol=1
    self.setNumCols(nCols+nExtraCols)
    self.setUpdatesEnabled(0)
    # set column headers
    hdr = self.horizontalHeader()
    hdr.setLabel(self._molCol,'Molecule')
    if self._includeCheckboxes:
      hdr.setLabel(self._checkboxCol,'  ')
      self.adjustColumn(self._checkboxCol)
    for i in range(nCols):
      hdr.setLabel(i+nExtraCols,colNames[i])

    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    if hasattr(mols,'__len__'):
      dlg = QProgressDialog('loading','',len(mols))
      dlg.setLabelText('Loading Compounds')
    else:
      dlg = None

    nFailed = 0
    row = 0
    for mol in mols:
      if mol:
        if processFunc:
          processFunc(mol)

        self.insertRows(row,1)
        if self._includeCheckboxes:
          itm = QCheckTableItem(self,'')
          itm.setChecked(1)
          self.setItem(row,self._checkboxCol,itm)
        smi = Chem.MolToSmiles(mol).strip()
        itm = MolTableItem(self,QTableItem.Never,smi)
        itm.setMol(mol)
        itm.setDrawTarget(drawTarget)
        self.setItem(row,self._molCol,itm)
        if not mol.HasProp('_Name'):
          d = mol.GetProp(colNames[0])
          if type(d) in self._displayableTypes:
            mol.SetProp('_Name',d)
        if not mol.HasProp('_Name') and type(d[0]) in self._displayableTypes:
          mol.SetProp('_Name',d[0])
        for col in range(nCols):
          colN = colNames[col]
          d = mol.GetProp(colN)
          if type(d) in self._displayableTypes:
            self.setText(row,col+nExtraCols,str(d))
          else:
            self.setText(row,col+nExtraCols,"<Undisplayable>")
        row += 1
        if dlg:
          dlg.setProgress(row+nFailed)
      else:
        nFailed += 1
    if nFailed>0:
      rowsToRemove = range(row,nRows)
      for entry in rowsToRemove:
        self.removeRow(self.numRows()-1)

    QApplication.restoreOverrideCursor()
    if adjustCols:
      for i in range(nCols):
        if i!=self._molCol:
          self.adjustColumn(i)
    self.setUpdatesEnabled(1)
    self.redraw()



  def invertChecks(self):
    if not self._includeCheckboxes:
      return
    for row in range(self.numRows()):
      itm = self.item(row,self._checkboxCol)
      itm.setChecked(not itm.isChecked())
  def delCheckedRows(self,invert=0):
    """

      invert=0 means delete checked rows
      invert=1 means delete unchecked rows

    """
    if not self._includeCheckboxes:
      return
    row = 0
    self.setUpdatesEnabled(0)
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    try:
      while row < self.numRows():
        itm = self.item(row,self._checkboxCol)
        checked = itm.isChecked()
        if (checked and not invert) or (not checked and invert):
          self.removeRow(row)
        else:
          row += 1
    except:
      QApplication.restoreOverrideCursor()
      qtUtils.warning('problems encountered deleting a row or rows',exc_info=True)
    QApplication.restoreOverrideCursor()
    self.setUpdatesEnabled(1)
    self.updateContents()

  def highlightDuplicates(self,highlightColor=Qt.cyan):
    if not hasattr(self,'_molCol') or \
       self._molCol < 0 or \
       self._molCol >= self.numCols():
      return
    seen = {}
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    for row in range(self.numRows()):
      itm = self.item(row,self._molCol)
      try:
        smi = str(itm.smiles())
      except:
        pass
      else:
        if seen.has_key(smi):
          # it's a duplicate, flag it and the other member of the pair:
          itm.setFlagged(1,highlightColor)
          itm2 = self.item(seen[smi],self._molCol)
          itm2.setFlagged(1,highlightColor)
        else:
          seen[smi] = row
          itm.setFlagged(0)
    QApplication.restoreOverrideCursor()
    self.updateContents()

  def toggleCheckboxes(self):
    self.setUpdatesEnabled(0)
    if self._includeCheckboxes:
      # remove the column
      self.removeColumn(self._checkboxCol)
      if self._molCol > self._checkboxCol:
        self._molCol -= 1
    else:
      # add the column
      self.insertColumns(self._checkboxCol,1)
      self.horizontalHeader().setLabel(self._checkboxCol,'  ')
      if self._molCol >= self._checkboxCol:
        self._molCol += 1
      for row in range(self.numRows()):
        self.setItem(row,self._checkboxCol,QCheckTableItem(self,''))
      self.adjustColumn(self._checkboxCol)
    self.setUpdatesEnabled(1)
    self.updateContents()
    self._includeCheckboxes = not self._includeCheckboxes

  def exportToMolWriter(self,writer,idCol=-1):
    if idCol < 0:
      idCol = self._molCol+1

    # each column will be used as a property
    header = self.horizontalHeader()
    props = []
    cols = []
    for col in range(header.count()):
      if col not in (self._molCol,self._checkboxCol,idCol):
        props.append(str(header.label(col)))
        cols.append(col)
    writer.SetProps(props)

    for row in range(self.numRows()):
      mol = self.item(row,self._molCol).mol()
      for i in range(len(props)):
        col = cols[i]
        mol.SetProp(props[i],str(self.text(row,col)))
      mol.SetProp('_Name',str(self.text(row,idCol)))
      writer.write(mol)
    writer.flush()  
        
  def toFile(self,fileN=None,idCol=-1,delim='\t'):
    import os.path
    #filters = [('SD File (*.sdf)',lambda x:Chem.SDWriter(x)),
    #           ('CSV File (.*.csv)',lambda x:Chem.Smiles
    #dlg = QFileDialog(
    #fName = QFileDialog.getSaveFileName(
  
  def toDatabase(self,conn=None,idCol=-1):
    pass

  def smilesToClipboard(self,row,col):
    if col != self._molCol:
      return
    itm = self.item(row,col)
    smi = itm.smiles()
    clip = qApp.clipboard()
    clip.setText(smi)

  
  #
  #  signals/slots
  #
  def currentChanged(self,row,col):
    for cb in self._changedCallbacks:
      cb(self,row,col)
    if col == self._molCol:
      itm = self.item(row,col)
      if itm:
        if isinstance(itm,MolTableItem):
          itm.draw()
    
  def contextMenuRequested(self,row,col,pos):
    """ arises from a right click """
    popup = QPopupMenu(self)
    itms = []
    if col == self._molCol:
      fn =lambda a:self.smilesToClipboard(row,col)
      popup.insertItem(self.trUtf8('&Copy to Clipboard'),fn)
      itms.append(fn)

    if self._allowSaves:
      popup.insertSeparator()
      fn =lambda a:self.toFile()
      popup.insertItem(self.trUtf8('Export to &File'),fn)
      itms.append(fn)
      fn =lambda a:self.toDatabase()
      popup.insertItem(self.trUtf8('Export to &Database'),fn)
      itms.append(fn)
      
      
    if self._includeCheckboxes and col == self._checkboxCol:
      popup.insertSeparator()
      fn =lambda a:self.invertChecks()
      popup.insertItem(self.trUtf8('&Invert Checks'),fn)
      itms.append(fn)
      fn =lambda a:self.delCheckedRows(0)
      popup.insertItem(self.trUtf8('Delete &Checked'),fn)
      itms.append(fn)
      fn = lambda a:self.delCheckedRows(1)
      popup.insertItem(self.trUtf8('Delete &Unchecked'),fn)
      itms.append(fn)
    elif col == self._molCol:
      popup.insertSeparator()
      fn =lambda a:self.highlightDuplicates()
      popup.insertItem(self.trUtf8('&Highlight Duplicates'),fn)
      itms.append(fn)
      fn =lambda a:self.toggleCheckboxes()
      itm = popup.insertItem(self.trUtf8('&Checkboxes'),fn)
      popup.setItemChecked(itm,self._includeCheckboxes)
      itms.append(fn)
    if itms:
      popup.exec_loop(pos)
    


def test1():
  import os
  # start out by grabbing some data:
  from Dbase.DbConnection import DbConnect
  from Chem.Suppliers.DbMolSupplier import RandomAccessDbMolSupplier
  from Chem.Suppliers.DbMolSupplier import ForwardDbMolSupplier
  dbName = os.path.join(RDConfig.RDCodeDir,'qtGui','GuiLib','demoData','data.gdb')
  conn = DbConnect(dbName,'some_mols_dupes')
  data = conn.GetData()
  if 1:
    mols = RandomAccessDbMolSupplier(data)
  else:
    mols = ForwardDbMolSupplier(data)


  # build the app and widget
  from qtGui import Gui
  app,widg = Gui.Launcher(MolTable,None)

  # and load up those molecules:
  widg.loadFromMolSupplier(mols,includeCheckboxes=0)

  w = Chem.SmilesWriter('foob.txt',delimiter='\t',nameHeader="Mol_ID")
  widg.exportToMolWriter(w)
  
  app.exec_loop()
  widg.destroy(1)
    
def test2():
  import os
  # start out by grabbing some data:
  from VLib.NodeLib import DbMolSupply
  dbName = os.path.join(RDConfig.RDCodeDir,'qtGui','GuiLib','demoData','data.gdb')
  node = DbMolSupply.GetNode(dbName,'some_mols_dupes')

  # build the app and widget
  from qtGui import Gui
  app,widg = Gui.Launcher(MolTable,None)

  # and load up those molecules:
  widg.loadFromVLib(node,includeCheckboxes=0)

  app.exec_loop()
  widg.destroy(1)
    
if __name__ == '__main__':
  #test1()
  test2()
    
