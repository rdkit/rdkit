# $Id$
#
#  Copyright (C) 2003-2005  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" implementation bits for SignatureViewDialogs, which are for working with
 2D pharmacophore signatures

"""    
import RDConfig
from qt import *
from qttable import *
from qtGui.GuiLib.forms.SignatureViewDialog import SignatureViewDialog as _Form
from qtGui.DbQueryWidgetImpl import insertQueryWidget
from qtGui.GuiLib import SignatureTable
from qtGui import qtUtils
import DataStructs
import cPickle
import Chem


class SignatureViewDialog(_Form):
  """ class used to handle viewing signatures
  

  """
  def __init__(self,parent=None,initDir=''):
    _Form.__init__(self,parent)
    self._dbDir = initDir
    self._colNames = None

    # add the dbQuery widget to the db page
    fn = lambda y=self.buttonRefresh:y.setEnabled(1)
    self.queryWidget = insertQueryWidget(self.dbPage,
                                         clickCallback=fn)
    self.queryWidget.show()
    self.tabWidget.setTabEnabled(self.sigPage,0)
    self._drawTarget = None
    self.sigTable = None
    
  def setDrawTarget(self,tgt):
    self._drawTarget = tgt
  def drawTarget(self):
    return self._drawTarget
  def loadInfoFromDb(self):
    """ #DOC

    """
    w = self.queryWidget
    conn = w.getConn()
    colNames=[x.upper() for x in conn.GetColumnNames(what=w.sqlWhat(),join=w.sqlJoin())]
    if 'SMI' in colNames:
      smiIdx = colNames.index('SMI')
    elif 'SMILES' in colNames:
      smiIdx = colNames.index('SMILES')
    else:
      smiIdx = -1

    if smiIdx < 0:
      qtUtils.error('no SMILES or SMI column found')
      return
      
    if 'SIG' in colNames:
      sigIdx = colNames.index('SIG')
    elif 'SIGNATURE' in colNames:
      sigIdx = colNames.index('SIGNATURE')
    else:
      sigIdx = -1
    if sigIdx < 0:
      qtUtils.error('no SIGNATURE or SIG column found')
      return

    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    data = conn.GetData(fields=w.sqlWhat(),where=w.sqlWhere(),join=w.sqlJoin())
    nRows = len(data)
    if not nRows:
      return 0

    nCols = len(colNames)
    if self.sigTable is None:
      self.sigTable = SignatureTable.insertSignatureTable(self.sigPage,addInPos=0,
                                                          clearExisting=0)
      #self.sigTable = SignatureTable.insertSignatureTable(self.sigPage)
    else:
      self.sigTable.setNumRows(0)
      self.sigTable.setNumCols(0)
      
    tbl = self.sigTable
    tbl.setNumRows(nRows)
    tbl.setNumCols(1)

    tbl.setUpdatesEnabled(0)
    # set column headers
    hdr = tbl.horizontalHeader()
    hdr.setLabel(0,'Molecule')

    dlg = QProgressDialog('loading','',nRows)
    dlg.setLabelText('Loading Signatures')
    for row in xrange(nRows):
      d = data[row]
      smi = d[smiIdx]
      pkl = str(d[sigIdx])
      try:
        sig = cPickle.loads(pkl)
      except cPickle.UnpicklingError:
        QApplication.restoreOverrideCursor()
        qtUtils.error("Could not depickle the signature.")
        break
      mol = Chem.MolFromSmiles(smi)
      itm = SignatureTable.SignatureTableItem(tbl,QTableItem.Never,smi)
      itm.setMol(mol)
      try:
        itm.setSig(sig)
      except ValueError:
        QApplication.restoreOverrideCursor()
        qtUtils.error("Bad signature type detected.\nThis widget is intended for 2D pharmacophore signatures.")
        break
        
      # FIX: update this
      itm.setDrawTarget(self.drawTarget())
      tbl.setItem(row,0,itm)

      dlg.setProgress(row)  
    tbl.setUpdatesEnabled(1)
    self.tabWidget.setTabEnabled(self.sigPage,1)
    QApplication.restoreOverrideCursor()
    return nRows

  def refreshIt(self):
    """ updates our current signatures

      #DOC
    
    """
    self.loadInfoFromDb()

  def addBit(self,bit):
    tbl = self.sigTable
    if not tbl: return

    tbl.setUpdatesEnabled(0)
    idx = tbl.numCols()
    tbl.insertColumns(idx)
    hdr = tbl.horizontalHeader()
    hdr.setLabel(idx,str(bit))

    for i in xrange(tbl.numRows()):
      item = tbl.item(i,0)
      sig = item.sig()
      mol = item.mol()
      label = str(sig[bit])
      newItem = SignatureTable.SignatureTableItem(tbl,QTableItem.Never,label)
      newItem.setSig(sig)
      newItem.setBit(bit)
      newItem.setMol(mol)
      newItem.setDrawTarget(self.drawTarget())
      tbl.setItem(i,idx,newItem)
    tbl.setUpdatesEnabled(1)
      
    
    

  #-------------------
  #
  #  Slots
  # 
  #-------------------
  def refreshClick(self):
    """ callback """
    self.refreshIt()
    
  def addBitClick(self):
    """ callback """
    idx = int(str(self.sig_bitLine.text()))
    self.addBit(idx)
    
  def loadBitsClick(self):
    """ callback """
    qtGui.logger.debug('load bits')
    

if __name__ == '__main__':
  from qtGui import Gui

  app,widg = Gui.Launcher(SignatureViewDialog,None,'SigView')
  widg.exec_loop()

