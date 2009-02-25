#
#  Copyright (C) 2003  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" implementation bits for Pubmed searching

"""    
import RDConfig
from qt import *
from qtGui.GuiLib.forms.PubmedRecordWidget import PubmedRecord as _Form
from qtGui import qtUtils
import os

class PubmedRecord(_Form):
  """ DOC
  

  """
  def __init__(self,*args,**kwargs):
    if kwargs.has_key('record'):
      rec = kwargs['record']
      del kwargs['record']
    else:
      rec = None
    _Form.__init__(self,*args,**kwargs)
    self._rec = rec
    if self._rec:
      self.update()
    self._dir = '.'
    
  def setRecord(self,rec):
    self._rec = rec
    self.update()
  def getRecord(self):
    return self._rec
  def update(self):
    rec = self._rec
    self.titleBox.setText(rec.Title)
    self.titleBox.setCursorPosition(0)
    self.authorBox.setText(rec.Authors)
    self.authorBox.setCursorPosition(0)
    self.sourceBox.setText(rec.Source)
    self.volumeBox.setText(rec.Volume)
    self.pageBox.setText(rec.Pages)
    self.yearBox.setText(rec.PubYear)
    self.idBox.setText(rec.PubMedId)
    if rec.Abstract:
      self.abstractBox.setEnabled(1)
      self.abstractViewer.setText(rec.Abstract)
    else:
      self.abstractBox.setEnabled(0)

    if rec.keywords:
      self.keywordBox.setEnabled(1)
      self.keywordViewer.setText('\n'.join(rec.keywords))
    else:
      self.keywordBox.setEnabled(0)
    if rec.chemicals:
      self.chemicalBox.setEnabled(1)
      self.chemicalViewer.setText('\n'.join(rec.chemicals))
    else:
      self.chemicalBox.setEnabled(0)

  def getXML(self):
    return self._rec.toXML()

  # ---------
  #  Signals and Slots
  #
  def saveButtonClicked(self,fileName=''):
    if not fileName:
      fileName = str(QFileDialog.getSaveFileName(self._dir,
                                                 'XML files (*.xml);;All files (*.*)'))
      if fileName:
        self._dir = os.path.dirname(fileName)
    if not fileName:
      return
    else:
      outF = open(fileName,'w+')
      outF.write(self.getXML())
      outF.close()
      
if __name__ == '__main__':
  import os.path
  from Dbase.Pubmed import Searches
  #fName = os.path.join(RDConfig.RDCodeDir,'Dbase','Pubmed','test_data','records.xml')
  #ids=['11960484']
  fName = os.path.join(RDConfig.RDCodeDir,'Dbase','Pubmed','test_data','encoding_author_record.xml')
  ids=['12653536']
  inF = open(fName,'r')
  recs = Searches.GetRecords(ids,conn=inF)
  from qtGui import Gui

  app,widg = Gui.Launcher(PubmedRecord,None,'Record')
  widg.setRecord(recs[0])
  app.exec_loop()

