# $Id: Molecules.py 5011 2006-02-22 15:24:33Z glandrum $
#
#  Copyright (C) 2002-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" mixin class for dealing with ChemDraw for Excel

"""
import time
from Excel.ExcelWrapper import ExcelWrapper
_commandBarName='ChemDraw'

# controls we use
_newWorksheetName='New ChemDraw Worksheet'
_convertName='Convert SMILES To Molecule'
_showPictureName='Show Picture'
_hidePictureName='Hide Picture'

def _CDNewWorksheet(self,name=''):
  bar = self.xlApp.CommandBars(_commandBarName)
  bar.Controls(_newWorksheetName).Execute()
  if name: self.xlApp.ActiveSheet.Name=name

def _CDRangeOp(self,op,row,col,endRow=-1,endCol=-1):
  whichTry=0
  while whichTry<5:
    # we've got this little loop here because Excel has shown
    #  a tendency to be somewhat cranky about selections when
    #  contents are fired in quickly
    try:
      if endRow < 0:
        self.xlApp.Cells(row,col).Select()
      else:
        self.xlApp.Range(self.xlApp.Cells(row,col),
                         self.xlApp.Cells(endRow,endCol)).Select()
    except:
      whichTry+=1
      time.sleep(.1)
    else:
      break
  bar = self.xlApp.CommandBars(_commandBarName)

  # carry out the operation
  bar.Controls(op).Execute()
  
def _CDConvertCellsToMols(self,row,col,endRow=-1,endCol=-1):
  _CDRangeOp(self,_convertName,row,col,endRow,endCol)
  
def _CDShowMols(self,row,col,endRow=-1,endCol=-1):
  _CDRangeOp(self,_showPictureName,row,col,endRow,endCol)
  
def _CDHideMols(self,row,col,endRow=-1,endCol=-1):
  _CDRangeOp(self,_hidePictureName,row,col,endRow,endCol)
  

# on import, install the methods
ExcelWrapper.ChemdrawNewWorksheet=_CDNewWorksheet
ExcelWrapper.ChemdrawConvertCellsToMols=_CDConvertCellsToMols
ExcelWrapper.ChemdrawShowPictures=_CDShowMols
ExcelWrapper.ChemdrawHidePictures=_CDHideMols

if __name__ == '__main__':
  smis = ['CCOC','c1ccccc1','CC(=O)O','C1C(=O)CC1']
  w = ExcelWrapper()
  w.Workbooks.Add()
  w.ChemdrawNewWorksheet()
  for i in range(len(smis)):
    w[i+1,1]=smis[i]
  w.ChemdrawConvertCellsToMols(1,1)
  w.ChemdrawConvertCellsToMols(2,1,len(smis),1)
  w.ChemdrawShowPictures(1,1,len(smis),1)
  w.ChemdrawHidePictures(2,1)
