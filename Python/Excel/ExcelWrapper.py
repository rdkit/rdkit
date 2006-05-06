# $Id$
#
#  Copyright (C) 2000-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" a wrapper for interacting with Excel

"""
from win32com.client import Dispatch, constants
import types

class ExcelWrapper(object):
  """ This is essentially a convenience class to wrap up Excel so that it's
   more friendly (from my point of view).  It is neither particularly pretty nor
   pythonic(TM), but it is VERY useful.

  """

  def GetHeadings(self,fromRow=1,startCol=1):
    """ gets a list of column names

      **Arguments**
              
        - fromRow: the row containing column names

        - startCol: the first column in which to look

      **Returns**

        a list of column names

    """
    sheet = self.xlApp.ActiveSheet
    lastCol = self.FindLastCol(startCol,fromRow)
    t = self[fromRow,startCol:lastCol]
    if lastCol-startCol>0:
      res = map(lambda x: str(x),t[0])
    else:
      res = [str(t)]
    return res


  def SelectCol(self,col,replace=1):
    """ selects a column in Excel

      **Arguments**

        - col: which column to select

        - replace: not currently used
     
    """
    sheet = self.xlApp.ActiveSheet
    colLab = self.NumToLetter(col)
    rowTop = 1
    rowBot = self.FindLastRow(rowTop,col)
    rangeStr='%s%d:%s%d'%(colLab,rowTop,colLab,rowBot)
    sheet.Range(rangeStr).Select()
    
  def ScatterPlotCol(self,col1,col2=None,name=None):
    """Constructs a scatter plot from one or two columns in the active sheet

      **Arguments**

        - col1: the first column to be used

        - col2: if present provides the second set of data for the scatter plot

        - name: an optional name for the plot  If this is not present, the
          heading of _col1_ will be used

      **Notes**

        - existing contents of chart _name_ are blown away

    """
    sheet = self.xlApp.ActiveSheet
    if name is None:
      name = sheet.Cells(1,col1).Value
    
    col1Name = self[1,col1]
    if col2 is not None:
      col2Name = self[1,col2]
    rowTop = 2
    rowBot = self.FindLastRow(rowTop,col1)


    rng = self.GetRange(rowTop,rowBot,col1,col1)
    if col2 is not None:
      rng2 = self.GetRange(rowTop,rowBot,col2,col2)
    else:
      rng2 = None

    charts = self.ActiveWorkbook.Charts
    self.DisplayAlerts=0
    try:
      charts(name).Delete()
    except:
      pass
    self.DisplayAlerts=1
    charts.Add()
    chart = self.ActiveChart
    chart.Name = name
    chart.ChartType = constants.xlXYScatter
    #chart.SetSourceData(Source=rng,PlotBy=constants.xlColumns)
    chart.SeriesCollection(1).Values=rng
    if rng2 is not None:
      chart.SeriesCollection(1).XValues=rng2
    chart.HasLegend=0
    chart.Axes(constants.xlValue,constants.xlPrimary).HasTitle=1
    chart.Axes(constants.xlValue,constants.xlPrimary).AxisTitle.Caption=col1Name
    chart.Axes(constants.xlCategory,constants.xlPrimary).HasTitle=1
    chart.Axes(constants.xlCategory,constants.xlPrimary).AxisTitle.Caption=col2Name
    
    chart.Activate()
    

  def FindLastRow(self,rowStart,colStart):
    """ Finds the last row in the active sheet that contains data

      **Arguments**

        - rowStart: the row in which the search should be started

        - colStart: the column in which to look

    """
    sheet = self.xlApp.ActiveSheet
      
    sheet.Cells.Item(rowStart,colStart).Activate()
    region = sheet.Cells.CurrentRegion
    lastRow = region.Rows.Count
    return lastRow

  def FindLastCol(self,colStart,rowStart):
    """ Finds the last row in the active sheet that contains data

      **Arguments**

        - colStart: the column in which the search should be started

        - rowStart: the row in which to look

    """
    sheet = self.xlApp.ActiveSheet

    sheet.Cells.Item(rowStart,colStart).Activate()
    region = sheet.Cells.CurrentRegion
    lastCol = region.Columns.Count
    return lastCol

  def NumToLetter(self,num):
    """ converts a number to a letter

      **Arguments**

        - num: an integer > 0

      **Returns**

        the corresponding letter(s)

      **Notes**
        - hi, I'm a gross hack and only deal with arrays up to 676 (26*26) wide...
          This makes me want to cry, but the 30 minutes which I spent applying
           half my brain to the problem didn't yield anything general.

        - it turns out that excel appears to be unhappy with tables more than
         256 wide anyway. :-)

        - FIX: it's possible to design this entire thing out by being smarter
          about the way ranges are selected... this should be done at some
          point.

    """    
    assert num>0,'Cell indexing should start from one'
    num = num -1
    val = ''
    if num < 26:
      val = chr(num+ord('a')) + val
      num = 0
    elif num <= 26*26:
      val = chr(num/26-1+ord('a'))+chr(num%26+ord('a')) + val
    else:
      raise ValueError,'value %d too large'%num
    return val

  def GetRange(self,rowStart,rowEnd,colStart,colEnd):
    """ Returns an Excel *Range* object for the range specified

      These can be used to set properties, values of ranges, etc.
      
      **Arguments**

         - rowStart: first row in the range

         - rowEnd: last row in the range

         - colStart: first col in the range

         - colEnd: last col in the range

      **Returns**

         an Excel *Range* object (this is a COM object)
         
    """
    c1 = self.xlApp.Cells(rowStart,colStart)
    c2 = self.xlApp.Cells(rowEnd,colEnd)
    return self.xlApp.ActiveSheet.Range(c1,c2)
    
  def Destroy(self):
    """ wipes out the wrapper

    """
    self.xlApp = None
    self.__dict__ = None
    
  def __getitem__(self,which):
    """ allows indexing

    """
    slicing = 0
    if type(which[0]) == types.SliceType:
      rowStart = which[0].start
      rowEnd = which[0].stop
      slicing = 1
    else:
      rowStart = which[0]
      rowEnd = which[0]
    if type(which[1]) == types.SliceType:
      colStart = which[1].start
      colEnd = which[1].stop
      slicing = 1
    else:
      colStart = which[1]
      colEnd = which[1]

    if slicing == 0:
      return self.xlApp.ActiveSheet.Cells(which[0],which[1]).Value
    else:
      if rowEnd == -1 or rowEnd is None:
        rowEnd = self.FindLastRow(rowStart,colStart)
      if colEnd == -1 or colEnd is None:
        colEnd = self.FindLastCol(colStart,rowStart)
      if rowStart is None:
        rowStart = 0
      if colStart is None:
        colStart = 0
      rangeStr = '%s%d:%s%d'%(self.NumToLetter(colStart),rowStart,
                              self.NumToLetter(colEnd),rowEnd)

      return self.xlApp.ActiveSheet.Range(rangeStr).Value

  def __setitem__(self,which,val):
    """ allows indexed setting of values

    """
    slicing = 0
    if type(which[0]) == types.SliceType:
      rowStart = which[0].start
      rowEnd = which[0].stop
      slicing = 1
    else:
      rowStart = which[0]
      rowEnd = which[0]
    if type(which[1]) == types.SliceType:
      colStart = which[1].start
      colEnd = which[1].stop
      slicing = 1
    else:
      colStart = which[1]
      colEnd = which[1]

    if slicing == 0:
      which = list(which)
      if which[0] == -1:
        which[0] = self.FindLastRow(1,1)+1
      if which[1] == -1:
        which[1] = self.FindLastCol(1,1)+1

      self.xlApp.ActiveSheet.Cells(which[0],which[1]).Value=val
    else:
      if rowEnd == -1:
        rowEnd = self.FindLastRow(rowStart,colStart)
      if colEnd == -1:
        colEnd = self.FindLastCol(colStart,rowStart)
      rangeStr = '%s%d:%s%d'%(self.NumToLetter(colStart),rowStart,
                              self.NumToLetter(colEnd),rowEnd)

      self.xlApp.ActiveSheet.Range(rangeStr).Value = val

  def __getattr__(self,attrName):
    """ this is a trick to allow attribute access for both this object
       and the Excel object it wraps.

     **Arguments**

       - attrName: the name of the attribute to get
       
     **Notes**

       -Excel takes priority, so ask it first, then do our own local thing

    """
    try:
      exec('res = self.__dict__["xlApp"].%s'%attrName)
    except KeyError:
      res = self.__dict__[attrName]
    return res

  def __setattr__(self,attrName,value):
    """ if the Excel application has the attribute, try to set it there
      before we set our local copy

      **Arguments**

        - attrName: the attribute to set

        - value: the new value

     **Notes**

       -Excel takes priority, so ask it first, then do our own local thing

    """
    try:
      exec('temp=self.__dict__["xlApp"].%s'%attrName)
    except AttributeError:
      self.__dict__[attrName] = value
    else:
      exec('self.__dict__["xlApp"].%s = value'%attrName)


  def __init__(self,app=None):
    """ Constructor

       **Arguments**

         - app: an Excel _Dispatch_ object 
    """
    if app is None:
      app = Dispatch('Excel.Application')
    self.__dict__["xlApp"] = app
    

if __name__ == '__main__':
  # test code
  foo = ExcelWrapper()
  foo.Visible = 1
  foo.Workbooks.Add()
  foo[1,1] = 11
  foo[2,1] = 21
  foo[3,1] = 31
  foo[1,2] = 12
  foo[2,2] = 22
  foo[1:3,1:2] = [[101,201],[102,202],[103,203]]
  print foo[1:3,1]
  print foo[1,1:-1]
  print foo[1:3,1:-1]
  # don't clutter up excel whilst we are testing
  #foo.ActiveWorkbook.Close(SaveChanges=0)
