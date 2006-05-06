#
#  Copyright (C) 2002,2003  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" general utilities for interacting with composites from GUIs

"""
from qt import *
from Dbase.DbConnection import DbConnect
from ML import ScreenComposite

#FIX: I'm not sure that this module should live in this directory

def ScreenCompositeFromDetails(details,composites,reportWin=None,
                               goodVotes=None,badVotes=None,noVotes=None):
  """ Screens a composite using a details structure

#DOC: list of composites

  **Arguments**

    - details: a _CompositeRun.CompositeRun_ structure containing the
      screening details

    - composite: the composite model to be screened

    - reportWin: (optional) if provided, this
      _qtGui.GuiTextViewer.GuiTextViewer_ window will be filled with
      the HTML screening details.
  
    - goodVotes,badVotes,noVotes: (optional)  if provided these should
      be lists (or anything supporting an _append()_ method) which
      will be used to pass the screening results back.

  **Notes**

    - heavy lifting is done by code in _ML.ScreenComposite_:

      - _ScreenFromDetails()_ does the screening itself

      - _ScreenToHtml()_ produces the HTML, if required
    
  """
  # get the database's column names
  conn = DbConnect(details.dbName,details.tableName,user=details.dbUser,password=details.dbPassword)
  names = conn.GetColumnNames(table=details.tableName,join=details.dbJoin,what=details.dbWhat)

  if type(composites) not in (type((1,)),type([])):
    composites = (composites,)
  for composite in composites:
    composite.ClearModelExamples()
    composite.SetInputOrder(names)

  nPts = conn.GetDataCount(table=details.tableName,join=details.dbJoin,where=details.dbWhere)
  QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
  dlg = QProgressDialog('Screen Progress','',0)
  dlg.setLabelText("Screening Progress")

  setup = dlg.setTotalSteps
  try:
    screenRes=ScreenComposite.ScreenFromDetails(composites,details,
                                                callback=dlg.setProgress,appendExamples=1,
                                                setup=setup,
                                                goodVotes=goodVotes,badVotes=badVotes,
                                                noVotes=noVotes)
  except:
    import traceback
    traceback.print_exc()
    screenRes=None
  else:
    # do the visual display:
    if reportWin is not None:
      nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = screenRes
      html = ScreenComposite.ScreenToHtml(nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl)
      reportWin.setText(html)
      reportWin.show()
  QApplication.restoreOverrideCursor()
  return screenRes
