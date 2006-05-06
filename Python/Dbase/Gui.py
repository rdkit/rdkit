# $Id$
#
#  Copyright (C) 2000-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" simple GUI tools for interacting with databases

"""
raise NotImplementedError,'wxPython no longer supported'
if 0:
  import RDConfig
  from GuiFramework import GuiHistoryField, GuiIDs, GuiHtml, GuiDialog
  from wxPython.wx import *
  from wxPython.lib.grids import wxGridSizer
  from wxPython.lib.filebrowsebutton import FileBrowseButton
  from Dbase import DbInfo
  import gvib
  import string
  
  defaultDBase=RDConfig.RDDataDir+"/Toxicity/BBB/BBB.GDB"
  defaultTable="original_data"

  class DbTableSelectDialog(GuiDialog.GuiDialog):
    def __init__(self,parent,id,name='Select Database and Table'):
      GuiDialog.GuiDialog.__init__(self,parent,id,name)
      self.SetBackgroundColour(wxWHITE)
      self._browseButton = FileBrowseButton(self,-1,
                                            fileMask='*.gdb|*.GDB',
                                            size=wxSize(200,-1),
                                            labelText='Database:',
                                            buttonText='Browse',
                                            changeCallback=lambda x,y=self:y._GetTableNames())
      self.AddWidget(self._browseButton)
      self._szr = wxBoxSizer(wxHORIZONTAL)
      self._tblText=wxStaticText(self,-1,'Table:')
      self._szr.Add(self._tblText,0,wxCENTER)
      self._choice = wxChoice(self,-1)
      self._choice.Enable(0)
      self._szr.Add(self._choice,1,wxLEFT|wxCENTER,30)
      self.AddWidget(self._szr)

      self._uSzr = wxBoxSizer(wxHORIZONTAL)
      self._uText=wxStaticText(self,-1,'User:')
      self._uSzr.Add(self._uText,0,wxCENTER)
      self._userEntry = wxTextCtrl(self,-1)
      self._userEntry.SetValue('sysdba')
      self._uSzr.Add(self._userEntry,1,wxLEFT|wxCENTER,30)
      self.AddWidget(self._uSzr)

      self._pSzr = wxBoxSizer(wxHORIZONTAL)
      self._pText=wxStaticText(self,-1,'Password:')
      self._pSzr.Add(self._pText,0,wxCENTER)
      self._passwordEntry = wxTextCtrl(self,-1,style=wxTE_PASSWORD)
      self._passwordEntry.SetValue('masterkey')
      self._pSzr.Add(self._passwordEntry,1,wxLEFT|wxCENTER,30)
      self.AddWidget(self._pSzr)

    def _GetTableNames(self):
      try:
        dbName = self._browseButton.GetValue()
      except AttributeError:
        return
      self._choice.Clear()
      self._choice.Enable(0)
      if not dbName:
        return
      tblNames = DbInfo.GetTableNames(dbName)
      for name in tblNames:
        self._choice.Append(name)
      self._choice.Enable(1)

    def GetDbName(self):
      return self._browseButton.GetValue()
    def GetTableName(self):
      return self._choice.GetStringSelection()
    def GetUserName(self):
      return self._userEntry.GetValue()
    def GetPassword(self):
      return self._passwordEntry.GetValue()
    def SetDbName(self,val):
      self._browseButton.SetValue(val)
    def SetTableName(self,val):
      return self._choice.SetStringSelection(val)


  class DBGui(wxFrame):
    """ a simple (basically demo level) GUI for querying databases
      and seeing the results

    """
    def ConstructHtmlDataTable(self,res,headers,binCols=[]):
      """ constructs an HTML page/table for the contents of _res_

        the results are displayed in the local _htmlWin_

        **Arguments**

          - res: a list of lists containing the query results

          - headers: names of the columns

          - binCols: a list of columns containing binary data (BLOBS)

      """         
      htmlStr = '<html><head><title>Query Results</title><head>\n'
      htmlStr = htmlStr + '<body><table border=2>\n'
      htmlStr = htmlStr + '<tr>'
      for i in xrange(len(headers)):
        if i not in binCols:
          htmlStr = htmlStr + '<th>%s</th>'%headers[i]
      htmlStr = htmlStr + '</tr>\n'
      lines = []
      for entry in res:
        lines.append('<tr>')
        for i in xrange(len(entry)):
          if i not in binCols:
            lines.append('<td>%s</td>'%entry[i])
        lines.append('</tr>')
      htmlStr = htmlStr + string.join(lines,'\n')
      htmlStr = htmlStr + '</table></body></html>\n'
      self.htmlWin.SetPage(htmlStr)
      self.htmlWin.Show(1)

    def Execute(self,user='sysdba',password='masterkey',showTable=0):
      """ executes the SQL command and calls _ConstructHtmlDataTable()_

        the query is pulled from the form data

        **Arguments**

          - user: the user name to use with the DB

          - password: the password to use with the DB

      """   
      db = self.dbBrowseButton.GetValue()
      table = self.dbTableChoice.GetValue()
      com = self.comCtrl.GetValue()
      join = self.joinCtrl.GetValue()
      where = self.whereCtrl.GetValue()
      command = 'select %s from %s'%(com,table)
      if join != '':
        command = command + ' join %s'%join
      if where != '':
        command = command + ' where %s'%where
      print 'command: ', command
      cn = gvib.connect(db,user,password)
      c = cn.cursor()
      c.execute(command)
      desc = c.description
      results = c.fetchall()
      binCols = []
      headers = []
      for i in xrange(len(desc)):
        item = desc[i]
        headers.append(item[0])
        if item[1] == gvib.gvibTypes.SQL_BLOB:
          binCols.append(i)

      if showTable:
        self.ConstructHtmlDataTable(results,headers,binCols)

      if self.callback:
        self.callback(headers,results)


    def BrowseColNames(self,dbName='',table='',user='sysdba',password='masterkey'):
      """  launches a dialog to allow columns to be selected

        **Arguments**

          - dbName: the DB to be queried

          - table: the table to use

          - user: the user name to use with the DB

          - password: the password to use with the DB


      """
      if dbName == '':
        dbName = self.dbBrowseButton.GetValue()
      if table == '':
        table = self.dbTableChoice.GetValue()
      join = self.joinCtrl.GetValue()
      colData = DbInfo.GetColumnNamesAndTypes(dbName,table,user,password,
                                               join=join)

      dlg = GuiDialog.GuiDialog(self,-1,'Select Column Names')
      colBox = wxListBox(dlg,-1,style=wxLB_EXTENDED|wxLB_NEEDED_SB)
      for colName in colData:
        colBox.Append(colName[0])
      dlg.AddWidget(colBox)
      dlg.Finalize()

      res = dlg.ShowModal()
      if res == wxID_OK:
        selNums = colBox.GetSelections()
        selNames = []
        for num in selNums:
          selNames.append(colData[num][0])
        self.comCtrl.SetValue(string.join(selNames,','))

    def GetTableNames(self,dbName='',user='sysdba',password='masterkey'):
      if dbName == '':
        dbName = self.dbBrowseButton.GetValue()
      cn = gvib.connect(dbName,user,password)
      c = cn.cursor()
      comm = """SELECT RDB$RELATION_NAME
        FROM RDB$RELATIONS
        WHERE ((RDB$SYSTEM_FLAG = 0) OR 
        (RDB$SYSTEM_FLAG IS NULL)) AND 
        (RDB$VIEW_SOURCE IS NULL)
        ORDER BY RDB$RELATION_NAME"""
      c.execute(comm)
      names = c.fetchall()
      for name in names:
        self.dbTableChoice.Append(name[0])
      self.sizer.Fit(self)
    def InitFrame(self):
      """ initialization

      """
      self.SetBackgroundColour(wxWHITE)
      self.sizer = wxBoxSizer(wxVERTICAL)
      toAdd = []
      self.dbBrowseButton = FileBrowseButton(self,-1,
                                             fileMask='*.gdb|*.GDB',
                                             labelText='Database:',
                                             buttonText='Browse')

      self.dbBrowseButton.SetValue(defaultDBase)
      self.sizer.Add(self.dbBrowseButton,0,wxEXPAND,10)


      self.choiceSizer=wxBoxSizer(wxHORIZONTAL)
      self.dbText = wxStaticText(self,-1,'Table:')
      self.choiceSizer.Add(self.dbText,0,wxCENTER)
      self.dbTableChoice = GuiHistoryField.GuiHistoryField(self,-1,size=(200,-1),choices=[defaultTable])
      self.dbTableChoice.SetValue(defaultTable)
      self.choiceSizer.Add(self.dbTableChoice,1,wxLEFT|wxCENTER,30)
      id = GuiIDs.GetNext()
      self.dbTableButton = wxButton(self,id,'Browse',size=(75,-1))
      EVT_BUTTON(self,id,lambda x,y=self:y.GetTableNames())
      self.choiceSizer.Add(self.dbTableButton,1,wxLEFT|wxCENTER,30)


      self.sizer.Add(self.choiceSizer,1,wxEXPAND|wxALL,10)


      self.comSizer=wxBoxSizer(wxHORIZONTAL)
      self.comText = wxStaticText(self,-1,'What:')
      self.comSizer.Add(self.comText,0,wxCENTER)
      self.comCtrl = GuiHistoryField.GuiHistoryField(self,-1,size=(200,-1),
                                                     choices=['*'])
      self.comCtrl.SetValue('*')
      self.comSizer.Add(self.comCtrl,1,wxLEFT|wxCENTER,30)
      id = GuiIDs.GetNext()
      self.comButton = wxButton(self,id,'Browse',size=(75,-1))
      EVT_BUTTON(self,id,lambda x,y=self:y.BrowseColNames())
      self.comSizer.Add(self.comButton,0,wxLEFT|wxCENTER,20)
      self.sizer.Add(self.comSizer,1,wxEXPAND|wxALL,30)

      self.joinSizer=wxBoxSizer(wxHORIZONTAL)
      self.joinText = wxStaticText(self,-1,'Join:')
      self.joinSizer.Add(self.joinText,0,wxCENTER)
      self.joinCtrl = GuiHistoryField.GuiHistoryField(self,-1,size=(200,-1))
      self.joinCtrl.SetValue('')
      self.joinSizer.Add(self.joinCtrl,1,wxLEFT|wxCENTER,30)
      self.sizer.Add(self.joinSizer,1,wxEXPAND|wxALL,30)


      self.whereSizer=wxBoxSizer(wxHORIZONTAL)
      self.whereText = wxStaticText(self,-1,'Where:')
      self.whereSizer.Add(self.whereText,0,wxCENTER)
      self.whereCtrl = GuiHistoryField.GuiHistoryField(self,-1,size=(200,-1))
      self.whereCtrl.SetValue('')
      self.whereSizer.Add(self.whereCtrl,1,wxLEFT|wxCENTER,30)
      self.sizer.Add(self.whereSizer,1,wxEXPAND|wxALL,30)


      self.sizer.Add(60,20,0,wxEXPAND)
      id = GuiIDs.GetNext()
      self.goButton = wxButton(self,id,'Execute')
      EVT_BUTTON(self,id,lambda x,y=self:y.Execute(showTable=y.showTable))
      self.sizer.Add(self.goButton,0,wxEXPAND,10)


      #self.sizer.AddMany(toAdd)
      self.sizer.Fit(self)
      self.SetAutoLayout(true)
      self.SetSizer(self.sizer)


    def OnClose(self,event):
      """ callback for when the window is closed

      """
      if self.htmlWin:
        self.htmlWin.Destroy()
      wxFrame.Destroy(self)

    def Destroy(self):
      if self.htmlWin:
        self.htmlWin.Destroy()
      wxFrame.Destroy(self)

    def __init__(self,*args,**kwargs):
      """ Constructor

        **Arguments**

          all arguments are passed to *wxFrame.__init__*, so they
            should be appropriate for that.

      """
      if kwargs.has_key('searchCallback'):
        self.callback = kwargs['searchCallback']
        del kwargs['searchCallback']
      else:
        self.callback = None
      if kwargs.has_key('showTable'):
        self.showTable = kwargs['showTable']
        del kwargs['showTable']
      else:
        self.showTable = 1
      apply(wxFrame.__init__,(self,)+args,kwargs)
      self.InitFrame()
      self.htmlWin = GuiHtml.HtmlFrame(title='Query Results')
      EVT_CLOSE(self,self.OnClose)


  if __name__=='__main__':
    class TestApp(wxApp):
      def OnInit(self):
        self.frame = DBGui(None,-1,'DBGui')
        return true
    app = TestApp(0)
    app.frame.Show(1)
    app.MainLoop()
  

