# $Id$
#
# Copyright (C) 2002-2006 Greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" tools for interacting with chemdraw

"""

import os
import tempfile
import time

try:
  import pythoncom
  import win32com.client.gencache
  from win32com.client import Dispatch, constants, gencache
  cdxModule = win32com.client.gencache.EnsureModule("{5F646AAB-3B56-48D2-904C-A68D7989C251}", 0, 7,
                                                    0)
except Exception:
  cdxModule = None
  _cdxVersion = 0
  raise ImportError("ChemDraw version (at least version 7) not found.")
else:
  _cdxVersion = 7

if cdxModule:
  from win32com.client import Dispatch
  import win32gui

import re

cdApp = None
theDoc = None
theObjs = None
selectItem = None
cleanItem = None
centerItem = None


def StartChemDraw(visible=True, openDoc=False, showDoc=False):
  """ launches chemdraw """
  global cdApp, theDoc, theObjs, selectItem, cleanItem, centerItem
  if cdApp is not None:
    # if called more than once, do a restart
    holder = None
    selectItem = None
    cleanItem = None
    centerItem = None
    theObjs = None
    theDoc = None
    cdApp = None

  cdApp = Dispatch('ChemDraw.Application')
  if openDoc:
    theDoc = cdApp.Documents.Add()
    theObjs = theDoc.Objects
  else:
    theDoc = None
  selectItem = cdApp.MenuBars(1).Menus(2).MenuItems(8)
  cleanItem = cdApp.MenuBars(1).Menus(5).MenuItems(6)
  if _cdxVersion == 6:
    centerItem = cdApp.MenuBars(1).Menus(4).MenuItems(1)
  else:
    centerItem = cdApp.MenuBars(1).Menus(4).MenuItems(7)
  if visible:
    cdApp.Visible = 1
    if theDoc and showDoc:
      theDoc.Activate()


def ReactivateChemDraw(openDoc=True, showDoc=True):
  global cdApp, theDoc, theObjs
  cdApp.Visible = 1
  if openDoc:
    theDoc = cdApp.Documents.Add()
    if theDoc and showDoc:
      theDoc.Activate()
    theObjs = theDoc.Objects


# ------------------------------------------------------------------
#  interactions with Chemdraw
# ------------------------------------------------------------------
def CDXConvert(inData, inFormat, outFormat):
  """converts the data passed in from one format to another

    inFormat should be one of the following:
       chemical/x-cdx                   chemical/cdx 
       chemical/x-daylight-smiles       chemical/daylight-smiles 
       chemical/x-mdl-isis              chemical/mdl-isis 
       chemical/x-mdl-molfile           chemical/mdl-molfile 
       chemical/x-mdl-rxn               chemical/mdl-rxn 
       chemical/x-mdl-tgf               chemical/mdl-tgf 
       chemical/x-questel-F1  
       chemical/x-questel-F1-query 

    outFormat should be one of the preceding or:
       image/x-png                      image/png 
       image/x-wmf                      image/wmf 
       image/tiff  
       application/postscript  
       image/gif  
  """
  global theObjs, theDoc
  if cdApp is None:
    StartChemDraw()
  if theObjs is None:
    if theDoc is None:
      theDoc = cdApp.Documents.Add()
    theObjs = theDoc.Objects
  theObjs.SetData(inFormat, inData, pythoncom.Missing)
  outD = theObjs.GetData(outFormat)
  theObjs.Clear()
  return outD


def CDXClean(inData, inFormat, outFormat):
  """calls the CDXLib Clean function on the data passed in.

    CDXLib_Clean attempts to clean (prettify) the data before
    doing an output conversion.  It can be thought of as CDXConvert++.

    CDXClean supports the same input and output specifiers as CDXConvert
    (see above)

  """
  global cdApp, theDoc, theObjs, selectItem, cleanItem
  if cdApp is None:
    StartChemDraw()
  if theObjs is None:
    if theDoc is None:
      theDoc = cdApp.Documents.Add()
    theObjs = theDoc.Objects
  theObjs.SetData(inFormat, inData, pythoncom.Missing)
  theObjs.Select()
  cleanItem.Execute()
  outD = theObjs.GetData(outFormat)
  theObjs.Clear()
  return outD


def CDXDisplay(inData, inFormat='chemical/cdx', clear=1):
  """ displays the data in Chemdraw """
  global cdApp, theDoc, theObjs, selectItem, cleanItem, centerItem
  if cdApp is None:
    StartChemDraw()
  try:
    theDoc.Activate()
  except Exception:
    ReactivateChemDraw()
    theObjs = theDoc.Objects
  if clear:
    theObjs.Clear()
  theObjs.SetData(inFormat, inData, pythoncom.Missing)
  return


def CDXGrab(outFormat='chemical/x-mdl-molfile'):
  """ returns the contents of the active chemdraw document

  """
  global cdApp, theDoc
  if cdApp is None:
    res = ""
  else:
    cdApp.Visible = 1
    if not cdApp.ActiveDocument:
      ReactivateChemDraw()
    try:
      res = cdApp.ActiveDocument.Objects.GetData(outFormat)
    except Exception:
      res = ""
  return res


def CloseChemdraw():
  """ shuts down chemdraw

  """
  global cdApp
  try:
    cdApp.Quit()
  except Exception:
    pass
  Exit()


def Exit():
  """ destroys our link to Chemdraw

  """
  global cdApp
  cdApp = None


def SaveChemDrawDoc(fileName='save.cdx'):
  """force chemdraw to save the active document

  NOTE: the extension of the filename will determine the format
   used to save the file.
  """
  d = cdApp.ActiveDocument
  d.SaveAs(fileName)


def CloseChemDrawDoc():
  """force chemdraw to save the active document

  NOTE: the extension of the filename will determine the format
   used to save the file.
  """
  d = cdApp.ActiveDocument
  d.Close()


def RaiseWindowNamed(nameRe):
  # start by getting a list of all the windows:
  cb = lambda x, y: y.append(x)
  wins = []
  win32gui.EnumWindows(cb, wins)

  # now check to see if any match our regexp:
  tgtWin = -1
  for win in wins:
    txt = win32gui.GetWindowText(win)
    if nameRe.match(txt):
      tgtWin = win
      break

  if tgtWin >= 0:
    win32gui.ShowWindow(tgtWin, 1)
    win32gui.BringWindowToTop(tgtWin)


def RaiseChemDraw():
  e = re.compile('^ChemDraw')
  RaiseWindowNamed(e)


try:
  from io import StringIO

  from PIL import Image

  def SmilesToPilImage(smilesStr):
    """takes a SMILES string and returns a PIL image using chemdraw

    """
    return MolToPilImage(smilesStr, inFormat='chemical/daylight-smiles', outFormat='image/gif')

  def MolToPilImage(dataStr, inFormat='chemical/daylight-smiles', outFormat='image/gif'):
    """takes a molecule string and returns a PIL image using chemdraw

    """
    # do the conversion...
    res = CDXConvert(dataStr, inFormat, outFormat)
    dataFile = StringIO(str(res))
    img = Image.open(dataFile).convert('RGB')
    return img
except ImportError:

  def SmilesToPilImage(smilesStr):
    print('You need to have PIL installed to use this functionality')
    return None

  def MolToPilImage(dataStr, inFormat='chemical/daylight-smiles', outFormat='image/gif'):
    print('You need to have PIL installed to use this functionality')
    return None


# ------------------------------------------------------------------
#  interactions with Chem3D
# ------------------------------------------------------------------
c3dApp = None


def StartChem3D(visible=0):
  """ launches Chem3D """
  global c3dApp
  c3dApp = Dispatch('Chem3D.Application')
  if not c3dApp.Visible:
    c3dApp.Visible = visible


def CloseChem3D():
  """ shuts down Chem3D """
  global c3dApp
  c3dApp.Quit()
  c3dApp = None


availChem3DProps = ('DipoleMoment', 'BendEnergy', 'Non14VDWEnergy', 'StericEnergy',
                    'StretchBendEnergy', 'StretchEnergy', 'TorsionEnergy', 'VDW14Energy')


def Add3DCoordsToMol(data, format, props={}):
  """ adds 3D coordinates to the data passed in using Chem3D

    **Arguments**

      - data: the molecular data

      - format: the format of _data_.  Should be something accepted by
        _CDXConvert_

      - props: (optional) a dictionary used to return calculated properties
  
  """
  global c3dApp
  if c3dApp is None:
    StartChem3D()
  if format != 'chemical/mdl-molfile':
    molData = CDXClean(data, format, 'chemical/mdl-molfile')
  else:
    molData = data
  with tempfile.NamedTemporaryFile(suffix='.mol', delete=False) as molF:
    molF.write(molData)
  doc = c3dApp.Documents.Open(molF.name)

  if not doc:
    print('cannot open molecule')
    raise ValueError('No Molecule')

  # set up the MM2 job
  job = Dispatch('Chem3D.MM2Job')
  job.Type = 1
  job.DisplayEveryIteration = 0
  job.RecordEveryIteration = 0

  # start the calculation...
  doc.MM2Compute(job)
  # and wait for it to finish
  while doc.ComputeStatus in [0x434f4d50, 0x50454e44]:
    pass
  #outFName = tempfile.mktemp('.mol')
  # this is horrible, but apparently Chem3D gets pissy with tempfiles:
  outFName = os.getcwd() + '/to3d.mol'
  doc.SaveAs(outFName)

  # generate the properties
  for prop in availChem3DProps:
    props[prop] = eval('doc.%s' % prop)

  doc.Close(0)

  os.unlink(molF.name)
  c3dData = open(outFName, 'r').read()
  gone = 0
  while not gone:
    try:
      os.unlink(outFName)
    except Exception:
      time.sleep(.5)
    else:
      gone = 1
  return c3dData


def OptimizeSDFile(inFileName, outFileName, problemFileName='problems.sdf', restartEvery=20):
  """  optimizes the structure of every molecule in the input SD file

    **Arguments**

      - inFileName: name of the input SD file

      - outFileName: name of the output SD file

      - problemFileName: (optional) name of the SD file used to store molecules which
        fail during the optimization process

      - restartEvery: (optional)  Chem3D will be shut down and restarted
        every _restartEvery_ molecules to try and keep core leaks under control

  """
  inFile = open(inFileName, 'r')
  outFile = open(outFileName, 'w+')
  problemFile = None
  props = {}
  lines = []
  nextLine = inFile.readline()
  skip = 0
  nDone = 0
  t1 = time.time()
  while nextLine != '':
    if nextLine.find('M  END') != -1:
      lines.append(nextLine)
      molBlock = ''.join(lines)

      try:
        newMolBlock = Add3DCoordsToMol(molBlock, 'chemical/mdl-molfile', props=props)
      except Exception:
        badBlock = molBlock
        skip = 1
        lines = []
      else:
        skip = 0
        lines = [newMolBlock]
    elif nextLine.find('$$$$') != -1:
      t2 = time.time()
      nDone += 1
      print('finished molecule %d in %f seconds' % (nDone, time.time() - t1))
      t1 = time.time()
      if nDone % restartEvery == 0:
        CloseChem3D()
        StartChem3D()
        outFile.close()
        outFile = open(outFileName, 'a')
      if not skip:
        for prop in props.keys():
          lines.append('> <%s>\n%f\n\n' % (prop, props[prop]))
        lines.append(nextLine)
        outFile.write(''.join(lines))
        lines = []
      else:
        skip = 0
        lines.append(nextLine)
        if problemFile is None:
          problemFile = open(problemFileName, 'w+')
        problemFile.write(badBlock)
        problemFile.write(''.join(lines))
        lines = []
    else:
      lines.append(nextLine)
    nextLine = inFile.readline()
  outFile.close()
  if problemFile is not None:
    problemFile.close()


if __name__ == '__main__':
  inStr = 'CCC(C=O)CCC'
  img = SmilesToPilImage(inStr)
  img.save('foo.jpg')
  convStr = CDXClean(inStr, 'chemical/x-daylight-smiles', 'chemical/x-daylight-smiles')
  print('in:', inStr)
  print('out:', convStr)
  convStr = CDXConvert(inStr, 'chemical/x-daylight-smiles', 'chemical/x-mdl-molfile')
  print('in:', inStr)
  print('out:', convStr)
  convStr2 = CDXClean(convStr, 'chemical/x-mdl-molfile', 'chemical/x-mdl-molfile')
  print('out2:', convStr2)

  inStr = 'COc1ccc(c2onc(c2C(=O)NCCc3ccc(F)cc3)c4ccc(F)cc4)c(OC)c1'
  convStr = CDXConvert(inStr, 'chemical/x-daylight-smiles', 'chemical/x-mdl-molfile')
  out = open('test.mol', 'w+')
  out.write(convStr)
  out.close()
