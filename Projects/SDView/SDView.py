# $Id$
#
#  Copyright (c) 2008, Novartis Institutes for BioMedical Research Inc.
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# version 2 as published by the Free Software Foundation
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details
#
# Created by Greg Landrum Sept 2006
#
_version = "0.7.0"

_usage="""
   SDView [optional arguments] sdfilename

   Note: using "-" as the sd filename will cause the program to read the SD data from stdin

"""
_welcomeMessage="This is SDView version %s"%(_version)
import sys,os,copy
from RDLogger import logger
logger = logger()
import Chem
from Chem import AllChem
from Chem.ChemUtils.AlignDepict import AlignDepict
from qt import *
from qtGui import qtUtils
from qtGui.GuiLib.MolTable import MolTable
from qtGui.GuiLib.MolCanvas import MolCanvasView
from qtGui.GuiBase import GuiBase

class SDViewTable(MolTable):
  def __init__(self,*args,**kwargs):
    MolTable.__init__(self,*args)
    tgt= MolCanvasView()
    tgt.initCanvas((300,300))
    tgt.setIcon(QPixmap(qtUtils.logoImageData))
    tgt.setCaption('Molecule')
    if kwargs.get('includeAtomNumbers'):
      tgt.includeAtomNumbers=True
    tgt.show()
    self.molCanvas=tgt
    self._alwaysDraw=True

class SDViewGui(GuiBase):
  def __init__(self,*args,**kwargs):
    GuiBase.__init__(self,*args,**kwargs)
    widg = SDViewTable(self,**kwargs)
    self.setCentralWidget(widg)
    widg.show()
    self.setCaption('SDView')
    self.resize(400,400)
    self._dir='.'
    self.sdvTable = widg
    
  def _initMenus(self):
    GuiBase._initMenus(self)
    self._fileMenu.insertItem(self.trUtf8('&Save As'),self.saveAs,Qt.CTRL+Qt.Key_S)
    self._fileMenu.insertItem(self.trUtf8('&Quit'),qApp,SLOT("quit()"),Qt.CTRL+Qt.Key_Q )
    self._replaceToggle = QCheckBox(self.trUtf8('&Replace Current'),self._viewMenu)
    self._viewMenu.insertItem(self._replaceToggle)
    self._replaceToggle.setChecked(1)
    self._zoomToggle = QCheckBox(self.trUtf8('&Zoom on Display'),self._viewMenu)
    self._viewMenu.insertItem(self._zoomToggle)
    self._zoomToggle.setChecked(1)

  def save(self,fileN):
    try:
      w= Chem.SDWriter(fileN)
    except:
      qtUtils.error('could not open file %s for writing'%fileN)
      return
    self.sdvTable.exportToMolWriter(w)
    
  def saveAs(self):
    fileN = str(QFileDialog.getSaveFileName(self._dir,'SD Files(*.sdf);;All files (*.*)'))
    if fileN:
      self._dir = os.sep.join(os.path.split(fileN)[0:-1])
      self.save(fileN)
    

def RunGui(suppl,redraw=False,do3D=False,showOnly=None,numberAtoms=False,
           alignCore=None,alignPatt=None,keepHs=False,matchConf=False):
  from qtGui import Gui
  app,win = Gui.Launcher(SDViewGui,None,includeLogo=False)
  widg = win.centralWidget()
  widg.molCanvas.includeAtomNumbers=numberAtoms
      
  processFunc=None
  if do3D:
    from Chem.PyMol import MolViewer
    try:
      pymol = MolViewer()
    except:
      print >>sys.stderr,'ERROR: Unable to connect to PyMol server.\nPlease run ~landrgr1/extern/PyMol/launch.sh to start it.'
      return
    def cb(tbl,row,col,viewer=pymol,owner=win):
      itm = tbl.item(row,tbl._molCol)
      if not itm: return
      mol = itm.mol()
      if hasattr(mol,'mol3d'):
        mol3d = mol.mol3d
      else:
        mol3d = mol
      confId=-1
      nm = mol.GetProp('_Name')
      if not nm:
        nm="Row_%d"%(row+1)
      else:
        nm=nm.replace(':','').replace(' ','_').replace('-','_')
      if owner._replaceToggle.isChecked():
        viewer.HideAll()
      viewer.ShowMol(mol3d,name=nm,confId=confId,showOnly=False,
                     zoom=owner._zoomToggle.isChecked())
    widg.addChangedCallback(cb)

    if redraw:
      def pf(mol,keepHs=keepHs,matchConf=matchConf,alignCore=alignCore):
        if keepHs:
          tmp = Chem.RemoveHs(mol)
        else:
          tmp = mol
          mol = copy.deepcopy(tmp)
          if tmp.HasProp('_Name'):
            mol.SetProp('_Name',tmp.GetProp('_Name'))
        try:
          if matchConf:
            AllChem.GenerateDepictionMatching3DStructure(tmp,tmp,canonOrient=False,clearConfs=False)
          elif alignCore:
            AlignDepict(tmp,alignCore,corePattern=alignPatt,acceptFailure=True)
          else:
            AllChem.Compute2DCoords(tmp,canonOrient=True,clearConfs=True)
        except:
          print Chem.MolToMolBlock(tmp)
          print '----'
          print Chem.MolToMolBlock(mol)
        tmp.mol3d = mol
        return tmp
      processFunc=pf
  else:
    if redraw or alignCore:
      if not alignCore:
        processFunc=lambda x:AllChem.Compute2DCoords(x,canonOrient=True,clearConfs=True)
      else:
        processFunc=lambda x:AlignDepict(x,alignCore,corePattern=alignPatt,acceptFailure=True)
  if showOnly:
    suppl = [x for x in suppl if (x.HasProp('_Name') and x.GetProp('_Name') in showOnly)]
    
  # and load up those molecules:
  widg.loadFromVLib(suppl,includeCheckboxes=0,processFunc=processFunc)
  
  app.exec_loop()
  win.destroy(1)


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  --- --- ---  --- --- ---  --- --- --- 
# set up the option parser
from optparse import OptionParser
import os,RDConfig
parser=OptionParser(_usage,version='%prog '+_version)
parser.add_option('--do3D','--3D','--3d','--do3d',default=False,action='store_true',
                  help='display the molecules in 3D using PyMol')
parser.add_option('--redraw',default=False,action='store_true',
                  help='generate new 2D coordinates for the molecules')
parser.add_option('--numberAtoms','-n',default=False,action='store_true',
                  help='show atom numbers in the 2D sketches')
parser.add_option('--showOnly',default=None,
                  help='enter a list of molecules to be displayed')
parser.add_option('--keepHs',dest="removeHs",default=True,action='store_false',
                  help="don't remove Hs from the molecules before displaying them")
parser.add_option('--alignPatt','--alignPattern',default='',
                  help="provide a smarts pattern to recognize the core")
parser.add_option('--alignCore',default='',
                  help="provide a molecule for the core. If this corresponds to a mol file it will be read in and used as is. If it's a SMILES, a default depiction will be generated.")
parser.add_option('--matchConf',default=False,action='store_true',
                  help="attempt to generate a 2D depiction similar to the 3D orientation (only makes sense if --redraw is also provided).")
parser.add_option('--smiles','--smi',default=False,action='store_true',
                  help="read smiles instead of an sd files")
parser.add_option('--delimiter','--delim',default='',
                  help="set the delimiter for the smiles input file (defaults to automatic)")
parser.add_option('--forceSDF','--SDF-','--sdf','--sd',
                  default=False,action='store_true',
                  help="force reading an SD file, regardless of the extension")

if __name__=='__main__':
  import csv
  options,args = parser.parse_args()
  if len(args)<1:
    parser.error('no SD file provided')

  if options.showOnly:
    val= options.showOnly
    val = val.replace('\n',',').replace(' ',',').replace('\t',',')
    options.showOnly = val.split(',')

  if options.alignCore and options.matchConf:
    logger.warning('The --alignCore and --matchConf options are not compatible. Only --matchConf will be used.')
    options.alignCore=''
  if options.alignCore and not options.redraw:
    logger.warning('The --alignCore option only makes sense in combination with --redraw (to do 3D alignments use SDAlign.sh). --alignCore option ignored.')
    options.alignCore=''
  if options.alignPatt and not options.alignCore:
    logger.warning('The --alignPattern option only makes sense in combination with --alignCore. --alignPattern ignored.')
    options.alignPatt=''
  if options.smiles and options.forceSDF:
    logger.error('the --smiles and --forceSDF options are mutually exclusive')

  if args[0]!='-' and not options.smiles and not options.forceSDF:
    nm,ext = os.path.splitext(args[0])
    ext = ext.lower()
    if ext not in ('.sdf','.mol'):
      options.smiles=True

  if options.smiles:
    options.redraw=True

  if not options.smiles:
    from Chem.FastSDMolSupplier import FastSDMolSupplier
    if args[0] != "-":
      suppl = FastSDMolSupplier(args[0],removeHs=options.removeHs)
    else:
      suppl = Chem.SDMolSupplier()
      d = sys.stdin.read()
      suppl.SetData(d,removeHs=options.removeHs)
      d = None
  else:
    if not options.delimiter:
      sniffer = csv.Sniffer()
    else:
      sniffer=None
    if args[0] != "-":
      if sniffer:
        d = file(args[0],'r').read(1000)
        dialect = sniffer.sniff(d)
        options.delimiter = dialect.delimiter
      suppl = Chem.SmilesMolSupplier(args[0],delimiter=options.delimiter)
    else:
      suppl = Chem.SmilesMolSupplier()
      d = sys.stdin.read()
      if sniffer:
        dialect = sniffer.sniff(d)
        options.delimiter = dialect.delimiter
      suppl.SetData(d,delimiter=options.delimiter)
      d = None
  if options.alignCore:
    import os
    core = None
    if os.path.exists(options.alignCore):
      try:
        core = Chem.MolFromMolFile(options.alignCore)
      except:
        core = None
    if not core:
      core = Chem.MolFromSmiles(options.alignCore)
      AllChem.Compute2DCoords(core,canonOrient=True)
    if not core:
      logger.error('could not construct a molecule from the core provided with --alignCore')
      sys.exit(1)
      
    options.alignCore = core
    if not options.alignPatt:
      options.alignPatt=options.alignCore
    else:
      patt = Chem.MolFromSmarts(options.alignPatt)
      if patt.GetNumAtoms()!=core.GetNumAtoms():
        logger.error('The core molecule (from --alignCore) and core pattern (from --alignPatt) do not have the same number of atoms.')
        sys.exit(1)
                     
      options.alignPatt=patt

  RunGui(suppl,options.redraw,options.do3D,
         showOnly=options.showOnly,
         numberAtoms=options.numberAtoms,
         alignCore=options.alignCore,
         alignPatt=options.alignPatt,
         keepHs=not options.removeHs,
         matchConf=options.matchConf)


  
      
  
