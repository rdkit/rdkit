# $Id: PyMol.py 5114 2006-04-06 15:28:45Z glandrum $
#
# Copyright (C) 2004-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" uses pymol to interact with molecules

"""
import Chem
import xmlrpclib


_server=None
class MolViewer(object):
  def __init__(self,host='localhost',port=9123,force=0,**kwargs):
    global _server
    if not force and _server is not None:
      self.server=_server
    else:
      _server=None
      serv = xmlrpclib.Server('http://%s:%d'%(host,port))
      serv.ping()
      _server = serv
      self.server=serv
    self.InitializePyMol()
    
  def InitializePyMol(self):
    """ does some initializations to set up PyMol according to our
    tastes

    """
    self.server.do('set valence,1')
    self.server.do('set stick_rad,0.15')
    self.server.do('set mouse_selection_mode,0')
    self.server.do('set line_width,2')
    self.server.do('set selection_width,10')
    self.server.do('set auto_zoom,0')
    

  def DeleteAll(self):
    self.server.deleteAll()

  def DeleteAllExcept(self,excludes):
    allNames = self.server.getNames('*',False)
    for nm in allNames:
      if nm not in excludes:
        self.server.deleteObject(nm)

  def LoadFile(self,filename,name,showOnly=False):
    if showOnly:
      self.DeleteAll()
    id = self.server.loadFile(filename,name)
    return id

  def ShowMol(self,mol,name='molecule',showOnly=True,highlightFeatures=[],
              molB="",confId=-1,zoom=True):
    """ special case for displaying a molecule or mol block """
  
    if not molB:
      molB = Chem.MolToMolBlock(mol,confId=confId)
    server = self.server
    if not zoom:
      self.server.do('view rdinterface,store')
    if showOnly:
      self.DeleteAll()
    id = server.loadMolBlock(molB,name)
    if highlightFeatures:
      nm = name+'-features'
      conf = mol.GetConformer(confId)
      for feat in highlightFeatures:
        pt = [0.0,0.0,0.0]
        for idx in feat:
          loc = conf.GetAtomPosition(idx)
          pt[0] += loc[0]/len(feat)
          pt[1] += loc[1]/len(feat)
          pt[2] += loc[2]/len(feat)
        server.sphere(pt,0.2,(1,1,1),nm)
    if zoom:
      server.zoom('visible')
    else:
      self.server.do('view rdinterface,recall')
    return id
  def GetSelectedAtoms(self,whichSelection=None):
    if not whichSelection:
      sels = self.server.getNames('selections')
      if sels:
        whichSelection = sels[-1]
      else:
        whichSelection=None
    if whichSelection:   
      items = self.server.index(whichSelection)
    else:
      items = []
    return items


  def SelectAtoms(self,itemId,atomIndices,selName='selection'):
    ids = '(id '
    ids += ','.join(['%d'%(x+1) for x in atomIndices])
    ids += ')'
    cmd = 'select %s,%s and %s'%(selName,ids,itemId)
    self.server.do(cmd)
  
  def HighlightAtoms(self,indices,where,extraHighlight=False):
    if extraHighlight:
      idxText = ','.join(['%s and (id %d)'%(where,x) for x in indices])
      self.server.do('edit %s'%idxText)
    else:
      idxText = ' or '.join(['id %d'%x for x in indices])
      self.server.do('select selection, %s and (%s)'%(where,idxText))

  def SetDisplayStyle(self,obj,style=''):
    self.server.do('hide everything,%s'%(obj,))
    if style:
      self.server.do('show %s,%s'%(style,obj))

  def SelectProteinNeighborhood(self,aroundObj,inObj,distance=5.0,
                                name='neighborhood',showSurface=False):
    self.server.do('select %(name)s,byres (%(aroundObj)s around %(distance)f) and %(inObj)s'%locals())

    if showSurface:
      self.server.do('show surface,%s'%name)
      self.server.do('disable %s'%name)


  def AddPharmacophore(self,locs,colors,label,sphereRad=0.5):
    self.server.do('view rdinterface,store')
    self.server.resetCGO(label)
    for i,loc in enumerate(locs):
      self.server.sphere(loc,sphereRad,colors[i],label,1)
    self.server.do('enable %s'%label)
    self.server.do('view rdinterface,recall')
    

  def SetDisplayUpdate(self,val):
    if not val:
      self.server.do('set defer_update,1')
    else:
      self.server.do('set defer_update,0')

      
  def GetAtomCoords(self,sels):
    res = {}
    for label,idx in sels:
      coords = self.server.getAtomCoords('(%s and id %d)'%(label,idx))
      res[(label,idx)] = coords
      print 'grab:',label,idx,coords
    return res

  def HideAll(self):
    self.server.do('disable all')
  def HideObject(self,objName):
    self.server.do('disable %s'%objName)
  def DisplayObject(self,objName):
    self.server.do('enable %s'%objName)

  def Redraw(self):
    self.server.do('refresh')
  def Zoom(self,objName):
    self.server.zoom(objName)

  def DisplayHBonds(self,objName,molName,proteinName,
                    molSelText='(%(molName)s)',
                    proteinSelText='(%(proteinName)s and not het)'):
    cmd = "delete %(objName)s;\n"
    cmd += "dist %(objName)s," + molSelText+","+proteinSelText+",mode=2;\n"
    cmd += "enable %(objName)s;"
    cmd = cmd%locals()

    self.server.do(cmd)

  def DisplayCollisions(self,objName,molName,proteinName,distCutoff=3.0,
                        color='red',
                        molSelText='(%(molName)s)',
                        proteinSelText='(%(proteinName)s and not het)'):
    cmd = "delete %(objName)s;\n"
    cmd += "dist %(objName)s," + molSelText+","+proteinSelText+",%(distCutoff)f,mode=0;\n"
    cmd += """enable %(objName)s
    color %(color)s, %(objName)s"""
    cmd = cmd%locals()
    self.server.do(cmd)

