#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
import VPItems
from qt import *
from qtcanvas import *
import math
from SupplierWindowImpl import SupplierWindow
from SinkWindowImpl import SinkWindow
from ReactionWindowImpl import ReactionWindow
from FilterWindowImpl import FilterWindow
import VPIcons
import VPUtils
from PIL import Image
from utils import PilTools
from qtGui import qtUtils

class VPCircleItem(QCanvasEllipse,VPItems.VPItemMixin):
  """

    attachment pt order:
      [left,right]

  """
  def __init__(self,canvas,width,height,start=0,end=360*16,id=None,
               graphParent=None):
    QCanvasEllipse.__init__(self,width,height,start,end,canvas)
    VPItems.VPItemMixin.__init__(self,id,graphParent=graphParent)
    self._parentClass = QCanvasEllipse


    p1 = VPItems.VPLinkPoint(canvas,(-width/2,0),self,0)
    self._vpLinkPts.append(p1)
    p1 = VPItems.VPLinkPoint(canvas,(width/2,0),self,1)
    self._vpLinkPts.append(p1)

    
class VPRectangleItem(QCanvasRectangle,VPItems.VPItemMixin):
  """

    attachment pt order:

      [left top,left bottom,right center]

  """
  def __init__(self,canvas,width,height,id=None,
               graphParent=None):
    QCanvasRectangle.__init__(self,0,0,width,height,canvas)
    VPItems.VPItemMixin.__init__(self,id,graphParent=graphParent)
    self._parentClass = QCanvasRectangle

    p1 = VPItems.VPLinkPoint(canvas,(0,height/3),self,0)
    self._vpLinkPts.append(p1)
    p1 = VPItems.VPLinkPoint(canvas,(0,2*height/3),self,1)
    self._vpLinkPts.append(p1)
    p1 = VPItems.VPLinkPoint(canvas,(width,height/2),self,2)
    self._vpLinkPts.append(p1)


class VPTriangleItem(QCanvasPolygon,VPItems.VPItemMixin):
  """ point to the right



  """
  heightScale = math.sqrt(3.)/2.
  def __init__(self,canvas,sideLen,id=None,
               graphParent=None):
    QCanvasPolygon.__init__(self,canvas)
    VPItems.VPItemMixin.__init__(self,id,graphParent=graphParent)
    self._parentClass = QCanvasPolygon
    p0 = [0,0]
    height = self.heightScale*sideLen
    p1 = [0,sideLen]
    p2 = [height,sideLen/2]
    ptArr = QPointArray(p0 + p1 + p2)
    self.setPoints(ptArr)
    
    p1 = VPItems.VPLinkPoint(canvas,p2,self,0)
    self._vpLinkPts.append(p1)

class VPSupplyItem(QCanvasRectangle,VPItems.VPItemMixin):
  """ point to the right

      attachment pt order:

  """
  def __init__(self,canvas,side=100,label="Supply",fillColor=(255,255,255),id=None,
               graphParent=None):
    QCanvasRectangle.__init__(self,0,0,side,side,canvas)
    self._canvas = canvas
    VPItems.VPItemMixin.__init__(self,id,graphParent=graphParent)
    self._parentClass = QCanvasRectangle
    self._width = side
    self._height = side
    self._label=label
    self._addLabel()
    self._addIcon()
    self._addLinkPts()
    self._inputIndices = []
    self._outputIndices = [0]
    self.moveBy(-self._width/2,-self._height/2)
    self.setBrush(QBrush(apply(QColor,fillColor)))
    


  def _addLinkPts(self):
    pt = [self._width,self._height/2]
    p1 = VPItems.VPLinkPoint(self._canvas,pt,self,0)
    self._vpLinkPts.append(p1)

  def showProperties(self):
    if self._propDlg is None:
      self._propDlg = SupplierWindow(self)
      
    self._propDlg.show()

  def _addIcon(self):
    self._pm = QCanvasPixmap(QPixmap(VPIcons.supplyIcon_xpmData),
                             QPoint(0,0))
    self._pma = QCanvasPixmapArray()
    self._pma.setImage(0,self._pm)
    self._icon = QCanvasSprite(self._pma,self._canvas)
    self._icon.setZ(self.z()+1)
    self._icon.show()
    self._icon.moveBy(self._width/2-self._pm.width()/2,5)


class VPReactionItem(QCanvasRectangle,VPItems.VPItemMixin):
  """

    attachment pt order:

      [ reactant1,...,reactantN,product1,...,productN ]

  """
  def __init__(self,canvas,width=100,height=100,rxnFile=None,rxnImg=None,
               nReactants=2,nProducts=1,
               #id=None,label="Rxn",fillColor=(250,180,250),
               id=None,label="Rxn",fillColor=(255,255,255),
               graphParent=None):
    QCanvasRectangle.__init__(self,0,0,width,height,canvas)
    VPItems.VPItemMixin.__init__(self,id,graphParent=graphParent)
    self._parentClass = QCanvasRectangle
    self._canvas = canvas

    self._width=width
    self._height=height
    #p0 = [0,0]
    #p1 = [self._width,0]
    #p2 = [self._width,self._height]
    #p3 = [0,self._height]
    #ptArr = QPointArray(p0+p1+p2+p3)
    #self.setPoints(ptArr)
    self._inputIndices = range(nReactants)
    self._outputIndices = range(nReactants,nReactants+nProducts)
    self._label = label
    self._addLabel()
    self._addIcon()
    if rxnFile is not None:
      nReactants,nProducts = VPUtils.quickCDXMLParse(rxnFile)
    self._rxnFile = rxnFile
    if rxnImg is not None:
      try:
        self._rxnImg = Image.open(rxnImg)
      except:
        qtUtils.warning('cannot open image file %s'%(rxnImg))
        self._rxnImg = None
    else:
      self._rxnImg = rxnImg
      
    self._nReactants = nReactants
    self._nProducts = nProducts
    self._addLinkPts()

    self.moveBy(-self._width/2,-self._height/2)
    self.setBrush(QBrush(apply(QColor,fillColor)))

  def _addIcon(self):
    self._pm = QCanvasPixmap(QPixmap(VPIcons.reactionIcon_xpmData),
                             QPoint(0,0))
    self._pma = QCanvasPixmapArray()
    self._pma.setImage(0,self._pm)
    self._icon = QCanvasSprite(self._pma,self._canvas)
    self._icon.setZ(self.z()+1)
    self._icon.show()
    self._icon.moveBy(self._width/2-self._pm.width()/2,5)

  def _addLinkPts(self):
    # FIX: destroy existing links first
    self._vpLinkPts = []
    # first reactants
    vertDivider = self._nReactants + 1
    nSoFar = 0
    for i in range(self._nReactants):
      p1 = VPItems.VPLinkPoint(self._canvas,
                               (0,(i+1)*self._height/vertDivider),
                               self,nSoFar)
      self._vpLinkPts.append(p1)
      nSoFar += 1
    # now products
    vertDivider = self._nProducts + 1
    for i in range(self._nProducts):
      p1 = VPItems.VPLinkPoint(self._canvas,
                               (self._width,(i+1)*self._height/vertDivider),
                               self,nSoFar)
      self._vpLinkPts.append(p1)
      nSoFar += 1


  def setNumReactants(self,num):
    """ NOTE: this blows out any existing links """
    self._nReactants = num
    self._addLinkPts()

  def setNumProducts(self,num):
    """ NOTE: this blows out any existing links """
    self._nProducts = num
    self._addLinkPts()
    
  def showProperties(self):
    if self._propDlg is None:
      self._propDlg = ReactionWindow(self,rxnFile=self._rxnFile,rxnImg=self._rxnImg)
    self._propDlg.show()


class VPFilterItem(QCanvasRectangle,VPItems.VPItemMixin):
  """ triangular

    attachment pt order:

      [ input(s),output1,output0 ]

  """
  def __init__(self,canvas,side=100,nInputs=1,
               includeOutput1=1,includeOutput2=0,
               id=None,label="Filter",fillColor=(255,255,255),
               graphParent=None):
    QCanvasRectangle.__init__(self,0,0,side,side,canvas)
    VPItems.VPItemMixin.__init__(self,id,graphParent=graphParent)
    self._parentClass = QCanvasRectangle
    self._canvas = canvas
    
    self._height = side
    self._width = side

    self._label = label
    self._addLabel()
    self._addIcon()

    self._nInputs = nInputs
    self._includeOutput1 = includeOutput1
    self._includeOutput2 = includeOutput2
    self._nOutputs = self._includeOutput1+self._includeOutput2
    self._addLinkPts()
    self._inputIndices = range(self._nInputs)
    self._outputIndices = [self._nInputs]
    self.moveBy(-self._width/2,-self._height/2)
    self.setBrush(QBrush(apply(QColor,fillColor)))

  def _addIcon(self):
    self._pm = QCanvasPixmap(QPixmap(VPIcons.filterIcon_xpmData),
                             QPoint(0,0))
    self._pma = QCanvasPixmapArray()
    self._pma.setImage(0,self._pm)
    self._icon = QCanvasSprite(self._pma,self._canvas)
    self._icon.setZ(self.z()+1)
    self._icon.show()
    self._icon.moveBy(self._width/2-self._pm.width()/2,5)

  def _addLinkPts(self):
    self._vpLinkPts = []
    # first reactants
    vertDivider = self._nInputs + 1
    nSoFar = 0
    for i in range(self._nInputs):
      p1 = VPItems.VPLinkPoint(self._canvas,
                               (0,(i+1)*self._height/vertDivider),
                               self,nSoFar)
      self._vpLinkPts.append(p1)
      nSoFar += 1
    # output1
    vertDivider = self._nOutputs+1
    if self._includeOutput1:
      p1 = VPItems.VPLinkPoint(self._canvas,
                               (self._width,self._height/vertDivider),
                               self,nSoFar)
      self._vpLinkPts.append(p1)
      nSoFar += 1
    if self._includeOutput2:
      p1 = VPItems.VPLinkPoint(self._canvas,
                               (self._width,2*self._height/vertDivider),
                               self,nSoFar)
      self._vpLinkPts.append(p1)
      nSoFar += 1
    
  def showProperties(self):
    if self._propDlg is None:
      self._propDlg = FilterWindow(self)
    self._propDlg.show()

class VPSinkItem(QCanvasRectangle,VPItems.VPItemMixin):
  """ point to the left

      attachment pt order:

  """
  def __init__(self,canvas,side=100,label="Sink",fillColor=(255,255,255),id=None,
               graphParent=None):
    QCanvasRectangle.__init__(self,0,0,side,side,canvas)
    self._canvas = canvas
    VPItems.VPItemMixin.__init__(self,id,graphParent=graphParent)
    self._parentClass = QCanvasRectangle
    self._width = side
    self._height = side
    self._label=label
    self._addLabel()
    self._addIcon()
    self._addLinkPts()
    self._inputIndices = []
    self._outputIndices = [0]
    self.moveBy(-self._width/2,-self._height/2)
    self.setBrush(QBrush(apply(QColor,fillColor)))


  def _addLinkPts(self):
    pt = [0,self._height/2]
    p1 = VPItems.VPLinkPoint(self._canvas,pt,self,0)
    self._vpLinkPts.append(p1)

  def showProperties(self):
    if self._propDlg is None:
      self._propDlg = SinkWindow(self)
      
    self._propDlg.show()

  def _addIcon(self):
    self._pm = QCanvasPixmap(QPixmap(VPIcons.sinkIcon_xpmData),
                             QPoint(0,0))
    self._pma = QCanvasPixmapArray()
    self._pma.setImage(0,self._pm)
    self._icon = QCanvasSprite(self._pma,self._canvas)
    self._icon.setZ(self.z()+1)
    self._icon.show()
    self._icon.moveBy(self._width/2-self._pm.width()/2,5)
