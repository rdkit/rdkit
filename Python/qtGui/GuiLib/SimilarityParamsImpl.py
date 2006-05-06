# $Id: SimilarityParamsImpl.py 4684 2005-05-25 21:50:46Z glandrum $
#
#  Copyright (C) 2003-2005  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" 
"""    
from qt import *
from forms.SimilarityParams import SimilarityParamsWidget as _Form
import DataStructs

class SimilarityParamsWidget(_Form):
  def __init__(self,*args,**kwargs):
    """

    """
    _Form.__init__(self,*args,**kwargs)

  def getMetric(self):
    res = None
    if self.tanimotoRadio.isChecked():
      res = DataStructs.TanimotoSimilarity
    elif self.diceRadio.isChecked():
      res = DataStructs.DiceSimilarity
    elif self.cosineRadio.isChecked():
      res = DataStructs.CosineSimilarity
    elif self.sokalRadio.isChecked():
      res = DataStructs.SokalSimilarity
    return res


  def updateFingerprinterDetails(self,details):
    """ fills the contents of the details structure passed in

    **Arguments**:

      - details: a _FingerprinterDetails_ instance

    """
    details.minPath = self.fragmentMinSize.value()
    details.maxPath = self.fragmentMaxSize.value()
    details.fpSize =  int(str(self.fragmentNumBits.text()))
    details.metric = self.getMetric()

    if self.topNRadio.isChecked():
      details.doThreshold = 0
      details.topN = int(str(self.topN.text()))
    elif self.greaterThanRadio.isChecked():
      details.doThreshold = 1
      details.screenThresh = float(str(self.greaterThan.text()))
      
