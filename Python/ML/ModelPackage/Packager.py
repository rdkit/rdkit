## Automatically adapted for numpy.oldnumeric Sep 23, 2006 by alter_code1.py

#
# Copyright (C) 2002 Greg Landrum and Rational Discovery LLC
# All rights are reserved.
#


import exceptions
class DescriptorCalculationError(exceptions.Exception):
  """ used to signal problems generating descriptor values """
  pass
class ClassificationError(exceptions.Exception):
  """ used to signal problems generating predictions """
  pass

class ModelPackage(object):
  """ a container class to package a composite model with a descriptor
  calculator so that objects needing predictions (compounds, molecules, etc.)
  can be passed directly in without worrying about generating descriptors

  """
  def __init__(self,descCalc=None,model=None,dataSet=None,notes=''):
    self._descCalc = descCalc
    self._model = model
    self._notes = notes
    self._dataSet = dataSet
    self._initialized = 0
    self._supplementalData = []
    
  def SetCalculator(self,calc):
    self._descCalc = calc
  def GetCalculator(self):
    return self._descCalc

  def SetModel(self,model):
    self._model = model
  def GetModel(self):
    return self._model

  def SetDataset(self,data):
    self._dataSet = data
  def GetDataset(self):
    return self._dataSet

  def SetNotes(self,notes):
    self._notes = notes
  def GetNotes(self):
    return self._notes

  def SetSupplementalData(self,suppD):
    self._supplementalData = suppD
  def GetSupplementalData(self):
    if not hasattr(self,'_supplementalData'):
      self._supplementalData = []
    return self._supplementalData
  def AddSupplementalData(self,data):
    if not hasattr(self,'_supplementalData'):
      self._supplementalData = []
    self._supplementalData.append(data)
  
  def Classify(self,obj,label='',threshold=0):
    if not self._initialized:
      self.Init()
    try:
      descs = self._descCalc.CalcDescriptors(obj)
    except:
      raise DescriptorCalculationError,'problems encountered generating descriptors'

    argVect = [label]+list(descs)+[0]
    try:
      res = self._model.ClassifyExample(argVect,threshold=threshold,appendExample=0)
    except:
      raise ClassificationError,'problems encountered generating prediction'

    return res

  def Init(self):
    if self._model is None or self._descCalc is None:
      return
    
    nms = self._model.GetDescriptorNames()
    lbl = nms[0]
    act = nms[-1]
    descs = self._descCalc.GetDescriptorNames()
    order = [lbl] + list(descs) + [act]
    self._model.SetInputOrder(order)
    
    self._initialized = 1

if __name__=='__main__':
  from Chem import *
  import cPickle
  from ML.ModelPackage import Packager
  
  calc = cPickle.load(open('test_data/Jan9_build3_calc.dsc','rb'))
  model = cPickle.load(open('test_data/Jan9_build3_model.pkl','rb'))
  pkg = Packager.ModelPackage(descCalc=calc,model=model)
  pkg.SetNotes('General purpose model built from PhysProp data')
  testD = [
    ('Fc1ccc(NC(=O)c2cccnc2Oc3cccc(c3)C(F)(F)F)c(F)c1',0,1.0 ),
    (r'CN/1(=C\C=C(/C=C1)\C\2=C\C=N(C)(Cl)\C=C2)Cl',0,0.70),
    (r'NS(=O)(=O)c1cc(ccc1Cl)C2(O)NC(=O)c3ccccc32',1,0.70),
    ]
    
  for smi,pred,conf in testD:
    m = MolFromSmi(smi)
    p,c = pkg.Classify(m)
    if pred!=p or conf!=c:
      raise ValueError,'Bad Prediction: %s'%(repr((smi,pred,conf,p,c)))
  cPickle.dump(pkg,open('test_data/Jan9_build3_pkg.pkl','wb+'))
  from Numeric import *
  import numpy.oldnumeric.random_array as RandomArray
  
  names = calc.GetDescriptorNames()
  perm = [names[x] for x in RandomArray.permutation(len(names))]
  calc.simpleList = perm
  calc.descriptorNames = perm
  pkg.Init()
  for smi,pred,conf in testD:
    m = MolFromSmi(smi)
    p,c = pkg.Classify(m)
    if pred!=p or conf!=c:
      raise ValueError,'Bad Prediction: %s'%(repr((smi,pred,conf,p,c)))
  
              
