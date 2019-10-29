#
# Copyright (C) 2002-2008 Greg Landrum and Rational Discovery LLC
# All rights are reserved.
#


class DescriptorCalculationError(Exception):
  """ used to signal problems generating descriptor values """
  pass


class ClassificationError(Exception):
  """ used to signal problems generating predictions """
  pass


class ModelPackage(object):
  """ a container class to package a composite model with a descriptor
  calculator so that objects needing predictions (compounds, molecules, etc.)
  can be passed directly in without worrying about generating descriptors

  """

  def __init__(self, descCalc=None, model=None, dataSet=None, notes=''):
    self._descCalc = descCalc
    self._model = model
    self._notes = notes
    self._dataSet = dataSet
    self._initialized = 0
    self._supplementalData = []

  def SetCalculator(self, calc):
    self._descCalc = calc

  def GetCalculator(self):
    return self._descCalc

  def SetModel(self, model):
    self._model = model

  def GetModel(self):
    return self._model

  def SetDataset(self, data):
    self._dataSet = data

  def GetDataset(self):
    return self._dataSet

  def SetNotes(self, notes):
    self._notes = notes

  def GetNotes(self):
    return self._notes

  def SetSupplementalData(self, suppD):
    self._supplementalData = suppD

  def GetSupplementalData(self):
    if not hasattr(self, '_supplementalData'):
      self.SetSupplementalData([])
    return self._supplementalData

  def AddSupplementalData(self, data):
    if not hasattr(self, '_supplementalData'):
      self.SetSupplementalData([])
    self._supplementalData.append(data)

  def Classify(self, obj, label='', threshold=0):
    if not self._initialized:
      self.Init()
    try:
      descs = self._descCalc.CalcDescriptors(obj)
    except Exception:
      raise DescriptorCalculationError('problems encountered generating descriptors')

    argVect = [label] + list(descs) + [0]
    try:
      res = self._model.ClassifyExample(argVect, threshold=threshold, appendExample=0)
    except Exception:
      import traceback
      traceback.print_exc()
      raise ClassificationError('problems encountered generating prediction')

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
