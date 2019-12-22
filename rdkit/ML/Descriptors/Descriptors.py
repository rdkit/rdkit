#
#  Copyright (C) 2001,2002  greg Landrum and Rational Discovery LLC
#
""" Various bits and pieces for calculating descriptors

"""

import pickle


class DescriptorCalculator:
  """ abstract base class for descriptor calculators

  """

  def __init__(self, *args, **kwargs):
    """ Constructor

    """
    self.simpleList = None
    self.descriptorNames = None
    self.compoundList = None

  # ------------
  #  methods used to calculate descriptors
  # ------------

  def ShowDescriptors(self):
    """ prints out a list of the descriptors

    """
    if self.simpleList is None:
      raise NotImplementedError('Need to have a simpleList defined')
    print('#---------')
    print('Simple:')
    for desc in self.simpleList:
      print(desc)
    if self.compoundList:
      print('#---------')
      print('Compound:')
      for desc in self.compoundList:
        print(desc)

  def GetDescriptorNames(self):
    """ returns a list of the names of the descriptors this calculator generates

    """
    raise NotImplementedError('abstract base class')

  def SaveState(self, fileName):
    """ Writes this calculator off to a file so that it can be easily loaded later

     **Arguments**

       - fileName: the name of the file to be written

    """
    try:
      f = open(fileName, 'wb+')
    except Exception:
      print('cannot open output file %s for writing' % (fileName))
      return
    pickle.dump(self, f)
    f.close()

  def CalcDescriptors(self, what, *args, **kwargs):
    raise NotImplementedError('abstract base class')
