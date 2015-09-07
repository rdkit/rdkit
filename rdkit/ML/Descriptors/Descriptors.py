#
#  Copyright (C) 2001,2002  greg Landrum and Rational Discovery LLC
#
""" Various bits and pieces for calculating descriptors

"""
from __future__ import print_function
from rdkit import RDConfig

class DescriptorCalculator:
  """ abstract base class for descriptor calculators

  """
  
  #------------
  #  methods used to calculate descriptors
  #------------

  def ShowDescriptors(self):
    """ prints out a list of the descriptors

    """
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
    pass

  def SaveState(self,fileName):
    """ Writes this calculator off to a file so that it can be easily loaded later

     **Arguments**

       - fileName: the name of the file to be written
       
    """
    from rdkit.six.moves import cPickle
    try:
      f = open(fileName,'wb+')
    except Exception:
      print('cannot open output file %s for writing'%(fileName))
      return
    cPickle.dump(self,f)
    f.close()

  def CalcDescriptors(self,what,*args,**kwargs):
    pass
    
  def __init__(self,*args,**kwargs):
    """ Constructor

    """
    self.simpleList = None
    self.descriptorNames = None
    self.compoundList = None
    
