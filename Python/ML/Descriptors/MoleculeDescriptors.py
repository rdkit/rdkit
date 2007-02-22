#
#  Copyright (C) 2002  greg Landrum and Rational Discovery LLC
#
""" Various bits and pieces for calculating Molecular descriptors

"""
import RDConfig
from ML.Descriptors import Descriptors
from Chem import AvailDescriptors
AvailDescriptors.Desensitize()
import re

class MolecularDescriptorCalculator(Descriptors.DescriptorCalculator):
  """ used for calculating descriptors for molecules

  """
  def __init__(self,simpleList,*args,**kwargs):
    """ Constructor

      **Arguments**

        - simpleList: list of simple descriptors to be calculated
              (see below for format) 

      **Note**

        - format of simpleList:

           a list of strings which are keys into _AvailDescriptors.descDict_

    """
    self.simpleList = simpleList[:]
    self.descriptorNames = self.simpleList[:]
    self.compoundList = None
  
  def SaveState(self,fileName):
    """ Writes this calculator off to a file so that it can be easily loaded later

     **Arguments**

       - fileName: the name of the file to be written
       
    """
    import cPickle
    try:
      f = open(fileName,'wb+')
    except:
      print 'cannot open output file %s for writing'%(fileName)
      return
    cPickle.dump(self,f)
    f.close()

  def CalcDescriptors(self,mol,*args,**kwargs):
    """ calculates all descriptors for a given molecule

      **Arguments**

        - mol: the molecule to be used

      **Returns**  
        a tuple of all descriptor values

    """
    res = [None]*len(self.simpleList)
    for i in range(len(self.simpleList)):
      nm = self.simpleList[i]
      fn = AvailDescriptors.descDict.get(nm,lambda x:777)
      #print '>',nm
      try:
        res[i] = fn(mol)
      except:
        import traceback
        traceback.print_exc()
        res[i]=-666
    return tuple(res)

  def GetDescriptorNames(self):
    """ returns a tuple of the names of the descriptors this calculator generates

    """
    return tuple(self.descriptorNames)

  def GetDescriptorSummaries(self):
    """ returns a tuple of summaries for the descriptors this calculator generates

    """
    res = []
    for nm in self.simpleList:
      fn = AvailDescriptors.descDict.get(nm,lambda x:777)
      if hasattr(fn,'__doc__') and fn.__doc__:
        doc = fn.__doc__.split('\n\n')[0].strip()
        doc = re.sub('\ *\n\ *',' ',doc)
      else:
        doc = 'N/A'
      res.append(doc)
    return res

  def GetDescriptorFuncs(self):
    """ returns a tuple of the functions used to generate this calculator's descriptors

    """
    res = []
    for nm in self.simpleList:
      fn = AvailDescriptors.descDict.get(nm,lambda x:777)
      res.append(fn)
    return tuple(res)  
    
if __name__ == '__main__':
  from Chem import *
  from ML.Descriptors import MoleculeDescriptors
  descs = ['MolLogP','IWd','Chi1v']
  smis = ['CCOC','CC=O','CCC(=O)O']
  calc = MoleculeDescriptors.MolecularDescriptorCalculator(descs)
  calc.SaveState('test_data/molcalc.dsc')
  print calc.GetDescriptorNames()

  for smi in smis:
    mol = MolFromSmi(smi)
    print smi,calc.CalcDescriptors(mol)
  
