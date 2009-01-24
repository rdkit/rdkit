#
#  Copyright (C) 2001 greg Landrum
#
""" the COM server for Descriptor Calculators

"""

from pyRDKit import RDConfig
from pyRDKit.ML.Descriptors import Descriptors,Parser
import cPickle
import sys
import winerror
from win32com.server import exception


class DescriptorServer:
  """ exposes an interface for descriptor calculator COM servers

    This interface does not support modifying the calculator, only
     classifying new examples

    **Public Methods**

      - LoadCalculator

      - CalcDescriptors

      - GetDescriptorNames

      - Close
      
    **Public Attributes**

      None
     
    **ProgID**
    
      RD.DescCalc

  """
  _public_methods_ = ['LoadCalculator','CalcDescriptors',
                      'GetDescriptorNames','Close']
  _public_attrs_ = []
  _reg_clsid_ = "{2DEC34F0-7FBA-4752-8BF7-65E4D81B46E4}"
  _reg_progid_ = "RD.DescCalc"

  def Close(self):
    """ Blows out the local calculator

      **Note**

        _LoadCalculator()_ must be called after this for further use of the
        calculator.

    """
    self._descCalculator = None
    self.propnames = None
    
  def LoadCalculator(self,fileName):
    """ Loads a (pickled) calculator from a file

     **Arguments**
     
       - fileName: the name of the file to load

    """
    try:
      f = open(str(fileName),'rb')
    except:
      raise exception.COMException('The file %s could not be opened'%(fileName),
                                   winerror.ERROR_FILE_NOT_FOUND)
    try:
      self._descCalculator = cPickle.load(f)
    except:
      raise exception.COMException('The calculator could not be loaded',
                                   winerror.DISP_E_EXCEPTION)
    return fileName

  def GetDescriptorNames(self):
    """ returns a list of the names of the descriptors this calculator generates

    """
    try:
      c = self._descCalculator
    except:
      return ['']
    else:
      return self._descCalculator.GetDescriptorNames()
    

  def CalcDescriptors(self,argVect,colNames):
    """ Calculates the descriptors for a new composition

      **Arguments**

        - argVect: a list of values

        - colNames: the names of the columns in argVect

      **Returns**

       a list of descriptor values

      **Note**

        - colNames should include the names of all composition descriptors
          returned by any compound descriptors in the calculator

    """
    argVect = list(argVect)
    argVect[0] = str(argVect[0])
    pDict = {}
    for i in xrange(len(colNames)):
      pDict[str(colNames[i])] = argVect[i]
    try:
      return self._descCalculator.CalcDescriptors(argVect,pDict)
    except:
      outF = open(RDConfig.RDCodeDir+'/ml/descriptors/log.txt','a+')
      outF.write('#------------------------------\n')
      outF.write('av: %s\n'%str(argVect))
      outF.write('cn: %s\n'%str(colNames))
      outF.write('pd: %s\n'%str(pDict))
      outF.close()

  
if __name__=='__main__':
    from win32com.server import register
    register.UseCommandLine(DescriptorServer)

