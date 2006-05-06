#
#  Copyright (C) 2001 greg Landrum
#
""" the COM server for Composite models

"""

import RDConfig, string
from ML.Composite import Composite
import cPickle
import sys
import winerror
from win32com.server import exception

class CompositeServer(object):
  """ exposes an interface for composite COM servers

    This interface does not support modifying the composite, only
     classifying new examples

    **Public Methods**

     - LoadComposite

     - ClassifyExample

     - ShowComposite

     - GetVoteDetails

     - SetInputOrder

     - GetDescriptorNames

     - Close

    **Public Attributes**

      None
     
    **ProgID**
    
      RD.Composite

  """
  _public_methods_ = ['LoadComposite','ClassifyExample','ShowComposite',
                      'GetVoteDetails','SetInputOrder','GetDescriptorNames',
                      'Close']
  _public_attrs_ = []
  _reg_clsid_ = '{9F62358E-9043-4BF9-93C7-47BED8BE522F}'
  _reg_progid_ = "RD.Composite"

  def Close(self):
    """ Blows out the local composite

      **Note**

        _LoadComposite()_ must be called after this for further use of the
        composite

    """
    self._composite = None

  def LoadComposite(self,fileName):
    """ Loads a (pickled) composite from a file

     **Arguments**
       - fileName: the name of the file to load

    """
    try:
      f = open(str(fileName),'rb')
    except:
      raise exception.COMException('The file %s could not be opened'%(fileName),
                                   winerror.ERROR_FILE_NOT_FOUND)
    try:
      self._composite = cPickle.load(f)
    except:
      raise exception.COMException('The composite could not be loaded',
                                   winerror.ERROR_INVALID_HANDLE)
      
    return fileName

  def ShowComposite(self):
    """ returns a string representation of the composite

    """
    try:
      c = self._composite
    except AttributeError:
      raise exception.COMException('The composite has not yet been loaded',
                                   winerror.ERROR_INVALID_HANDLE)
    else:
      return str(c)
    
  def GetVoteDetails(self):
    """ returns a list of the results of the last vote

    """
    try:
      c = self._composite
    except AttributeError:
      raise exception.COMException('The composite has not yet been loaded',
                                   winerror.ERROR_FILE_NOT_FOUND)
    else:
      try:
        return c.GetVoteDetails()
      except:
        raise exception.COMException('VoteDetails unavaliable',winerror.ERROR_INVALID_HANDLE)
  
  def SetInputOrder(self,colNames):
    """ Sets the input order for the composite

      this is used so that the composite can remap inputs which are not
      in an order it expects

      **Arguments**

        - colNames: a list of the names of the columns in the data which will be
           passed in to the composite

      **Notes**

        - there must be a column name which matches each of the composite's
          descriptors (accessible via _GetDescriptorNames()_)

    """
    try:
      c = self._composite
    except AttributeError:
      raise exception.COMException('The composite has not yet been loaded',
                                   winerror.ERROR_INVALID_HANDLE)

    try:
      colNames = map(str,colNames)
      c.SetInputOrder(colNames)
    except:
      raise exception.COMException('SetInputOrder failed.',winerror.ERROR_INVALID_HANDLE)
    return 0 
      
  def GetDescriptorNames(self):
    """ returns a list of the descriptor names the composite expects

    """
    try:
      c = self._composite
    except AttributeError:
      raise exception.COMException('The composite has not yet been loaded',
                                   winerror.ERROR_INVALID_HANDLE)
    else:
      try:
        return c.GetDescriptorNames()
      except:
        raise exception.COMException('GetDescriptorNames failed.',winerror.ERROR_INVALID_HANDLE)

  def ClassifyExample(self,example,threshold=0):
    """ classifies a new example

      **Arguments**

        - example: a list containing the example to be classified

        - threshold: the threshold to be used for high-confidence predictions

      **Returns**

        a list two elements long:

           1) the classification (this will be -1 if the confidence of the prediction
              was below _threshold_

           2) the confidence
      
    """       
    try:
      c = self._composite
    except AttributeError:
      raise exception.COMException('The composite has not yet been loaded',
                                   winerror.ERROR_INVALID_HANDLE)
    if 0:
      try:
        res = c.ClassifyExample(example,threshold=threshold)
      except:
        import traceback,sys
        outF = open(RDConfig.RDCodeDir+'/ml/composite/log.txt','a+')
        outF.write('#------------------------------\n')
        outF.write('ex: %s\n'%str(example))
        traceback.print_tb(sys.exc_info()[2],file=outF)
        outF.close()
        res = []
    else:
      try:
        res = c.ClassifyExample(example,threshold=threshold)
      except:
        raise exception.COMException('ClassifyExample failed.',winerror.ERROR_INVALID_HANDLE)
    return list(res)
  
if __name__=='__main__':
    from win32com.server import register
    register.UseCommandLine(CompositeServer)

