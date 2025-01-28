# $Id$
#
#  Copyright (C) 2002-2008  greg Landrum and Rational Discovery LLC
#
""" Various bits and pieces for calculating Molecular descriptors

"""

import pickle
import re

from rdkit.Chem import Descriptors as DescriptorsMod
from rdkit.ML.Descriptors import Descriptors
from rdkit.RDLogger import logger

logger = logger()


class MolecularDescriptorCalculator(Descriptors.DescriptorCalculator):
  """ used for calculating descriptors for molecules

  """

  def __init__(self, simpleList, *args, **kwargs):
    """ Constructor

      **Arguments**

        - simpleList: list of simple descriptors to be calculated
              (see below for format)

      **Note**

        - format of simpleList:

           a list of strings which are functions in the rdkit.Chem.Descriptors module

    """
    self.simpleList = tuple(simpleList)
    self.descriptorNames = tuple(self.simpleList)
    self.compoundList = None
    self._findVersions()

  def _findVersions(self):
    """ returns a tuple of the versions of the descriptor calculators

    """
    self.descriptorVersions = []
    for nm in self.simpleList:
      vers = 'N/A'
      if hasattr(DescriptorsMod, nm):
        fn = getattr(DescriptorsMod, nm)
        if hasattr(fn, 'version'):
          vers = fn.version
      self.descriptorVersions.append(vers)

  def SaveState(self, fileName):
    """ Writes this calculator off to a file so that it can be easily loaded later

     **Arguments**

       - fileName: the name of the file to be written

    """
    try:
      f = open(fileName, 'wb+')
    except Exception:
      logger.error('cannot open output file %s for writing' % (fileName))
      return
    pickle.dump(self, f)
    f.close()

  def CalcDescriptors(self, mol, *args, **kwargs):
    """ calculates all descriptors for a given molecule

      **Arguments**

        - mol: the molecule to be used

      **Returns**
        a tuple of all descriptor values

    """
    res = [-666] * len(self.simpleList)
    for i, nm in enumerate(self.simpleList):
      fn = getattr(DescriptorsMod, nm, lambda x: 777)
      try:
        res[i] = fn(mol)
      except Exception:
        import traceback
        traceback.print_exc()
    return tuple(res)

  def GetDescriptorNames(self):
    """ returns a tuple of the names of the descriptors this calculator generates

    """
    return self.descriptorNames

  def GetDescriptorSummaries(self):
    """ returns a tuple of summaries for the descriptors this calculator generates

    """
    res = []
    for nm in self.simpleList:
      fn = getattr(DescriptorsMod, nm, lambda x: 777)
      if hasattr(fn, '__doc__') and fn.__doc__:
        doc = fn.__doc__.split('\n\n')[0].strip()
        doc = re.sub(' *\n *', ' ', doc)
      else:
        doc = 'N/A'
      res.append(doc)
    return res

  def GetDescriptorFuncs(self):
    """ returns a tuple of the functions used to generate this calculator's descriptors

    """
    res = []
    for nm in self.simpleList:
      fn = getattr(DescriptorsMod, nm, lambda x: 777)
      res.append(fn)
    return tuple(res)

  def GetDescriptorVersions(self):
    """ returns a tuple of the versions of the descriptor calculators

    """
    return tuple(self.descriptorVersions)
