#
# Copyright (C) 2003-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Supplies an abstract class for working with sequences of molecules

"""


class MolSupplier(object):
  """ we must, at minimum, support forward iteration

  """

  def __init__(self):
    raise ValueError('cannot instantiate MolSuppliers')

  def Reset(self):
    pass

  def __iter__(self):
    self.Reset()
    return self

  def next(self):
    res = self.NextMol()
    if res is not None:
      return res
    else:
      raise StopIteration

  def NextMol(self):
    """   Must be implemented in child class

    """
    pass

  __next__ = next  # PY3
