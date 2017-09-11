#
# Copyright (C) 2017 Peter Gedeck
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
'''
Collection of utilities to be used with descriptors

'''


def setDescriptorVersion(version='1.0.0'):
  """ Set the version on the descriptor function.

  Use as a decorator """
  def wrapper(func):
    func.version = version
    return func
  return wrapper
