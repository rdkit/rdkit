# $Id$
#
#  Copyright (C) 2003-2006  Greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" making random numbers consistent so we get good regressions

"""

import sys

import random as _random
if sys.hexversion >= 0x20303f0:
  _randGen = _random.WichmannHill()
  random = _randGen.random
  randrange = _randGen.randrange
  def seed(val):
    global _randGen,random,randrange
    _randGen = _random.WichmannHill()
    _randGen.whseed(val)
    random = _randGen.random
    randrange = _randGen.randrange
else:
  random = _random.random
  randrange = _random.randrange
  seed = _random.whseed
