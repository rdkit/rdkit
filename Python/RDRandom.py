# $Id$
#
#  Copyright (C) 2003-2006  Greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" making random numbers consistent so we get good regressions

"""

import sys

import random as _random
if sys.hexversion >= 0x20303f0:
  _randGen = _random.WichmannHill()
  random = _randGen.random
  randrange = _randGen.randrange
  shuffle = _randGen.shuffle
  def seed(val):
    global _randGen,random,randrange,shuffle
    _randGen = _random.WichmannHill()
    _randGen.whseed(val)
    random = _randGen.random
    randrange = _randGen.randrange
    shuffle = _randGen.shuffle
else:
  random = _random.random
  randrange = _random.randrange
  seed = _random.whseed
  def shuffle(val):
    raise NotImplementedError,'shuffle not implemented for older python versions'
