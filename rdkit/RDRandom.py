#
#  Copyright (C) 2003-2023  Greg Landrum and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" making random numbers consistent so we get good regressions

"""

import random as _random
import sys

random = _random.random
randrange = _random.randrange
seed = _random.seed

if sys.hexversion >= 0x30b0000:
  # Python 3.11 deprecated the random argument to random.shuffle, which changed
  # the behavior and made it impossible to exactly reproduce old results. this is
  # a slightly adapted form of random.shuffle from python 3.10:
  #  https://github.com/python/cpython/blob/b05352e4c2f25b292fb7de0ab927e74415bc2dd8/Lib/random.py#LL380-L404C40
  def shuffle(x, random=None):
    """Shuffle list x in place, and return None.
        Optional argument random is a 0-argument function returning a
        random float in [0.0, 1.0); if it is the default None, the
        standard random.random will be used.
        """

    if random is None:
      randbelow = _random._randbelow
      for i in reversed(range(1, len(x))):
        # pick an element in x[:i+1] with which to exchange x[i]
        j = randbelow(i + 1)
        x[i], x[j] = x[j], x[i]
    else:
      floor = _random._floor
      for i in reversed(range(1, len(x))):
        # pick an element in x[:i+1] with which to exchange x[i]
        j = floor(random() * (i + 1))
        x[i], x[j] = x[j], x[i]
else:
  shuffle = _random.shuffle
