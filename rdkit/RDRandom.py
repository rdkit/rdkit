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
random = _random.random
randrange = _random.randrange
seed = _random.seed
