# $Id$
#
#  Copyright (C) 2007-2009 Greg Landrum
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
try:
  import sqlalchemy
except ImportError:
  from Loader_orig import *
else:
  from Loader_sa import *
