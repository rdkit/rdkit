# $Id$
#
#  Copyright (C) 2007-2009 Greg Landrum
#   @@ All Rights Reserved @@
#
try:
  import sqlalchemy
except ImportError:
  from Loader_orig import *
else:
  from Loader_sa import *

