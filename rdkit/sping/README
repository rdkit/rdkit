Basic ideas for a package format for piddle:

1. Name change
--------------
Provisionally, I'm using the name "sping" for Simple
Platform Independent Graphics .  If someone has a better idea, I'm happy to
use it.

2. Package Structure
--------------------
A.  Canvas template

piddle.py -> sping.pid   # "pid" will stand for Plug In Drawing
          -> colors.py   # colors have been separated although sping.pid does a
                           "from colors import *"

B.  Backend renaming
piddleXY.py ->  sping/XY/pidXY.py  # backends are all capitals "XY"
                OR
                sping/XY.py

# It might be a bad idea to have two different ways to specify backends.  Comments?


For users, you can always get the XYCanvas by doing:

    > from sping.XY import XYCanvas # should also work

In addition 

    > from sping.XY import * # should expose a XYCanvas


C.  Backend helper files (e.g. pdfdoc.py)
    Are placed under the backend subdirectory, for example:

    piddle/pdfdoc.py   ->  sping/PDF/pdfdoc.py


D.  Other directories:

sping/*.py  -> core pid.py and any other core classes/utilties exposed to user

sping/lib  -> contains code useful to all backends.  Examples: bezier geometry code,
         coordinate transformation code 

sping/util ->  tools and other useful classes

sping/exapmles -> instructive examples for using various backends
                  (perhaps should be in separate doc distribution?)

sping/tests ->  place all verification tests here 


sping/contrib -> a stub directory where people can add their own code and contributions from
                 others that may not yet belong in the core yet

---------------
PIL usage
--------------

PIL is now used as a package.  I have added a patch to PIL to enable to it
to support all.  I am making a win32 (known to work on Windows NT 4)
available from the piddle web site.
