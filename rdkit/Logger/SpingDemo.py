# $Id$
#
#       Copyright (c) 2001-2006, Greg Landrum and Rational Discovery LLC,
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" demo code for the Logger class

"""

from sping.SVG import pidSVG
from sping.PIL import pidPIL
from sping import pid
import Logger

# create a logged canvas and draw on it
sz = (300,300)
c1 = Logger.Logger(pidSVG.SVGCanvas,sz,'foo.svg',loggerFlushCommand='clear')
c1.drawPolygon([(100,100),(100,200),(200,200),(200,100)],fillColor=pid.Color(0,0,1))
c1.drawLines([(100,100,200,200),(100,200,200,100)],color=pid.Color(0,1,0),width=2)

# because the log has been instantiated with clear() as the loggerFlushCommand, 
# this will blow out the log as well as the contents of the canvas.
c1.clear()

# draw some more stuff
c1.drawPolygon([(100,100),(100,200),(200,200),(200,100)],fillColor=pid.Color(1,0,0))
c1.drawLines([(100,100,200,200),(100,200,200,100)],color=pid.Color(0,0,0),width=2)
# and write the resulting file.
c1.save()

# save the log by pickling it.
from rdkit.six.moves import cPickle
cPickle.dump(c1._LoggerGetLog(),open('foo.pkl','wb+'))

# create a new canvas
c2 = pidPIL.PILCanvas(sz,'foo.png')
# read the pickled log back in 
t = cPickle.load(open('foo.pkl','rb'))
# and play the log on the new canvas
Logger.replay(t,c2)
# there should now be a file 'foo.png' with the image


