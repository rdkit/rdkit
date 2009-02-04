#!/usr/bin/env python
import pidtest
import sping.PS.pidPS

canvas = sping.PS.pidPS.PSCanvas(name="testallps.ps", size=(400,500))

pidtest.drawStrings(canvas)

canvas.flush()
canvas.clear() # resets the canvas


pidtest.drawRotstring(canvas)
canvas.flush()
canvas.clear()

x1 = 100
y1 = 100
radius = 20
canvas.drawString( "Test",100,250 )
canvas.drawEllipse(x1,y1-radius,x1+2*radius,y1+radius,
                   edgeWidth=1)
canvas.flush()
canvas.clear()

pidtest.drawBasics(canvas)
canvas.flush()
canvas.clear()

import pstests
pstests.testLatin1Chars(canvas)

canvas.save()
