import math

from numpy import *

from rdkit.sping import pid


def DrawSpiral(canvas, startColor, endColor, startRadius, endRadius, nLoops, degsPerSlice=70,
               degsPerStep=1, startAngle=0, centerPos=None, dir=1):
  if centerPos is None:
    centerPos = (canvas.size[0] / 2, canvas.size[1] / 2)
  nSlices = int(math.ceil(360 * nLoops / degsPerSlice))
  radPerStep = math.pi * degsPerStep / 180.
  stepsPerSlice = degsPerSlice / degsPerStep
  radiusStep = float(endRadius - startRadius) / (stepsPerSlice * nSlices)
  colorStep = (array(endColor, float) - array(startColor, float)) / nSlices
  print('INFO:', nSlices, radPerStep, stepsPerSlice, radiusStep, colorStep)

  angle = math.pi * startAngle / 180.
  radius = startRadius
  color = array(startColor, float)

  for i in range(nSlices):
    pts = [(centerPos[0], centerPos[1])]
    for j in range(stepsPerSlice):
      xPos = centerPos[0] + radius * math.cos(angle)
      yPos = centerPos[1] + radius * math.sin(angle)
      pts.append((xPos, yPos))

      angle += dir * radPerStep
      radius += radiusStep
    xPos = centerPos[0] + radius * math.cos(angle)
    yPos = centerPos[1] + radius * math.sin(angle)
    pts.append((xPos, yPos))
    canvas.drawPolygon(pts, edgeColor=pid.transparent,
                       fillColor=pid.Color(color[0], color[1], color[2]), closed=1)
    angle -= dir * radPerStep
    color += colorStep


if __name__ == '__main__':
  #from sping.PIL.pidPIL import PILCanvas
  #canv = PILCanvas(size=(600,600),name='test.png')
  from rdkit.sping.SVG.pidSVG import SVGCanvas

  #from rdkit.sping.PDF.pidPDf import PDFCanvas
  canv = SVGCanvas(size=(600, 600), name='test.svg')
  #canv = PDFCanvas(size=(600,600),name='test.pdf')
  DrawSpiral(canv, (.2, .2, 1), (.9, .9, 1.), 200, 50, 8, startAngle=-135, degsPerSlice=11, dir=-1)
  canv.save()
