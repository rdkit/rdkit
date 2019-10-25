# $Id$
#
#  Copyright (C)2003-2006  Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" functionality for drawing hierarchical catalogs on sping
    canvases

"""
from sping import pid as piddle


class VisOpts(object):
  circRad = 10
  minCircRad = 4
  maxCircRad = 16
  circColor = piddle.Color(0.6, 0.6, 0.9)
  terminalEmptyColor = piddle.Color(.8, .8, .2)
  terminalOnColor = piddle.Color(0.8, 0.8, 0.8)
  terminalOffColor = piddle.Color(0.2, 0.2, 0.2)
  outlineColor = piddle.transparent
  lineColor = piddle.Color(0, 0, 0)
  lineWidth = 1
  horizOffset = 5
  vertOffset = 75
  topMargin = 20
  labelFont = piddle.Font(face='helvetica', size=10)
  highlightColor = piddle.Color(1., 1., .4)
  highlightWidth = 2


visOpts = VisOpts()


def GetMinCanvasSize(adjList, levelList):
  maxAcross = -1
  for k in levelList.keys():
    nHere = len(levelList[k])
    maxAcross = max(maxAcross, nHere)
  nLevs = len(levelList.keys())
  minSize = (maxAcross * (visOpts.minCircRad * 2 + visOpts.horizOffset),
             visOpts.topMargin + nLevs * visOpts.vertOffset)
  return minSize


def DrawHierarchy(adjList, levelList, canvas, entryColors=None, bitIds=None, minLevel=-1,
                  maxLevel=1e8):
  """

  Arguments:

   - adjList: adjacency list representation of the hierarchy to be drawn

   - levelList: dictionary mapping level -> list of ids

  """
  if bitIds is None:
    bitIds = []
  if entryColors is None:
    entryColors = {}

  levelLengths = levelList.keys()
  levelLengths.sort()
  minLevel = max(minLevel, levelLengths[0])
  maxLevel = min(maxLevel, levelLengths[-1])

  dims = canvas.size
  drawLocs = {}
  # start at the bottom of the hierarchy and work up:
  for levelLen in range(maxLevel, minLevel - 1, -1):
    nLevelsDown = levelLen - minLevel
    pos = [0, visOpts.vertOffset * nLevelsDown + visOpts.topMargin]

    ids = levelList.get(levelLen, [])

    # FIX: we'll eventually want to figure out some kind of sorting here:
    nHere = len(ids)
    canvas.defaultFont = visOpts.labelFont
    if nHere:
      # figure the size of each node at this level:
      spacePerNode = float(dims[0]) / nHere
      spacePerNode -= visOpts.horizOffset
      nodeRad = max(spacePerNode / 2, visOpts.minCircRad)
      nodeRad = min(nodeRad, visOpts.maxCircRad)
      spacePerNode = nodeRad * 2 + visOpts.horizOffset
      # start in the middle of the canvas:
      pos[0] = dims[0] / 2.
      # maybe we need to offset a little:
      if nHere % 2:
        pos[0] -= spacePerNode / 2

      # move to the left by half the number of nodes:
      pos[0] -= (nHere // 2 - .5) * spacePerNode

      # Find the locations and draw connectors:
      for ID in ids:
        if not bitIds or ID in bitIds:
          # first do lines down to the next level:
          if levelLen != maxLevel:
            for neighbor in adjList[ID]:
              if neighbor in drawLocs:
                p2 = drawLocs[neighbor][0]
                canvas.drawLine(pos[0], pos[1], p2[0], p2[1], visOpts.lineColor, visOpts.lineWidth)
          drawLocs[ID] = tuple(pos), nodeRad
          pos[0] += spacePerNode

  for ID in drawLocs.keys():
    pos, nodeRad = drawLocs[ID]
    x1, y1 = pos[0] - nodeRad, pos[1] - nodeRad
    x2, y2 = pos[0] + nodeRad, pos[1] + nodeRad
    drawColor = entryColors.get(ID, visOpts.circColor)
    canvas.drawEllipse(x1, y1, x2, y2, visOpts.outlineColor, 0, drawColor)
    label = str(ID)
    # txtLoc = ( pos[0]-canvas.stringWidth(label)/2,
    #           pos[1]+canvas.fontHeight()/4 )
    txtLoc = (pos[0] + canvas.fontHeight() / 4, pos[1] + canvas.stringWidth(label) / 2)
    canvas.drawString(label, txtLoc[0], txtLoc[1], angle=90)

  return drawLocs
