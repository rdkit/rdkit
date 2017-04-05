# $Id$
#
#  Copyright (C) 2002,2003  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" functionality for drawing trees on sping canvases

"""
import math

from rdkit.sping import pid as piddle


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
  lineWidth = 2
  horizOffset = 10
  vertOffset = 50
  labelFont = piddle.Font(face='helvetica', size=10)
  highlightColor = piddle.Color(1., 1., .4)
  highlightWidth = 2


visOpts = VisOpts()


def CalcTreeNodeSizes(node):
  """Recursively calculate the total number of nodes under us.

    results are set in node.totNChildren for this node and
    everything underneath it.
  """
  children = node.GetChildren()
  if len(children) > 0:
    nHere = 0
    nBelow = 0
    for child in children:
      CalcTreeNodeSizes(child)
      nHere = nHere + child.totNChildren
      if child.nLevelsBelow > nBelow:
        nBelow = child.nLevelsBelow
  else:
    nBelow = 0
    nHere = 1

  node.nExamples = len(node.GetExamples())
  node.totNChildren = nHere
  node.nLevelsBelow = nBelow + 1


def _ExampleCounter(node, min, max):
  if node.GetTerminal():
    cnt = node.nExamples
    if cnt < min:
      min = cnt
    if cnt > max:
      max = cnt
  else:
    for child in node.GetChildren():
      provMin, provMax = _ExampleCounter(child, min, max)
      if provMin < min:
        min = provMin
      if provMax > max:
        max = provMax
  return min, max


def _ApplyNodeScales(node, min, max):
  if node.GetTerminal():
    if max != min:
      loc = float(node.nExamples - min) / (max - min)
    else:
      loc = .5
    node._scaleLoc = loc
  else:
    for child in node.GetChildren():
      _ApplyNodeScales(child, min, max)


def SetNodeScales(node):
  min, max = 1e8, -1e8
  min, max = _ExampleCounter(node, min, max)
  node._scales = min, max
  _ApplyNodeScales(node, min, max)


def DrawTreeNode(node, loc, canvas, nRes=2, scaleLeaves=False, showPurity=False):
  """Recursively displays the given tree node and all its children on the canvas
  """
  try:
    nChildren = node.totNChildren
  except AttributeError:
    nChildren = None
  if nChildren is None:
    CalcTreeNodeSizes(node)

  if not scaleLeaves or not node.GetTerminal():
    rad = visOpts.circRad
  else:
    scaleLoc = getattr(node, "_scaleLoc", 0.5)

    rad = visOpts.minCircRad + node._scaleLoc * (visOpts.maxCircRad - visOpts.minCircRad)

  x1 = loc[0] - rad
  y1 = loc[1] - rad
  x2 = loc[0] + rad
  y2 = loc[1] + rad

  if showPurity and node.GetTerminal():
    examples = node.GetExamples()
    nEx = len(examples)
    if nEx:
      tgtVal = int(node.GetLabel())
      purity = 0.0
      for ex in examples:
        if int(ex[-1]) == tgtVal:
          purity += 1. / len(examples)
    else:
      purity = 1.0

    deg = purity * math.pi
    xFact = rad * math.sin(deg)
    yFact = rad * math.cos(deg)
    pureX = loc[0] + xFact
    pureY = loc[1] + yFact

  children = node.GetChildren()
  # just move down one level
  childY = loc[1] + visOpts.vertOffset
  # this is the left-hand side of the leftmost span
  childX = loc[0] - ((visOpts.horizOffset + visOpts.circRad) * node.totNChildren) / 2
  for i in range(len(children)):
    # center on this child's space
    child = children[i]
    halfWidth = ((visOpts.horizOffset + visOpts.circRad) * child.totNChildren) / 2

    childX = childX + halfWidth
    childLoc = [childX, childY]
    canvas.drawLine(loc[0], loc[1], childLoc[0], childLoc[1], visOpts.lineColor, visOpts.lineWidth)
    DrawTreeNode(child, childLoc, canvas, nRes=nRes, scaleLeaves=scaleLeaves, showPurity=showPurity)

    # and move over to the leftmost point of the next child
    childX = childX + halfWidth

  if node.GetTerminal():
    lab = node.GetLabel()
    cFac = float(lab) / float(nRes - 1)
    if hasattr(node, 'GetExamples') and node.GetExamples():
      theColor = (1. - cFac) * visOpts.terminalOffColor + cFac * visOpts.terminalOnColor
      outlColor = visOpts.outlineColor
    else:
      theColor = (1. - cFac) * visOpts.terminalOffColor + cFac * visOpts.terminalOnColor
      outlColor = visOpts.terminalEmptyColor
    canvas.drawEllipse(x1, y1, x2, y2, outlColor, visOpts.lineWidth, theColor)
    if showPurity:
      canvas.drawLine(loc[0], loc[1], pureX, pureY, piddle.Color(1, 1, 1), 2)
  else:
    theColor = visOpts.circColor
    canvas.drawEllipse(x1, y1, x2, y2, visOpts.outlineColor, visOpts.lineWidth, theColor)

    # this does not need to be done every time
    canvas.defaultFont = visOpts.labelFont

    labelStr = str(node.GetLabel())
    strLoc = (loc[0] - canvas.stringWidth(labelStr) / 2, loc[1] + canvas.fontHeight() / 4)

    canvas.drawString(labelStr, strLoc[0], strLoc[1])
  node._bBox = (x1, y1, x2, y2)


def CalcTreeWidth(tree):
  try:
    tree.totNChildren
  except AttributeError:
    CalcTreeNodeSizes(tree)
  totWidth = tree.totNChildren * (visOpts.circRad + visOpts.horizOffset)
  return totWidth


def DrawTree(tree, canvas, nRes=2, scaleLeaves=False, allowShrink=True, showPurity=False):
  dims = canvas.size
  loc = (dims[0] / 2, visOpts.vertOffset)
  if scaleLeaves:
    # try:
    #   l = tree._scales
    # except AttributeError:
    #   l = None
    # if l is None:
    SetNodeScales(tree)
  if allowShrink:
    treeWid = CalcTreeWidth(tree)
    while treeWid > dims[0]:
      visOpts.circRad /= 2
      visOpts.horizOffset /= 2
      treeWid = CalcTreeWidth(tree)
  DrawTreeNode(tree, loc, canvas, nRes, scaleLeaves=scaleLeaves, showPurity=showPurity)


def ResetTree(tree):
  tree._scales = None
  tree.totNChildren = None
  for child in tree.GetChildren():
    ResetTree(child)


def _simpleTest(canv):
  from .Tree import TreeNode as Node
  root = Node(None, 'r', label='r')
  c1 = root.AddChild('l1_1', label='l1_1')
  c2 = root.AddChild('l1_2', isTerminal=1, label=1)
  c3 = c1.AddChild('l2_1', isTerminal=1, label=0)
  c4 = c1.AddChild('l2_2', isTerminal=1, label=1)

  DrawTreeNode(root, (150, visOpts.vertOffset), canv)


if __name__ == '__main__':
  from rdkit.sping.PIL.pidPIL import PILCanvas
  canv = PILCanvas(size=(300, 300), name='test.png')
  _simpleTest(canv)
  canv.save()
