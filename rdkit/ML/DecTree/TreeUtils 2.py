#
#  Copyright (C) 2001-2004  greg Landrum and Rational Discovery LLC
#  All Rights Reserved
#
""" Utilities for working with trees

"""


def CollectLabelLevels(tree, levels, level=0, maxDepth=1e8):
  if level < maxDepth:
    if not tree.GetTerminal():
      l = tree.GetLabel()
      currLevel = levels.get(l, 1e8)
      if level < currLevel:
        levels[l] = level
      for child in tree.GetChildren():
        CollectLabelLevels(child, levels, level + 1, maxDepth)
  return levels


def CollectDescriptorNames(tree, names, level=0, maxDepth=1e8):
  if level < maxDepth:
    if not tree.GetTerminal():
      names[tree.GetLabel()] = tree.GetName()
      for child in tree.GetChildren():
        CollectDescriptorNames(child, names, level + 1, maxDepth)
  return names
