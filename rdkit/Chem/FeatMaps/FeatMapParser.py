# $Id$
#
# Copyright (C) 2006 Greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import re

from rdkit import Geometry
from rdkit.Chem.FeatMaps import FeatMapPoint, FeatMaps

"""

ScoreMode=All
DirScoreMode=Ignore

BeginParams
  family=Aromatic radius=2.5 width=1.0 profile=Gaussian
  family=Acceptor radius=1.5
EndParams

# optional
BeginPoints
  family=Acceptor pos=(1.0, 0.0, 5.0) weight=1.25 dir=(1, 1, 0)
  family=Aromatic pos=(0.0,1.0,0.0) weight=2.0 dir=(0,0,1) dir=(0,0,-1)
  family=Acceptor pos=(1.0,1.0,2.0) weight=1.25
EndPoints

"""


class FeatMapParseError(ValueError):
  pass


class FeatMapParser(object):
  data = None

  def __init__(self, file=None, data=None):
    if file:
      self.data = file.readlines()
    elif data:
      self.SetData(data)
    self._lineNum = 0

  def SetData(self, data):
    if isinstance(data, str):
      self.data = data.split('\n')
    else:
      self.data = data
    self._lineNum = 0

  def _NextLine(self):
    txt = ''
    while 1:
      try:
        l = self.data[self._lineNum].split('#')[0].strip()
      except IndexError:
        break
      self._lineNum += 1
      if l:
        txt += l
        if l[-1] != '\\':
          break
    return txt

  def Parse(self, featMap=None):
    if featMap is None:
      featMap = FeatMaps.FeatMap()

    l = self._NextLine().strip()
    while l:
      splitL = l.split('=')
      if len(splitL) == 1:
        keyword = splitL[0].strip().lower()
        if keyword == 'beginpoints':
          pts = self.ParseFeatPointBlock()
          for pt in pts:
            featMap.AddFeatPoint(pt)
        elif keyword == 'beginparams':
          featMap.params = self.ParseParamBlock()
        else:
          raise FeatMapParseError('Unrecognized keyword %s on line %d' % (keyword, self._lineNum))
      else:
        keyword = splitL[0].strip().lower()
        val = splitL[1].strip()
        if keyword == 'scoremode':
          try:
            featMap.scoreMode = getattr(FeatMaps.FeatMapScoreMode, val)
          except AttributeError:
            raise FeatMapParseError('ScoreMode %s not recognized on line %d' % (val, self._lineNum))
        elif keyword == 'dirscoremode':
          try:
            featMap.dirScoreMode = getattr(FeatMaps.FeatDirScoreMode, val)
          except AttributeError:
            raise FeatMapParseError('DirScoreMode %s not recognized on line %d' %
                                    (val, self._lineNum))
        else:
          raise FeatMapParseError('Unrecognized keyword %s on line %d' % (keyword, self._lineNum))
      l = self._NextLine().strip()
    return featMap

  def ParseParamBlock(self):
    paramLineSplitter = re.compile(r'([a-zA-Z]+) *= *(\S+)')
    params = {}

    l = self._NextLine()
    while l and l != 'EndParams':
      param = FeatMaps.FeatMapParams()
      vals = paramLineSplitter.findall(l)
      for name, val in vals:
        name = name.lower()
        if name == 'family':
          family = val
        elif name == 'radius':
          param.radius = float(val)
        elif name == 'width':
          param.width = float(val)
        elif name == 'profile':
          try:
            param.featProfile = getattr(param.FeatProfile, val)
          except AttributeError:
            raise FeatMapParseError('Profile %s not recognized on line %d' % (val, self._lineNum))
        else:
          raise FeatMapParseError('FeatMapParam option %s not recognized on line %d' %
                                  (name, self._lineNum))
      params[family] = param
      l = self._NextLine()

    if l != 'EndParams':
      raise FeatMapParseError('EndParams line not found')

    return params

  def _parsePoint(self, txt):
    txt = txt.strip()
    startP = 0
    endP = len(txt)
    if txt[0] == '(':
      startP += 1
    if txt[-1] == ')':
      endP -= 1
    txt = txt[startP:endP]
    splitL = txt.split(',')
    if len(splitL) != 3:
      raise ValueError('Bad location string')
    return Geometry.Point3D(float(splitL[0]), float(splitL[1]), float(splitL[2]))

  def ParseFeatPointBlock(self):
    featLineSplitter = re.compile(r'([a-zA-Z]+) *= *')
    feats = []

    l = self._NextLine()
    while l and l != 'EndPoints':
      vals = featLineSplitter.split(l)
      while vals.count(''):
        vals.remove('')
      p = FeatMapPoint.FeatMapPoint()
      for i in range(0, len(vals), 2):
        name = vals[i].lower()
        value = vals[i + 1]
        if name == 'family':
          p.SetFamily(value.strip())
        elif name == 'weight':
          p.weight = float(value)
        elif name == 'pos':
          p.SetPos(self._parsePoint(value))
        elif name == 'dir':
          p.featDirs.append(self._parsePoint(value))
        else:
          raise FeatMapParseError(f'FeatPoint option {name} not recognized on line {self._lineNum}')

      feats.append(p)
      l = self._NextLine()
    return feats
