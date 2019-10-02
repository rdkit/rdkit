#
#  Copyright (C) 2008 Greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import os
import re

from aggdraw import Brush, Pen
from aggdraw import Draw
from aggdraw import Font
from rdkit import RDConfig
from rdkit.Chem.Draw.canvasbase import CanvasBase

faceMap = {'sans': os.path.join(RDConfig.RDCodeDir, 'Chem', 'Draw', 'FreeSans.ttf')}


def convertColor(color):
  color = (int(color[0] * 255), int(color[1] * 255), int(color[2] * 255))
  return color


class Canvas(CanvasBase):
  # fonts appear smaller in aggdraw than with cairo
  # fix that here:
  fontScale = 1.2

  def __init__(self,
               img=None,
               imageType=None,  # determines file type
               fileName=None,  # if set determines output file name
               size=None, ):
    if img is None:
      try:
        import Image
      except ImportError:
        from PIL import Image
      if size is None:
        raise ValueError('please provide either an image or a size')
      img = Image.new('RGBA', size, "white")
    self.image = img
    self.draw = Draw(img)
    self.draw.setantialias(True)
    if size is None:
      self.size = self.draw.size
    else:
      self.size = size
    if imageType and imageType not in ('png', 'jpg'):
      raise ValueError('unsupported image type for agg canvas')
    self.drawType = imageType
    self.fileName = fileName

  def _doLine(self, p1, p2, pen, **kwargs):
    if kwargs.get('dash', (0, 0)) == (0, 0):
      self.draw.line((p1[0], p1[1], p2[0], p2[1]), pen)
    else:
      dash = kwargs['dash']
      pts = self._getLinePoints(p1, p2, dash)

      currDash = 0
      dashOn = True
      while currDash < (len(pts) - 1):
        if dashOn:
          p1 = pts[currDash]
          p2 = pts[currDash + 1]
          self.draw.line((p1[0], p1[1], p2[0], p2[1]), pen)
        currDash += 1
        dashOn = not dashOn

  def addCanvasLine(self, p1, p2, color=(0, 0, 0), color2=None, **kwargs):
    if color2 and color2 != color:
      mp = (p1[0] + p2[0]) / 2., (p1[1] + p2[1]) / 2.
      color = convertColor(color)
      self._doLine(p1, mp, Pen(color, kwargs.get('linewidth', 1)), **kwargs)
      color2 = convertColor(color2)
      self._doLine(mp, p2, Pen(color2, kwargs.get('linewidth', 1)), **kwargs)
    else:
      color = convertColor(color)
      self._doLine(p1, p2, Pen(color, kwargs.get('linewidth', 1)), **kwargs)

  def addCanvasText(self, text, pos, font, color=(0, 0, 0), **kwargs):
    orientation = kwargs.get('orientation', 'E')
    color = convertColor(color)
    aggFont = Font(color, faceMap[font.face], size=font.size * self.fontScale)

    blocks = list(re.finditer(r'\<(.+?)\>(.+?)\</\1\>', text))
    w, h = 0, 0
    supH = 0
    subH = 0
    if not len(blocks):
      w, h = self.draw.textsize(text, aggFont)
      tw, th = w, h
      offset = w * pos[2]
      dPos = pos[0] - w / 2. + offset, pos[1] - h / 2.
      self.draw.text(dPos, text, aggFont)
    else:
      dblocks = []
      idx = 0
      for block in blocks:
        blockStart, blockEnd = block.span(0)
        if blockStart != idx:
          # untagged text:
          tblock = text[idx:blockStart]
          tw, th = self.draw.textsize(tblock, aggFont)
          w += tw
          h = max(h, th)
          dblocks.append((tblock, '', tw, th))
        fmt = block.groups()[0]
        tblock = block.groups()[1]
        if fmt in ('sub', 'sup'):
          lFont = Font(color, faceMap[font.face], size=0.8 * font.size * self.fontScale)
        else:
          lFont = aggFont
        tw, th = self.draw.textsize(tblock, lFont)
        w += tw
        if fmt == 'sub':
          subH = max(subH, th)
        elif fmt == 'sup':
          supH = max(supH, th)
        else:
          h = max(h, th)
        dblocks.append((tblock, fmt, tw, th))
        idx = blockEnd
      if idx != len(text):
        # untagged text:
        tblock = text[idx:]
        tw, th = self.draw.textsize(tblock, aggFont)
        w += tw
        h = max(h, th)
        dblocks.append((tblock, '', tw, th))

      supH *= 0.5
      subH *= 0.5
      h += supH + subH
      offset = w * pos[2]
      dPos = [pos[0] - w / 2. + offset, pos[1] - h / 2.]
      if orientation == 'W':
        dPos = [pos[0] - w + offset, pos[1] - h / 2.]
      elif orientation == 'E':
        dPos = [pos[0] + offset, pos[1] - h / 2.]
      else:
        dPos = [pos[0] - w / 2 + offset, pos[1] - h / 2.]

      if supH:
        dPos[1] += supH
      for txt, fmt, tw, th in dblocks:
        tPos = dPos[:]
        if fmt == 'sub':
          tPos[1] += subH
        elif fmt == 'sup':
          tPos[1] -= supH
        if fmt in ('sub', 'sup'):
          lFont = Font(color, faceMap[font.face], size=0.8 * font.size * self.fontScale)
        else:
          lFont = aggFont
        self.draw.text(tPos, txt, lFont)
        dPos[0] += tw
    return (tw + th * .4, th + th * .4, offset)

  def addCanvasPolygon(self, ps, color=(0, 0, 0), fill=True, stroke=False, **kwargs):
    if not fill and not stroke:
      return
    dps = []
    for p in ps:
      dps.extend(p)
    color = convertColor(color)
    brush = None
    pen = None
    if fill:
      brush = Brush(color)
    if stroke:
      pen = Pen(color)
    self.draw.polygon(dps, pen, brush)

  def addCanvasDashedWedge(self, p1, p2, p3, dash=(2, 2), color=(0, 0, 0), color2=None, **kwargs):
    pen = Pen(color, kwargs.get('linewidth', 1))
    dash = (3, 3)
    pts1 = self._getLinePoints(p1, p2, dash)
    pts2 = self._getLinePoints(p1, p3, dash)

    if len(pts2) < len(pts1):
      pts2, pts1 = pts1, pts2

    for i in range(len(pts1)):
      self.draw.line((pts1[i][0], pts1[i][1], pts2[i][0], pts2[i][1]), pen)

  def flush(self):
    self.draw.flush()
    if self.fileName:
      self.image.save(self.fileName)
