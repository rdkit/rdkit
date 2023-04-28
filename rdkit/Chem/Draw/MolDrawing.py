# $Id$
#
#  Copyright (C) 2008-2011 Greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import copy
import functools
import math

import numpy

from rdkit import Chem


def cmp(t1, t2):
  return (t1 < t2) * -1 or (t1 > t2) * 1


periodicTable = Chem.GetPeriodicTable()


class Font(object):

  def __init__(self, face=None, size=None, name=None, weight=None):
    self.face = face or 'sans'
    self.size = size or '12'
    self.weight = weight or 'normal'
    self.name = name


class DrawingOptions(object):
  dotsPerAngstrom = 30
  useFraction = 0.85

  atomLabelFontFace = "sans"
  atomLabelFontSize = 12
  atomLabelMinFontSize = 7
  atomLabelDeuteriumTritium = False

  bondLineWidth = 1.2
  dblBondOffset = .25
  dblBondLengthFrac = .8

  defaultColor = (1, 0, 0)
  selectColor = (1, 0, 0)
  bgColor = (1, 1, 1)

  colorBonds = True
  noCarbonSymbols = True
  includeAtomNumbers = False
  atomNumberOffset = 0
  radicalSymbol = u'\u2219'

  dash = (4, 4)

  wedgeDashedBonds = True
  showUnknownDoubleBonds = True

  # used to adjust overall scaling for molecules that have been laid out with non-standard
  # bond lengths
  coordScale = 1.0

  elemDict = {
    1: (0.55, 0.55, 0.55),
    7: (0, 0, 1),
    8: (1, 0, 0),
    9: (.2, .8, .8),
    15: (1, .5, 0),
    16: (.8, .8, 0),
    17: (0, .8, 0),
    35: (.5, .3, .1),
    53: (.63, .12, .94),
    0: (.5, .5, .5),
  }


class MolDrawing(object):

  def __init__(self, canvas=None, drawingOptions=None):
    self.canvas = canvas
    self.canvasSize = None
    if canvas:
      self.canvasSize = canvas.size
    self.drawingOptions = drawingOptions or DrawingOptions()

    self.atomPs = {}
    self.boundingBoxes = {}

    if self.drawingOptions.bgColor is not None:
      self.canvas.addCanvasPolygon(
        ((0, 0), (canvas.size[0], 0), (canvas.size[0], canvas.size[1]), (0, canvas.size[1])),
        color=self.drawingOptions.bgColor, fill=True, stroke=False)

  def transformPoint(self, pos):
    res = [0, 0]
    res[0] = (pos[0] + self.molTrans[0]
              ) * self.currDotsPerAngstrom * self.drawingOptions.useFraction + self.drawingTrans[0]
    res[1] = self.canvasSize[1] - (
      (pos[1] + self.molTrans[1]) * self.currDotsPerAngstrom * self.drawingOptions.useFraction +
      self.drawingTrans[1])
    return res

  def _getBondOffset(self, p1, p2):
    # get the vector between the points:
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]

    # figure out the angle and the perpendicular:
    ang = math.atan2(dy, dx)
    perp = ang + math.pi / 2.

    # here's the offset for the parallel bond:
    offsetX = math.cos(perp) * self.drawingOptions.dblBondOffset * self.currDotsPerAngstrom
    offsetY = math.sin(perp) * self.drawingOptions.dblBondOffset * self.currDotsPerAngstrom

    return perp, offsetX, offsetY

  def _getOffsetBondPts(self, p1, p2, offsetX, offsetY, lenFrac=None):
    lenFrac = lenFrac or self.drawingOptions.dblBondLengthFrac

    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    # ----
    # now figure out where to start and end it:

    # offset the start point:
    fracP1 = p1[0] + offsetX, p1[1] + offsetY

    # now move a portion of the way along the line to the neighbor:
    frac = (1. - lenFrac) / 2
    fracP1 = fracP1[0] + dx * frac, fracP1[1] + dy * frac
    fracP2 = fracP1[0] + dx * lenFrac, fracP1[1] + dy * lenFrac
    return fracP1, fracP2

  def _offsetDblBond(self, p1, p2, bond, a1, a2, conf, direction=1, lenFrac=None):
    perp, offsetX, offsetY = self._getBondOffset(p1, p2)
    offsetX = offsetX * direction
    offsetY = offsetY * direction

    # if we're a ring bond, we may need to flip over to the other side:
    if bond.IsInRing():
      bondIdx = bond.GetIdx()
      a2Idx = a2.GetIdx()
      # find a ring bond from a1 to an atom other than a2:
      for otherBond in a1.GetBonds():
        if otherBond.GetIdx() != bondIdx and otherBond.IsInRing():
          sharedRing = False
          for ring in self.bondRings:
            if bondIdx in ring and otherBond.GetIdx() in ring:
              sharedRing = True
              break
          if not sharedRing:
            continue
          a3 = otherBond.GetOtherAtom(a1)
          if a3.GetIdx() != a2Idx:
            p3 = self.transformPoint(
              conf.GetAtomPosition(a3.GetIdx()) * self.drawingOptions.coordScale)
            dx2 = p3[0] - p1[0]
            dy2 = p3[1] - p1[1]
            dotP = dx2 * offsetX + dy2 * offsetY
            if dotP < 0:
              perp += math.pi
              offsetX = math.cos(
                perp) * self.drawingOptions.dblBondOffset * self.currDotsPerAngstrom
              offsetY = math.sin(
                perp) * self.drawingOptions.dblBondOffset * self.currDotsPerAngstrom

    fracP1, fracP2 = self._getOffsetBondPts(p1, p2, offsetX, offsetY, lenFrac=lenFrac)
    return fracP1, fracP2

  def _getBondAttachmentCoordinates(self, p1, p2, labelSize):
    newpos = [None, None]
    if labelSize is not None:
      labelSizeOffset = [
        labelSize[0][0] / 2 + (cmp(p2[0], p1[0]) * labelSize[0][2]), labelSize[0][1] / 2
      ]
      if p1[1] == p2[1]:
        newpos[0] = p1[0] + cmp(p2[0], p1[0]) * labelSizeOffset[0]
      else:
        if abs(labelSizeOffset[1] * (p2[0] - p1[0]) / (p2[1] - p1[1])) < labelSizeOffset[0]:
          newpos[0] = p1[0] + cmp(p2[0], p1[0]) * abs(labelSizeOffset[1] * (p2[0] - p1[0]) /
                                                      (p2[1] - p1[1]))
        else:
          newpos[0] = p1[0] + cmp(p2[0], p1[0]) * labelSizeOffset[0]
      if p1[0] == p2[0]:
        newpos[1] = p1[1] + cmp(p2[1], p1[1]) * labelSizeOffset[1]
      else:
        if abs(labelSizeOffset[0] * (p1[1] - p2[1]) / (p2[0] - p1[0])) < labelSizeOffset[1]:
          newpos[1] = p1[1] + cmp(p2[1], p1[1]) * abs(labelSizeOffset[0] * (p1[1] - p2[1]) /
                                                      (p2[0] - p1[0]))
        else:
          newpos[1] = p1[1] + cmp(p2[1], p1[1]) * labelSizeOffset[1]
    else:
      newpos = copy.deepcopy(p1)
    return newpos

  def _drawWedgedBond(self, bond, pos, nbrPos, width=None, color=None, dash=None):
    width = width or self.drawingOptions.bondLineWidth
    color = color or self.drawingOptions.defaultColor
    _, offsetX, offsetY = self._getBondOffset(pos, nbrPos)
    offsetX *= .75
    offsetY *= .75
    poly = ((pos[0], pos[1]), (nbrPos[0] + offsetX, nbrPos[1] + offsetY), (nbrPos[0] - offsetX,
                                                                           nbrPos[1] - offsetY))
    # canvas.drawPolygon(poly,edgeColor=color,edgeWidth=1,fillColor=color,closed=1)
    if not dash:
      self.canvas.addCanvasPolygon(poly, color=color)
    elif self.drawingOptions.wedgeDashedBonds and self.canvas.addCanvasDashedWedge:
      self.canvas.addCanvasDashedWedge(poly[0], poly[1], poly[2], color=color)
    else:
      self.canvas.addCanvasLine(pos, nbrPos, linewidth=width * 2, color=color, dashes=dash)

  def _drawBond(self, bond, atom, nbr, pos, nbrPos, conf, width=None, color=None, color2=None,
                labelSize1=None, labelSize2=None):
    width = width or self.drawingOptions.bondLineWidth
    color = color or self.drawingOptions.defaultColor
    color2 = color2 or self.drawingOptions.defaultColor
    p1_raw = copy.deepcopy(pos)
    p2_raw = copy.deepcopy(nbrPos)
    newpos = self._getBondAttachmentCoordinates(p1_raw, p2_raw, labelSize1)
    newnbrPos = self._getBondAttachmentCoordinates(p2_raw, p1_raw, labelSize2)
    addDefaultLine = functools.partial(self.canvas.addCanvasLine, linewidth=width, color=color,
                                       color2=color2)
    bType = bond.GetBondType()
    if bType == Chem.BondType.SINGLE:
      bDir = bond.GetBondDir()
      if bDir in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH):
        # if the bond is "backwards", change the drawing direction:
        if bond.GetBeginAtom().GetChiralTag() in (Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                                                  Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
          p1, p2 = newpos, newnbrPos
          wcolor = color
        else:
          p2, p1 = newpos, newnbrPos
          wcolor = color2
        if bDir == Chem.BondDir.BEGINWEDGE:
          self._drawWedgedBond(bond, p1, p2, color=wcolor, width=width)
        elif bDir == Chem.BondDir.BEGINDASH:
          self._drawWedgedBond(bond, p1, p2, color=wcolor, width=width,
                               dash=self.drawingOptions.dash)
      else:
        addDefaultLine(newpos, newnbrPos)
    elif bType == Chem.BondType.DOUBLE:
      crossBond = (self.drawingOptions.showUnknownDoubleBonds
                   and bond.GetStereo() == Chem.BondStereo.STEREOANY)
      if (not crossBond and (bond.IsInRing() or
                             (atom.GetDegree() != 1 and bond.GetOtherAtom(atom).GetDegree() != 1))):
        addDefaultLine(newpos, newnbrPos)
        fp1, fp2 = self._offsetDblBond(newpos, newnbrPos, bond, atom, nbr, conf)
        addDefaultLine(fp1, fp2)
      else:
        fp1, fp2 = self._offsetDblBond(newpos, newnbrPos, bond, atom, nbr, conf, direction=.5,
                                       lenFrac=1.0)
        fp3, fp4 = self._offsetDblBond(newpos, newnbrPos, bond, atom, nbr, conf, direction=-.5,
                                       lenFrac=1.0)
        if crossBond:
          fp2, fp4 = fp4, fp2
        addDefaultLine(fp1, fp2)
        addDefaultLine(fp3, fp4)

    elif bType == Chem.BondType.AROMATIC:
      addDefaultLine(newpos, newnbrPos)
      fp1, fp2 = self._offsetDblBond(newpos, newnbrPos, bond, atom, nbr, conf)
      addDefaultLine(fp1, fp2, dash=self.drawingOptions.dash)
    elif bType == Chem.BondType.TRIPLE:
      addDefaultLine(newpos, newnbrPos)
      fp1, fp2 = self._offsetDblBond(newpos, newnbrPos, bond, atom, nbr, conf)
      addDefaultLine(fp1, fp2)
      fp1, fp2 = self._offsetDblBond(newpos, newnbrPos, bond, atom, nbr, conf, direction=-1)
      addDefaultLine(fp1, fp2)
    else:
      addDefaultLine(newpos, newnbrPos, dash=(1, 2))

  def scaleAndCenter(self, mol, conf, coordCenter=False, canvasSize=None, ignoreHs=False):
    canvasSize = canvasSize or self.canvasSize
    xAccum = 0
    yAccum = 0
    minX = 1e8
    minY = 1e8
    maxX = -1e8
    maxY = -1e8

    nAts = mol.GetNumAtoms()
    for i in range(nAts):
      if ignoreHs and mol.GetAtomWithIdx(i).GetAtomicNum() == 1:
        continue
      pos = conf.GetAtomPosition(i) * self.drawingOptions.coordScale
      xAccum += pos[0]
      yAccum += pos[1]
      minX = min(minX, pos[0])
      minY = min(minY, pos[1])
      maxX = max(maxX, pos[0])
      maxY = max(maxY, pos[1])

    dx = abs(maxX - minX)
    dy = abs(maxY - minY)
    xSize = dx * self.currDotsPerAngstrom
    ySize = dy * self.currDotsPerAngstrom

    if coordCenter:
      molTrans = -xAccum / nAts, -yAccum / nAts
    else:
      molTrans = -(minX + (maxX - minX) / 2), -(minY + (maxY - minY) / 2)
    self.molTrans = molTrans

    if xSize >= .95 * canvasSize[0]:
      scale = .9 * canvasSize[0] / xSize
      xSize *= scale
      ySize *= scale
      self.currDotsPerAngstrom *= scale
      self.currAtomLabelFontSize = max(self.currAtomLabelFontSize * scale,
                                       self.drawingOptions.atomLabelMinFontSize)
    if ySize >= .95 * canvasSize[1]:
      scale = .9 * canvasSize[1] / ySize
      xSize *= scale
      ySize *= scale
      self.currDotsPerAngstrom *= scale
      self.currAtomLabelFontSize = max(self.currAtomLabelFontSize * scale,
                                       self.drawingOptions.atomLabelMinFontSize)
    drawingTrans = canvasSize[0] / 2, canvasSize[1] / 2
    self.drawingTrans = drawingTrans

  def _drawLabel(self, label, pos, baseOffset, font, color=None, **kwargs):
    color = color or self.drawingOptions.defaultColor
    x1 = pos[0]
    y1 = pos[1]
    labelSize = self.canvas.addCanvasText(label, (x1, y1, baseOffset), font, color, **kwargs)
    return labelSize

  def AddMol(self, mol, centerIt=True, molTrans=None, drawingTrans=None, highlightAtoms=[],
             confId=-1, flagCloseContactsDist=2, highlightMap=None, ignoreHs=False,
             highlightBonds=[], **kwargs):
    """Set the molecule to be drawn.

    Parameters:
      hightlightAtoms -- list of atoms to highlight (default [])
      highlightMap -- dictionary of (atom, color) pairs (default None)

    Notes:
      - specifying centerIt will cause molTrans and drawingTrans to be ignored
    """
    conf = mol.GetConformer(confId)
    if 'coordScale' in kwargs:
      self.drawingOptions.coordScale = kwargs['coordScale']

    self.currDotsPerAngstrom = self.drawingOptions.dotsPerAngstrom
    self.currAtomLabelFontSize = self.drawingOptions.atomLabelFontSize
    if centerIt:
      self.scaleAndCenter(mol, conf, ignoreHs=ignoreHs)
    else:
      self.molTrans = molTrans or (0, 0)
      self.drawingTrans = drawingTrans or (0, 0)

    font = Font(face=self.drawingOptions.atomLabelFontFace, size=self.currAtomLabelFontSize)

    obds = None
    if not mol.HasProp('_drawingBondsWedged'):
      # this is going to modify the molecule, get ready to undo that
      obds = [x.GetBondDir() for x in mol.GetBonds()]
      Chem.WedgeMolBonds(mol, conf)

    includeAtomNumbers = kwargs.get('includeAtomNumbers', self.drawingOptions.includeAtomNumbers)
    self.atomPs[mol] = {}
    self.boundingBoxes[mol] = [0] * 4
    self.activeMol = mol
    self.bondRings = mol.GetRingInfo().BondRings()
    labelSizes = {}
    for atom in mol.GetAtoms():
      labelSizes[atom.GetIdx()] = None
      if ignoreHs and atom.GetAtomicNum() == 1:
        drawAtom = False
      else:
        drawAtom = True
      idx = atom.GetIdx()
      pos = self.atomPs[mol].get(idx, None)
      if pos is None:
        pos = self.transformPoint(conf.GetAtomPosition(idx) * self.drawingOptions.coordScale)
        self.atomPs[mol][idx] = pos
        if drawAtom:
          self.boundingBoxes[mol][0] = min(self.boundingBoxes[mol][0], pos[0])
          self.boundingBoxes[mol][1] = min(self.boundingBoxes[mol][1], pos[1])
          self.boundingBoxes[mol][2] = max(self.boundingBoxes[mol][2], pos[0])
          self.boundingBoxes[mol][3] = max(self.boundingBoxes[mol][3], pos[1])

      if not drawAtom:
        continue
      nbrSum = [0, 0]
      for bond in atom.GetBonds():
        nbr = bond.GetOtherAtom(atom)
        if ignoreHs and nbr.GetAtomicNum() == 1:
          continue
        nbrIdx = nbr.GetIdx()
        if nbrIdx > idx:
          nbrPos = self.atomPs[mol].get(nbrIdx, None)
          if nbrPos is None:
            nbrPos = self.transformPoint(
              conf.GetAtomPosition(nbrIdx) * self.drawingOptions.coordScale)
            self.atomPs[mol][nbrIdx] = nbrPos
            self.boundingBoxes[mol][0] = min(self.boundingBoxes[mol][0], nbrPos[0])
            self.boundingBoxes[mol][1] = min(self.boundingBoxes[mol][1], nbrPos[1])
            self.boundingBoxes[mol][2] = max(self.boundingBoxes[mol][2], nbrPos[0])
            self.boundingBoxes[mol][3] = max(self.boundingBoxes[mol][3], nbrPos[1])

        else:
          nbrPos = self.atomPs[mol][nbrIdx]
        nbrSum[0] += nbrPos[0] - pos[0]
        nbrSum[1] += nbrPos[1] - pos[1]

      iso = atom.GetIsotope()
      labelIt = (not self.drawingOptions.noCarbonSymbols or iso or atom.GetAtomicNum() != 6
                 or atom.GetFormalCharge() != 0 or atom.GetNumRadicalElectrons()
                 or includeAtomNumbers or atom.HasProp('molAtomMapNumber') or atom.GetDegree() == 0)
      orient = ''
      if labelIt:
        baseOffset = 0
        if includeAtomNumbers:
          symbol = str(atom.GetIdx())
          symbolLength = len(symbol)
        else:
          base = atom.GetSymbol()
          if base == 'H' and (iso == 2
                              or iso == 3) and self.drawingOptions.atomLabelDeuteriumTritium:
            if iso == 2:
              base = 'D'
            else:
              base = 'T'
            iso = 0
          symbolLength = len(base)

          nHs = 0
          if not atom.HasQuery():
            nHs = atom.GetTotalNumHs()
          hs = ''
          if nHs == 1:
            hs = 'H'
            symbolLength += 1
          elif nHs > 1:
            hs = 'H<sub>%d</sub>' % nHs
            symbolLength += 1 + len(str(nHs))

          chg = atom.GetFormalCharge()
          if chg == 0:
            chg = ''
          elif chg == 1:
            chg = '+'
          elif chg == -1:
            chg = '-'
          else:
            chg = '%+d' % chg
          symbolLength += len(chg)
          if chg:
            chg = '<sup>%s</sup>' % chg

          if atom.GetNumRadicalElectrons():
            rad = self.drawingOptions.radicalSymbol * atom.GetNumRadicalElectrons()
            rad = '<sup>%s</sup>' % rad
            symbolLength += atom.GetNumRadicalElectrons()
          else:
            rad = ''

          isotope = ''
          isotopeLength = 0
          if iso:
            isotope = '<sup>%d</sup>' % atom.GetIsotope()
            isotopeLength = len(str(atom.GetIsotope()))
            symbolLength += isotopeLength
          mapNum = ''
          mapNumLength = 0
          if atom.HasProp('molAtomMapNumber'):
            mapNum = ':' + atom.GetProp('molAtomMapNumber')
            mapNumLength = 1 + len(str(atom.GetProp('molAtomMapNumber')))
            symbolLength += mapNumLength
          deg = atom.GetDegree()
          # This should be done in a better way in the future:
          # 'baseOffset' should be determined by getting the size of 'isotope' and
          # the size of 'base', or the size of 'mapNum' and the size of 'base'
          # (depending on 'deg' and 'nbrSum[0]') in order to determine the exact
          # position of the base

          if deg == 0:
            tSym = periodicTable.GetElementSymbol(atom.GetAtomicNum())
            if tSym in ('O', 'S', 'Se', 'Te', 'F', 'Cl', 'Br', 'I', 'At'):
              symbol = '%s%s%s%s%s%s' % (hs, isotope, base, chg, rad, mapNum)
            else:
              symbol = '%s%s%s%s%s%s' % (isotope, base, hs, chg, rad, mapNum)
          elif deg > 1 or nbrSum[0] < 1:
            symbol = '%s%s%s%s%s%s' % (isotope, base, hs, chg, rad, mapNum)
            baseOffset = 0.5 - (isotopeLength + len(base) / 2.) / symbolLength
          else:
            symbol = '%s%s%s%s%s%s' % (rad, chg, hs, isotope, base, mapNum)
            baseOffset = -0.5 + (mapNumLength + len(base) / 2.) / symbolLength

          if deg == 1:
            if abs(nbrSum[1]) > 1:
              islope = nbrSum[0] / abs(nbrSum[1])
            else:
              islope = nbrSum[0]
            if abs(islope) > .3:
              if islope > 0:
                orient = 'W'
              else:
                orient = 'E'
            elif abs(nbrSum[1]) > 10:
              if nbrSum[1] > 0:
                orient = 'N'
              else:
                orient = 'S'
          else:
            orient = 'C'

        if highlightMap and idx in highlightMap:
          color = highlightMap[idx]
        elif highlightAtoms and idx in highlightAtoms:
          color = self.drawingOptions.selectColor
        else:
          color = self.drawingOptions.elemDict.get(atom.GetAtomicNum(), (0, 0, 0))
        labelSize = self._drawLabel(symbol, pos, baseOffset, font, color=color, orientation=orient)
        labelSizes[atom.GetIdx()] = [labelSize, orient]

    for bond in mol.GetBonds():
      atom, idx = bond.GetBeginAtom(), bond.GetBeginAtomIdx()
      nbr, nbrIdx = bond.GetEndAtom(), bond.GetEndAtomIdx()
      pos = self.atomPs[mol].get(idx, None)
      nbrPos = self.atomPs[mol].get(nbrIdx, None)
      if highlightBonds and bond.GetIdx() in highlightBonds:
        width = 2.0 * self.drawingOptions.bondLineWidth
        color = self.drawingOptions.selectColor
        color2 = self.drawingOptions.selectColor
      elif highlightAtoms and idx in highlightAtoms and nbrIdx in highlightAtoms:
        width = 2.0 * self.drawingOptions.bondLineWidth
        color = self.drawingOptions.selectColor
        color2 = self.drawingOptions.selectColor
      elif highlightMap is not None and idx in highlightMap and nbrIdx in highlightMap:
        width = 2.0 * self.drawingOptions.bondLineWidth
        color = highlightMap[idx]
        color2 = highlightMap[nbrIdx]
      else:
        width = self.drawingOptions.bondLineWidth
        if self.drawingOptions.colorBonds:
          color = self.drawingOptions.elemDict.get(atom.GetAtomicNum(), (0, 0, 0))
          color2 = self.drawingOptions.elemDict.get(nbr.GetAtomicNum(), (0, 0, 0))
        else:
          color = self.drawingOptions.defaultColor
          color2 = color
      self._drawBond(bond, atom, nbr, pos, nbrPos, conf, color=color, width=width, color2=color2,
                     labelSize1=labelSizes[idx], labelSize2=labelSizes[nbrIdx])

  # if we modified the bond wedging state, undo those changes now
    if obds:
      for i, d in enumerate(obds):
        mol.GetBondWithIdx(i).SetBondDir(d)

    if flagCloseContactsDist > 0:
      tol = flagCloseContactsDist * flagCloseContactsDist
      for i, _ in enumerate(mol.GetAtoms()):
        pi = numpy.array(self.atomPs[mol][i])
        for j in range(i + 1, mol.GetNumAtoms()):
          pj = numpy.array(self.atomPs[mol][j])
          d = pj - pi
          dist2 = d[0] * d[0] + d[1] * d[1]
          if dist2 <= tol:
            self.canvas.addCanvasPolygon(
              ((pi[0] - 2 * flagCloseContactsDist, pi[1] - 2 * flagCloseContactsDist),
               (pi[0] + 2 * flagCloseContactsDist, pi[1] - 2 * flagCloseContactsDist),
               (pi[0] + 2 * flagCloseContactsDist, pi[1] + 2 * flagCloseContactsDist),
               (pi[0] - 2 * flagCloseContactsDist, pi[1] + 2 * flagCloseContactsDist)),
              color=(1., 0, 0), fill=False, stroke=True)
