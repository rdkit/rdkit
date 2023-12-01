#
#  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Created by Greg Landrum, October 2006
#
import os
import re
import weakref

from rdkit import Chem, RDConfig


class FGHierarchyNode(object):
  children = None
  name = ""
  label = ""
  pattern = None
  smarts = ""
  rxnSmarts = ""
  parent = None
  removalReaction = None

  def __init__(self, name, patt, smarts="", label="", rxnSmarts="", parent=None):
    self.name = name
    self.pattern = patt
    if parent:
      self.parent = weakref.ref(parent)
    self.label = label
    self.smarts = smarts
    self.children = []
    self.rxnSmarts = rxnSmarts

  def __len__(self):
    res = 1
    for child in self.children:
      res += len(child)
    return res


class FuncGroupFileParseError(ValueError):
  pass


groupDefns = {}
hierarchy = None
lastData = None
lastFilename = None


def BuildFuncGroupHierarchy(fileNm=None, data=None, force=False):
  global groupDefns, hierarchy, lastData, lastFilename
  if (not force and hierarchy and (not data or data == lastData)
      and (not fileNm or fileNm == lastFilename)):
    return hierarchy[:]
  lastData = data
  splitter = re.compile('\t+')

  if not fileNm and not data:
    fileNm = os.path.join(RDConfig.RDDataDir, 'Functional_Group_Hierarchy.txt')

  if fileNm:
    with open(fileNm, 'r') as inF:
      data = inF.readlines()
    lastFilename = fileNm
  elif data:
    data = data.splitlines()
  else:
    raise ValueError("need data or filename")

  groupDefns = {}
  res = []
  for lineNo, line in enumerate(data, 1):
    line = line.strip()
    line = line.split('//')[0]
    if not line:
      continue
    splitL = splitter.split(line)
    if len(splitL) < 3:
      raise FuncGroupFileParseError("Input line %d (%s) is not long enough." % (lineNo, repr(line)))
    label = splitL[0].strip()
    if label in groupDefns:
      raise FuncGroupFileParseError("Duplicate label on line %d." % lineNo)
    labelHierarchy = label.split('.')
    if len(labelHierarchy) > 1:
      for i in range(len(labelHierarchy) - 1):
        tmp = '.'.join(labelHierarchy[:i + 1])
        if tmp not in groupDefns:
          raise FuncGroupFileParseError("Hierarchy member %s (line %d) not found." % (tmp, lineNo))
      parent = groupDefns['.'.join(labelHierarchy[:-1])]
    else:
      parent = None
    smarts = splitL[1]
    patt = Chem.MolFromSmarts(smarts)
    if not patt:
      raise FuncGroupFileParseError('Smarts "%s" (line %d) could not be parsed.' % (smarts, lineNo))

    name = splitL[2].strip()

    rxnSmarts = ''
    if len(splitL) > 3:
      rxnSmarts = splitL[3]

    node = FGHierarchyNode(name, patt, smarts=smarts, label=label, parent=parent,
                           rxnSmarts=rxnSmarts)
    if parent:
      parent.children.append(node)
    else:
      res.append(node)
    groupDefns[label] = node
  hierarchy = res[:]
  return res


def _SetNodeBits(mol, node, res, idx):
  ms = mol.GetSubstructMatches(node.pattern)
  count = 0
  seen = {}
  for m in ms:
    if m[0] not in seen:
      count += 1
      seen[m[0]] = 1
  if count:
    res[idx] = count
    idx += 1
    for child in node.children:
      idx = _SetNodeBits(mol, child, res, idx)
  else:
    idx += len(node)
  return idx


def CreateMolFingerprint(mol, hierarchy):
  totL = 0
  for entry in hierarchy:
    totL += len(entry)
  res = [0] * totL
  idx = 0
  for entry in hierarchy:
    idx = _SetNodeBits(mol, entry, res, idx)
  return res
