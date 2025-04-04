#
#  Copyright (c) 2017, Novartis Institutes for BioMedical Research Inc.
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
# Created by Nadine Schneider & Peter Ertl, July 2017

import re
from collections import Counter, defaultdict, namedtuple

import numpy as np
import seaborn as sns
from numpy import linalg

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D


# build an svg grid image to print
def _svgsToGrid(svgs, labels, svgsPerRow=4, molSize=(250, 150), fontSize=12):

  matcher = re.compile(r'^(<.*>\n)(<rect .*</rect>\n)(.*)</svg>', re.DOTALL)
  hdr = ''
  ftr = '</svg>'
  rect = ''
  nRows = len(svgs) // svgsPerRow
  if len(svgs) % svgsPerRow:
    nRows += 1
  blocks = [''] * (nRows * svgsPerRow)
  labelSizeDist = fontSize * 5
  fullSize = (svgsPerRow * (molSize[0] + molSize[0] / 10.0), nRows * (molSize[1] + labelSizeDist))

  count = 0
  for svg, name in zip(svgs, labels):
    h, r, b = matcher.match(svg).groups()
    if hdr == '':
      hdr = h.replace("width='{}px'".format(molSize[0]), "width='{}px'".format(fullSize[0]))
      hdr = hdr.replace("height='{}px'".format(molSize[1]), "height='{}px'".format(fullSize[1]))
    if rect == '':
      rect = r

    tspanFmt = '<tspan x="{0}" y="{1}">{2}</tspan>'
    names = name.split('|')
    legend = []
    legend.append(
      '<text font-family="sans-serif" font-size="{}px" text-anchor="middle" fill="black">'.format(
        fontSize))
    legend.append(tspanFmt.format(molSize[0] / 2., molSize[1] + fontSize * 2, names[0]))
    if len(names) > 1:
      legend.append(tspanFmt.format(molSize[0] / 2., molSize[1] + fontSize * 3.5, names[1]))
    legend.append('</text>')
    legend = '\n'.join(legend)

    blocks[count] = b + legend
    count += 1

  for i, elem in enumerate(blocks):
    row = i // svgsPerRow
    col = i % svgsPerRow
    elem = rect + elem
    blocks[i] = '<g transform="translate(%d,%d)" >%s</g>' % (col *
                                                             (molSize[0] + molSize[0] / 10.0), row *
                                                             (molSize[1] + labelSizeDist), elem)
  res = hdr + '\n'.join(blocks) + ftr
  return res


def determineAtomSubstituents(atomID, mol, distanceMatrix, verbose=False):
  atomPaths = distanceMatrix[atomID]
  # determine the direct neighbors of the atom
  neighbors = [n for n, i in enumerate(atomPaths) if i == 1]
  # store the ids of the neighbors (substituents)
  subs = defaultdict(list)
  # track in how many substituents an atom is involved (can happen in rings)
  sharedNeighbors = defaultdict(int)
  # determine the max path length for each substituent
  maxShell = defaultdict(int)
  for n in neighbors:
    subs[n].append(n)
    sharedNeighbors[n] += 1
    maxShell[n] = 0
  # second shell of neighbors
  mindist = 2
  # max distance from atom
  maxdist = int(np.max(atomPaths))
  for d in range(mindist, maxdist + 1):
    if verbose:
      print("Shell: ", d)
    newShell = [n for n, i in enumerate(atomPaths) if i == d]
    for aidx in newShell:
      if verbose:
        print("Atom ", aidx, " in shell ", d)
      atom = mol.GetAtomWithIdx(aidx)
      # find neighbors of the current atom that are part of the substituent already
      for n in atom.GetNeighbors():
        nidx = n.GetIdx()
        for k, v in subs.items():
          # is the neighbor in the substituent and is not in the same shell as the current atom
          # and we haven't added the current atom already then put it in the correct substituent list
          if nidx in v and nidx not in newShell and aidx not in v:
            subs[k].append(aidx)
            sharedNeighbors[aidx] += 1
            maxShell[k] = d
            if verbose:
              print("Atom ", aidx, " assigned to ", nidx)
  if verbose:
    print(subs)
    print(sharedNeighbors)

  return subs, sharedNeighbors, maxShell


def _getSizeOfSubstituents(sub, sharedNeighbors, weighdownShared=True):
  if weighdownShared:
    return sum(1.0 / sharedNeighbors[a] for a in sub)
  else:
    return len(sub)


def getBondsSubstituent(mol, atoms):
  bonds = []
  for b in mol.GetBonds():
    a1 = b.GetBeginAtomIdx()
    a2 = b.GetEndAtomIdx()
    if a1 in atoms and a2 in atoms:
      bonds.append(b.GetIdx())
  return bonds


def getNumAromaticBondsSubstituent(mol, subAtoms):
  return sum(1 for b in getBondsSubstituent(mol, subAtoms) if mol.GetBondWithIdx(b).GetIsAromatic())


def getNumRotatableBondsSubstituent(mol, subAtoms):
  rotatableBond = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')
  matches = mol.GetSubstructMatches(rotatableBond)
  numRotBonds = 0
  for a1, a2 in matches:
    if a1 in subAtoms and a2 in subAtoms:
      numRotBonds += 1
  return numRotBonds


substituentDescriptor = namedtuple('substituentDescriptor', [
  'size', 'relSize', 'numNO', 'relNumNO', 'relNumNO_2', 'pathLength', 'relPathLength',
  'relPathLength_2', 'sharedNeighbors', 'numRotBonds', 'numAroBonds'
])


def calcSizeSubstituents(mol, subs, sharedNeighbors, maxShell):
  sizeDict = defaultdict()
  numAtoms = mol.GetNumAtoms()
  for sidx, sub in sorted(subs.items(), key=lambda x: len(x[1])):
    size = _getSizeOfSubstituents(sub, sharedNeighbors)
    numNOs = 0
    numShared = 0
    # determine the number of oxygen and nitrogen atoms
    for i in sub:
      if mol.GetAtomWithIdx(i).GetAtomicNum() in [7, 8]:
        numNOs += 1.0 / sharedNeighbors[i]
      if sharedNeighbors[i] > 1:
        numShared += 1
    numRotBs = getNumRotatableBondsSubstituent(mol, set(sub))
    aroBonds = getNumAromaticBondsSubstituent(mol, set(sub))
    # fill the substituentDescriptor tuple
    sizeDict[sidx] = substituentDescriptor(
      size=size, relSize=size / numAtoms, numNO=numNOs, relNumNO=numNOs / numAtoms,
      relNumNO_2=numNOs / size, pathLength=maxShell[sidx], relPathLength=maxShell[sidx] / numAtoms,
      relPathLength_2=maxShell[sidx] / size, sharedNeighbors=numShared, numRotBonds=numRotBs,
      numAroBonds=aroBonds)
  # if we have less then 4 substituents the missing ones need to be an hydrogen atoms
  if len(sizeDict) < 4:
    for i in range(4 - len(sizeDict)):
      sizeDict['H' + str(i)] = substituentDescriptor(size=0, relSize=0, numNO=0, relNumNO=0,
                                                     relNumNO_2=0, pathLength=0, relPathLength=0,
                                                     relPathLength_2=0, sharedNeighbors=0,
                                                     numRotBonds=0, numAroBonds=0)
  return sizeDict


# Visualization of the substituents
def visualizeSubstituentsGrid(
    mol,
    aIdx,
    molSize=(300, 150),
    kekulize=True,
):
  dists = Chem.GetDistanceMatrix(mol)
  idxChiral = Chem.FindMolChiralCenters(mol)[0][0]
  subs, sharedNeighbors, maxShell = determineAtomSubstituents(aIdx, mol, dists, False)

  colors = sns.husl_palette(len(subs), s=.6)
  mc = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=kekulize)
  count = 0
  svgs = []
  labels = []
  for sub in sorted(subs.values(), key=lambda x: _getSizeOfSubstituents(x, sharedNeighbors)):
    color = tuple(colors[count])
    count += 1
    atColors = {atom: color for atom in sub}

    bonds = getBondsSubstituent(mol, set(sub))
    bnColors = {bond: color for bond in bonds}

    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0], molSize[1])
    drawer.DrawMolecule(mc, highlightAtoms=atColors.keys(), highlightAtomColors=atColors,
                        highlightBonds=bonds, highlightBondColors=bnColors)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    svgs.append(svg.replace('svg:', ''))
    labels.append("Substituent " + str(count) + " (#atoms: " + str(len(sub)) + ", size normed: " +
                  str(_getSizeOfSubstituents(sub, sharedNeighbors)) + ")")
  return _svgsToGrid(svgs, labels, svgsPerRow=len(svgs), molSize=molSize, fontSize=12)


def visualizeChiralSubstituentsGrid(mol):
  idxChiral = Chem.FindMolChiralCenters(mol)[0][0]
  return visualizeSubstituentsGrid(mol, idxChiral)


# Chiral moment descriptor
def calcSP3CarbonSubstituentMoment(subSizes):

  if len(subSizes) != 4:
    raise ValueError(
      'Function "calcSP3CarbonSubstituentMoment" expects an array of size 4 as parameter')

  # tetrahedron unit vectors
  x1 = np.array([1, 1, 1])
  x2 = np.array([-1, 1, -1])
  x3 = np.array([1, -1, -1])
  x4 = np.array([-1, -1, 1])

  substituentMoment = linalg.norm((subSizes[0] * x1) + (subSizes[1] * x2) + (subSizes[2] * x3) +
                                  (subSizes[3] * x4))
  return substituentMoment


def calculateChiralDescriptors(mol, idxChiral, dists, verbose=False):

  desc = {}
  subs, sharedNeighbors, maxShell = determineAtomSubstituents(idxChiral, mol, dists, verbose)
  sizes = calcSizeSubstituents(mol, subs, sharedNeighbors, maxShell)
  paths = dists[idxChiral]
  # set some basic descriptors
  desc['numAtoms'] = mol.GetNumAtoms()
  desc['numBonds'] = mol.GetNumBonds()
  desc['numRotBonds'] = AllChem.CalcNumRotatableBonds(mol)
  desc['ringChiralCenter'] = int(mol.GetAtomWithIdx(idxChiral).IsInRing())
  # determine the max path length in the molecule and the mean pairwise distance of all atom pairs
  desc['meanDist'] = np.sum(dists) / ((desc['numAtoms'] - 1) * (desc['numAtoms']))
  desc['maxDist'] = int(np.max(dists))
  # determine the max path length from the chiral center and the mean pairwise distance of
  # all atom pairs from the chiral center
  desc['meanDistFromCC'] = np.sum(paths) / (desc['numAtoms'] - 1)
  desc['maxDistfromCC'] = int(np.max(paths))
  # determine the number of neighbors per shell/distance level
  nlevels = Counter(paths.astype(int))
  # consider the levels until a path length of 10
  for i in range(1, 11):
    desc['nLevel' + str(i)] = nlevels[i]
  # determine the number of nitrogen and oxygen atoms in a certain level around the chiral center
  for i in range(1, 4):
    desc['phLevel' + str(i)] = len(
      [n for n, j in enumerate(paths) if j == i and mol.GetAtomWithIdx(n).GetAtomicNum() in [7, 8]])
  # determine the number of aromatic atoms in a certain level around the chiral center
  for i in range(1, 4):
    desc['arLevel' + str(i)] = len(
      [n for n, j in enumerate(paths) if j == i and mol.GetAtomWithIdx(n).GetIsAromatic()])
  # set the size descriptors for each substituent, sort them from smallest to largest
  for n, v in enumerate(sorted(sizes.values(), key=lambda x: x.size), 1):
    sn = 's' + str(n)
    desc[sn + '_size'] = v.size
    desc[sn + '_relSize'] = v.relSize
    desc[sn + '_phSize'] = v.numNO
    desc[sn + '_phRelSize'] = v.relNumNO
    desc[sn + '_phRelSize_2'] = v.relNumNO_2
    desc[sn + '_pathLength'] = v.pathLength
    desc[sn + '_relPathLength'] = v.relPathLength
    desc[sn + '_relPathLength_2'] = v.relPathLength_2
    desc[sn + '_numSharedNeighbors'] = v.sharedNeighbors
    desc[sn + '_numRotBonds'] = v.numRotBonds
    desc[sn + '_numAroBonds'] = v.numAroBonds
  # some combination of substituent sizes
  desc['s34_size'] = desc['s3_size'] + desc['s4_size']
  desc['s34_phSize'] = desc['s3_phSize'] + desc['s4_phSize']
  desc['s34_relSize'] = desc['s3_relSize'] + desc['s4_relSize']
  desc['s34_phRelSize'] = desc['s3_phRelSize'] + desc['s4_phRelSize']
  # calculate the chiral moment --> kind of 3D descriptor
  desc['chiralMoment'] = calcSP3CarbonSubstituentMoment(
    [desc['s1_size'], desc['s2_size'], desc['s3_size'], desc['s4_size']])
  desc['chiralPhMoment'] = calcSP3CarbonSubstituentMoment(
    [desc['s1_phSize'], desc['s2_phSize'], desc['s3_phSize'], desc['s4_phSize']])
  return desc


def generateChiralDescriptorsForAllCenters(mol, verbose=False):
  """
    Generates descriptors for all chiral centers in the molecule.
    Details of these descriptors are described in: 
    Schneider et al., Chiral Cliffs: Investigating the Influence of Chirality on Binding Affinity
    https://doi.org/10.1002/cmdc.201700798. 
    >>> # test molecules are taken from the publication above (see Figure 3 and Figure 8)
    >>> testmols = {
    ...   "CHEMBL319180" : 'CCCN1C(=O)[C@@H](NC(=O)Nc2cccc(C)c2)N=C(N3CCN(C)CC3)c4ccccc14',
    ...   }
    >>> mol = Chem.MolFromSmiles(testmols['CHEMBL319180'])
    >>> desc = generateChiralDescriptorsForAllCenters(mol)
    >>> desc.keys()
    dict_keys([6])
    >>> desc[6]['arLevel2']
    0
    >>> desc[6]['s4_pathLength']
    7
    >>> desc[6]['maxDist']
    14
    >>> desc[6]['maxDistfromCC']
    7
    """

  desc = {}
  dists = Chem.GetDistanceMatrix(mol)
  for idxChiral, _ in Chem.FindMolChiralCenters(mol):
    desc[idxChiral] = calculateChiralDescriptors(mol, idxChiral, dists, verbose=False)
  return desc


def generateChiralDescriptors(mol, verbose=False):
  """
    Generates descriptors for the 'first' chiral centers in the molecule.
    Details of these descriptors are described in: 
    Schneider et al., Chiral Cliffs: Investigating the Influence of Chirality on Binding Affinity
    https://doi.org/10.1002/cmdc.201700798. 
    >>> # test molecules are taken from the publication above (see Figure 3 and Figure 8)
    >>> testmols = {
    ...   "CHEMBL319180" : 'CCCN1C(=O)[C@@H](NC(=O)Nc2cccc(C)c2)N=C(N3CCN(C)CC3)c4ccccc14',
    ...   "CHEMBL3350250" : 'CC(C)[C@@]1(CCc2ccc(O)cc2)CC(=O)C(=C(O)O1)Sc3cc(C)c(NS(=O)(=O)c4ccc(cn4)C(F)(F)F)cc3C(C)(C)C',
    ...   "CHEMBL3698720" : 'N[C@@H]1CCN(C1)c2cnc(Nc3ncc4c5ccncc5n(C6CCCC6)c4n3)cn2'
    ...   }
    >>> mol = Chem.MolFromSmiles(testmols['CHEMBL319180'])
    >>> desc = generateChiralDescriptors(mol)
    >>> desc['arLevel2']
    0
    >>> desc['s4_pathLength']
    7
    >>> desc['maxDist']
    14
    >>> desc['maxDistfromCC']
    7
    >>> mol = Chem.MolFromSmiles(testmols['CHEMBL3350250'])
    >>> desc = generateChiralDescriptors(mol)
    >>> desc['nLevel8']
    5
    >>> desc['s1_pathLength']
    2
    >>> desc['s2_size']
    9.0
    >>> desc['numRotBonds']
    9
    >>> mol = Chem.MolFromSmiles(testmols['CHEMBL3698720'])
    >>> desc = generateChiralDescriptors(mol)
    >>> desc['s3_numRotBonds']
    0
    >>> desc['s34_size']
    29.0
    >>> desc['phLevel2']
    1
    >>> desc['s2_numSharedNeighbors']
    0
    
    """

  dists = Chem.GetDistanceMatrix(mol)
  idxChiral = Chem.FindMolChiralCenters(mol)[0][0]
  return calculateChiralDescriptors(mol, idxChiral, dists, verbose=False)


#------------------------------------
#
#  doctest boilerplate
#
def _test():
  import doctest
  import sys
  return doctest.testmod(sys.modules["__main__"])


if __name__ == '__main__':
  import sys
  failed, tried = _test()
  sys.exit(failed)
