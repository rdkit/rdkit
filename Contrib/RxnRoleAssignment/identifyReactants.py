#
#  Copyright (c) 2016, Novartis Institutes for BioMedical Research Inc.
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
# Created by Nadine Schneider, July 2016

import itertools
from collections import Counter, defaultdict

import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem, rdqueries

from . import utils


class MoleculeDetails(object):

  __slots__ = [
    'detailFP', 'scaffoldFP', 'bitInfoDetailFP', 'bitInfoScaffoldFP', 'reactivity', 'bitReactivity',
    'molecule'
  ]

  def _atomDetailInvariant(self, mol):
    mol.UpdatePropertyCache(False)
    num_atoms = mol.GetNumAtoms()
    Chem.GetSSSR(mol)
    rinfo = mol.GetRingInfo()
    invariants = [0] * num_atoms
    for i, a in enumerate(mol.GetAtoms()):
      descriptors = []
      descriptors.append(a.GetAtomicNum())
      descriptors.append(a.GetTotalDegree())
      descriptors.append(a.GetTotalNumHs())
      descriptors.append(rinfo.IsAtomInRingOfSize(a.GetIdx(), 6))
      descriptors.append(rinfo.IsAtomInRingOfSize(a.GetIdx(), 5))
      descriptors.append(a.IsInRing())
      descriptors.append(a.GetIsAromatic())
      invariants[i] = hash(tuple(descriptors)) & 0xffffffff
    return invariants

  def _atomScaffoldInvariant(self, mol):
    num_atoms = mol.GetNumAtoms()
    invariants = [0] * num_atoms
    for i, a in enumerate(mol.GetAtoms()):
      descriptors = []
      descriptors.append(a.GetAtomicNum())
      invariants[i] = hash(tuple(descriptors)) & 0xffffffff
    return invariants

  def _createFP(self, mol, invariant, bitinfo, useBondTypes=True, radius=1):
    return AllChem.GetMorganFingerprint(mol=mol, radius=radius, invariants=invariant,
                                        useBondTypes=useBondTypes, bitInfo=bitinfo)

  def _isHeteroAtom(self, a):
    return a.GetAtomicNum() not in (6, 1)

  def _isSp3OrAromaticCarbon(self, a):
    if a.GetAtomicNum() != 6:
      return False
    if a.GetIsAromatic():
      return True
    for b in a.GetBonds():
      if b.GetBondTypeAsDouble() > 1.5:
        return False
    return True

  def _calcReactivityAtom(self, a):
    # exclude sp3 carbons or uncharged single heavy atoms such as water molecules
    if self._isSp3OrAromaticCarbon(a) or (len(a.GetNeighbors()) == 0 and a.GetFormalCharge() == 0):
      return 0
    # all other atoms have at least a reactivity of one
    reactivity = 1
    b = a.GetBonds()
    # if it is a heteroatom or has an H (we already know it's not SP3 or aromatic) increase the reactivity
    if self._isHeteroAtom(a) or a.GetTotalNumHs() > 0:
      reactivity += 1
    # slightly increase reactivity for atoms in aromatic rings compared to aliphatic rings
    if a.IsInRing():
      if a.GetIsAromatic():
        reactivity += 0.5
    # but prefer non-ring atoms
    else:
      reactivity += 1
    # increase reactivity of charged atoms
    if a.GetFormalCharge():
      reactivity += 2
    for bo in b:
      # look at the direct neighbors of the atom
      ni = bo.GetOtherAtom(a)
      # for non-single bonds increase the reactivity
      if bo.GetBondTypeAsDouble() > 1.5:
        reactivity += 1
        # if there are hydrogens attached, increase the reactivity
        if ni.GetTotalNumHs() > 0:
          reactivity += 1
      # if it is a bond to a hetero atom further increase the reactivity
      if self._isHeteroAtom(ni):
        reactivity += 1
        # bonds between nitrogens and oxygen or between oxygen and oxygen or between nitrogen and nitrogen are more reactive
        if a.GetAtomicNum() in (7, 8) and ni.GetAtomicNum() in (7, 8):
          reactivity += 2
        # if the neighbor is a Mg, Si, P, Pd, or Sn atom increase the reactivity
        elif ni.GetAtomicNum() in (12, 14, 15, 46, 50):
          reactivity += 1
    return reactivity

  def _calcReactivityMolecule(self, mol):
    reactivityAtoms = [self._calcReactivityAtom(a) for a in mol.GetAtoms()]
    return reactivityAtoms

  def __init__(self, molecule, verbose=0):
    self.molecule = molecule
    self.bitInfoDetailFP = {}
    self.detailFP = self._createFP(molecule, self._atomDetailInvariant(molecule),
                                   self.bitInfoDetailFP)
    self.bitInfoScaffoldFP = {}
    self.scaffoldFP = self._createFP(molecule, self._atomScaffoldInvariant(molecule),
                                     self.bitInfoScaffoldFP, useBondTypes=False)
    reactivityAtoms = self._calcReactivityMolecule(molecule)
    reactivity = sum(reactivityAtoms)
    if Chem.MolToSmiles(molecule) in frequentReagents:
      reactivity *= 0.8
    self.reactivity = reactivity


def _calcScore(reactantFP, productFP, bitInfoProd=None, output=False):
  if output:
    print("--- _calcScore ---")
  score = 0
  dFP = productFP - reactantFP
  numRBits = float(utils.getNumPositiveCounts(reactantFP))
  if output > 2:
    print("num RBits: ", numRBits)
  numPBits = float(utils.getNumPositiveCounts(productFP))
  if output > 2:
    print("num PBits: ", numPBits)
  numUnmappedPBits = float(utils.getNumPositiveCounts(dFP))
  if output > 2:
    print("num UnmappedPBits: ", numUnmappedPBits)
  numUnmappedRBits = float(utils.getNumNegativeCounts(dFP))
  if output > 2:
    print("num UnmappedRBits: ", numUnmappedRBits)

  numUnmappedPAtoms = -1
  bitsUnmappedPAtoms = -1
  if bitInfoProd is not None:
    numUnmappedPAtoms, bitsUnmappedPAtoms = utils.getNumPositiveBitCountsOfRadius0(dFP, bitInfoProd)
    if output > 2:
      print("num UnmappedPAtoms: ", numUnmappedPAtoms)
  ratioMappedPBits = 1 - (numUnmappedPBits / numPBits)
  ratioUnmappedRBits = numUnmappedRBits / numRBits
  score = max(ratioMappedPBits - ratioUnmappedRBits * ratioUnmappedRBits, 0)

  if output > 1:
    print("score: ", score, "(", ratioMappedPBits, ",", ratioUnmappedRBits * ratioUnmappedRBits,
          ",", ratioUnmappedRBits, ")")

  return [score, numUnmappedPBits, numUnmappedPAtoms, bitsUnmappedPAtoms]


# Set of frequent reagents derived from all patent reactions
frequentReagents = set([
  'CCN(CC)CC', '[Li+]', '[Na+]', 'O=C(O)CC(O)(CC(=O)O)C(=O)O', 'O=S(=O)(O)O', 'CN1CCCC1=O',
  'CCN(C(C)C)C(C)C', 'c1ccncc1', '[K]', 'CC(C)(C)O', 'CCO', 'Cc1ccc(S(=O)(=O)O)cc1',
  'ClC(Cl)(Cl)Cl', '[Na]', 'CC(C)(C)[O-]', 'O=C([O-])O', 'COCCOC', '[NH4+]', 'CC(C)OC(C)C',
  'O=C([O-])[O-]', 'CC(=O)OC(C)=O', 'O=C=O', '[Cl-]', 'c1ccc(P(c2ccccc2)c2ccccc2)cc1', '[H-]',
  'N#N', 'CN1CCOCC1', 'C1COCCO1', 'c1ccccc1', '[Cs+]', '[K+]', '[OH-]', 'CCCCCC', 'CCCCC',
  'CN(C)C=O', 'C[O-]', 'Cc1ccccc1', 'C1CCC2=NCCCN2CC1', 'CO', 'CCCCO', 'O=C(O)C(F)(F)F',
  'O=P([O-])([O-])[O-]', 'CCOC(C)=O', '[Mg+2]', 'C1CCCCC1', 'O', 'N', 'II', 'O=CO', 'CC(=O)N(C)C',
  'CC(=O)O', 'CCOCC', 'CC(C)O', 'C[Si](C)(C)Cl', 'Cc1ccccc1C', 'CC(C)=O', 'CS(=O)(=O)O',
  'CN(C)c1ccncc1', 'Cl', 'ClCCCl', 'O=S(Cl)Cl', 'ClC(Cl)Cl', '[Li]CCCC', '[Pd]', '[H][H]', '[Br-]',
  'CS(C)=O', 'COC(C)(C)C', 'O=S(=O)([O-])[O-]', 'CC(Cl)Cl', 'CC(=O)[O-]',
  'CCCC[N+](CCCC)(CCCC)CCCC', 'ClCCl', 'CC#N', 'C1CCOC1', 'CCCCCCC'
])


def _getBestCombination(rfps, pfps, output=False):

  if output:
    print("--- _getBestCombination ---")

  tests = []
  numReactants = len(rfps)
  # generate first all reactant combinations
  for i in range(1, numReactants + 1):
    for x in itertools.combinations(range(numReactants), i):
      temp = []
      for j in x:
        # don't include frequent reagents
        if not rfps[j][1]:
          numAtms = rfps[j][0].molecule.GetNumAtoms()
          # not test single ions
          if numAtms > 1:
            # store the number of reactant atoms for later
            temp.append((rfps[j][0].molecule.GetNumAtoms(), j))
        else:
          if output > 3:
            print("Frequent reagent found: ", j)
      if temp not in tests:
        tests.append(temp)
  # initialisation of the results
  maxScore = 0
  maxDetailScore = 0
  finalReacts = [[]]
  # get the product fingerprints
  productsDetailFP = utils.getSumFps([i.detailFP for i in pfps])
  productsScaffoldFP = utils.getSumFps([i.scaffoldFP for i in pfps])
  # get the number of atoms for the product
  numProductAtoms = 0
  for i in pfps:
    numProductAtoms += i.molecule.GetNumAtoms()
  # get the bitinfo for the product FP
  productsDetailFPBitInfo = {}
  productsScaffoldFPBitInfo = {}
  for i in pfps:
    productsDetailFPBitInfo.update(i.bitInfoDetailFP)
    productsScaffoldFPBitInfo.update(i.bitInfoScaffoldFP)
  # set some initial values
  numUnmappedPAtoms, bitsUnmappedPAtoms = utils.getNumPositiveBitCountsOfRadius0(
    productsScaffoldFP, productsScaffoldFPBitInfo)
  finalNumUnmappedProdAtoms = [[
    len(productsDetailFP.GetNonzeroElements()),
    len(productsScaffoldFP.GetNonzeroElements()), numUnmappedPAtoms, bitsUnmappedPAtoms
  ]]

  for test in tests:
    if len(test) < 1:
      continue
    # get the number of involved reactant atoms
    numReactantAtoms = np.array(test)[:, 0].sum()
    # ignore combinations including too many or too few atoms
    if numReactantAtoms > 5 * numProductAtoms or numReactantAtoms < numProductAtoms * 0.8:
      continue

    if output > 0:
      print("Combination: ", test)

    #build the combined reactant FPs
    reactantsDetailFP = utils.getSumFps([rfps[i[1]][0].detailFP for i in test])
    reactantsScaffoldFP = utils.getSumFps([rfps[i[1]][0].scaffoldFP for i in test])

    # get the scores for both FPs
    detailFPScore = _calcScore(reactantsDetailFP, productsDetailFP,
                               bitInfoProd=productsDetailFPBitInfo, output=output)
    scaffoldFPScore = _calcScore(reactantsScaffoldFP, productsScaffoldFP,
                                 bitInfoProd=productsScaffoldFPBitInfo, output=output)
    # final score
    score = detailFPScore[0] + scaffoldFPScore[0]

    if output > 0:
      print(">>>> score: ", score)
      print(">>>> scores (detail, scaffold): ", detailFPScore[0], scaffoldFPScore[0])
      print(">>>> num unmapped productFP bits: ", detailFPScore[1], scaffoldFPScore[1],
            detailFPScore[2], scaffoldFPScore[2])

    if score > maxScore:
      maxScore = score
      maxDetailScore = detailFPScore[0]
      del finalReacts[:]
      del finalNumUnmappedProdAtoms[:]
      # set the final reactants
      finalReacts.append([i[1] for i in test])
      # for tracking the mapping of the product atoms include the number of unmapped detailedFP bits, the number of unmapped
      # atoms based on the scaffold FP,  the number of unmapped scaffoldFP bits, and the unmapped scaffoldFP bits
      finalNumUnmappedProdAtoms.append(
        [detailFPScore[1], scaffoldFPScore[2], scaffoldFPScore[1], scaffoldFPScore[-1]])
      if output > 0:
        print(" >> maxScore: ", maxScore)
        print(" >> Final reactants: ", finalReacts)
      # test for almost perfect matchings (e.g. oxidations, reduction etc.)
      if scaffoldFPScore[0] > 0.9999 and detailFPScore[0] > 0.8:
        return finalReacts, finalNumUnmappedProdAtoms
      # test for number of mapped product atoms e.g. to capture deprotections earlier
      if len(finalNumUnmappedProdAtoms) > 0 and len(test) == 1:
        if finalNumUnmappedProdAtoms[0][1] == 0 and finalNumUnmappedProdAtoms[0][0] <= 3:
          return finalReacts, finalNumUnmappedProdAtoms
    # include alternative solutions
    elif abs(score - maxScore) < 0.0000001 and score > 0.0:
      finalReacts.append([i[1] for i in test])
      finalNumUnmappedProdAtoms.append(
        [detailFPScore[1], scaffoldFPScore[2], scaffoldFPScore[1], scaffoldFPScore[-1]])
      if output > 0:
        print(" >> Added alternative result")
        print(" >> Final reactants: ", finalReacts)

  return finalReacts, finalNumUnmappedProdAtoms


def _findMissingReactiveReactants(rfps, pfps, currentReactants, unmappedPAtoms, output=False):
  if output:
    print("--- _findMissingReactiveReactants ---")
  if not len(unmappedPAtoms):
    return currentReactants
  # if there are unmapped product bits find possible reactants for those
  else:
    finalReactants = []
    numReactants = len(rfps)
    # investigate all possible solutions of the scoring before
    for reacts, umPA in zip(currentReactants, unmappedPAtoms):
      # if there are unmapped product atoms find possible reactants for those
      finalReactants.append(reacts)
      if umPA[1] > 0:
        remainingReactants = set(range(numReactants)).difference(set(reacts))
        # sort the possible reactants by the reactivity
        remainingReactants = sorted(
          remainingReactants,
          key=lambda x: rfps[x].reactivity / float(rfps[x].molecule.GetNumAtoms()), reverse=True)
        missingPAtoms = []
        # get the missing atoms and counts
        for bit, c in umPA[-1]:
          for pbi in range(len(pfps)):
            if bit in pfps[pbi].bitInfoScaffoldFP:
              a = pfps[pbi].bitInfoScaffoldFP[bit][0]
              missingPAtoms.extend([pfps[pbi].molecule.GetAtomWithIdx(a[0]).GetAtomicNum()] * c)
        missingPAtoms = Counter(missingPAtoms)
        if output > 0:
          print(missingPAtoms)
        # build queries for the missing atoms
        queries = [(rdqueries.AtomNumEqualsQueryAtom(a), a) for a in missingPAtoms]
        maxFullfilledQueries = 0
        maxReactivity = -1
        addReactants = []
        # search for the most reactive reactants capturing all/most of the unmapped product atoms
        for r in remainingReactants:
          if output > 0:
            print(" >> Reactant", r, rfps[r].reactivity / float(rfps[r].molecule.GetNumAtoms()))
          countFullfilledQueries = 0
          for q, a in queries:
            if len(rfps[r].molecule.GetAtomsMatchingQuery(q)) >= missingPAtoms[a]:
              countFullfilledQueries += 1
          if output > 0:
            print(" Max reactivity", maxReactivity)
            print(" Max fulfilled queries", maxFullfilledQueries)
          if countFullfilledQueries > maxFullfilledQueries:
            maxFullfilledQueries = countFullfilledQueries
            maxReactivity = rfps[r].reactivity / float(rfps[r].molecule.GetNumAtoms())
            addReactants = [r]
          elif maxFullfilledQueries and countFullfilledQueries == maxFullfilledQueries and \
               rfps[r].reactivity/float(rfps[r].molecule.GetNumAtoms()) >= maxReactivity:
            maxFullfilledQueries = countFullfilledQueries
            addReactants.append(r)
          if output > 0:
            print(" Added reactants", addReactants)
        finalReactants[-1].extend(addReactants)
  if output > 0:
    print(" >> Final reactants", finalReactants)
  return finalReactants


def _detectObviousReagents(reactants, products):
  unchangedReacts = set()
  unchangedProds = set()
  for i, r in enumerate(reactants):
    for j, p in enumerate(products):
      if r == p:
        unchangedReacts.add(i)
        unchangedProds.add(j)
  return unchangedReacts, unchangedProds


def identifyReactants(reaction, output=False):
  rxn = AllChem.ChemicalReaction(reaction)
  AllChem.RemoveMappingNumbersFromReactions(rxn)
  if output:
    print("--- identifyReactants ---")
  reactants = rxn.GetReactants()
  products = rxn.GetProducts()
  ### Preprocessing
  uniqueReactants, reactantSmiles = utils.uniqueMolecules(reactants)
  uniqueProducts, productSmiles = utils.uniqueMolecules(products)
  # find molecules which do not change in the rxn
  unmodifiedReactants, unmodifiedProducts = _detectObviousReagents(reactantSmiles, productSmiles)
  if output:
    print("  >>> Found reagents in reactants:", unmodifiedReactants)
    print("  >>> Found reagents in products:", unmodifiedProducts)
  if len(products) == len(unmodifiedProducts):
    unmodifiedProducts = set()
  uniquePotentialReactants = [r for r in sorted(set(uniqueReactants.values()))]
  uniquePotentialProducts = [
    p for p in sorted(set(uniqueProducts.values())) if p not in unmodifiedProducts
  ]

  ### Find the most probable reactants
  # only generate moleculeDetail objects for unique, potential reactants and products
  rfps = [MoleculeDetails(reactants[r]) for r in uniquePotentialReactants]
  pfps = [MoleculeDetails(products[p]) for p in uniquePotentialProducts]

  rfpsPrep = [(MoleculeDetails(reactants[r]), reactantSmiles[r] in frequentReagents)
              for r in uniquePotentialReactants]

  reacts, unmappedProdAtoms = _getBestCombination(rfpsPrep, pfps, output=output)
  # no reactants where found try again including the frequent reagents
  if np.array(reacts).shape == (1, 0):
    rfpsPrep = [(MoleculeDetails(reactants[r]), 0) for r in uniquePotentialReactants]
    reacts, unmappedProdAtoms = _getBestCombination(rfpsPrep, pfps, output=output)

  ### Postprocessing
  # identify missing reactants
  reacts = _findMissingReactiveReactants(rfps, pfps, reacts, unmappedProdAtoms, output=output)
  finalreacts = []
  for i in reacts:
    temp = [uniquePotentialReactants[j] for j in i]
    finalreacts.append(set(temp))

  return finalreacts, unmodifiedReactants, unmodifiedProducts


# reassign the reaction roles of a reaction
def reassignRXNRoles(rxn):
  utils.transferAgentsToReactants(rxn)
  reacts, rAgents, pAgents = identifyReactants(rxn)
  if len(reacts) < 1:
    return None
  new_rxn = AllChem.ChemicalReaction()
  for i in range(rxn.GetNumProductTemplates()):
    new_rxn.AddProductTemplate(rxn.GetProductTemplate(i))
  for i in range(rxn.GetNumReactantTemplates()):
    if i in reacts[0]:
      new_rxn.AddReactantTemplate(rxn.GetReactantTemplate(i))
    else:
      new_rxn.AddAgentTemplate(rxn.GetReactantTemplate(i))
  return new_rxn


# clean-up the reaction smiles
def reassignReactionRoles(smi):
  rxn = AllChem.ReactionFromSmarts(smi, useSmiles=True)
  new_rxn = reassignRXNRoles(rxn)
  if new_rxn is None:
    return ''
  smi_new = AllChem.ReactionToSmiles(new_rxn)
  return smi_new
