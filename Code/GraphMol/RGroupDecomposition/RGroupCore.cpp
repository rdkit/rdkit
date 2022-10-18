//
//  Copyright (c) 2017-2021, Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "RGroupCore.h"
#include <GraphMol/Substruct/SubstructUtils.h>

namespace RDKit {

namespace {
// From answer 12 in
// https://stackoverflow.com/questions/5279051/how-can-i-create-cartesian-product-of-vector-of-vectors
// by anumi
static std::vector<std::vector<int>> cartesianProduct(
    const std::vector<std::vector<int>> &v) {
  std::vector<std::vector<int>> s = {{}};
  for (const auto &u : v) {
    std::vector<std::vector<int>> r;
    for (const auto &x : s) {
      for (const auto y : u) {
        r.push_back(x);
        r.back().push_back(y);
      }
    }
    s = std::move(r);
  }
  return s;
}
}  // namespace

// move this to constructor if the create new core path can be removed from
// RGroupDecomposition::add
void RCore::init() {
  findIndicesWithRLabel();
  countUserRGroups();
  buildMatchingMol();
}

void RCore::findIndicesWithRLabel() {
  // Find all the core atoms that have user
  // label and set their indices to 1 in core_atoms_with_user_labels
  core_atoms_with_user_labels.resize(core->getNumAtoms());
  for (const auto atom : core->atoms()) {
    int label;
    if (atom->getPropIfPresent(RLABEL, label) && label > 0) {
      core_atoms_with_user_labels.set(atom->getIdx());
    }
  }
}

// Return a copy of core where dummy atoms are replaced by
// the respective matching atom in mol, while other atoms have
// their aromatic flag and formal charge copied from
// the respective matching atom in mol
ROMOL_SPTR RCore::replaceCoreAtomsWithMolMatches(
    bool &hasCoreDummies, const ROMol &mol, const MatchVectType &match) const {
  auto coreReplacedAtoms = boost::make_shared<RWMol>(*core);
  hasCoreDummies = false;
  for (const auto &p : match) {
    auto atom = coreReplacedAtoms->getAtomWithIdx(p.first);
    if (atom->getAtomicNum() == 0) {
      hasCoreDummies = true;
    }
    if (isAtomWithMultipleNeighborsOrNotDummyRGroupAttachment(*atom)) {
      auto molAtom = mol.getAtomWithIdx(p.second);
      replaceCoreAtom(*coreReplacedAtoms, *atom, *molAtom);
    }
  }

  // look for dummy rgroup attachment points in core that are not in mapping
  // there is no r group for these guys.  I've marked them as missing.
  // It might be better to delete them coreReplaceAtoms, but I'm not sure of the
  // implications of that.
  for (auto atom : coreReplacedAtoms->atoms()) {
    if (atom->getAtomicNum() == 0 && atom->hasProp(RLABEL)) {
      auto missing = std::find_if(match.cbegin(), match.cend(),
                                  [atom](const std::pair<int, int> &p) {
                                    return static_cast<unsigned int>(p.first) ==
                                           atom->getIdx();
                                  }) == match.end();
      if (missing) {
        atom->setProp<bool>(MISSING_RGROUP, true);
      }
    }
  }

  std::map<int, int> matchLookup(match.cbegin(), match.cend());
  for (auto bond : coreReplacedAtoms->bonds()) {
    if (bond->hasQuery()) {
      hasCoreDummies = true;
      const auto molBond =
          mol.getBondBetweenAtoms(matchLookup[bond->getBeginAtomIdx()],
                                  matchLookup[bond->getEndAtomIdx()]);
      if (molBond == nullptr) {
        // this can happen if we have a user-defined R group that is not
        // matched in the query
        CHECK_INVARIANT(bond->getBeginAtom()->getAtomicNum() == 0 ||
                            bond->getEndAtom()->getAtomicNum() == 0,
                        "Failed to find core bond in molecule");
      } else {
        Bond newBond(molBond->getBondType());
        newBond.setIsAromatic(molBond->getIsAromatic());
        coreReplacedAtoms->replaceBond(bond->getIdx(), &newBond, true);
      }
    }
  }

#ifdef VERBOSE
  std::cerr << "Original core smarts  " << MolToSmarts(*core) << std::endl;
  std::cerr << "Dummy replaced core smarts  " << MolToSmarts(*coreReplacedAtoms)
            << std::endl;
#endif

  if (mol.getNumConformers() > 0) {
    // if the input structure has coordinates copy them to the core
    if (!coreReplacedAtoms->getNumConformers()) {
      coreReplacedAtoms->addConformer(
          new Conformer(coreReplacedAtoms->getNumAtoms()));
    }
    auto &replacedConformer = coreReplacedAtoms->getConformer();
    const auto &molConformer = mol.getConformer();

    for (const auto &p : match) {
      auto molPoint = molConformer.getAtomPos(p.second);
      replacedConformer.setAtomPos(p.first, molPoint);
    }
  } else {
    // otherwise, delete all core coordinates from the replaced core
    coreReplacedAtoms->clearConformers();
  }

  return coreReplacedAtoms;
}

void RCore::replaceCoreAtom(RWMol &mol, Atom &atom, const Atom &other) const {
  auto atomicNumber = other.getAtomicNum();
  auto targetAtom = &atom;
  bool wasDummy = (atom.getAtomicNum() == 0);
  if (wasDummy || atom.hasQuery()) {
    if (atom.hasQuery()) {
      Atom newAtom(atomicNumber);
      auto atomIdx = atom.getIdx();
      mol.replaceAtom(atomIdx, &newAtom, false, true);
      targetAtom = mol.getAtomWithIdx(atomIdx);
    } else {
      atom.setAtomicNum(atomicNumber);
    }
  }
  targetAtom->setIsAromatic(other.getIsAromatic());
  targetAtom->setFormalCharge(other.getFormalCharge());
  if (wasDummy) {
    targetAtom->setNoImplicit(true);
    unsigned int numHs = 0;
    const auto &otherMol = other.getOwningMol();
    for (const auto &nbri :
         boost::make_iterator_range(otherMol.getAtomNeighbors(&other))) {
      const auto nbrAtom = otherMol[nbri];
      if (nbrAtom->getAtomicNum() == 1) {
        ++numHs;
      }
    }
    targetAtom->setNumExplicitHs(numHs + other.getTotalNumHs());
    targetAtom->updatePropertyCache(false);
  }
}

// Final core returned to user with dummy atoms and bonds set to those in the
// match
RWMOL_SPTR RCore::coreWithMatches(const ROMol &coreReplacedAtoms) const {
  auto finalCore = boost::make_shared<RWMol>(*labelledCore);
  for (size_t atomIdx = 0; atomIdx < coreReplacedAtoms.getNumAtoms();
       ++atomIdx) {
    auto coreAtom = finalCore->getAtomWithIdx(atomIdx);
    auto templateAtom = coreReplacedAtoms.getAtomWithIdx(atomIdx);
    auto unlabelledCoreAtom = core->getAtomWithIdx(atomIdx);
    if (templateAtom->getAtomicNum() > 0 &&
        isAtomWithMultipleNeighborsOrNotDummyRGroupAttachment(
            *unlabelledCoreAtom)) {
      replaceCoreAtom(*finalCore, *coreAtom, *templateAtom);
    }
  }

  for (size_t bondIdx = 0; bondIdx < coreReplacedAtoms.getNumBonds();
       ++bondIdx) {
    auto coreBond = finalCore->getBondWithIdx(bondIdx);
    if (coreBond->hasQuery()) {
      auto templateBond = coreReplacedAtoms.getBondWithIdx(bondIdx);
      Bond newBond(templateBond->getBondType());
      newBond.setIsAromatic(templateBond->getIsAromatic());
      finalCore->replaceBond(bondIdx, &newBond, true);
    }
  }

  // update core coordinates to input structures
  if (coreReplacedAtoms.getNumConformers() && finalCore->getNumConformers()) {
    const auto &replacedConformer = coreReplacedAtoms.getConformer();
    auto &finalConformer = finalCore->getConformer();

    size_t atomIdx = 0;
    for (; atomIdx < coreReplacedAtoms.getNumAtoms(); ++atomIdx) {
      auto molPoint = replacedConformer.getAtomPos(atomIdx);
      finalConformer.setAtomPos(atomIdx, molPoint);
    }
    // Dummy atom coordinates are not in input structure, so calculate them here
    for (; atomIdx < finalCore->getNumAtoms(); atomIdx++) {
      const auto atom = finalCore->getAtomWithIdx(atomIdx);
      const int neighborIdx = *finalCore->getAtomNeighbors(atom).first;
      MolOps::setTerminalAtomCoords(*finalCore, atomIdx, neighborIdx);
    }
  }

  // Remove unmapped dummies
  std::vector<Atom *> atomsToRemove;
  for (const auto atom : coreReplacedAtoms.atoms()) {
    if (atom->getAtomicNum() == 0 && atom->hasProp(MISSING_RGROUP)) {
      atomsToRemove.push_back(finalCore->getAtomWithIdx(atom->getIdx()));
    }
  }
  finalCore->beginBatchEdit();
  for (auto atom : atomsToRemove) {
    finalCore->removeAtom(atom);
  }
  finalCore->commitBatchEdit();

  finalCore->updatePropertyCache(false);
  return finalCore;
}

// matching the core to the target molecule is a two step process
// First match to a reduced representation (the core minus terminal user
// R-groups). Next, match the R-groups. We do this as the core may not be a
// substructure match for the molecule if a single molecule atom matches 2
// terminal user defined RGroup attachments (see
// https://github.com/rdkit/rdkit/pull/4002) buildMatchingMol() creates the
// reduced representation from the core and matchTerminalUserRGroups() adds in
// the terminal R Groups to the match

// Builds a matching molecule which is the core with terminal user R groups
// removed Also creates the data structures used in matching the R groups
void RCore::buildMatchingMol() {
  matchingMol = boost::make_shared<RWMol>(*core);
  terminalRGroupDummyAtoms.clear();
  terminalRGroupAtomToNeighbor.clear();
  RWMol::ATOM_PTR_VECT atomsToRemove;
  for (auto atom : matchingMol->atoms()) {
    // keep track of the original core index in the matching molecule atom
    atom->setProp<int>(RLABEL_CORE_INDEX, atom->getIdx());
    if (atom->getAtomicNum() == 0 && atom->getDegree() == 1 &&
        isDummyRGroupAttachment(*atom)) {
      // remove terminal user R groups and save core atom index and mapping to
      // heavy neighbor
      atomsToRemove.push_back(atom);
      terminalRGroupDummyAtoms.insert(atom->getIdx());
      const int neighborIdx = *matchingMol->getAtomNeighbors(atom).first;
      terminalRGroupAtomToNeighbor.emplace(atom->getIdx(), neighborIdx);
    }
  }

  for (auto atom : atomsToRemove) {
    matchingMol->removeAtom(atom);
  }
}

// Given a matching molecule substructure match to a target molecule, return
// core matches with terminal user R groups matched
std::vector<MatchVectType> RCore::matchTerminalUserRGroups(
    const RWMol &target, MatchVectType match,
    const SubstructMatchParameters &sssParams) const {
  // Transform match indexed by matching molecule atoms to a map
  // indexed by core atoms
  std::transform(match.begin(), match.end(), match.begin(),
                 [this](const std::pair<int, int> &mapping) {
                   auto queryIdx =
                       this->matchingIndexToCoreIndex(mapping.first);
                   std::pair<int, int> newMapping(queryIdx, mapping.second);
                   return newMapping;
                 });
  std::map<int, int> matchMap(match.cbegin(), match.cend());

  std::vector<MatchVectType> allMappings;
  if (terminalRGroupDummyAtoms.empty()) {
    allMappings.push_back(match);
    return allMappings;
  }

  // build a set of target atoms currently mapped
  std::set<int> mappedTargetIdx;
  std::transform(
      match.cbegin(), match.cend(),
      std::inserter(mappedTargetIdx, mappedTargetIdx.begin()),
      [](const std::pair<int, int> &mapping) { return mapping.second; });

  // Dummy atoms/r group attachments that cannot be mapped to target atoms
  std::vector<int> missingDummies;
  // A map of terminal dummies/R group attachment points to a list of possible
  // target atoms that the R group can map to
  std::map<int, std::vector<int>> availableMappingsForDummyMap;

  // Loop over all the user terminal R groups and see if they can be included
  // in the match and, if so, record the target atoms that the R group can
  // be mapped to.
  for (const auto dummyIdx : terminalRGroupDummyAtoms) {
    // find the heavy atom in the core attached to the dummy
    const int neighborIdx = terminalRGroupAtomToNeighbor.find(dummyIdx)->second;
    const auto coreBond = core->getBondBetweenAtoms(dummyIdx, neighborIdx);
    // find the atom in the target mapped to the neighbor in the core
    const int targetIdx = matchMap[neighborIdx];
    const auto targetAtom = target.getAtomWithIdx(targetIdx);
    ROMol::ADJ_ITER nbrIter, endNbrs;
    std::vector<int> available;
    // now look for neighbors of that target atom that are not mapped to a core
    // atom- the dummy atom can potentially be mapped to each of those
    boost::tie(nbrIter, endNbrs) = target.getAtomNeighbors(targetAtom);
    while (nbrIter != endNbrs) {
      if (mappedTargetIdx.find(*nbrIter) == mappedTargetIdx.end()) {
        const auto targetBond = target.getBondBetweenAtoms(targetIdx, *nbrIter);
        // check for bond compatibility
        if (bondCompat(coreBond, targetBond, sssParams)) {
          available.push_back(*nbrIter);
        }
      }
      ++nbrIter;
    }

    if (available.size()) {
      availableMappingsForDummyMap[dummyIdx] = available;
    } else {
      // We can't map this dummy to a target atom.  That is OK if it is an
      // unlabelled core attachment atom or a user R label connected to a
      // query or wildcard atom
      const auto dummy = core->getAtomWithIdx(dummyIdx);
      const auto neighbor = core->getAtomWithIdx(neighborIdx);
      if (dummy->hasProp(UNLABELLED_CORE_ATTACHMENT)) {
        missingDummies.push_back(dummyIdx);
      } else if (isUserRLabel(*dummy) &&
                 (neighbor->getAtomicNum() == 0 || neighbor->hasQuery())) {
        // https://github.com/rdkit/rdkit/issues/4505
        missingDummies.push_back(dummyIdx);
      } else {
        return allMappings;
      }
    }
  }
  if (availableMappingsForDummyMap.empty()) {
    allMappings.push_back(match);
    return allMappings;
  }

  // find singleton duplicate mappings. This can happen if there are 2 user
  // specified R groups attached to a query atom and only one of the R groups
  // can match (see second case in https://github.com/rdkit/rdkit/issues/4505)
  // We remove the R group with the highest label
  for (auto mapping = availableMappingsForDummyMap.begin();
       mapping != availableMappingsForDummyMap.end();) {
    auto values = mapping->second;
    bool erased = false;
    if (values.size() == 1) {
      const auto dummy = core->getAtomWithIdx(mapping->first);
      auto [nbrIter, endNbrs] = core->getAtomNeighbors(dummy);
      auto heavyNeighbor = core->getAtomWithIdx(*nbrIter);
      bool userQueryDummy =
          isUserRLabel(*dummy) &&
          (heavyNeighbor->getAtomicNum() == 0 || heavyNeighbor->hasQuery());
      if (userQueryDummy) {
        boost::tie(nbrIter, endNbrs) = core->getAtomNeighbors(heavyNeighbor);
        while (nbrIter != endNbrs) {
          int otherIdx = *nbrIter;
          if (otherIdx != mapping->first) {
            auto mappingIter = availableMappingsForDummyMap.find(otherIdx);
            if (mappingIter != availableMappingsForDummyMap.end() &&
                mappingIter->second.size() == 1 &&
                mappingIter->second[0] == values[0]) {
              auto otherDummy = core->getAtomWithIdx(otherIdx);
              bool remove = isUserRLabel(*otherDummy) &&
                            dummy->getProp<int>(RLABEL) >
                                otherDummy->getProp<int>(RLABEL);
              if (remove) {
                missingDummies.push_back(mapping->first);
                availableMappingsForDummyMap.erase(mapping++);
                erased = true;
                break;
              }
            }
          }
          nbrIter++;
        }
      }
    }
    if (!erased) {
      ++mapping;
    }
  }

  std::vector<int> dummiesWithMapping;
  std::vector<std::vector<int>> availableMappingsForDummy;
  for (const auto &mapping : availableMappingsForDummyMap) {
    dummiesWithMapping.push_back(mapping.first);
    availableMappingsForDummy.push_back(mapping.second);
  }

  // enumerate over all available atoms using a cartesian product.
  // this is not the most ideal way to do things as the straightforward product
  // will allow duplicates.  In the unlikely evert that there is a performance
  // issue a custom enumerator will be needed
  const auto allAvailableMappings = cartesianProduct(availableMappingsForDummy);
  // the size of the final mapping
  size_t size = allAvailableMappings[0].size() + match.size();
  // these indices are needed for the whole molecule match check functor

  std::unique_ptr<RWMol> checkCore;
  std::map<size_t, size_t> coreToCheck;
  const std::string indexProp("__core_index__");
  bool hasMissing = !missingDummies.empty();
  if (hasMissing) {
    // if there are dummies that we can map these need to be removed from the
    // query before atom-by-atom matching.  Create a copy of the query for that
    // and use properties to map atoms back to the core
    for (auto atom : core->atoms()) {
      atom->setProp(indexProp, atom->getIdx());
    }
    checkCore = std::make_unique<RWMol>(*core);
    std::sort(missingDummies.begin(), missingDummies.end(),
              std::greater<int>());
    for (int index : missingDummies) {
      auto [nbrIdx, endNbrs] =
          checkCore->getAtomNeighbors(checkCore->getAtomWithIdx(index));
      auto neighborAtom = checkCore->getAtomWithIdx(*nbrIdx);
      checkCore->removeAtom(index);
      neighborAtom->updatePropertyCache(false);
    }
    size_t index = 0U;
    for (const auto atom : checkCore->atoms()) {
      auto coreIndex = atom->getProp<int>(indexProp);
      coreToCheck[coreIndex] = index++;
    }
    for (auto atom : core->atoms()) {
      atom->clearProp(indexProp);
    }
  }

  auto queryIndices = new std::uint32_t[size];
  auto targetIndices = new std::uint32_t[size];
  for (size_t position = 0; position < match.size(); position++) {
    const auto &pair = match[position];
    auto queryIndex = hasMissing ? coreToCheck[pair.first] : pair.first;
    queryIndices[position] = queryIndex;
    targetIndices[position] = pair.second;
  }

  auto queryMatchingMol = hasMissing ? checkCore.get() : core.get();
  MolMatchFinalCheckFunctor molMatchFunctor(*queryMatchingMol, target,
                                            sssParams);
  boost::dynamic_bitset<> targetBondsPresent(target.getNumBonds());

  // Filter all available mappings removing those that contain duplicates,
  // or violate chirality
  for (const auto &dummyMapping : allAvailableMappings) {
    CHECK_INVARIANT(match.size() + dummyMapping.size() == size,
                    "Size error in dummy mapping");
    auto duplicateBonds = false;
    targetBondsPresent.reset();
    for (size_t i = 0; i < dummyMapping.size(); i++) {
      size_t position = match.size() + i;
      auto queryIndex = hasMissing ? coreToCheck[dummiesWithMapping[i]]
                                   : dummiesWithMapping[i];
      queryIndices[position] = queryIndex;
      targetIndices[position] = dummyMapping[i];
      const int neighborIdx =
          terminalRGroupAtomToNeighbor.find(dummiesWithMapping[i])->second;
      const int targetNeighborIdx = matchMap[neighborIdx];
      const auto targetBond =
          target.getBondBetweenAtoms(dummyMapping[i], targetNeighborIdx);
      CHECK_INVARIANT(targetBond != nullptr, "Matching target bond not found");
      const auto targetBondIdx = targetBond->getIdx();
      // check for duplicates- some of this could also be handled within a
      // modified cartesian product.
      if (targetBondsPresent[targetBondIdx]) {
        duplicateBonds = true;
        break;
      }
      targetBondsPresent[targetBondIdx] = 1;
    }
    // use MolMatchFinalCheckFunctor to check this match works with chirality
    if (!duplicateBonds && molMatchFunctor(queryIndices, targetIndices)) {
      MatchVectType matchWithDummy(match);
      for (size_t i = 0; i < dummyMapping.size(); i++) {
        matchWithDummy.emplace_back(dummiesWithMapping[i], dummyMapping[i]);
      }
      allMappings.push_back(matchWithDummy);
    }
  }

  delete[] queryIndices;
  delete[] targetIndices;
  return allMappings;
}

// This function checks the bond environment is valid when onlyMatchAtRGroups
// is set.  When this function is called attachmentIdx is the index of an atom
// in the target that is mapped to an R group.  Validates that all core bonds to
// the attachment point are present - in certain circumstances there may be two
// core bonds to a target attachment point and when onlyMatchAtRGroups is set
// both bonds should be present.
bool RCore::checkAllBondsToAttachmentPointPresent(
    const ROMol &mol, const int attachmentIdx,
    const MatchVectType &mapping) const {
  const auto atom = mol.getAtomWithIdx(attachmentIdx);
  std::set<int> coreNeighbors;
  for (const auto &nbri :
       boost::make_iterator_range(mol.getAtomNeighbors(atom))) {
    const auto &nbr = mol[nbri];
    // could a neighbor to an r group attachment match another r group
    // attachment?  I don't think so.
    if (nbr->getAtomicNum() >= 1) {
      const auto match = std::find_if(
          mapping.cbegin(), mapping.cend(), [nbri](std::pair<int, int> p) {
            return p.second == static_cast<int>(nbri);
          });
      if (match != mapping.end()) {
        auto coreAtom = core->getAtomWithIdx(match->first);
        // don't need to match a non terminal user R group
        // if (!(coreAtom->getDegree() > 1 && isUserRLabel(*coreAtom))) {
        if (!(coreAtom->getAtomicNum() == 0 && isUserRLabel(*coreAtom))) {
          coreNeighbors.insert(match->first);
        }
      }
    }
  }

  CHECK_INVARIANT(
      coreNeighbors.size() >= 1,
      "Unable to find target atom(s) matching core for attachment point");
  if (coreNeighbors.size() == 1) {
    // currently this routine is only called when we know the attachment to
    // one core atom exists.
    return true;
  }

  // at this point we know the target atom is connected to two or more core
  // atoms.  Now check we have attachment points for each.
  //  There should be a map entry for a different attachment point for each of
  //  the core atoms.
  std::vector<int> coreDummyIndices;
  std::for_each(mapping.cbegin(), mapping.cend(),
                [&coreDummyIndices, attachmentIdx](std::pair<int, int> p) {
                  if (p.second == attachmentIdx) {
                    coreDummyIndices.push_back(p.first);
                  }
                });
  if (coreDummyIndices.size() != coreNeighbors.size()) {
    return false;
  }
  // check that there is a bond to an attachment point for every neighbor
  for (const auto neighborIdx : coreNeighbors) {
    bool foundBond = false;
    for (const auto dummyIdx : coreDummyIndices) {
      if (core->getBondBetweenAtoms(neighborIdx, dummyIdx) != nullptr) {
        foundBond = true;
        break;
      }
    }
    if (!foundBond) {
      return false;
    }
  }
  return true;
}

// Convert a matching molecule index to a core index
int RCore::matchingIndexToCoreIndex(int matchingIndex) const {
  auto atom = matchingMol->getAtomWithIdx(matchingIndex);
  CHECK_INVARIANT(atom->hasProp(RLABEL_CORE_INDEX),
                  "Matched atom missing core index");
  return atom->getProp<int>(RLABEL_CORE_INDEX);
}

}  // namespace RDKit
