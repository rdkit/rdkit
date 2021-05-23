
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

std::vector<MatchVectType> RCore::matchTerminalUserRGroups(
    const RWMol &target, MatchVectType match) const {
  std::transform(match.begin(), match.end(), match.begin(),
                 [this](const std::pair<int, int> &mapping) {
                   auto queryIdx =
                       this->matchingIndexToCoreIndex(mapping.first);
                   std::pair<int, int> newMapping(queryIdx, mapping.second);
                   return newMapping;
                 });
  std::map<int, int> matchMap(match.cbegin(), match.cend());

  std::vector<MatchVectType> allMappings;
  if (terminalRGroupAtomsWithUserLabels.size() == 0) {
    allMappings.push_back(match);
    return allMappings;
  }

  std::set<int> mappedTargetIdx;
  std::transform(
      match.cbegin(), match.cend(),
      std::inserter(mappedTargetIdx, mappedTargetIdx.begin()),
      [](const std::pair<int, int> &mapping) { return mapping.second; });

  std::vector<int> dummiesWithMapping;
  std::vector<std::vector<int>> availableMappingsForDummy;

  SubstructMatchParameters ssParameters;
  ssParameters.useChirality = true;
  for (const auto dummyIdx : terminalRGroupAtomsWithUserLabels) {
    const int neighborIdx = terminalRGroupAtomToNeighbor.find(dummyIdx)->second;
    const auto coreBond = core->getBondBetweenAtoms(dummyIdx, neighborIdx);
    const int targetIdx = matchMap[neighborIdx];
    const auto targetAtom = target.getAtomWithIdx(targetIdx);
    ROMol::ADJ_ITER nbrIter, endNbrs;
    std::vector<int> available;
    boost::tie(nbrIter, endNbrs) = target.getAtomNeighbors(targetAtom);
    while (nbrIter != endNbrs) {
      if (mappedTargetIdx.find(*nbrIter) == mappedTargetIdx.end()) {
        const auto targetBond = target.getBondBetweenAtoms(targetIdx, *nbrIter);
        if (bondCompat(coreBond, targetBond, ssParameters)) {
          available.push_back(*nbrIter);
        }
      }
      ++nbrIter;
    }

    if (available.size()) {
      dummiesWithMapping.push_back(dummyIdx);
      availableMappingsForDummy.push_back(available);
    } else {
      // We could continue here and allow a userR group to be unmapped
      // However, the code below will fail at the molMatchFunctor check as all
      // query atoms are not matched to the target- the query would need to be
      // edited before passing to the molMatchFunctor.
      return allMappings;
    }
  }
  if (availableMappingsForDummy.size() == 0) {
    allMappings.push_back(match);
    return allMappings;
  }

  // enumerate
  const auto allAvailableMappings = cartesianProduct(availableMappingsForDummy);
  size_t size = allAvailableMappings[0].size() + match.size();
  auto queryIndices = new std::uint32_t[size];
  auto targetIndices = new std::uint32_t[size];
  for (size_t position = 0; position < match.size(); position++) {
    const auto &pair = match[position];
    queryIndices[position] = pair.first;
    targetIndices[position] = pair.second;
  }

  MolMatchFinalCheckFunctor molMatchFunctor(*core, target, ssParameters);
  boost::dynamic_bitset<> targetBondsPresent(target.getNumBonds());

  for (const auto &dummyMapping : allAvailableMappings) {
    CHECK_INVARIANT(match.size() + dummyMapping.size() == size,
                    "Size error in dummy mapping");
    auto duplicateBonds = false;
    targetBondsPresent.reset();
    for (size_t i = 0; i < dummyMapping.size(); i++) {
      size_t position = match.size() + i;
      queryIndices[position] = dummiesWithMapping[i];
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

  delete queryIndices;
  delete targetIndices;
  return allMappings;
}

bool RCore::checkAllBondsToAttachmentPointPresent(
    const ROMol &mol, const int attachmentIdx,
    const MatchVectType &mapping) const {
  const auto atom = mol.getAtomWithIdx(attachmentIdx);
  std::set<int> coreNeighbors;
  for (const auto &nbri :
       boost::make_iterator_range(mol.getAtomNeighbors(atom))) {
    const auto &nbr = mol[nbri];
    // could a neighbor to an r group attachment match another r group
    // attachment?
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

int RCore::matchingIndexToCoreIndex(int matchingIndex) const {
  auto atom = matchingMol->getAtomWithIdx(matchingIndex);
  CHECK_INVARIANT(atom->hasProp(RLABEL_CORE_INDEX),
                  "Matched atom missing core index");
  return atom->getProp<int>(RLABEL_CORE_INDEX);
}

void RCore::buildMatchingMol() {
  matchingMol = boost::make_shared<RWMol>(*core);
  terminalRGroupAtomsWithUserLabels.clear();
  terminalRGroupAtomToNeighbor.clear();
  RWMol::ATOM_PTR_VECT atomsToRemove;
  for (auto atom : matchingMol->atoms()) {
    atom->setProp<int>(RLABEL_CORE_INDEX, atom->getIdx());
    if (atom->getAtomicNum() == 0 && atom->getDegree() == 1 &&
        isUserRLabel(*atom)) {
      atomsToRemove.push_back(atom);
      terminalRGroupAtomsWithUserLabels.insert(atom->getIdx());
      const int neighborIdx = *matchingMol->getAtomNeighbors(atom).first;
      terminalRGroupAtomToNeighbor.emplace(atom->getIdx(), neighborIdx);
    }
  }

  for (auto atom : atomsToRemove) {
    matchingMol->removeAtom(atom);
  }
}

}  // namespace RDKit