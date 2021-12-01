//
//  Copyright (C) 2021 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "GenericGroups.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>
#include <algorithm>

namespace RDKit {
class ROMol;

namespace GenericGroups {

namespace Matchers {

using AtomMatcherFunc = std::function<bool(const Atom &)>;
using BondMatcherFunc = std::function<bool(const Bond &)>;

bool AllAtomsMatch(const ROMol &mol, const Atom &atom,
                   boost::dynamic_bitset<> ignore, AtomMatcherFunc matcher,
                   BondMatcherFunc bondMatcher = nullptr,
                   AtomMatcherFunc atLeastOneAtom = nullptr,
                   BondMatcherFunc atLeastOneBond = nullptr) {
  PRECONDITION(&atom.getOwningMol() == &mol, "atom not owned by molecule");
  if (matcher && !matcher(atom)) {
    return false;
  }
  bool atomAtLeast = atLeastOneAtom == nullptr;
  bool bondAtLeast = atLeastOneBond == nullptr;

  std::deque<const Atom *> nbrs;
  nbrs.push_back(&atom);
  while (!nbrs.empty()) {
    const auto atm = nbrs.front();
    if (!atomAtLeast && atLeastOneAtom(*atm)) {
      atomAtLeast = true;
    }
    nbrs.pop_front();
    ignore.set(atm->getIdx());
    for (const auto nbr : mol.atomNeighbors(atm)) {
      if (ignore[nbr->getIdx()]) {
        continue;
      }
      if (matcher && !matcher(*nbr)) {
        return false;
      }
      if (bondMatcher || !bondAtLeast) {
        const auto bnd = mol.getBondBetweenAtoms(atm->getIdx(), nbr->getIdx());
        if (bondMatcher && !(bondMatcher)(*bnd)) {
          return false;
        }
        if (!bondAtLeast && atLeastOneBond(*bnd)) {
          bondAtLeast = true;
        }
      }
      nbrs.push_back(nbr);
    }
  }
  return atomAtLeast && bondAtLeast;
}

bool AlkylAtomMatcher(const ROMol &mol, const Atom &atom,
                      boost::dynamic_bitset<> ignore) {
  if (!mol.getRingInfo() || !mol.getRingInfo()->isInitialized()) {
    MolOps::fastFindRings(mol);
  }
  auto atomMatcher = [](const Atom &at) -> bool {
    return !at.getIsAromatic() &&
           (at.getAtomicNum() == 6 || at.getAtomicNum() == 1);
  };
  auto atLeastMatcher = [](const Atom &at) -> bool {
    return at.getAtomicNum() == 6;
  };
  auto bondMatcher = [](const Bond &bnd) -> bool {
    return bnd.getBondType() == Bond::BondType::SINGLE &&
           !bnd.getIsAromatic() && !queryIsBondInRing(&bnd);
  };
  return AllAtomsMatch(mol, atom, ignore, atomMatcher, bondMatcher,
                       atLeastMatcher);
}

bool AcyclicAtomMatcher(const ROMol &mol, const Atom &atom,
                        boost::dynamic_bitset<> ignore) {
  if (!mol.getRingInfo() || !mol.getRingInfo()->isInitialized()) {
    MolOps::fastFindRings(mol);
  }
  auto atomMatcher = [](const Atom &at) -> bool {
    return at.getOwningMol().getRingInfo()->numAtomRings(at.getIdx()) == 0;
  };
  return AllAtomsMatch(mol, atom, ignore, atomMatcher);
}

bool CarboacyclicAtomMatcher(const ROMol &mol, const Atom &atom,
                             boost::dynamic_bitset<> ignore) {
  if (!mol.getRingInfo() || !mol.getRingInfo()->isInitialized()) {
    MolOps::fastFindRings(mol);
  }
  auto atomMatcher = [](const Atom &at) -> bool {
    return at.getAtomicNum() == 6 &&
           at.getOwningMol().getRingInfo()->numAtomRings(at.getIdx()) == 0;
  };
  return AllAtomsMatch(mol, atom, ignore, atomMatcher);
}

bool HeteroacyclicAtomMatcher(const ROMol &mol, const Atom &atom,
                              boost::dynamic_bitset<> ignore) {
  if (!mol.getRingInfo() || !mol.getRingInfo()->isInitialized()) {
    MolOps::fastFindRings(mol);
  }
  auto atomMatcher = [](const Atom &at) -> bool {
    return at.getOwningMol().getRingInfo()->numAtomRings(at.getIdx()) == 0;
  };
  auto atLeastOne = [](const Atom &at) -> bool {
    return at.getAtomicNum() != 6 && at.getAtomicNum() != 1;
  };
  BondMatcherFunc bondMatcher = nullptr;

  return AllAtomsMatch(mol, atom, ignore, atomMatcher, bondMatcher, atLeastOne);
}

bool AlkoxyacyclicAtomMatcher(const ROMol &mol, const Atom &atom,
                              boost::dynamic_bitset<> ignore) {
  if (!mol.getRingInfo() || !mol.getRingInfo()->isInitialized()) {
    MolOps::fastFindRings(mol);
  }
  if (atom.getDegree() != 2 || atom.getAtomicNum() != 8) {
    return false;
  }
  const Atom *nnbr = nullptr;
  for (const auto *nbr : mol.atomNeighbors(&atom)) {
    if (!ignore[nbr->getIdx()]) {
      nnbr = nbr;
      break;
    }
  }
  if (nnbr == nullptr) {
    return false;
  }
  auto atomMatcher = [](const Atom &at) -> bool {
    return !at.getIsAromatic() &&
           (at.getAtomicNum() == 6 || at.getAtomicNum() == 1);
  };
  auto atLeastMatcher = [](const Atom &at) -> bool {
    return at.getAtomicNum() == 6;
  };
  auto bondMatcher = [](const Bond &bnd) -> bool {
    return bnd.getBondType() == Bond::BondType::SINGLE &&
           !bnd.getIsAromatic() && !queryIsBondInRing(&bnd);
  };
  return AllAtomsMatch(mol, *nnbr, ignore, atomMatcher, bondMatcher,
                       atLeastMatcher);
}

namespace {
bool UnsatAlkXAtomMatcher(const ROMol &mol, const Atom &atom,
                          boost::dynamic_bitset<> ignore,
                          Bond::BondType extraBondType) {
  // nominally requires at least two Cs, but since it can only
  // contain Cs and Hs and since a multiple bond is required, that condition is
  // redundant
  if (!mol.getRingInfo() || !mol.getRingInfo()->isInitialized()) {
    MolOps::fastFindRings(mol);
  }
  auto atomMatcher = [](const Atom &at) -> bool {
    return !at.getIsAromatic() &&
           (at.getAtomicNum() == 6 || at.getAtomicNum() == 1);
  };
  auto bondMatcher = [extraBondType](const Bond &bnd) -> bool {
    return (bnd.getBondType() == Bond::BondType::SINGLE ||
            bnd.getBondType() == extraBondType) &&
           !bnd.getIsAromatic() && !queryIsBondInRing(&bnd);
  };
  auto atLeastMatcher = [extraBondType](const Bond &bnd) -> bool {
    return bnd.getBondType() == extraBondType;
  };
  AtomMatcherFunc atomAtLeast = nullptr;
  return AllAtomsMatch(mol, atom, ignore, atomMatcher, bondMatcher, atomAtLeast,
                       atLeastMatcher);
}
}  // namespace

bool AlkenylAtomMatcher(const ROMol &mol, const Atom &atom,
                        boost::dynamic_bitset<> ignore) {
  return UnsatAlkXAtomMatcher(mol, atom, ignore, Bond::BondType::DOUBLE);
}

bool AlkynylAtomMatcher(const ROMol &mol, const Atom &atom,
                        boost::dynamic_bitset<> ignore) {
  return UnsatAlkXAtomMatcher(mol, atom, ignore, Bond::BondType::TRIPLE);
}

namespace {
bool checkAtomRing(const ROMol &mol, const Atom &atom,
                   const boost::dynamic_bitset<> &ignore,
                   const std::vector<int> &ring, AtomMatcherFunc matcher,
                   AtomMatcherFunc atLeastOne) {
  bool atLeast = atLeastOne == nullptr;
  for (auto aidx : ring) {
    if (aidx != static_cast<int>(atom.getIdx()) &&
        (ignore[aidx] || (matcher && !matcher(*mol.getAtomWithIdx(aidx))))) {
      return false;
    }
    if (!atLeast && atLeastOne(*mol.getAtomWithIdx(aidx))) {
      atLeast = true;
    }
  }
  return atLeast;
}
bool checkBondRing(const ROMol &mol, const std::vector<int> &bring,
                   BondMatcherFunc matcher, BondMatcherFunc atLeastOne) {
  bool atLeast = atLeastOne == nullptr;
  for (auto bidx : bring) {
    if (matcher && !matcher(*mol.getBondWithIdx(bidx))) {
      return false;
    }
    if (!atLeast && atLeastOne(*mol.getBondWithIdx(bidx))) {
      atLeast = true;
    }
  }
  return atLeast;
}

bool FusedRingMatch(const ROMol &mol, const Atom &atom,
                    boost::dynamic_bitset<> ignore,
                    AtomMatcherFunc atomMatcher = nullptr,
                    BondMatcherFunc bondMatcher = nullptr,
                    AtomMatcherFunc atLeastOneAtomPerRing = nullptr,
                    BondMatcherFunc atLeastOneBondPerRing = nullptr,
                    AtomMatcherFunc atLeastOneAtom = nullptr) {
  PRECONDITION(&atom.getOwningMol() == &mol, "atom not owned by molecule");
  if (atomMatcher && !atomMatcher(atom)) {
    return false;
  }

  if (!mol.getRingInfo() || !mol.getRingInfo()->isInitialized()) {
    MolOps::findSSSR(mol);
  }
  if (!mol.getRingInfo()->numAtomRings(atom.getIdx())) {
    return false;
  }

  // we're going to traverse over the entire fused ring system involving this
  // atom start by finding the first ring:
  std::set<int> ringAtoms;
  for (auto i = 0u; i < mol.getRingInfo()->numRings(); ++i) {
    const auto &ring = mol.getRingInfo()->atomRings()[i];
    if (std::find(ring.begin(), ring.end(), atom.getIdx()) != ring.end()) {
      if (!checkAtomRing(mol, atom, ignore, ring, atomMatcher,
                         atLeastOneAtomPerRing)) {
        return false;
      }
      if (!checkBondRing(mol, mol.getRingInfo()->bondRings()[i], bondMatcher,
                         atLeastOneBondPerRing)) {
        return false;
      }

      ringAtoms.insert(ring.begin(), ring.end());
      break;
    }
  }
  // now loop over all rings and find the ones which share at least two atoms
  // with what we've seen so far
  for (auto i = 0u; i < mol.getRingInfo()->numRings(); ++i) {
    const auto &ring = mol.getRingInfo()->atomRings()[i];
    // check overlap of this ring with what we've seen so far and make sure that
    // the new atoms we are adding all pass the test
    std::set<int> sring(ring.begin(), ring.end());
    std::vector<int> diff(sring.size());
    auto dit =
        std::set_difference(sring.begin(), sring.end(), ringAtoms.begin(),
                            ringAtoms.end(), diff.begin());
    auto numNewAtoms = dit - diff.begin();
    if (!numNewAtoms || sring.size() - numNewAtoms < 2) {
      // we don't overlap by at least two atoms
      continue;
    }
    if (!checkAtomRing(mol, atom, ignore, ring, atomMatcher,
                       atLeastOneAtomPerRing)) {
      return false;
    }

    if (!checkBondRing(mol, mol.getRingInfo()->bondRings()[i], bondMatcher,
                       atLeastOneBondPerRing)) {
      return false;
    }
    ringAtoms.insert(diff.begin(), dit);
  }

  if (atLeastOneAtom) {
    return std::find_if(ringAtoms.begin(), ringAtoms.end(),
                        [&mol, atLeastOneAtom](auto idx) -> bool {
                          return atLeastOneAtom(*mol.getAtomWithIdx(idx));
                        }) != ringAtoms.end();
  }

  return true;
}
}  // namespace

bool CarbocycloalkylAtomMatcher(const ROMol &mol, const Atom &atom,
                                boost::dynamic_bitset<> ignore) {
  auto atomMatcher = [](const Atom &at) -> bool {
    return !at.getIsAromatic() && at.getAtomicNum() == 6;
  };
  auto bondMatcher = [](const Bond &bnd) -> bool {
    return !bnd.getIsAromatic() && bnd.getBondType() == Bond::BondType::SINGLE;
  };
  return FusedRingMatch(mol, atom, ignore, atomMatcher, bondMatcher);
}

bool CarbocycloalkenylAtomMatcher(const ROMol &mol, const Atom &atom,
                                  boost::dynamic_bitset<> ignore) {
  auto atomMatcher = [](const Atom &at) -> bool {
    return at.getAtomicNum() == 6;
  };
  auto atLeastOneBond = [](const Bond &bnd) -> bool {
    return bnd.getIsAromatic() || bnd.getBondType() == Bond::BondType::DOUBLE ||
           bnd.getBondType() == Bond::BondType::AROMATIC;
  };
  AtomMatcherFunc atLeastOne = nullptr;
  BondMatcherFunc bondMatcher = nullptr;
  return FusedRingMatch(mol, atom, ignore, atomMatcher, bondMatcher, atLeastOne,
                        atLeastOneBond);
}

bool CarboarylAtomMatcher(const ROMol &mol, const Atom &atom,
                          boost::dynamic_bitset<> ignore) {
  auto atomMatcher = [](const Atom &at) -> bool {
    return at.getIsAromatic() && at.getAtomicNum() == 6;
  };
  auto bondMatcher = [](const Bond &bnd) -> bool {
    return bnd.getIsAromatic() || bnd.getBondType() == Bond::BondType::AROMATIC;
  };
  return FusedRingMatch(mol, atom, ignore, atomMatcher, bondMatcher);
}

bool CarbocyclicAtomMatcher(const ROMol &mol, const Atom &atom,
                            boost::dynamic_bitset<> ignore) {
  auto atomMatcher = [](const Atom &at) -> bool {
    return at.getAtomicNum() == 6;
  };
  return FusedRingMatch(mol, atom, ignore, atomMatcher);
}

bool NoCarbonRingAtomMatcher(const ROMol &mol, const Atom &atom,
                             boost::dynamic_bitset<> ignore) {
  auto atomMatcher = [](const Atom &at) -> bool {
    return at.getAtomicNum() != 6;
  };
  return FusedRingMatch(mol, atom, ignore, atomMatcher);
}

bool HeterocyclicAtomMatcher(const ROMol &mol, const Atom &atom,
                             boost::dynamic_bitset<> ignore) {
  auto atLeastOne = [](const Atom &at) -> bool {
    return at.getAtomicNum() != 6 && at.getAtomicNum() != 1;
  };
  AtomMatcherFunc atomMatcher = nullptr;
  AtomMatcherFunc oneAtomPerRing = nullptr;
  BondMatcherFunc bondMatcher = nullptr;
  BondMatcherFunc oneBondPerRing = nullptr;
  return FusedRingMatch(mol, atom, ignore, atomMatcher, bondMatcher,
                        oneAtomPerRing, oneBondPerRing, atLeastOne);
}

bool HeteroarylAtomMatcher(const ROMol &mol, const Atom &atom,
                           boost::dynamic_bitset<> ignore) {
  auto atomMatcher = [](const Atom &at) -> bool { return at.getIsAromatic(); };
  auto bondMatcher = [](const Bond &bnd) -> bool {
    return bnd.getIsAromatic() || bnd.getBondType() == Bond::BondType::AROMATIC;
  };
  auto atLeastOne = [](const Atom &at) -> bool {
    return at.getAtomicNum() != 6 && at.getAtomicNum() != 1;
  };
  AtomMatcherFunc oneAtomPerRing = nullptr;
  BondMatcherFunc oneBondPerRing = nullptr;
  return FusedRingMatch(mol, atom, ignore, atomMatcher, bondMatcher,
                        oneAtomPerRing, oneBondPerRing, atLeastOne);
}

bool CyclicAtomMatcher(const ROMol &mol, const Atom &atom,
                       boost::dynamic_bitset<> ignore) {
  if (!mol.getRingInfo() || !mol.getRingInfo()->isInitialized()) {
    MolOps::fastFindRings(mol);
  }
  auto atomMatcher = [](const Atom &at) -> bool {
    return at.getOwningMol().getRingInfo()->numAtomRings(at.getIdx()) > 0;
  };
  return FusedRingMatch(mol, atom, ignore, atomMatcher);
}

}  // namespace Matchers

bool genericAtomMatcher(const ROMol &mol, const ROMol &query,
                        const std::vector<unsigned int> &match) {
  boost::dynamic_bitset<> ignore(mol.getNumAtoms());
  for (const auto idx : match) {
    ignore.set(idx);
  }

  for (const auto atom : query.atoms()) {
    if (atom->getDegree() != 1) {
      continue;
    }
    std::string genericLabel;
    if (atom->getPropIfPresent(common_properties::_QueryAtomGenericLabel,
                               genericLabel)) {
      auto found = genericMatchers.find(genericLabel);
      if (found != genericMatchers.end() &&
          !found->second(mol, *mol.getAtomWithIdx(match[atom->getIdx()]),
                         ignore)) {
        return false;
      }
    }
  }
  return true;
}

void convertGenericQueriesToSubstanceGroups(ROMol &mol) {
  for (const auto atom : mol.atoms()) {
    std::string label;
    if (atom->getPropIfPresent(common_properties::_QueryAtomGenericLabel,
                               label)) {
      SubstanceGroup sg(&mol, "SUP");
      sg.setProp("LABEL", label);
      sg.addAtomWithIdx(atom->getIdx());
      addSubstanceGroup(mol, sg);
      atom->clearProp(common_properties::_QueryAtomGenericLabel);
    }
  }
}

void setGenericQueriesFromProperties(ROMol &mol, bool useAtomLabels,
                                     bool useSGroups) {
  if (useAtomLabels) {
    for (const auto atom : mol.atoms()) {
      std::string label;
      if (atom->getPropIfPresent(common_properties::atomLabel, label)) {
        // pseudoatom labels from CXSMILES end with "_p"... strip that if
        // present
        if (label.size() > 4 && label.compare(label.size() - 2, 2, "_p") == 0) {
          label = label.substr(0, label.size() - 2);
        }
        if (genericMatchers.find(label) != genericMatchers.end()) {
          atom->setProp(common_properties::_QueryAtomGenericLabel, label);
          atom->clearProp(common_properties::atomLabel);
        }
      }
    }
  }
  if (useSGroups) {
    auto &sgs = getSubstanceGroups(mol);
    auto iter = sgs.begin();
    while (iter != sgs.end()) {
      const auto &sgroup = *iter;
      if (sgroup.getProp<std::string>("TYPE") == "SUP") {
        std::string label;
        if (sgroup.getPropIfPresent("LABEL", label) &&
            genericMatchers.find(label) != genericMatchers.end()) {
          for (auto aidx : sgroup.getAtoms()) {
            mol.getAtomWithIdx(aidx)->setProp(
                common_properties::_QueryAtomGenericLabel, label);
          }
          iter = sgs.erase(iter);
        } else {
          ++iter;
        }
      } else {
        ++iter;
      }
    }
  }
}
}  // namespace GenericGroups
}  // namespace RDKit
