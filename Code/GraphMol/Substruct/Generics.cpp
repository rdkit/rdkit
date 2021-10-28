//
//  Copyright (C) 2021 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Generics.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>
#include <algorithm>

namespace RDKit {
class ROMol;

namespace SubstructSearch {

namespace Generics {
bool AllAtomsMatch(const ROMol &mol, const Atom &atom,
                   boost::dynamic_bitset<> ignore,
                   std::function<bool(const Atom &)> matcher,
                   std::function<bool(const Bond &)> bondMatcher,
                   std::function<bool(const Atom &)> atLeastOneAtom,
                   std::function<bool(const Bond &)> atLeastOneBond) {
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
                       atLeastMatcher, nullptr);
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
  return AllAtomsMatch(mol, atom, ignore, atomMatcher, bondMatcher, nullptr,
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

bool FusedRingMatch(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore,
    std::function<bool(const Atom &)> atomMatcher = nullptr,
    std::function<bool(const Bond &)> bondMatcher = nullptr,
    std::function<bool(const Atom &)> atLeastOneAtom = nullptr,
    std::function<bool(const Bond &)> atLeastOneBond = nullptr) {
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
      bool atomAtLeast = atLeastOneAtom == nullptr;
      for (auto aidx : ring) {
        if (aidx != static_cast<int>(atom.getIdx()) &&
            (ignore[aidx] ||
             (atomMatcher && !atomMatcher(*mol.getAtomWithIdx(aidx))))) {
          return false;
        }
        if (!atomAtLeast && atLeastOneAtom(*mol.getAtomWithIdx(aidx))) {
          atomAtLeast = true;
        }
      }
      if (!atomAtLeast) {
        return false;
      }

      bool bondAtLeast = atLeastOneBond == nullptr;
      for (auto bidx : mol.getRingInfo()->bondRings()[i]) {
        const auto bond = mol.getBondWithIdx(bidx);
        if (bondMatcher && !bondMatcher(*bond)) {
          return false;
        }
        if (!bondAtLeast && atLeastOneBond(*bond)) {
          bondAtLeast = true;
        }
      }
      if (!bondAtLeast) {
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
    // check the atoms
    bool atomAtLeast = atLeastOneAtom == nullptr;
    for (auto aidx : ring) {
      if ((aidx != static_cast<int>(atom.getIdx()) && ignore[aidx]) ||
          (atomMatcher && !atomMatcher(*mol.getAtomWithIdx(aidx)))) {
        return false;
      }
      if (!atomAtLeast && atLeastOneAtom(*mol.getAtomWithIdx(aidx))) {
        atomAtLeast = true;
      }
    }
    if (!atomAtLeast) {
      return false;
    }
    // check the bonds:
    bool bondAtLeast = atLeastOneBond == nullptr;
    for (auto bidx : mol.getRingInfo()->bondRings()[i]) {
      const auto bond = mol.getBondWithIdx(bidx);
      if (bondMatcher && !bondMatcher(*bond)) {
        return false;
      }
      if (!bondAtLeast && atLeastOneBond(*bond)) {
        bondAtLeast = true;
      }
    }
    if (!bondAtLeast) {
      return false;
    }
    ringAtoms.insert(diff.begin(), dit);
  }

  return true;
}

bool CycloalkylAtomMatcher(const ROMol &mol, const Atom &atom,
                           boost::dynamic_bitset<> ignore) {
  auto atomMatcher = [](const Atom &at) -> bool {
    return !at.getIsAromatic() && at.getAtomicNum() == 6;
  };
  auto bondMatcher = [](const Bond &bnd) -> bool {
    return !bnd.getIsAromatic() && bnd.getBondType() == Bond::BondType::SINGLE;
  };
  return FusedRingMatch(mol, atom, ignore, atomMatcher, bondMatcher);
}

bool CycloalkenylAtomMatcher(const ROMol &mol, const Atom &atom,
                             boost::dynamic_bitset<> ignore) {
  auto atomMatcher = [](const Atom &at) -> bool {
    return at.getAtomicNum() == 6;
  };
  auto atLeastOneBond = [](const Bond &bnd) -> bool {
    return bnd.getIsAromatic() || bnd.getBondType() == Bond::BondType::DOUBLE ||
           bnd.getBondType() == Bond::BondType::AROMATIC;
  };
  return FusedRingMatch(mol, atom, ignore, atomMatcher, nullptr, nullptr,
                        atLeastOneBond);
  // if (!atomMatcher(atom)) {
  //   return false;
  // }
  // if (!mol.getRingInfo() || !mol.getRingInfo()->isInitialized()) {
  //   MolOps::findSSSR(mol);
  // }
  // if (!mol.getRingInfo()->numAtomRings(atom.getIdx())) {
  //   return false;
  // }
  // bool seenDoubleOrAromaticBond = false;
  // for (const auto &ring : mol.getRingInfo()->atomRings()) {
  //   if (std::find(ring.begin(), ring.end(), atom.getIdx()) != ring.end()) {
  //     for (auto aidx : ring) {
  //       if (aidx != static_cast<int>(atom.getIdx()) &&
  //           (ignore[aidx] || !atomMatcher(*mol.getAtomWithIdx(aidx)))) {
  //         return false;
  //       }
  //     }
  //     if (!seenDoubleOrAromaticBond) {
  //       for (unsigned int i = 1; i < ring.size(); ++i) {
  //         auto bnd = mol.getBondBetweenAtoms(ring[i - 1], ring[i]);
  //         ASSERT_INVARIANT(bnd, "expected bond not found");
  //         if (bnd->getIsAromatic() ||
  //             bnd->getBondType() == Bond::BondType::DOUBLE ||
  //             bnd->getBondType() == Bond::BondType::AROMATIC) {
  //           seenDoubleOrAromaticBond = true;
  //           break;
  //         }
  //       }
  //       if (!seenDoubleOrAromaticBond) {
  //         auto bnd = mol.getBondBetweenAtoms(ring[ring.size() - 1], ring[0]);
  //         ASSERT_INVARIANT(bnd, "expected bond not found");
  //         if (bnd->getIsAromatic() ||
  //             bnd->getBondType() == Bond::BondType::DOUBLE ||
  //             bnd->getBondType() == Bond::BondType::AROMATIC) {
  //           seenDoubleOrAromaticBond = true;
  //         }
  //       }
  //     }
  //   }
  // }
  // return seenDoubleOrAromaticBond;
}

}  // namespace Generics

bool GenericAtomMatcher(const ROMol &mol, const ROMol &query,
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
      auto found = Generics::genericMatchers.find(genericLabel);
      if (found != Generics::genericMatchers.end() &&
          !found->second(mol, *mol.getAtomWithIdx(match[atom->getIdx()]),
                         ignore)) {
        return false;
      }
    }
  }
  return true;
}
}  // namespace SubstructSearch

}  // namespace RDKit
