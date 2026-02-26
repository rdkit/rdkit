//
//  Copyright (C) 2001-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Canon.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/new_canon.h>

#include <GraphMol/SmilesParse/SmilesParseOps.h>
#include <GraphMol/RDKitQueries.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/hash/hash.hpp>
#include <RDGeneral/utils.h>

#include <ranges>
#include <queue>
#include <algorithm>

namespace RDKit {
namespace Canon {
namespace {
static constexpr Bond::BondDir flipStereoBondDir(Bond::BondDir bondDir) {
  switch (bondDir) {
    case Bond::ENDUPRIGHT:
      return Bond::ENDDOWNRIGHT;
    case Bond::ENDDOWNRIGHT:
      return Bond::ENDUPRIGHT;
    default:
      return bondDir;
  }
}

void setDirectionFromNeighboringBond(Bond &sourceBond, bool isSourceBondFlipped,
                                     Bond &targetBond,
                                     bool isTargetBondFlipped) {
  auto dir = sourceBond.getBondDir();

  // By default, both bonds on the same side of the double bond
  // should have opposite directions, but this can change if
  // one (and only one) of the bonds is flipped (if both were
  // flipped, the flips in the direction would cancel out).
  if (isSourceBondFlipped == isTargetBondFlipped) {
    dir = flipStereoBondDir(dir);
  }

  targetBond.setBondDir(dir);
}

Bond::BondDir getReferenceDirection(const Bond &dblBond, const Atom &refAtom,
                                    const Atom &targetAtom,
                                    const Bond &refControllingBond,
                                    bool refIsFlipped, const Bond &targetBond,
                                    bool targetIsFlipped) {
  Bond::BondDir dir = Bond::NONE;
  if (dblBond.getStereo() == Bond::STEREOE ||
      dblBond.getStereo() == Bond::STEREOTRANS) {
    dir = refControllingBond.getBondDir();
  } else if (dblBond.getStereo() == Bond::STEREOZ ||
             dblBond.getStereo() == Bond::STEREOCIS) {
    dir = flipStereoBondDir(refControllingBond.getBondDir());
  }
  CHECK_INVARIANT(dir != Bond::NONE, "stereo not set");

  // If we're not looking at the bonds used to determine the
  // stereochemistry, we need to flip the setting on the other bond:
  const INT_VECT &stereoAtoms = dblBond.getStereoAtoms();

  if (refAtom.getDegree() == 3 &&
      std::ranges::find(stereoAtoms,
                        static_cast<int>(refControllingBond.getOtherAtomIdx(
                            refAtom.getIdx()))) == stereoAtoms.end()) {
    dir = flipStereoBondDir(dir);
  }
  if (targetAtom.getDegree() == 3 &&
      std::ranges::find(stereoAtoms,
                        static_cast<int>(targetBond.getOtherAtomIdx(
                            targetAtom.getIdx()))) == stereoAtoms.end()) {
    dir = flipStereoBondDir(dir);
  }

  // XOR: flips will cancel out if we do both
  if (refIsFlipped != targetIsFlipped) {
    dir = flipStereoBondDir(dir);
  }

  return dir;
}

void checkSameSideDirsAreCompatible(const Bond &firstBond,
                                    const Bond &secondBond,
                                    bool isFirstBondFlipped,
                                    bool isSecondBondFlipped) {
  auto dirsShouldMatch = isFirstBondFlipped != isSecondBondFlipped;
  auto dirsMatch = firstBond.getBondDir() == secondBond.getBondDir();

  CHECK_INVARIANT(dirsMatch == dirsShouldMatch,
                  "inconsistent bond direction state");
}

}  // namespace

namespace details {
bool isUnsaturated(const Atom *atom, const ROMol &mol) {
  for (auto bond : mol.atomBonds(atom)) {
    // can't just check for single bonds, because dative bonds also have an
    // order of 1
    if (bond->getBondTypeAsDouble() > 1) {
      return true;
    }
  }
  return false;
}

bool hasSingleHQuery(const Atom::QUERYATOM_QUERY *q) {
  // list queries are series of nested ors of AtomAtomicNum queries
  PRECONDITION(q, "bad query");
  bool res = false;
  const auto &descr = q->getDescription();
  if (descr == "AtomAnd") {
    for (auto cIt = q->beginChildren(); cIt != q->endChildren(); ++cIt) {
      const auto &cDescr = (*cIt)->getDescription();
      if (cDescr == "AtomHCount") {
        return !(*cIt)->getNegation() &&
               ((ATOM_EQUALS_QUERY *)(*cIt).get())->getVal() == 1;
      } else if (cDescr == "AtomAnd") {
        res = hasSingleHQuery((*cIt).get());
        if (res) {
          return true;
        }
      }
    }
  }
  return res;
}

bool atomHasFourthValence(const Atom *atom) {
  if (atom->getNumExplicitHs() == 1 ||
      (!atom->needsUpdatePropertyCache() &&
       atom->getValence(Atom::ValenceType::IMPLICIT) == 1)) {
    return true;
  }
  if (atom->hasQuery()) {
    // the SMARTS [C@@H] produces an atom with a H query, but we also
    // need to treat this like an explicit H for chirality purposes
    // This was Github #1489
    return hasSingleHQuery(atom->getQuery());
  }
  return false;
}
}  // namespace details

bool chiralAtomNeedsTagInversion(const RDKit::ROMol &mol,
                                 const RDKit::Atom *atom, bool isAtomFirst,
                                 size_t numClosures) {
  PRECONDITION(atom, "bad atom");
  return atom->getDegree() == 3 &&
         ((isAtomFirst && atom->getNumExplicitHs() == 1) ||
          (!details::atomHasFourthValence(atom) && numClosures == 1 &&
           !details::isUnsaturated(atom, mol)));
}

auto _possibleCompare = [](const PossibleType &arg1, const PossibleType &arg2) {
  return (std::get<0>(arg1) < std::get<0>(arg2));
};

// FIX: this may only be of interest from the SmilesWriter, should we
// move it there?
//
//
void canonicalizeDoubleBond(Bond *dblBond, const UINT_VECT &bondVisitOrders,
                            const UINT_VECT &atomVisitOrders,
                            std::vector<int8_t> &bondDirCounts,
                            std::vector<int8_t> &atomDirCounts) {
  PRECONDITION(dblBond, "bad bond");
  PRECONDITION(dblBond->getBondType() == Bond::DOUBLE, "bad bond order");
  PRECONDITION(dblBond->getStereo() > Bond::STEREOANY, "bad bond stereo");
  PRECONDITION(dblBond->getStereoAtoms().size() >= 2, "bad bond stereo atoms");
  PRECONDITION(atomVisitOrders[dblBond->getBeginAtomIdx()] > 0 ||
                   atomVisitOrders[dblBond->getEndAtomIdx()] > 0,
               "neither end atom traversed");

  Atom *atom1 = dblBond->getBeginAtom();
  Atom *atom2 = dblBond->getEndAtom();
  // we only worry about double bonds that begin and end at atoms
  // of degree 2 or 3:
  if ((atom1->getDegree() != 2 && atom1->getDegree() != 3) ||
      (atom2->getDegree() != 2 && atom2->getDegree() != 3)) {
    return;
  }
  // ensure that atom1 is the lower numbered atom of the double bond (the one
  // traversed first)
  if (atomVisitOrders[dblBond->getBeginAtomIdx()] >=
      atomVisitOrders[dblBond->getEndAtomIdx()]) {
    std::swap(atom1, atom2);
  }

  Bond *firstFromAtom1 = nullptr, *secondFromAtom1 = nullptr;
  Bond *firstFromAtom2 = nullptr, *secondFromAtom2 = nullptr;

  ROMol &mol = dblBond->getOwningMol();

  // -------------------------------------------------------
  // find the lowest visit order bonds from each end and determine
  // if anything is already constraining our choice of directions:
  bool dir1Set = false, dir2Set = false;

  auto findNeighborBonds = [&mol, &dblBond, &bondDirCounts, &bondVisitOrders](
                               auto atom, auto &firstNeighborBond,
                               auto &secondNeighborBond, auto &dirSet) {
    auto firstVisitOrder = mol.getNumBonds() + 1;
    for (const auto bond : mol.atomBonds(atom)) {
      if (bond == dblBond || !canSetDoubleBondStereo(*bond)) {
        continue;
      }

      auto bondIdx = bond->getIdx();
      if (bondDirCounts[bondIdx] > 0) {
        dirSet = true;
      }
      if (!firstNeighborBond || bondVisitOrders[bondIdx] < firstVisitOrder) {
        if (firstNeighborBond) {
          secondNeighborBond = firstNeighborBond;
        }
        firstNeighborBond = bond;
        firstVisitOrder = bondVisitOrders[bondIdx];
      } else {
        secondNeighborBond = bond;
      }
    }
  };

  findNeighborBonds(atom1, firstFromAtom1, secondFromAtom1, dir1Set);
  findNeighborBonds(atom2, firstFromAtom2, secondFromAtom2, dir2Set);

  // Make sure we found everything we need to find.
  //   This really shouldn't be a problem, but molecules can end up in odd
  //   states; for example, allenes can end up here. Instead of checking for
  //   them explicitly, exit early in any such possible state.
  if (!firstFromAtom1 || !firstFromAtom2) {
    return;
  }

  // We interpret double bonds like this (this is a TRANS bond,
  // a CIS one would be similar, but both anchors would be either
  // above or below the double bond):
  //
  //   anchor1
  //          |
  //           atom1 === atom2
  //                          |
  //                           anchor2
  //
  // When parsing a SMILES, we expect anchor1 to come before atom1,
  // and atom2 before anchor2. But that is not always the case.
  // Inverting the order has effects on the bond directions, so
  // we define the variables below to help us find the right direction
  // for the bond directions. These are interpreted as follows:
  //
  // - firstFromAtom1 and secondFromAtom1 are considered "flipped" if
  // the anchor would come after the double bond start atom1, e.g.
  // in the SMILES C(\N)(/O)=C(/N)\O they are both flipped. Note
  // that, with this definition secondFromAtom1 is usually "flipped".
  //
  bool isFirstFromAtom1Flipped = [&]() {
    auto anchorIdx = firstFromAtom1->getOtherAtom(atom1)->getIdx();
    return (atomVisitOrders[atom1->getIdx()] < atomVisitOrders[anchorIdx]) !=
           firstFromAtom1->hasProp(
               common_properties::_TraversalRingClosureBond);
  }();

  bool isSecondFromAtom1Flipped = [&]() {
    if (secondFromAtom1 == nullptr) {
      return false;
    }
    auto anchorIdx = secondFromAtom1->getOtherAtom(atom1)->getIdx();
    return (atomVisitOrders[atom1->getIdx()] < atomVisitOrders[anchorIdx]) !=
           secondFromAtom1->hasProp(
               common_properties::_TraversalRingClosureBond);
  }();

  // - firstFromAtom2 and secondFromAtom2 are considered "flipped" if
  // the anchor atom comes before the double bond end atom2 (I think
  // this always requires rings to be involved). An example of firstFromAtom2
  //  being flipped would be the C atom in "C\N=2" in C=c1s/c2n(c1=O)CCCCC\N=2

  bool isFirstFromAtom2Flipped = [&]() {
    auto anchorIdx = firstFromAtom2->getOtherAtom(atom2)->getIdx();
    return (atomVisitOrders[anchorIdx] < atomVisitOrders[atom2->getIdx()]) !=
           firstFromAtom2->hasProp(
               common_properties::_TraversalRingClosureBond);
  }();

  bool isSecondFromAtom2Flipped = [&]() {
    if (secondFromAtom2 == nullptr) {
      return false;
    }
    auto anchorIdx = secondFromAtom2->getOtherAtom(atom2)->getIdx();
    return (atomVisitOrders[anchorIdx] < atomVisitOrders[atom2->getIdx()]) !=
           secondFromAtom2->hasProp(
               common_properties::_TraversalRingClosureBond);
  }();

  // Both directions are already set. Update accounting
  // and check if both directions on each side are set.
  // We hit this in cases with cycles like CO/C1=C/C=C\C=C/C=N\1.
  if (dir1Set && dir2Set) {
    // Check that directions on atom1 side are present and consistent
    if (secondFromAtom1) {
      if (!bondDirCounts[firstFromAtom1->getIdx()]) {
        setDirectionFromNeighboringBond(
            *secondFromAtom1, isSecondFromAtom1Flipped, *firstFromAtom1,
            isFirstFromAtom1Flipped);
      } else if (!bondDirCounts[secondFromAtom1->getIdx()]) {
        setDirectionFromNeighboringBond(
            *firstFromAtom1, isFirstFromAtom1Flipped, *secondFromAtom1,
            isSecondFromAtom1Flipped);
      } else {
        checkSameSideDirsAreCompatible(*firstFromAtom1, *secondFromAtom1,
                                       isFirstFromAtom1Flipped,
                                       isSecondFromAtom1Flipped);
      }

      bondDirCounts[secondFromAtom1->getIdx()] += 1;
      atomDirCounts[atom1->getIdx()] += 1;
    }
    bondDirCounts[firstFromAtom1->getIdx()] += 1;
    atomDirCounts[atom1->getIdx()] += 1;

    // Check that directions on atom2 side are present and consistent
    if (secondFromAtom2) {
      if (!bondDirCounts[firstFromAtom2->getIdx()]) {
        setDirectionFromNeighboringBond(
            *secondFromAtom2, isSecondFromAtom2Flipped, *firstFromAtom2,
            isFirstFromAtom2Flipped);
      } else if (!bondDirCounts[secondFromAtom2->getIdx()]) {
        setDirectionFromNeighboringBond(
            *firstFromAtom2, isFirstFromAtom2Flipped, *secondFromAtom2,
            isSecondFromAtom2Flipped);
      } else {
        checkSameSideDirsAreCompatible(*firstFromAtom2, *secondFromAtom2,
                                       isFirstFromAtom2Flipped,
                                       isSecondFromAtom2Flipped);
      }

      bondDirCounts[secondFromAtom2->getIdx()] += 1;
      atomDirCounts[atom2->getIdx()] += 1;
    }
    bondDirCounts[firstFromAtom2->getIdx()] += 1;
    atomDirCounts[atom2->getIdx()] += 1;

    // Finally, check that directions across the double bond are consistent
    auto expectedFirstFromAtom2Dir = getReferenceDirection(
        *dblBond, *atom1, *atom2, *firstFromAtom1, isFirstFromAtom1Flipped,
        *firstFromAtom2, isFirstFromAtom2Flipped);
    CHECK_INVARIANT(expectedFirstFromAtom2Dir == firstFromAtom2->getBondDir(),
                    "inconsistent bond direction state");

    return;
  }

  bool setFromBond1 = true;
  Bond *atom1ControllingBond = firstFromAtom1;
  Bond *atom2ControllingBond = firstFromAtom2;
  if (!dir1Set && !dir2Set) {
    // ----------------------------------
    // nothing has touched our bonds so far, so set the
    // directions to "arbitrary" values:

    // the bond we came in on becomes ENDUPRIGHT:
    auto atom1Dir = Bond::ENDUPRIGHT;
    firstFromAtom1->setBondDir(atom1Dir);

    bondDirCounts[firstFromAtom1->getIdx()] += 1;
    atomDirCounts[atom1->getIdx()] += 1;

  } else if (!dir2Set) {
    // at least one of the bonds on atom1 has its directionality set already:
    if (bondDirCounts[firstFromAtom1->getIdx()] > 0) {
      // The first bond's direction has been set at some earlier point:
      bondDirCounts[firstFromAtom1->getIdx()] += 1;
      atomDirCounts[atom1->getIdx()] += 1;

      if (secondFromAtom1 && bondDirCounts[secondFromAtom1->getIdx()]) {
        // both bonds have their directionalities set, check if
        // they are compatible.
        checkSameSideDirsAreCompatible(*firstFromAtom1, *secondFromAtom1,
                                       isFirstFromAtom1Flipped,
                                       isSecondFromAtom1Flipped);

        bondDirCounts[secondFromAtom1->getIdx()] += 1;
        atomDirCounts[atom1->getIdx()] += 1;
      }
    } else {
      // the second bond must be present and setting the direction:
      CHECK_INVARIANT(secondFromAtom1, "inconsistent state");
      CHECK_INVARIANT(bondDirCounts[secondFromAtom1->getIdx()] > 0,
                      "inconsistent state");

      setDirectionFromNeighboringBond(*secondFromAtom1,
                                      isSecondFromAtom1Flipped, *firstFromAtom1,
                                      isFirstFromAtom1Flipped);

      // acknowledge that secondFromAtom1 is relevant for this bond,
      // and prevent removeRedundantBondDirSpecs from removing this
      // direction.
      bondDirCounts[secondFromAtom1->getIdx()] += 1;

      bondDirCounts[firstFromAtom1->getIdx()] += 1;
      atomDirCounts[atom1->getIdx()] += 2;
      atom1ControllingBond = secondFromAtom1;
    }
  } else {
    // dir2 has been set, and dir1 hasn't: we're dealing with a stereochem
    // specification on a ring double bond:
    setFromBond1 = false;
    // at least one of the bonds on atom2 has its directionality set already:
    if (bondDirCounts[firstFromAtom2->getIdx()] > 0) {
      // The second bond's direction has been set at some earlier point:
      bondDirCounts[firstFromAtom2->getIdx()] += 1;
      atomDirCounts[atom2->getIdx()] += 1;

      if (secondFromAtom2 && bondDirCounts[secondFromAtom2->getIdx()]) {
        // both bonds have their directionalities set, check if
        // they are compatible.
        checkSameSideDirsAreCompatible(*firstFromAtom2, *secondFromAtom2,
                                       isFirstFromAtom2Flipped,
                                       isSecondFromAtom2Flipped);

        bondDirCounts[secondFromAtom2->getIdx()] += 1;
        atomDirCounts[atom2->getIdx()] += 1;
      }

    } else {
      // the second bond must be present and setting the direction:
      CHECK_INVARIANT(secondFromAtom2, "inconsistent state");
      CHECK_INVARIANT(bondDirCounts[secondFromAtom2->getIdx()] > 0,
                      "inconsistent state");

      setDirectionFromNeighboringBond(*secondFromAtom2,
                                      isSecondFromAtom2Flipped, *firstFromAtom2,
                                      isFirstFromAtom2Flipped);

      // acknowledge that secondFromAtom2 is relevant for this bond,
      // and prevent removeRedundantBondDirSpecs from removing this
      // direction.
      bondDirCounts[secondFromAtom2->getIdx()] += 1;

      bondDirCounts[firstFromAtom2->getIdx()] += 1;
      atomDirCounts[atom2->getIdx()] += 2;
      atom2ControllingBond = secondFromAtom2;
    }
  }  // end of the ring stereochemistry if

  // now set the directionality on the other side:
  if (setFromBond1) {
    auto isControllingAtomFlipped =
        (atom1ControllingBond == firstFromAtom1 ? isFirstFromAtom1Flipped
                                                : isSecondFromAtom1Flipped);

    auto atom2Dir = getReferenceDirection(
        *dblBond, *atom1, *atom2, *atom1ControllingBond,
        isControllingAtomFlipped, *firstFromAtom2, isFirstFromAtom2Flipped);

    firstFromAtom2->setBondDir(atom2Dir);

    bondDirCounts[firstFromAtom2->getIdx()] += 1;
    atomDirCounts[atom2->getIdx()] += 1;
  } else {
    auto isControllingAtomFlipped =
        (atom2ControllingBond == firstFromAtom2 ? isFirstFromAtom2Flipped
                                                : isSecondFromAtom2Flipped);
    auto atom1Dir = getReferenceDirection(
        *dblBond, *atom2, *atom1, *atom2ControllingBond,
        isControllingAtomFlipped, *firstFromAtom1, isFirstFromAtom1Flipped);

    firstFromAtom1->setBondDir(atom1Dir);

    bondDirCounts[firstFromAtom1->getIdx()] += 1;
    atomDirCounts[atom1->getIdx()] += 1;
  }

  // -----------------------------------
  //
  // Check if there are other bonds from atoms 1 and 2 that need
  // to have their directionalities set:
  ///
  if (atom1->getDegree() == 3 && secondFromAtom1) {
    if (!bondDirCounts[secondFromAtom1->getIdx()]) {
      setDirectionFromNeighboringBond(*firstFromAtom1, isFirstFromAtom1Flipped,
                                      *secondFromAtom1,
                                      isSecondFromAtom1Flipped);
    }
    bondDirCounts[secondFromAtom1->getIdx()] += 1;
    atomDirCounts[atom1->getIdx()] += 1;
  }

  if (atom2->getDegree() == 3 && secondFromAtom2) {
    if (!bondDirCounts[secondFromAtom2->getIdx()]) {
      setDirectionFromNeighboringBond(*firstFromAtom2, isFirstFromAtom2Flipped,
                                      *secondFromAtom2,
                                      isSecondFromAtom2Flipped);
    }
    bondDirCounts[secondFromAtom2->getIdx()] += 1;
    atomDirCounts[atom2->getIdx()] += 1;
  }
}

void canonicalizeDoubleBonds(ROMol &mol, const UINT_VECT &bondVisitOrders,
                             const UINT_VECT &atomVisitOrders,
                             std::vector<int8_t> &bondDirCounts,
                             std::vector<int8_t> &atomDirCounts,
                             const MolStack &molStack) {
  // start by removing the current directions on single bonds
  // around double bonds. At the same time, we build a prioritized
  // queue to decide the order in which we will canonicalize bonds.

  // We want to start with bonds with the most neighboring stereo
  // bonds, and in case of ties, start with the bond that has
  // the lowest position in the molStack

  auto getNeighboringStereoBond = [&mol](const Atom *dblBndAtom,
                                         const Bond *nbrBnd) -> Bond * {
    auto otherAtom = nbrBnd->getOtherAtom(dblBndAtom);
    for (const auto bond : mol.atomBonds(otherAtom)) {
      if (bond != nbrBnd && bond->getBondType() == Bond::DOUBLE &&
          bond->getStereo() > Bond::STEREOANY) {
        return bond;
      }
    }
    return nullptr;
  };

  std::greater<const unsigned int &> molStackComparer;
  std::less<const unsigned int &> numStereoNbrsComparer;

  std::unordered_map<const Bond *, std::vector<Bond *>> stereoBondNbrs;
  auto compareBondPriority = [&stereoBondNbrs, &bondVisitOrders,
                              &molStackComparer, &numStereoNbrsComparer](
                                 const Bond *aBnd, const Bond *bBnd) {
    const auto aNumStereoNbrs = stereoBondNbrs[aBnd].size();
    const auto bNumStereoNbrs = stereoBondNbrs[bBnd].size();

    if (aNumStereoNbrs == bNumStereoNbrs) {
      return molStackComparer(bondVisitOrders[aBnd->getIdx()],
                              bondVisitOrders[bBnd->getIdx()]);
    }
    return numStereoNbrsComparer(aNumStereoNbrs, bNumStereoNbrs);
  };

  std::priority_queue<Bond *, std::vector<Bond *>,
                      decltype(compareBondPriority)>
      q{compareBondPriority};

  for (auto &msI : molStack) {
    if (msI.type != MOL_STACK_BOND) {
      // not a bond, skip it
      continue;
    }

    auto bond = msI.obj.bond;
    Bond::BondDir dir = bond->getBondDir();
    if (dir == Bond::ENDDOWNRIGHT || dir == Bond::ENDUPRIGHT) {
      bond->setBondDir(Bond::NONE);
    }

    if (bond->getBondType() != Bond::DOUBLE ||
        bond->getStereo() <= Bond::STEREOANY ||
        bond->getStereoAtoms().size() < 2) {
      // not a bond that can have stereo or that needs canonicalization
      bond->setStereo(Bond::STEREONONE);
      continue;
    }

    auto &currentNbrs = stereoBondNbrs[bond];
    for (const auto *dblBondAtom : {bond->getBeginAtom(), bond->getEndAtom()}) {
      for (const auto *nbrBond : mol.atomBonds(dblBondAtom)) {
        if (!canHaveDirection(*nbrBond)) {
          continue;
        }
        auto nbrDblBnd = getNeighboringStereoBond(dblBondAtom, nbrBond);
        if (nbrDblBnd != nullptr) {
          currentNbrs.push_back(nbrDblBnd);
        }
      }
    }
    std::ranges::sort(currentNbrs, [&molStackComparer, &bondVisitOrders](
                                       const Bond *aBnd, const Bond *bBnd) {
      // Reversing the bonds is intentional: molStackComparer
      // is a std::greater comparer (priority queue returns
      // the highest element), but here we want to sort in
      // increasing order, so we want a std::less comparer,
      // which can be achieved by reversing the std::greater
      // because we can have no ties here
      return molStackComparer(bondVisitOrders[bBnd->getIdx()],
                              bondVisitOrders[aBnd->getIdx()]);
    });

    q.emplace(bond);
  }

  // Now that we have bonds in the order we want to handle them,
  // do the canonicalization
  std::vector<bool> seen_bonds(mol.getNumBonds());
  while (!q.empty()) {
    const auto bond = q.top();
    q.pop();
    if (seen_bonds[bond->getIdx()]) {
      continue;
    }

    std::queue<Bond *> connectedBondsQ;
    connectedBondsQ.push(bond);

    while (!connectedBondsQ.empty()) {
      const auto currentBond = connectedBondsQ.front();
      connectedBondsQ.pop();
      if (seen_bonds[currentBond->getIdx()]) {
        continue;
      }

      Canon::canonicalizeDoubleBond(currentBond, bondVisitOrders,
                                    atomVisitOrders, bondDirCounts,
                                    atomDirCounts);
      seen_bonds[currentBond->getIdx()] = true;
      for (auto nbrStereoBnd : stereoBondNbrs[currentBond]) {
        if (!seen_bonds[nbrStereoBnd->getIdx()]) {
          connectedBondsQ.push(nbrStereoBnd);
        }
      }
    }
  }
}

// finds cycles
void dfsFindCycles(ROMol &mol, int atomIdx, int inBondIdx,
                   std::vector<AtomColors> &colors, const UINT_VECT &ranks,
                   VECT_INT_VECT &atomRingClosures,
                   const boost::dynamic_bitset<> *bondsInPlay,
                   const std::vector<std::string> *bondSymbols, bool doRandom) {
  Atom *atom = mol.getAtomWithIdx(atomIdx);

  colors[atomIdx] = GREY_NODE;

  // ---------------------
  //
  //  Build the list of possible destinations from here
  //
  // ---------------------
  std::vector<PossibleType> possibles;
  auto bondsPair = mol.getAtomBonds(atom);
  possibles.reserve(bondsPair.second - bondsPair.first);

  while (bondsPair.first != bondsPair.second) {
    Bond *theBond = mol[*(bondsPair.first)];
    ++bondsPair.first;
    if (bondsInPlay && !(*bondsInPlay)[theBond->getIdx()]) {
      continue;
    }
    if (inBondIdx < 0 ||
        theBond->getIdx() != static_cast<unsigned int>(inBondIdx)) {
      int otherIdx = theBond->getOtherAtomIdx(atomIdx);
      auto rank = ranks[otherIdx];
      // ---------------------
      //
      // things are a bit more complicated if we are sitting on a
      // ring atom. we would like to traverse first to the
      // ring-closure atoms, then to atoms outside the ring first,
      // then to atoms in the ring that haven't already been visited
      // (non-ring-closure atoms).
      //
      //  Here's how the black magic works:
      //   - non-ring atom neighbors have their original ranks
      //   - ring atom neighbors have this added to their ranks:
      //       (MAX_BONDTYPE - bondOrder)*MAX_NATOMS*MAX_NATOMS
      //   - ring-closure neighbors lose a factor of:
      //       (MAX_BONDTYPE+1)*MAX_NATOMS*MAX_NATOMS
      //
      //  This tactic biases us to traverse to non-ring neighbors first,
      //  original ordering if bond orders are all equal... crafty, neh?
      //
      // ---------------------
      if (!doRandom) {
        if (colors[otherIdx] == GREY_NODE) {
          rank -= static_cast<int>(MAX_BONDTYPE + 1) * MAX_NATOMS * MAX_NATOMS;
          if (!bondSymbols) {
            rank += static_cast<int>(MAX_BONDTYPE - theBond->getBondType()) *
                    MAX_NATOMS;
          } else {
            const std::string &symb = (*bondSymbols)[theBond->getIdx()];
            std::uint32_t hsh = gboost::hash_range(symb.begin(), symb.end());
            rank += (hsh % MAX_NATOMS) * MAX_NATOMS;
          }
        } else if (theBond->getOwningMol().getRingInfo()->numBondRings(
                       theBond->getIdx())) {
          if (!bondSymbols) {
            rank += static_cast<int>(MAX_BONDTYPE - theBond->getBondType()) *
                    MAX_NATOMS * MAX_NATOMS;
          } else {
            const std::string &symb = (*bondSymbols)[theBond->getIdx()];
            std::uint32_t hsh = gboost::hash_range(symb.begin(), symb.end());
            rank += (hsh % MAX_NATOMS) * MAX_NATOMS * MAX_NATOMS;
          }
        }
      } else {
        // randomize the rank
        rank = getRandomGenerator()();
      }
      // std::cerr << "            " << atomIdx << ": " << otherIdx << " " <<
      // rank
      //           << std::endl;
      // std::cerr<<"aIdx: "<< atomIdx <<"   p: "<<otherIdx<<" Rank:
      // "<<ranks[otherIdx] <<" "<<colors[otherIdx]<<"
      // "<<theBond->getBondType()<<" "<<rank<<std::endl;
      possibles.emplace_back(rank, otherIdx, theBond);
    }
  }

  // ---------------------
  //
  //  Sort on ranks
  //
  // ---------------------
  std::sort(possibles.begin(), possibles.end(), _possibleCompare);
  // if (possibles.size())
  //   std::cerr << " aIdx1: " << atomIdx
  //             << " first: " << possibles.front()std:std::get<0>() << " "
  //             << possibles.front()std:std::get<1>() << std::endl;
  // // ---------------------
  //
  //  Now work the children
  //
  // ---------------------
  for (auto &possible : possibles) {
    int possibleIdx = std::get<1>(possible);
    Bond *bond = std::get<2>(possible);
    switch (colors[possibleIdx]) {
      case WHITE_NODE:
        // -----
        // we haven't seen this node at all before, traverse
        // -----
        dfsFindCycles(mol, possibleIdx, bond->getIdx(), colors, ranks,
                      atomRingClosures, bondsInPlay, bondSymbols, doRandom);
        break;
      case GREY_NODE:
        // -----
        // we've seen this, but haven't finished it (we're finishing a ring)
        // -----
        atomRingClosures[possibleIdx].push_back(bond->getIdx());
        atomRingClosures[atomIdx].push_back(bond->getIdx());
        break;
      default:
        // -----
        // this node has been finished. don't do anything.
        // -----
        break;
    }
  }
  colors[atomIdx] = BLACK_NODE;
}  // namespace Canon

void dfsBuildStack(ROMol &mol, int atomIdx, int inBondIdx,
                   std::vector<AtomColors> &colors, const UINT_VECT &ranks,
                   boost::dynamic_bitset<> &cyclesAvailable, MolStack &molStack,
                   VECT_INT_VECT &atomRingClosures,
                   std::vector<INT_LIST> &atomTraversalBondOrder,
                   const boost::dynamic_bitset<> *bondsInPlay,
                   const std::vector<std::string> *bondSymbols, bool doRandom) {
  Atom *atom = mol.getAtomWithIdx(atomIdx);
  boost::dynamic_bitset<> seenFromHere(mol.getNumAtoms());

  seenFromHere.set(atomIdx);
  molStack.push_back(MolStackElem(atom));
  colors[atomIdx] = GREY_NODE;

  INT_LIST travList;
  if (inBondIdx >= 0) {
    travList.push_back(inBondIdx);
  }

  // ---------------------
  //
  //  Add any ring closures
  //
  // ---------------------
  if (!atomRingClosures[atomIdx].empty()) {
    std::vector<unsigned int> ringsClosed;
    for (auto bIdx : atomRingClosures[atomIdx]) {
      travList.push_back(bIdx);
      Bond *bond = mol.getBondWithIdx(bIdx);
      seenFromHere.set(bond->getOtherAtomIdx(atomIdx));
      unsigned int ringIdx = std::numeric_limits<unsigned int>::max();
      if (bond->getPropIfPresent(common_properties::_TraversalRingClosureBond,
                                 ringIdx)) {
        // this is end of the ring closure
        // we can just pull the ring index from the bond itself:
        molStack.push_back(MolStackElem(bond, atomIdx));
        molStack.push_back(MolStackElem(ringIdx));
        // don't make the ring digit immediately available again: we don't want
        // to have the same
        // ring digit opening and closing rings on an atom.
        ringsClosed.push_back(ringIdx - 1);
      } else {
        // this is the beginning of the ring closure, we need to come up with a
        // ring index:
        auto lowestRingIdx = cyclesAvailable.find_first();
        if (lowestRingIdx == boost::dynamic_bitset<>::npos) {
          throw ValueErrorException(
              "Too many rings open at once. SMILES cannot be generated.");
        }
        cyclesAvailable.set(lowestRingIdx, false);
        ++lowestRingIdx;
        bond->setProp(common_properties::_TraversalRingClosureBond,
                      static_cast<unsigned int>(lowestRingIdx));
        molStack.push_back(MolStackElem(lowestRingIdx));
      }
    }
    for (auto ringIdx : ringsClosed) {
      cyclesAvailable.set(ringIdx);
    }
  }

  // ---------------------
  //
  //  Build the list of possible destinations from here
  //
  // ---------------------
  std::vector<PossibleType> possibles;
  possibles.reserve(atom->getDegree());
  for (auto theBond : mol.atomBonds(atom)) {
    if (bondsInPlay && !(*bondsInPlay)[theBond->getIdx()]) {
      continue;
    }
    if (inBondIdx < 0 ||
        theBond->getIdx() != static_cast<unsigned int>(inBondIdx)) {
      int otherIdx = theBond->getOtherAtomIdx(atomIdx);
      // ---------------------
      //
      // This time we skip the ring-closure atoms (we did them
      // above); we want to traverse first to atoms outside the ring
      // then to atoms in the ring that haven't already been visited
      // (non-ring-closure atoms).
      //
      // otherwise it's the same ranking logic as above
      //
      // ---------------------
      if (colors[otherIdx] != WHITE_NODE || seenFromHere[otherIdx]) {
        // ring closure or finished atom... skip it.
        continue;
      }
      auto rank = ranks[otherIdx];
      if (!doRandom) {
        if (theBond->getOwningMol().getRingInfo()->numBondRings(
                theBond->getIdx())) {
          if (!bondSymbols) {
            rank += static_cast<int>(MAX_BONDTYPE - theBond->getBondType()) *
                    MAX_NATOMS * MAX_NATOMS;
          } else {
            const std::string &symb = (*bondSymbols)[theBond->getIdx()];
            std::uint32_t hsh = gboost::hash_range(symb.begin(), symb.end());
            rank += (hsh % MAX_NATOMS) * MAX_NATOMS * MAX_NATOMS;
          }
        }
      } else {
        // randomize the rank
        rank = getRandomGenerator()();
      }

      possibles.emplace_back(rank, otherIdx, theBond);
    }
  }

  // ---------------------
  //
  //  Sort on ranks
  //
  // ---------------------
  std::sort(possibles.begin(), possibles.end(), _possibleCompare);
  // if (possibles.size())
  //   std::cerr << " aIdx2: " << atomIdx
  //             << " first: " << possibles.front()std:std::get<0>() << " "
  //             << possibles.front()std:std::get<1>() << std::endl;

  // ---------------------
  //
  //  Now work the children
  //
  // ---------------------
  for (auto possiblesIt = possibles.begin(); possiblesIt != possibles.end();
       ++possiblesIt) {
    int possibleIdx = std::get<1>(*possiblesIt);
    if (colors[possibleIdx] != WHITE_NODE) {
      // we're either done or it's a ring-closure, which we already processed...
      // this test isn't strictly required, because we only added WHITE notes to
      // the possibles list, but it seems logical to document it
      continue;
    }
    Bond *bond = std::get<2>(*possiblesIt);
    Atom *otherAtom = mol.getAtomWithIdx(possibleIdx);
    // ww might have some residual data from earlier calls, clean that up:
    otherAtom->clearProp(common_properties::_TraversalBondIndexOrder);
    travList.push_back(bond->getIdx());
    if (possiblesIt + 1 != possibles.end()) {
      // we're branching
      molStack.push_back(
          MolStackElem("(", rdcast<int>(possiblesIt - possibles.begin())));
    }
    molStack.push_back(MolStackElem(bond, atomIdx));
    dfsBuildStack(mol, possibleIdx, bond->getIdx(), colors, ranks,
                  cyclesAvailable, molStack, atomRingClosures,
                  atomTraversalBondOrder, bondsInPlay, bondSymbols, doRandom);
    if (possiblesIt + 1 != possibles.end()) {
      molStack.push_back(
          MolStackElem(")", rdcast<int>(possiblesIt - possibles.begin())));
    }
  }

  atomTraversalBondOrder[atom->getIdx()] = travList;
  colors[atomIdx] = BLACK_NODE;
}

void canonicalDFSTraversal(ROMol &mol, int atomIdx, int inBondIdx,
                           std::vector<AtomColors> &colors,
                           const UINT_VECT &ranks, MolStack &molStack,
                           VECT_INT_VECT &atomRingClosures,
                           std::vector<INT_LIST> &atomTraversalBondOrder,
                           const boost::dynamic_bitset<> *bondsInPlay,
                           const std::vector<std::string> *bondSymbols,
                           bool doRandom) {
  PRECONDITION(colors.size() >= mol.getNumAtoms(), "vector too small");
  PRECONDITION(ranks.size() >= mol.getNumAtoms(), "vector too small");
  PRECONDITION(atomRingClosures.size() >= mol.getNumAtoms(),
               "vector too small");
  PRECONDITION(atomTraversalBondOrder.size() >= mol.getNumAtoms(),
               "vector too small");
  PRECONDITION(!bondsInPlay || bondsInPlay->size() >= mol.getNumBonds(),
               "bondsInPlay too small");
  PRECONDITION(!bondSymbols || bondSymbols->size() >= mol.getNumBonds(),
               "bondSymbols too small");

  std::vector<AtomColors> tcolors(colors.begin(), colors.end());
  dfsFindCycles(mol, atomIdx, inBondIdx, tcolors, ranks, atomRingClosures,
                bondsInPlay, bondSymbols, doRandom);

  boost::dynamic_bitset<> cyclesAvailable(MAX_CYCLES);
  cyclesAvailable.set();
  dfsBuildStack(mol, atomIdx, inBondIdx, colors, ranks, cyclesAvailable,
                molStack, atomRingClosures, atomTraversalBondOrder, bondsInPlay,
                bondSymbols, doRandom);
}

void clearBondDirs(ROMol &mol, Bond *refBond, const Atom *fromAtom,
                   std::vector<int8_t> &bondDirCounts,
                   std::vector<int8_t> &atomDirCounts) {
  PRECONDITION(bondDirCounts.size() >= mol.getNumBonds(), "bad dirCount size");
  PRECONDITION(refBond, "bad bond");
  PRECONDITION(&refBond->getOwningMol() == &mol, "bad bond");
  PRECONDITION(fromAtom, "bad atom");
  PRECONDITION(&fromAtom->getOwningMol() == &mol, "bad bond");

  auto clearDirection = [&atomDirCounts, &bondDirCounts](Bond *bond) {
    --bondDirCounts[bond->getIdx()];
    if (!bondDirCounts[bond->getIdx()]) {
      bond->setBondDir(Bond::NONE);
      --atomDirCounts[bond->getBeginAtomIdx()];
      --atomDirCounts[bond->getEndAtomIdx()];
    }
  };

  for (auto oBond : mol.atomBonds(fromAtom)) {
    if (oBond != refBond && canHaveDirection(*oBond)) {
      if ((bondDirCounts[oBond->getIdx()] >=
           bondDirCounts[refBond->getIdx()]) &&
          atomDirCounts[oBond->getBeginAtomIdx()] != 1 &&
          atomDirCounts[oBond->getEndAtomIdx()] != 1) {
        clearDirection(oBond);
      } else if (atomDirCounts[refBond->getBeginAtomIdx()] != 1 &&
                 atomDirCounts[refBond->getEndAtomIdx()] != 1) {
        // we found a neighbor that could have directionality set,
        // but it had a lower bondDirCount than us, so we must
        // need to be adjusted:
        clearDirection(refBond);
      }
      break;
    }
  }
}

// CanonicalizeDoubleBonds tries to add as many directions as possible
// to stereo double bonds, but some of these may coerce STEREONONE or
// STEREOANY into stereo just because they are "in the wrong place",
// in the middle of direction bonds of neighboring stereo bonds. Here
// we try to fix some of them (we probably can't fix all) before
// removing redundant ones in removeRedundantBondDirSpecs.
void removeUnwantedBondDirSpecs(ROMol &mol, MolStack &molStack,
                                std::vector<int8_t> &bondDirCounts,
                                std::vector<int8_t> &atomDirCounts,
                                std::vector<unsigned int> &bondVisitOrders) {
  PRECONDITION(bondDirCounts.size() >= mol.getNumBonds(), "bad dirCount size");

  for (auto &msI : molStack) {
    if (msI.type != MOL_STACK_BOND) {
      continue;
    }

    if (msI.obj.bond->getBondType() != Bond::DOUBLE ||
        msI.obj.bond->getStereo() > Bond::STEREOANY) {
      continue;
    }

    auto firstAtom = msI.obj.bond->getBeginAtom();
    auto secondAtom = msI.obj.bond->getEndAtom();
    if (firstAtom->getDegree() == 1 || secondAtom->getDegree() == 1) {
      // One side of the bond does not have any neighbors. There's no way for
      // this double bond to have stereo!
      continue;
    }

    std::vector<Bond *> removalCandidates;

    // Look at the first side of the non-stereo double bond

    for (auto bond : mol.atomBonds(firstAtom)) {
      if (bondDirCounts[bond->getIdx()]) {
        removalCandidates.push_back(bond);
      }
    }
    if (removalCandidates.empty()) {
      // No bonds with direction on this side, so this non-stereo
      // bond won't be coerced into stereo.
      continue;
    }

    if (atomDirCounts[firstAtom->getIdx()]) {
      // We only keep atomDirCounts for atoms at the end of a stereo double
      // bond. This means that if an end of this non-stereo double bond has
      // a dir count, then both bonds have this atom in common (like the
      // two bonds in S=P(=N\C)/C have P in common), and we can't remove
      // the direction from that side, as it will remove stereo from the
      // stereo bond too.
      removalCandidates.clear();
    }

    // Now look at the other side

    uint8_t candidatesOnSecondEnd = 0;
    for (auto bond : mol.atomBonds(secondAtom)) {
      if (bondDirCounts[bond->getIdx()]) {
        removalCandidates.push_back(bond);
        ++candidatesOnSecondEnd;
      }
    }

    if (candidatesOnSecondEnd == 0) {
      // No bonds with direction on this side, so this non-stereo
      // bond won't be coerced into stereo.
      continue;
    }

    if (atomDirCounts[secondAtom->getIdx()]) {
      // If we got here, and can't remove bonds on this side, this probably
      // means there's nothing we can do, and this bond will be coerced
      // into stereo.
      continue;
    }

    // Sort by position in the molStack, prefer the bond closest to the start
    std::ranges::sort(
        removalCandidates, [&bondVisitOrders](const auto &a, const auto &b) {
          return bondVisitOrders[a->getIdx()] < bondVisitOrders[b->getIdx()];
        });

    for (auto candidateBond : removalCandidates) {
      Atom *otherAtom = nullptr;
      if (candidateBond->getBeginAtom() == firstAtom ||
          candidateBond->getEndAtom() == firstAtom) {
        otherAtom = candidateBond->getOtherAtom(firstAtom);
      } else if (candidateBond->getBeginAtom() == secondAtom ||
                 candidateBond->getEndAtom() == secondAtom) {
        otherAtom = candidateBond->getOtherAtom(secondAtom);
      } else {
        CHECK_INVARIANT(false, "inconsistent bond ends");
      }

      // to be able to remove the bond, the "other end", the atom that
      // is part of a stereo double bond, must have 2 directions, so that
      // that bond keeps stereo even if we remove one of the directions.
      if (atomDirCounts[otherAtom->getIdx()] == 2) {
        bondDirCounts[candidateBond->getIdx()] = 0;
        candidateBond->setBondDir(Bond::NONE);
        atomDirCounts[otherAtom->getIdx()] -= 1;
        break;
      }
    }
  }
}

void removeRedundantBondDirSpecs(ROMol &mol, MolStack &molStack,
                                 std::vector<int8_t> &bondDirCounts,
                                 std::vector<int8_t> &atomDirCounts) {
  PRECONDITION(bondDirCounts.size() >= mol.getNumBonds(), "bad dirCount size");

  auto clearBondDirsFromAtom = [&mol, &bondDirCounts, &atomDirCounts](
                                   Bond *tBond, const Atom *atom) {
    for (auto bond : mol.atomBonds(atom)) {
      if (bond != tBond && bond->getBondType() == Bond::DOUBLE &&
          bond->getStereo() > Bond::STEREOANY) {
        clearBondDirs(mol, tBond, atom, bondDirCounts, atomDirCounts);
        return;
      }
    }
  };

  // find bonds that have directions indicated that are redundant:
  for (auto &msI : molStack) {
    if (msI.type != MOL_STACK_BOND) {
      continue;
    }
    Bond *tBond = msI.obj.bond;
    const Atom *canonBeginAtom = mol.getAtomWithIdx(msI.number);
    const Atom *canonEndAtom =
        mol.getAtomWithIdx(tBond->getOtherAtomIdx(msI.number));
    if (canHaveDirection(*tBond) && bondDirCounts[tBond->getIdx()]) {
      clearBondDirsFromAtom(tBond, canonBeginAtom);
      clearBondDirsFromAtom(tBond, canonEndAtom);
    } else if (tBond->getBondDir() != Bond::NONE) {
      // we aren't supposed to have a direction set, but we do:
      tBond->setBondDir(Bond::NONE);
    }
  }
}

// insert (-1) for hydrogens or missing ligands, where these are placed
// depends on if it is the first atom or not
static void insertImplicitNbors(INT_LIST &bonds, const Atom::ChiralType tag,
                                const bool firstAtom) {
  unsigned int ref_max = Chirality::getMaxNbors(tag);
  if (bonds.size() < ref_max) {
    if (firstAtom) {
      bonds.insert(bonds.begin(), ref_max - bonds.size(), -1);
    } else {
      bonds.insert(++bonds.begin(), ref_max - bonds.size(), -1);
    }
  }
}

void canonicalizeFragment(ROMol &mol, int atomIdx,
                          std::vector<AtomColors> &colors,
                          const UINT_VECT &ranks, MolStack &molStack,
                          const boost::dynamic_bitset<> *bondsInPlay,
                          const std::vector<std::string> *bondSymbols,
                          bool doIsomericSmiles, bool doRandom,
                          bool doChiralInversions) {
  boost::dynamic_bitset<> atomsInPlay(mol.getNumAtoms());
  if (!bondsInPlay) {
    // if we weren't given a bondsInPlay, then all bonds are in play, so we need
    // to set both those and the atomsInPlay here:
    atomsInPlay.set();
  } else {
    for (const auto bnd : mol.bonds()) {
      if ((*bondsInPlay)[bnd->getIdx()]) {
        atomsInPlay.set(bnd->getBeginAtomIdx());
        atomsInPlay.set(bnd->getEndAtomIdx());
      }
    }
  }
  canonicalizeFragment(mol, atomIdx, colors, ranks, molStack, &atomsInPlay,
                       bondsInPlay, bondSymbols, doIsomericSmiles, doRandom,
                       doChiralInversions);
}
RDKIT_GRAPHMOL_EXPORT void canonicalizeFragment(
    ROMol &mol, int atomIdx, std::vector<AtomColors> &colors,
    const std::vector<unsigned int> &ranks, MolStack &molStack,
    const boost::dynamic_bitset<> *atomsInPlay,
    const boost::dynamic_bitset<> *bondsInPlay,
    const std::vector<std::string> *bondSymbols, bool doIsomericSmiles,
    bool doRandom, bool doChiralInversions) {
  PRECONDITION(colors.size() >= mol.getNumAtoms(), "vector too small");
  PRECONDITION(ranks.size() >= mol.getNumAtoms(), "vector too small");
  PRECONDITION(!atomsInPlay || atomsInPlay->size() >= mol.getNumAtoms(),
               "atomsInPlay too small");
  PRECONDITION(!bondsInPlay || bondsInPlay->size() >= mol.getNumBonds(),
               "bondsInPlay too small");
  PRECONDITION(!bondSymbols || bondSymbols->size() >= mol.getNumBonds(),
               "bondSymbols too small");
  unsigned int nAtoms = mol.getNumAtoms();

  // make sure that we've done the stereo perception:
  if (!mol.hasProp(common_properties::_StereochemDone)) {
    MolOps::assignStereochemistry(mol, false);
  }

  // we need ring information; make sure findSSSR has been called before
  // if not call now
  // NOTE: if called from the SMARTS code, the ring info will be set to SSSR,
  // but no ring infor in actually set
  if (!mol.getRingInfo()->isSymmSssr()) {
    MolOps::findSSSR(mol);
  }
  mol.getAtomWithIdx(atomIdx)->setProp(common_properties::_TraversalStartPoint,
                                       true);

  VECT_INT_VECT atomRingClosures(nAtoms);
  std::vector<INT_LIST> atomTraversalBondOrder(nAtoms);
  Canon::canonicalDFSTraversal(mol, atomIdx, -1, colors, ranks, molStack,
                               atomRingClosures, atomTraversalBondOrder,
                               bondsInPlay, bondSymbols, doRandom);

  CHECK_INVARIANT(!molStack.empty(), "Empty stack.");
  CHECK_INVARIANT(molStack.begin()->type == MOL_STACK_ATOM,
                  "Corrupted stack. First element should be an atom.");

  // collect some information about traversal order on chiral atoms
  boost::dynamic_bitset<> numSwapsChiralAtoms(nAtoms);
  std::vector<int> atomPermutationIndices(nAtoms, 0);
  if (doIsomericSmiles) {
    for (const auto atom : mol.atoms()) {
      if (atomsInPlay && !(*atomsInPlay)[atom->getIdx()]) {
        continue;
      }
      if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
        // check if all of this atom's bonds are in play
        for (const auto bnd : mol.atomBonds(atom)) {
          if (bondsInPlay && !(*bondsInPlay)[bnd->getIdx()]) {
            atom->setProp(common_properties::_brokenChirality, true);
            break;
          }
        }
        if (atom->hasProp(common_properties::_brokenChirality)) {
          continue;
        }

        // Extra check needed if/when @AL1/@AL2 supported
        if (Chirality::detail::isAtomPotentialTetrahedralCenter(atom) ||
            Chirality::hasNonTetrahedralStereo(atom)) {
          int perm = 0;
          if (Chirality::hasNonTetrahedralStereo(atom)) {
            atom->getPropIfPresent(common_properties::_chiralPermutation, perm);
          }

          const unsigned int firstIdx = molStack.begin()->obj.atom->getIdx();
          const bool firstInPart = atom->getIdx() == firstIdx;

          // Check if the atom can be chiral, and if chirality needs inversion
          const INT_LIST &trueOrder = atomTraversalBondOrder[atom->getIdx()];

          // We have to make sure that trueOrder contains all the
          // bonds, even if they won't be written to the SMILES
          int nSwaps = 0;
          if (trueOrder.size() < atom->getDegree()) {
            INT_LIST tOrder = trueOrder;
            for (const auto bnd : mol.atomBonds(atom)) {
              int bndIdx = bnd->getIdx();
              if (std::find(trueOrder.begin(), trueOrder.end(), bndIdx) ==
                  trueOrder.end()) {
                tOrder.push_back(bndIdx);
              }
            }
            if (!perm) {
              nSwaps = atom->getPerturbationOrder(tOrder);
            } else {
              insertImplicitNbors(tOrder, atom->getChiralTag(), firstInPart);
              perm = Chirality::getChiralPermutation(atom, tOrder);
            }
          } else {
            if (!perm) {
              nSwaps = atom->getPerturbationOrder(trueOrder);
            } else {
              INT_LIST tOrder = trueOrder;
              insertImplicitNbors(tOrder, atom->getChiralTag(), firstInPart);
              perm = Chirality::getChiralPermutation(atom, tOrder);
            }
          }

          // in future this should be moved up and simplified, there should not
          // be an option to not do chiral inversions
          if (doChiralInversions &&
              chiralAtomNeedsTagInversion(
                  mol, atom, firstInPart,
                  atomRingClosures[atom->getIdx()].size())) {
            // This is a special case. Here's an example:
            //   Our internal representation of a chiral center is equivalent
            //   to:
            //     [C@](F)(O)(C)[H]
            //   we'll be dumping it without the H, which entails a
            //   reordering:
            //     [C@@H](F)(O)C
            ++nSwaps;
          }
          if (nSwaps % 2) {
            numSwapsChiralAtoms.set(atom->getIdx());
          }
          atomPermutationIndices[atom->getIdx()] = perm;
        }
      }
    }
  }

  std::vector<unsigned int> atomVisitOrders(mol.getNumAtoms());
  std::vector<unsigned int> bondVisitOrders(mol.getNumBonds());

  unsigned int pos = 0;
  for (const auto &msI : molStack) {
    if (msI.type == MOL_STACK_ATOM) {
      atomVisitOrders[msI.obj.atom->getIdx()] = pos;
    } else if (msI.type == MOL_STACK_BOND) {
      bondVisitOrders[msI.obj.bond->getIdx()] = pos;
      auto dir = msI.obj.bond->getBondDir();
      if (dir == Bond::ENDDOWNRIGHT || dir == Bond::ENDUPRIGHT) {
        msI.obj.bond->setBondDir(Bond::NONE);
      }
    }
    ++pos;
  }

  std::vector<int8_t> bondDirCounts(mol.getNumBonds(), 0);
  std::vector<int8_t> atomDirCounts(nAtoms, 0);
  canonicalizeDoubleBonds(mol, bondVisitOrders, atomVisitOrders, bondDirCounts,
                          atomDirCounts, molStack);

  // traverse the stack and canonicalize atoms with (ring) stereochemistry
  if (doIsomericSmiles) {
    boost::dynamic_bitset<> ringStereoChemAdjusted(nAtoms);
    for (auto &msI : molStack) {
      if (msI.type == MOL_STACK_ATOM &&
          msI.obj.atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
          !msI.obj.atom->hasProp(common_properties::_brokenChirality)) {
        if (msI.obj.atom->hasProp(common_properties::_ringStereoAtoms)) {
          // FIX: handle stereogroups here too
          if (!ringStereoChemAdjusted[msI.obj.atom->getIdx()]) {
            msI.obj.atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
            ringStereoChemAdjusted.set(msI.obj.atom->getIdx());
          }
          const INT_VECT &ringStereoAtoms = msI.obj.atom->getProp<INT_VECT>(
              common_properties::_ringStereoAtoms);
          for (auto nbrV : ringStereoAtoms) {
            int nbrIdx = abs(nbrV) - 1;
            // Adjust the chirality flag of the ring stereo atoms according to
            // the first one
            if (!ringStereoChemAdjusted[nbrIdx] &&
                atomVisitOrders[nbrIdx] >
                    atomVisitOrders[msI.obj.atom->getIdx()]) {
              mol.getAtomWithIdx(nbrIdx)->setChiralTag(
                  msI.obj.atom->getChiralTag());
              if (nbrV < 0) {
                mol.getAtomWithIdx(nbrIdx)->invertChirality();
              }
              // Odd number of swaps for first chiral ring atom --> needs to be
              // swapped but we want to retain chirality
              if (numSwapsChiralAtoms[msI.obj.atom->getIdx()]) {
                // Odd number of swaps for chiral ring neighbor --> needs to be
                // swapped but we want to retain chirality
                if (!numSwapsChiralAtoms[nbrIdx]) {
                  mol.getAtomWithIdx(nbrIdx)->invertChirality();
                }
              }
              // Even number of swaps for first chiral ring atom --> don't need
              // to be swapped
              else {
                // Odd number of swaps for chiral ring neighbor --> needs to be
                // swapped
                if (numSwapsChiralAtoms[nbrIdx]) {
                  mol.getAtomWithIdx(nbrIdx)->invertChirality();
                }
              }
              ringStereoChemAdjusted.set(nbrIdx);
            }
          }
        } else if (size_t sgidx;
                   msI.obj.atom->getPropIfPresent("_stereoGroup", sgidx) &&
                   mol.getStereoGroups().size() > sgidx) {
          // make sure that the reference atom in the stereogroup is CCW
          auto &sg = mol.getStereoGroups()[sgidx];
          bool swapIt =
              msI.obj.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW;
          if (swapIt) {
            msI.obj.atom->invertChirality();
          }
          if (swapIt || numSwapsChiralAtoms[msI.obj.atom->getIdx()]) {
            for (auto at : sg.getAtoms()) {
              if (at == msI.obj.atom) {
                continue;
              }
              at->invertChirality();
            }
          }

        } else {
          if (msI.obj.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
              msI.obj.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
            if ((numSwapsChiralAtoms[msI.obj.atom->getIdx()])) {
              msI.obj.atom->invertChirality();
            }
          } else if (atomPermutationIndices[msI.obj.atom->getIdx()]) {
            msI.obj.atom->setProp(
                common_properties::_chiralPermutation,
                atomPermutationIndices[msI.obj.atom->getIdx()]);
          }
        }
      }
    }
  }
  Canon::removeUnwantedBondDirSpecs(mol, molStack, bondDirCounts, atomDirCounts,
                                    bondVisitOrders);

  Canon::removeRedundantBondDirSpecs(mol, molStack, bondDirCounts,
                                     atomDirCounts);
}

void canonicalizeEnhancedStereo(ROMol &mol,
                                const std::vector<unsigned int> *atomRanks) {
  const auto &sgs = mol.getStereoGroups();
  if (sgs.empty()) {
    return;
  }

  std::vector<unsigned int> lranks;
  if (!atomRanks) {
    bool breakTies = true;
    rankMolAtoms(mol, lranks, breakTies);
    atomRanks = &lranks;
  }
  // one thing that makes this all easier is that the stereogroups are
  // independent of each other
  std::vector<StereoGroup> newSgs;
  for (auto &sg : sgs) {
    // we don't do anything to ABS groups
    if (sg.getGroupType() == StereoGroupType::STEREO_ABSOLUTE) {
      newSgs.push_back(sg);
      continue;
    }

    // sort the atoms by rank:
    auto getAtomRank = [&atomRanks](const Atom *at1, const Atom *at2) {
      return atomRanks->at(at1->getIdx()) < atomRanks->at(at2->getIdx());
    };
    auto sgAtoms = sg.getAtoms();
    std::sort(sgAtoms.begin(), sgAtoms.end(), getAtomRank);

    // sort the bonds by atom rank:
    auto getBondRank = [&atomRanks](const Bond *bd1, const Bond *bd2) {
      unsigned int bd1at1 = atomRanks->at(bd1->getBeginAtomIdx());
      unsigned int bd1at2 = atomRanks->at(bd1->getEndAtomIdx());
      unsigned int bd2at1 = atomRanks->at(bd2->getBeginAtomIdx());
      unsigned int bd2at2 = atomRanks->at(bd2->getEndAtomIdx());
      if (bd1at1 < bd1at2) {
        std::swap(bd1at1, bd1at2);
      }
      if (bd2at1 < bd2at2) {
        std::swap(bd2at1, bd2at2);
      }
      if (bd1at1 != bd2at1) {
        return bd1at1 < bd2at1;
      }
      return bd1at2 < bd2at2;
    };
    auto sgBonds = sg.getBonds();
    std::sort(sgBonds.begin(), sgBonds.end(), getBondRank);

    // find the reference (lowest-ranked) atom (or lowest-ranked bond)

    Atom::ChiralType foundRefState = Atom::ChiralType::CHI_TETRAHEDRAL_CCW;
    if (sgAtoms.size() > 0) {
      foundRefState = sgAtoms.front()->getChiralTag();
    } else if (sgBonds.size() > 0) {
      if (sgBonds.front()->getStereo() == Bond::BondStereo::STEREOATROPCCW) {
        foundRefState =
            Atom::ChiralType::CHI_TETRAHEDRAL_CCW;  // convert atropisomer CCW
                                                    // to atom CCW
      } else {
        foundRefState =
            Atom::ChiralType::CHI_TETRAHEDRAL_CW;  // convert atropisomer CW
                                                   // to atom CW
      }
    }
    // we will use CCW as the "canonical" state for chirality, so if the
    // referenceAtom is already CCW then we don't need to do anything more
    // with this stereogroup
    auto refState = Atom::ChiralType::CHI_TETRAHEDRAL_CCW;
    if (foundRefState != refState) {
      // we need to flip everyone... so loop over the other atoms and bonds
      // and flip them all:

      for (auto atom : sgAtoms) {
        atom->invertChirality();
      }
      for (auto bond : sgBonds) {
        bond->invertChirality();
      }
    }
    newSgs.emplace_back(
        StereoGroup(sg.getGroupType(), std::move(sgAtoms), std::move(sgBonds)));

    // note that we do not forward the Group Ids: this is intentional, so that
    // the Ids are reassigned based on the canonicalized order.
    if (sgAtoms.size() > 0) {
      sgAtoms.front()->setProp("_stereoGroup", newSgs.size() - 1, true);
    }
  }
  mol.setStereoGroups(newSgs);
}

void addSingleAbsGroup(ROMol &mol) {
  // all chiral centers are added to an abs group
  // if there are not chiral centers, no group is added

  std::vector<StereoGroup> sgs;
  std::vector<Atom *> chiralAtoms;
  std::vector<Bond *> chiralBonds;
  for (auto &atom : mol.atoms()) {
    if (atom->getChiralTag() == Atom::ChiralType::CHI_TETRAHEDRAL_CCW ||
        atom->getChiralTag() == Atom::ChiralType::CHI_TETRAHEDRAL_CW) {
      chiralAtoms.push_back(atom);
    }
  }
  for (auto &bond : mol.bonds()) {
    if (bond->getStereo() == Bond::BondStereo::STEREOATROPCW ||
        bond->getStereo() == Bond::BondStereo::STEREOATROPCCW) {
      chiralBonds.push_back(bond);
    }
  }

  if (!chiralAtoms.empty() || !chiralBonds.empty()) {
    sgs.emplace_back(StereoGroupType::STEREO_ABSOLUTE, chiralAtoms,
                     chiralBonds);
  }
  mol.setStereoGroups(sgs);  // could be empty, or have one abs group
}

void clearStereoGroups(ROMol &mol) {
  // all chiral centers are added to an abs group
  // if there are not chiral centers, no group is added
  std::vector<StereoGroup> sgs;
  mol.setStereoGroups(sgs);
}

}  // namespace Canon

}  // namespace RDKit
