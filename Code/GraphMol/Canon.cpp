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

Bond::BondDir getReferenceDirection(const Bond *dblBond, const Atom *refAtom,
                                    const Atom *targetAtom,
                                    const Bond *refControllingBond,
                                    bool refIsFlipped, const Bond *targetBond,
                                    bool targetIsFlipped) {
  Bond::BondDir dir = Bond::NONE;
  if (dblBond->getStereo() == Bond::STEREOE ||
      dblBond->getStereo() == Bond::STEREOTRANS) {
    dir = refControllingBond->getBondDir();
  } else if (dblBond->getStereo() == Bond::STEREOZ ||
             dblBond->getStereo() == Bond::STEREOCIS) {
    dir = flipStereoBondDir(refControllingBond->getBondDir());
  }
  CHECK_INVARIANT(dir != Bond::NONE, "stereo not set");

  // If we're not looking at the bonds used to determine the
  // stereochemistry, we need to flip the setting on the other bond:
  const INT_VECT &stereoAtoms = dblBond->getStereoAtoms();

  if (refAtom->getDegree() == 3 &&
      std::ranges::find(stereoAtoms,
                        static_cast<int>(refControllingBond->getOtherAtomIdx(
                            refAtom->getIdx()))) == stereoAtoms.end()) {
    dir = flipStereoBondDir(dir);
  }
  if (targetAtom->getDegree() == 3 &&
      std::ranges::find(stereoAtoms,
                        static_cast<int>(targetBond->getOtherAtomIdx(
                            targetAtom->getIdx()))) == stereoAtoms.end()) {
    dir = flipStereoBondDir(dir);
  }

  // XOR: flips will cancel out if we do both
  if (refIsFlipped != targetIsFlipped) {
    dir = flipStereoBondDir(dir);
  }

  return dir;
}
}  // namespace

namespace details {
bool isUnsaturated(const Atom *atom, const ROMol &mol) {
  for (const auto &bndItr :
       boost::make_iterator_range(mol.getAtomBonds(atom))) {
    // can't just check for single bonds, because dative bonds also have an
    // order of 1
    if (mol[bndItr]->getBondTypeAsDouble() > 1) {
      return true;
    }
  }
  return false;
}

bool hasSingleHQuery(const Atom::QUERYATOM_QUERY *q) {
  // list queries are series of nested ors of AtomAtomicNum queries
  PRECONDITION(q, "bad query");
  bool res = false;
  std::string descr = q->getDescription();
  if (descr == "AtomAnd") {
    for (auto cIt = q->beginChildren(); cIt != q->endChildren(); ++cIt) {
      auto cDescr = (*cIt)->getDescription();
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
                            UINT_VECT &bondDirCounts,
                            UINT_VECT &atomDirCounts) {
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

  if (dir1Set && dir2Set) {
    // Both directions are already set. Nothing to do.

    // To do: check that the directions are consistent with each other.

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
      if (secondFromAtom1) {
        // both bonds have their directionalities set, make sure
        // they are compatible:
        if (firstFromAtom1->getBondDir() == secondFromAtom1->getBondDir() &&
            bondDirCounts[firstFromAtom2->getIdx()]) {
          CHECK_INVARIANT(
              ((firstFromAtom1->getBeginAtomIdx() == atom1->getIdx()) ^
               (secondFromAtom1->getBeginAtomIdx() == atom1->getIdx())),
              "inconsistent state");
        }
      }
    } else {
      // the second bond must be present and setting the direction:
      CHECK_INVARIANT(secondFromAtom1, "inconsistent state");
      CHECK_INVARIANT(bondDirCounts[secondFromAtom1->getIdx()] > 0,
                      "inconsistent state");
      auto atom1Dir = secondFromAtom1->getBondDir();

      // By default, both bonds on the same side of the double bond
      // should have opposite directions, but this can change if
      // one (and only one) of the bonds is flipped (if both were
      // flipped, the flips in the direction would cancel out).
      if (isSecondFromAtom1Flipped == isFirstFromAtom1Flipped) {
        atom1Dir = flipStereoBondDir(atom1Dir);
      }

      firstFromAtom1->setBondDir(atom1Dir);

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
    } else {
      // the second bond must be present and setting the direction:
      CHECK_INVARIANT(secondFromAtom2, "inconsistent state");
      CHECK_INVARIANT(bondDirCounts[secondFromAtom2->getIdx()] > 0,
                      "inconsistent state");

      auto atom2Dir = secondFromAtom2->getBondDir();

      // By default, both bonds on the same side of the double bond
      // should have opposite directions, but this can change if
      // one (and only one) of the bonds is flipped (if both were
      // flipped, the flips in the direction would cancel out).
      if (isSecondFromAtom2Flipped == isFirstFromAtom2Flipped) {
        atom2Dir = flipStereoBondDir(atom2Dir);
      }

      firstFromAtom2->setBondDir(atom2Dir);

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
        dblBond, atom1, atom2, atom1ControllingBond, isControllingAtomFlipped,
        firstFromAtom2, isFirstFromAtom2Flipped);

    firstFromAtom2->setBondDir(atom2Dir);

    bondDirCounts[firstFromAtom2->getIdx()] += 1;
    atomDirCounts[atom2->getIdx()] += 1;
  } else {
    auto isControllingAtomFlipped =
        (atom2ControllingBond == firstFromAtom2 ? isFirstFromAtom2Flipped
                                                : isSecondFromAtom2Flipped);
    auto atom1Dir = getReferenceDirection(
        dblBond, atom2, atom1, atom2ControllingBond, isControllingAtomFlipped,
        firstFromAtom1, isFirstFromAtom1Flipped);

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
      auto otherDir = firstFromAtom1->getBondDir();

      // Again, bonds on the same side of the double bond default
      // to have opposite directions, but we need to consider if
      // they are flipped.
      if (isFirstFromAtom1Flipped == isSecondFromAtom1Flipped) {
        otherDir = flipStereoBondDir(otherDir);
      }

      secondFromAtom1->setBondDir(otherDir);
    }
    bondDirCounts[secondFromAtom1->getIdx()] += 1;
    atomDirCounts[atom1->getIdx()] += 1;
  }

  if (atom2->getDegree() == 3 && secondFromAtom2) {
    if (!bondDirCounts[secondFromAtom2->getIdx()]) {
      auto otherDir = firstFromAtom2->getBondDir();

      // Again, bonds on the same side of the double bond default
      // to have opposite directions, but we need to consider if
      // they are flipped.
      if (isFirstFromAtom2Flipped == isSecondFromAtom2Flipped) {
        otherDir = flipStereoBondDir(otherDir);
      }

      secondFromAtom2->setBondDir(otherDir);
    }
    bondDirCounts[secondFromAtom2->getIdx()] += 1;
    atomDirCounts[atom2->getIdx()] += 1;
  }
}

void canonicalizeDoubleBonds(ROMol &mol, const UINT_VECT &bondVisitOrders,
                             const UINT_VECT &atomVisitOrders,
                             UINT_VECT &bondDirCounts, UINT_VECT &atomDirCounts,
                             const MolStack &molStack) {
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

  // start by removing the current directions on single bonds
  // around double bonds. At the same time, we build a prioritized
  // queue to decide the order in which we will canonicalize bonds.

  // We want to start with bonds with the most neighboring stereo
  // bonds, and in case of ties, start with the bond that has
  // the lowest position in the molStack
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
      // not a bond that can have stereo nor needs canonicalization
      bond->setStereo(Bond::STEREONONE);
      continue;
    }

    auto &currentNbrs = stereoBondNbrs[bond];
    for (const auto *dblBondAtom :
         std::array<Atom *, 2>{bond->getBeginAtom(), bond->getEndAtom()}) {
      for (const auto *bond : mol.atomBonds(dblBondAtom)) {
        if (!canHaveDirection(*bond)) {
          continue;
        }
        auto nbrDblBnd = getNeighboringStereoBond(dblBondAtom, bond);
        if (nbrDblBnd != nullptr) {
          currentNbrs.push_back(nbrDblBnd);
        }
      }
    }
    std::ranges::sort(currentNbrs,
                      [&molStackComparer](const Bond *aBnd, const Bond *bBnd) {
                        // Reversing the bonds is intentional: molStackComparer
                        // is a std::greater comparer (priority queue returns
                        // the highest element), but here we want to sort in
                        // increasing order, so we want a std::less comparer,
                        // which can be achieved by reversing the std::greater
                        // because we can have no ties here
                        return molStackComparer(bBnd->getIdx(), aBnd->getIdx());
                      });

    q.emplace(bond);
  }

  // Now that we have bonds in the order we want to handle them,
  // do the canonicalization
  boost::dynamic_bitset<> seen_bonds(mol.getNumBonds());
  while (!q.empty()) {
    const auto bond = q.top();
    q.pop();
    if (seen_bonds.test(bond->getIdx())) {
      continue;
    }

    std::queue<Bond *> connectedBondsQ;
    connectedBondsQ.push(bond);

    while (!connectedBondsQ.empty()) {
      const auto currentBond = connectedBondsQ.front();
      connectedBondsQ.pop();

      Canon::canonicalizeDoubleBond(currentBond, bondVisitOrders,
                                    atomVisitOrders, bondDirCounts,
                                    atomDirCounts);
      seen_bonds.set(currentBond->getIdx());
      for (auto nbrStereoBnd : stereoBondNbrs[currentBond]) {
        if (!seen_bonds.test(nbrStereoBnd->getIdx())) {
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
  possibles.resize(0);
  ROMol::OBOND_ITER_PAIR bondsPair = mol.getAtomBonds(atom);
  possibles.reserve(bondsPair.second - bondsPair.first);

  while (bondsPair.first != bondsPair.second) {
    Bond *theBond = mol[*(bondsPair.first)];
    bondsPair.first++;
    if (bondsInPlay && !(*bondsInPlay)[theBond->getIdx()]) {
      continue;
    }
    if (inBondIdx < 0 ||
        theBond->getIdx() != static_cast<unsigned int>(inBondIdx)) {
      int otherIdx = theBond->getOtherAtomIdx(atomIdx);
      long rank = ranks[otherIdx];
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
                   std::vector<AtomColors> &colors, VECT_INT_VECT &cycles,
                   const UINT_VECT &ranks, UINT_VECT &cyclesAvailable,
                   MolStack &molStack, VECT_INT_VECT &atomRingClosures,
                   std::vector<INT_LIST> &atomTraversalBondOrder,
                   const boost::dynamic_bitset<> *bondsInPlay,
                   const std::vector<std::string> *bondSymbols, bool doRandom) {
  Atom *atom = mol.getAtomWithIdx(atomIdx);
  INT_LIST directTravList, cycleEndList;
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
  if (atomRingClosures[atomIdx].size()) {
    std::vector<unsigned int> ringsClosed;
    for (auto bIdx : atomRingClosures[atomIdx]) {
      travList.push_back(bIdx);
      Bond *bond = mol.getBondWithIdx(bIdx);
      seenFromHere.set(bond->getOtherAtomIdx(atomIdx));
      unsigned int ringIdx;
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
        auto cAIt =
            std::find(cyclesAvailable.begin(), cyclesAvailable.end(), 1);
        if (cAIt == cyclesAvailable.end()) {
          throw ValueErrorException(
              "Too many rings open at once. SMILES cannot be generated.");
        }
        unsigned int lowestRingIdx = cAIt - cyclesAvailable.begin();
        cyclesAvailable[lowestRingIdx] = 0;
        ++lowestRingIdx;
        bond->setProp(common_properties::_TraversalRingClosureBond,
                      lowestRingIdx);
        molStack.push_back(MolStackElem(lowestRingIdx));
      }
    }
    for (auto ringIdx : ringsClosed) {
      cyclesAvailable[ringIdx] = 1;
    }
  }

  // ---------------------
  //
  //  Build the list of possible destinations from here
  //
  // ---------------------
  std::vector<PossibleType> possibles;
  possibles.resize(0);
  ROMol::OBOND_ITER_PAIR bondsPair = mol.getAtomBonds(atom);
  possibles.reserve(bondsPair.second - bondsPair.first);

  while (bondsPair.first != bondsPair.second) {
    Bond *theBond = mol[*(bondsPair.first)];
    bondsPair.first++;
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
      unsigned long rank = ranks[otherIdx];
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
       possiblesIt++) {
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
    dfsBuildStack(mol, possibleIdx, bond->getIdx(), colors, cycles, ranks,
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
                           VECT_INT_VECT &cycles, const UINT_VECT &ranks,
                           UINT_VECT &cyclesAvailable, MolStack &molStack,
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

  std::vector<AtomColors> tcolors;
  tcolors.resize(colors.size());
  std::copy(colors.begin(), colors.end(), tcolors.begin());
  dfsFindCycles(mol, atomIdx, inBondIdx, tcolors, ranks, atomRingClosures,
                bondsInPlay, bondSymbols, doRandom);
  dfsBuildStack(mol, atomIdx, inBondIdx, colors, cycles, ranks, cyclesAvailable,
                molStack, atomRingClosures, atomTraversalBondOrder, bondsInPlay,
                bondSymbols, doRandom);
}

void clearBondDirs(ROMol &mol, Bond *refBond, const Atom *fromAtom,
                   UINT_VECT &bondDirCounts, UINT_VECT &atomDirCounts,
                   const UINT_VECT &) {
  PRECONDITION(bondDirCounts.size() >= mol.getNumBonds(), "bad dirCount size");
  PRECONDITION(refBond, "bad bond");
  PRECONDITION(&refBond->getOwningMol() == &mol, "bad bond");
  PRECONDITION(fromAtom, "bad atom");
  PRECONDITION(&fromAtom->getOwningMol() == &mol, "bad bond");

  ROMol::OEDGE_ITER beg, end;
  boost::tie(beg, end) = mol.getAtomBonds(fromAtom);
  bool nbrPossible = false, adjusted = false;
  while (beg != end) {
    Bond *oBond = mol[*beg];
    if (oBond != refBond && canHaveDirection(*oBond)) {
      nbrPossible = true;
      if ((bondDirCounts[oBond->getIdx()] >=
           bondDirCounts[refBond->getIdx()]) &&
          atomDirCounts[oBond->getBeginAtomIdx()] != 1 &&
          atomDirCounts[oBond->getEndAtomIdx()] != 1) {
        adjusted = true;
        bondDirCounts[oBond->getIdx()] -= 1;
        if (!bondDirCounts[oBond->getIdx()]) {
          // no one is setting the direction here:
          oBond->setBondDir(Bond::NONE);
          atomDirCounts[oBond->getBeginAtomIdx()] -= 1;
          atomDirCounts[oBond->getEndAtomIdx()] -= 1;
          // std::cerr<<"ob:"<<oBond->getIdx()<<" ";
        }
      }
    }
    beg++;
  }
  if (nbrPossible && !adjusted &&
      atomDirCounts[refBond->getBeginAtomIdx()] != 1 &&
      atomDirCounts[refBond->getEndAtomIdx()] != 1) {
    // we found a neighbor that could have directionality set,
    // but it had a lower bondDirCount than us, so we must
    // need to be adjusted:
    bondDirCounts[refBond->getIdx()] -= 1;
    if (!bondDirCounts[refBond->getIdx()]) {
      refBond->setBondDir(Bond::NONE);
      atomDirCounts[refBond->getBeginAtomIdx()] -= 1;
      atomDirCounts[refBond->getEndAtomIdx()] -= 1;
    }
  }
}

void removeRedundantBondDirSpecs(ROMol &mol, MolStack &molStack,
                                 UINT_VECT &bondDirCounts,
                                 UINT_VECT &atomDirCounts,
                                 const UINT_VECT &bondVisitOrders) {
  PRECONDITION(bondDirCounts.size() >= mol.getNumBonds(), "bad dirCount size");
  // find bonds that have directions indicated that are redundant:
  for (auto &msI : molStack) {
    if (msI.type == MOL_STACK_BOND) {
      Bond *tBond = msI.obj.bond;
      const Atom *canonBeginAtom = mol.getAtomWithIdx(msI.number);
      const Atom *canonEndAtom =
          mol.getAtomWithIdx(tBond->getOtherAtomIdx(msI.number));
      if (canHaveDirection(*tBond) && bondDirCounts[tBond->getIdx()] >= 1) {
        // start by finding the double bond that sets tBond's direction:
        const Atom *dblBondAtom = nullptr;
        ROMol::OEDGE_ITER beg, end;
        boost::tie(beg, end) = mol.getAtomBonds(canonBeginAtom);
        while (beg != end) {
          if (mol[*beg] != tBond && mol[*beg]->getBondType() == Bond::DOUBLE &&
              mol[*beg]->getStereo() > Bond::STEREOANY) {
            dblBondAtom =
                canonBeginAtom;  // tBond->getOtherAtom(canonBeginAtom);
            break;
          }
          beg++;
        }
        if (dblBondAtom != nullptr) {
          clearBondDirs(mol, tBond, dblBondAtom, bondDirCounts, atomDirCounts,
                        bondVisitOrders);
        }
        dblBondAtom = nullptr;
        boost::tie(beg, end) = mol.getAtomBonds(canonEndAtom);
        while (beg != end) {
          if (mol[*beg] != tBond && mol[*beg]->getBondType() == Bond::DOUBLE &&
              mol[*beg]->getStereo() > Bond::STEREOANY) {
            dblBondAtom = canonEndAtom;  // tBond->getOtherAtom(canonEndAtom);
            break;
          }
          beg++;
        }
        if (dblBondAtom != nullptr) {
          clearBondDirs(mol, tBond, dblBondAtom, bondDirCounts, atomDirCounts,
                        bondVisitOrders);
        }
      } else if (tBond->getBondDir() != Bond::NONE) {
        // we aren't supposed to have a direction set, but we do:
        tBond->setBondDir(Bond::NONE);
      }
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
  PRECONDITION(colors.size() >= mol.getNumAtoms(), "vector too small");
  PRECONDITION(ranks.size() >= mol.getNumAtoms(), "vector too small");
  PRECONDITION(!bondsInPlay || bondsInPlay->size() >= mol.getNumBonds(),
               "bondsInPlay too small");
  PRECONDITION(!bondSymbols || bondSymbols->size() >= mol.getNumBonds(),
               "bondSymbols too small");
  unsigned int nAtoms = mol.getNumAtoms();

  UINT_VECT bondDirCounts(mol.getNumBonds(), 0);
  UINT_VECT atomDirCounts(nAtoms, 0);
  UINT_VECT cyclesAvailable(MAX_CYCLES, 1);
  VECT_INT_VECT cycles(nAtoms);

  boost::dynamic_bitset<> ringStereoChemAdjusted(nAtoms);

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
  Canon::canonicalDFSTraversal(mol, atomIdx, -1, colors, cycles, ranks,
                               cyclesAvailable, molStack, atomRingClosures,
                               atomTraversalBondOrder, bondsInPlay, bondSymbols,
                               doRandom);

  CHECK_INVARIANT(!molStack.empty(), "Empty stack.");
  CHECK_INVARIANT(molStack.begin()->type == MOL_STACK_ATOM,
                  "Corrupted stack. First element should be an atom.");

  // collect some information about traversal order on chiral atoms
  boost::dynamic_bitset<> numSwapsChiralAtoms(nAtoms);
  std::vector<int> atomPermutationIndices(nAtoms, 0);
  if (doIsomericSmiles) {
    for (const auto atom : mol.atoms()) {
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
          int nSwaps = 0;
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
          if (trueOrder.size() < atom->getDegree()) {
            INT_LIST tOrder = trueOrder;
            for (const auto bnd : mol.atomBonds(atom)) {
              int bndIdx = bnd->getIdx();
              if (std::find(trueOrder.begin(), trueOrder.end(), bndIdx) ==
                  trueOrder.end()) {
                tOrder.push_back(bndIdx);
                break;
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

  canonicalizeDoubleBonds(mol, bondVisitOrders, atomVisitOrders, bondDirCounts,
                          atomDirCounts, molStack);

  // traverse the stack and canonicalize atoms with (ring) stereochemistry
  if (doIsomericSmiles) {
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
  Canon::removeRedundantBondDirSpecs(mol, molStack, bondDirCounts,
                                     atomDirCounts, bondVisitOrders);
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

};  // namespace Canon

}  // namespace RDKit
