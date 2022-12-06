//
//  Copyright (C) 2001-2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/Canon.h>
#include <GraphMol/Chirality.h>
#include "SmilesParse.h"
#include "SmilesParseOps.h"
#include <list>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <boost/format.hpp>
#include <RDGeneral/RDLog.h>

namespace SmilesParseOps {
using namespace RDKit;

void ClearAtomChemicalProps(RDKit::Atom *atom) {
  TEST_ASSERT(atom);
  atom->setIsotope(0);
  atom->setFormalCharge(0);
  atom->setNumExplicitHs(0);
}

void CheckRingClosureBranchStatus(RDKit::Atom *atom, RDKit::RWMol *mp) {
  // github #786 and #1652: if the ring closure comes after a branch,
  // the stereochem is wrong.
  // This function is called while closing a branch during construction of
  // the molecule from SMILES and corrects for what happens when parsing odd
  // (and arguably wrong) SMILES constructs like:
  //   1) [C@@](F)1(C)CCO1
  //   2) C1CN[C@](O)(N)1
  //   3) [C@](Cl)(F)1CC[C@H](F)CC1
  // In the first two cases the stereochemistry at the chiral atom
  // needs to be reversed. In the third case the stereochemistry should be
  // reversed when the Cl is added, but left alone when the F is added.
  // We recognize these situations using the index of the chiral atom
  // and the degree of that chiral atom at the time the ring closure
  // digit is encountered during parsing.
  // ----------
  // github #1972 adds these examples:
  //   1) [C@@]1(Cl)(F)I.Br1    (ok)
  //   2) [C@@](Cl)1(F)I.Br1    (reverse)
  //   3) [C@@](Cl)(F)1I.Br1    (ok)
  //   4) [C@@](Cl)(F)(I)1.Br1  (reverse)
  PRECONDITION(atom, "bad atom");
  PRECONDITION(mp, "bad mol");
  if (atom->getIdx() != mp->getNumAtoms(true) - 1 &&
      (atom->getDegree() == 1 ||
       (atom->getDegree() == 2 && atom->getIdx() != 0) ||
       (atom->getDegree() == 3 && atom->getIdx() == 0)) &&
      (atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
       atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW)) {
    // std::cerr << "crcbs: " << atom->getIdx() << std::endl;
    atom->invertChirality();
  }
}

void ReportParseError(const char *message, bool throwIt) {
  PRECONDITION(message, "bad message");
  if (!throwIt) {
    BOOST_LOG(rdErrorLog) << "SMILES Parse Error: " << message << std::endl;
  } else {
    throw SmilesParseException(message);
  }
}

void CleanupAfterParseError(RWMol *mol) {
  PRECONDITION(mol, "no molecule");
  // blow out any partial bonds:
  for (auto markI : *mol->getBondBookmarks()) {
    RWMol::BOND_PTR_LIST &bonds = markI.second;
    for (auto &bond : bonds) {
      delete bond;
    }
  }
}

namespace {
bool couldBeRingClosure(int val) { return val < 100000 && val >= 0; }
}  // namespace
//
// set bondOrder to Bond::IONIC to skip the formation of a bond
//  between the fragment and the molecule
//
void AddFragToMol(RWMol *mol, RWMol *frag, Bond::BondType bondOrder,
                  Bond::BondDir bondDir) {
  PRECONDITION(mol, "no molecule");
  PRECONDITION(frag, "no fragment");
  PRECONDITION(mol->getActiveAtom(), "no active atom");
  Atom *lastAt = mol->getActiveAtom();
  int nOrigAtoms = mol->getNumAtoms();
  int nOrigBonds = mol->getNumBonds();

  //
  // Add the fragment's atoms and bonds to the molecule:
  //
  mol->insertMol(*frag);

  //
  // update ring-closure order information on the added atoms:
  //
  for (const auto atom : frag->atoms()) {
    INT_VECT tmpVect;
    if (atom->getPropIfPresent(common_properties::_RingClosures, tmpVect)) {
      for (auto &v : tmpVect) {
        // if the ring closure is not already a bond, don't touch it:
        if (v >= 0) {
          v += nOrigBonds;
        }
      }
      auto newAtom = mol->getAtomWithIdx(nOrigAtoms + atom->getIdx());
      newAtom->setProp(common_properties::_RingClosures, tmpVect);
    }
  }

  //
  //  ses up the bond between the mol and the branch
  //
  if (bondOrder != Bond::IONIC) {
    // FIX: this is not so much with the elegance...
    auto firstAt = mol->getAtomWithIdx(nOrigAtoms);
    int atomIdx1 = firstAt->getIdx();
    int atomIdx2 = lastAt->getIdx();
    if (frag->hasBondBookmark(ci_LEADING_BOND)) {
      // std::cout << "found it" << std::endl;
      const ROMol::BOND_PTR_LIST &leadingBonds =
          frag->getAllBondsWithBookmark(ci_LEADING_BOND);
      for (auto leadingBond : leadingBonds) {
        // we've already got a bond, so just set its local info
        // and then add it to the molecule intact (no sense doing
        // any extra work).
        leadingBond->setOwningMol(mol);
        leadingBond->setEndAtomIdx(leadingBond->getBeginAtomIdx() + nOrigAtoms);
        leadingBond->setBeginAtomIdx(atomIdx2);
        mol->addBond(leadingBond, true);
      }
      mol->clearBondBookmark(ci_LEADING_BOND);
    } else {
      // SMARTS semantics: unspecified bonds can be single or aromatic
      if (bondOrder == Bond::UNSPECIFIED) {
        auto *newB = new QueryBond(Bond::SINGLE);
        newB->setQuery(makeSingleOrAromaticBondQuery());
        newB->setOwningMol(mol);
        newB->setBeginAtomIdx(atomIdx1);
        newB->setEndAtomIdx(atomIdx2);
        newB->setProp(RDKit::common_properties::_unspecifiedOrder, 1);
        mol->addBond(newB);
        delete newB;
      } else {
        Bond::BondType bo = bondOrder;
        if (bo == Bond::DATIVEL) {
          std::swap(atomIdx1, atomIdx2);
          bo = Bond::DATIVE;
        } else if (bo == Bond::DATIVER) {
          bo = Bond::DATIVE;
        }
        int idx = mol->addBond(atomIdx2, atomIdx1, bo) - 1;
        mol->getBondWithIdx(idx)->setBondDir(bondDir);
      }
    }
  }

  //
  // okay, the next thing we have to worry about is the possibility
  // that there might be ring opening/closing in the fragment we just
  // dealt with e.g. for things like C1C(C1) and C1C.C1
  // We deal with this by copying in the bookmarks and partial bonds
  // that exist in the fragment
  //
  for (auto atIt : *frag->getAtomBookmarks()) {
    // don't bother even considering bookmarks outside
    // the range used for cycles
    if (couldBeRingClosure(atIt.first)) {
      for (auto at2 : atIt.second) {
        int newIdx = at2->getIdx() + nOrigAtoms;
        mol->setAtomBookmark(mol->getAtomWithIdx(newIdx), atIt.first);
        while (frag->hasBondBookmark(atIt.first)) {
          Bond *b = frag->getBondWithBookmark(atIt.first);
          int atomIdx1 = b->getBeginAtomIdx() + nOrigAtoms;
          b->setOwningMol(mol);
          b->setBeginAtomIdx(atomIdx1);
          mol->setBondBookmark(b, atIt.first);
          frag->clearBondBookmark(atIt.first, b);
        }
      }
    }
  }

  frag->clearAllAtomBookmarks();
  frag->clearAllBondBookmarks();
};

typedef std::pair<size_t, int> SIZET_PAIR;
typedef std::pair<int, int> INT_PAIR;
template <typename T>
bool operator<(const std::pair<T, T> &p1, const std::pair<T, T> &p2) {
  return p1.first < p2.first;
}

void AdjustAtomChiralityFlags(RWMol *mol) {
  PRECONDITION(mol, "no molecule");
  for (auto atom : mol->atoms()) {
    Atom::ChiralType chiralType = atom->getChiralTag();
    if (chiralType == Atom::CHI_TETRAHEDRAL_CW ||
        chiralType == Atom::CHI_TETRAHEDRAL_CCW) {
      //
      // The atom is marked as chiral, set the SMILES-order of the
      // atom's bonds.  This is easy for non-ring-closure bonds,
      // because the SMILES order is determined solely by the atom
      // indices.  Things are trickier for ring-closure bonds, which we
      // need to insert into the list in a particular order
      //
      INT_VECT ringClosures;
      atom->getPropIfPresent(common_properties::_RingClosures, ringClosures);

#if 0
      std::cerr << "CLOSURES: ";
      std::copy(ringClosures.begin(), ringClosures.end(),
                std::ostream_iterator<int>(std::cerr, " "));
      std::cerr << std::endl;
#endif
      std::list<SIZET_PAIR> neighbors;
      // push this atom onto the list of neighbors (we'll use this
      // to find our place later):
      neighbors.emplace_back(atom->getIdx(), -1);
      std::list<size_t> bondOrder;
      for (auto nbrIdx :
           boost::make_iterator_range(mol->getAtomNeighbors(atom))) {
        Bond *nbrBond = mol->getBondBetweenAtoms(atom->getIdx(), nbrIdx);
        if (std::find(ringClosures.begin(), ringClosures.end(),
                      static_cast<int>(nbrBond->getIdx())) ==
            ringClosures.end()) {
          neighbors.emplace_back(nbrIdx, nbrBond->getIdx());
        }
      }
      // sort the list of non-ring-closure bonds:
      neighbors.sort();

      // find the location of this atom.  it pretty much has to be
      // first in the list, e.g for smiles like [C@](F)(Cl)(Br)I, or
      // second (everything else).
      auto selfPos = neighbors.begin();
      if (selfPos->first != atom->getIdx()) {
        ++selfPos;
      }
      CHECK_INVARIANT(selfPos->first == atom->getIdx(), "weird atom ordering");

      // copy over the bond ids:
      INT_LIST bondOrdering;
      for (auto neighborIt = neighbors.begin(); neighborIt != neighbors.end();
           ++neighborIt) {
        if (neighborIt != selfPos) {
          bondOrdering.push_back(rdcast<int>(neighborIt->second));
        } else {
          // we are not going to add the atom itself, but we will push on
          // ring closure bonds at this point (if required):
          bondOrdering.insert(bondOrdering.end(), ringClosures.begin(),
                              ringClosures.end());
        }
      }

      // ok, we now have the SMILES ordering of the bonds, figure out the
      // permutation order.
      //
      //  This whole thing is necessary because the ring-closure bonds
      //  in the SMILES come before the bonds to the other neighbors, but
      //  they come after the neighbors in the molecule we build.
      //  A crude example:
      //   in F[C@](Cl)(Br)I the C-Cl bond is index 1 in both SMILES
      //         and as built
      //   in F[C@]1(Br)I.Cl1 the C-Cl bond is index 1 in the SMILES
      //         and index 3 as built.
      //
      int nSwaps = atom->getPerturbationOrder(bondOrdering);
      // FIX: explain this one:
      // At least part of what's going on here for degree 3 atoms:
      //   - The first part: if we're at the beginning of the SMILES and have
      //      an explicit H, we need to add a swap.
      //      This is to reflect that [C@](Cl)(F)C is equivalent to Cl[C@@](F)C
      //      but [C@H](Cl)(F)C is fine as-is (The H-C bond is the one you look
      //      down).
      //   - The second part is more complicated and deals with situations like
      //      F[C@]1CCO1. In this case we otherwise end up looking like we need
      //      to invert the chirality, which is bogus. The chirality here needs
      //      to remain @ just as it does in F[C@](Cl)CCO1
      //   - We have to be careful with the second part to not sweep things like
      //      C[S@]2(=O).Cl2 into the same bin (was github #760). We detect
      //      those cases by looking for unsaturated atoms
      //
      if (Canon::chiralAtomNeedsTagInversion(
              *mol, atom, atom->hasProp(common_properties::_SmilesStart),
              ringClosures.size())) {
        ++nSwaps;
      }
      // std::cerr << "nswaps " << atom->getIdx() << " " << nSwaps
      //           << std::endl;
      // std::copy(bondOrdering.begin(), bondOrdering.end(),
      //           std::ostream_iterator<int>(std::cerr, ", "));
      // std::cerr << std::endl;
      if (nSwaps % 2) {
        atom->invertChirality();
      }
    }
  }
}  // namespace SmilesParseOps

Bond::BondType GetUnspecifiedBondType(const RWMol *mol, const Atom *atom1,
                                      const Atom *atom2) {
  PRECONDITION(mol, "no molecule");
  PRECONDITION(atom1, "no atom1");
  PRECONDITION(atom2, "no atom2");
  Bond::BondType res;
  if (atom1->getIsAromatic() && atom2->getIsAromatic()) {
    res = Bond::AROMATIC;
  } else {
    res = Bond::SINGLE;
  }
  return res;
}
void SetUnspecifiedBondTypes(RWMol *mol) {
  PRECONDITION(mol, "no molecule");
  for (auto bond : mol->bonds()) {
    if (bond->hasProp(RDKit::common_properties::_unspecifiedOrder)) {
      bond->setBondType(GetUnspecifiedBondType(mol, bond->getBeginAtom(),
                                               bond->getEndAtom()));
      if (bond->getBondType() == Bond::AROMATIC) {
        bond->setIsAromatic(true);
      } else {
        bond->setIsAromatic(false);
      }
    }
  }
}

namespace {
void swapBondDirIfNeeded(Bond *bond1, const Bond *bond2) {
  PRECONDITION(bond1, "bad bond1");
  PRECONDITION(bond2, "bad bond2");
  if (bond1->getBondDir() == Bond::NONE && bond2->getBondDir() != Bond::NONE) {
    bond1->setBondDir(bond2->getBondDir());
    if (bond1->getBeginAtom() != bond2->getBeginAtom()) {
      switch (bond1->getBondDir()) {
        case Bond::ENDDOWNRIGHT:
          bond1->setBondDir(Bond::ENDUPRIGHT);
          break;
        case Bond::ENDUPRIGHT:
          bond1->setBondDir(Bond::ENDDOWNRIGHT);
          break;
        default:
          break;
      }
    }
  }
}
}  // namespace

static const std::map<int, int> permutationLimits = {
    {RDKit::Atom::ChiralType::CHI_TETRAHEDRAL, 2},
    {RDKit::Atom::ChiralType::CHI_ALLENE, 2},
    {RDKit::Atom::ChiralType::CHI_SQUAREPLANAR, 3},
    {RDKit::Atom::ChiralType::CHI_OCTAHEDRAL, 30},
    {RDKit::Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL, 20}};

bool checkChiralPermutation(int chiralTag, int permutation) {
  if (chiralTag > RDKit::Atom::ChiralType::CHI_OTHER &&
      permutationLimits.find(chiralTag) != permutationLimits.end() &&
      (permutation < 0 || permutation > permutationLimits.at(chiralTag))) {
    return false;
  }
  return true;
}

void CheckChiralitySpecifications(RDKit::RWMol *mol, bool strict) {
  PRECONDITION(mol, "no molecule");
  for (const auto atom : mol->atoms()) {
    int permutation;
    if (atom->getChiralTag() > RDKit::Atom::ChiralType::CHI_OTHER &&
        permutationLimits.find(atom->getChiralTag()) !=
            permutationLimits.end() &&
        atom->getPropIfPresent(common_properties::_chiralPermutation,
                               permutation)) {
      if (!checkChiralPermutation(atom->getChiralTag(), permutation)) {
        std::string error =
            (boost::format("Invalid chiral specification on atom %d") %
             atom->getIdx())
                .str();
        BOOST_LOG(rdWarningLog) << error << std::endl;
        if (strict) {
          throw SmilesParseException(error);
        }
      }
      // directly convert @TH1 -> @ and @TH2 -> @@
      if (atom->getChiralTag() == RDKit::Atom::ChiralType::CHI_TETRAHEDRAL) {
        if (permutation == 0 || permutation == 1) {
          atom->setChiralTag(RDKit::Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
          atom->clearProp(common_properties::_chiralPermutation);
        } else if (permutation == 2) {
          atom->setChiralTag(RDKit::Atom::ChiralType::CHI_TETRAHEDRAL_CW);
          atom->clearProp(common_properties::_chiralPermutation);
        }
      }
    }
  }
}

void CloseMolRings(RWMol *mol, bool toleratePartials) {
  //  Here's what we want to do here:
  //    loop through the molecule's atom bookmarks
  //    for each bookmark:
  //       connect pairs of atoms sharing that bookmark
  //          left to right (in the order in which they were
  //          inserted into the molecule).
  //       whilst doing this, we have to be cognizant of the fact that
  //          there may well be partial bonds in the molecule which need
  //          to be tied in as well.  WOO HOO! IT'S A BIG MESS!
  PRECONDITION(mol, "no molecule");

  auto bookmarkIt = mol->getAtomBookmarks()->begin();
  while (bookmarkIt != mol->getAtomBookmarks()->end()) {
    auto &bookmark = *bookmarkIt;
    // don't bother even considering bookmarks outside
    // the range used for rings
    if (couldBeRingClosure(bookmark.first)) {
      RWMol::ATOM_PTR_LIST bookmarkedAtomsToRemove;
      auto atomIt = bookmark.second.begin();
      auto atomsEnd = bookmark.second.end();
      while (atomIt != atomsEnd) {
        Atom *atom1 = *atomIt;
        ++atomIt;
        if (!toleratePartials && atomIt == atomsEnd) {
          ReportParseError("unclosed ring");
        } else if (atomIt != atomsEnd && *atomIt == atom1) {
          // make sure we don't try to connect an atom to itself
          // this was github #1925
          auto fmt =
              boost::format{
                  "duplicated ring closure %1% bonds atom %2% to itself"} %
              bookmark.first % atom1->getIdx();
          std::string msg = fmt.str();
          ReportParseError(msg.c_str(), true);
        } else if (mol->getBondBetweenAtoms(atom1->getIdx(),
                                            (*atomIt)->getIdx()) != nullptr) {
          auto fmt =
              boost::format{
                  "ring closure %1% duplicates bond between atom %2% and atom "
                  "%3%"} %
              bookmark.first % atom1->getIdx() % (*atomIt)->getIdx();
          std::string msg = fmt.str();
          ReportParseError(msg.c_str(), true);
        } else if (atomIt != atomsEnd) {
          // we actually found an atom, so connect it to the first
          Atom *atom2 = *atomIt;
          ++atomIt;

          int bondIdx = -1;
          // We're guaranteed two partial bonds, one for each time
          // the ring index was used.  We give the first specification
          // priority.
          CHECK_INVARIANT(mol->hasBondBookmark(bookmark.first),
                          "Missing bond bookmark");

          // now use the info from the partial bond:
          // The partial bond itself will have a proper order and
          // directionality (with a minor caveat documented below) and will
          // have its beginning atom set already:
          RWMol::BOND_PTR_LIST bonds =
              mol->getAllBondsWithBookmark(bookmark.first);
          auto bondIt = bonds.begin();
          CHECK_INVARIANT(bonds.size() >= 2, "Missing bond");

          // get pointers to the two bonds:
          Bond *bond1 = *bondIt;
          ++bondIt;
          Bond *bond2 = *bondIt;

          // remove those bonds from the bookmarks:
          mol->clearBondBookmark(bookmark.first, bond1);
          mol->clearBondBookmark(bookmark.first, bond2);

          // Make sure the bonds have the correct starting atoms:
          CHECK_INVARIANT(bond1->getBeginAtomIdx() == atom1->getIdx(),
                          "bad begin atom");
          CHECK_INVARIANT(bond2->getBeginAtomIdx() == atom2->getIdx(),
                          "bad begin atom");

          // we use the _cxsmilesBondIdx value from the second one, if it's
          // there
          if (bond2->hasProp("_cxsmilesBondIdx")) {
            bond1->setProp("_cxsmilesBondIdx",
                           bond2->getProp<unsigned int>("_cxsmilesBondIdx"));
          }

          Bond *matchedBond;

          // figure out which (if either) bond has a specified type, we'll
          // keep that one.  We also need to update the end atom index to
          // match FIX: daylight barfs when you give it multiple specs for the
          // closure
          //   bond, we'll just take the first one and ignore others
          //   NOTE: we used to do this the other way (take the last
          //   specification),
          //   but that turned out to be troublesome in odd cases like
          //   C1CC11CC1.
          // std::cerr << ">-------------" << std::endl;
          // std::cerr << atom1->getIdx() << "-" << atom2->getIdx() << ": "
          //           << bond1->getBondType() << "("
          //           << bond1->hasProp(common_properties::_unspecifiedOrder)
          //           << "):" << bond1->getBondDir() << " "
          //           << bond2->getBondType() << "("
          //           << bond2->hasProp(common_properties::_unspecifiedOrder)
          //           << "):" << bond2->getBondDir() << std::endl;
          if (!bond1->hasProp(common_properties::_unspecifiedOrder)) {
            matchedBond = bond1;
            if (matchedBond->getBondType() == Bond::DATIVEL) {
              matchedBond->setBeginAtomIdx(atom2->getIdx());
              matchedBond->setEndAtomIdx(atom1->getIdx());
              matchedBond->setBondType(Bond::DATIVE);
            } else if (matchedBond->getBondType() == Bond::DATIVER) {
              matchedBond->setEndAtomIdx(atom2->getIdx());
              matchedBond->setBondType(Bond::DATIVE);
            } else {
              matchedBond->setEndAtomIdx(atom2->getIdx());
            }
            swapBondDirIfNeeded(bond1, bond2);
            delete bond2;
          } else {
            matchedBond = bond2;
            if (matchedBond->getBondType() == Bond::DATIVEL) {
              matchedBond->setBeginAtomIdx(atom1->getIdx());
              matchedBond->setEndAtomIdx(atom2->getIdx());
              matchedBond->setBondType(Bond::DATIVE);
            } else if (matchedBond->getBondType() == Bond::DATIVER) {
              matchedBond->setEndAtomIdx(atom1->getIdx());
              matchedBond->setBondType(Bond::DATIVE);
            } else {
              matchedBond->setEndAtomIdx(atom1->getIdx());
            }
            swapBondDirIfNeeded(bond2, bond1);
            delete bond1;
          }
          if (matchedBond->getBondType() == Bond::UNSPECIFIED &&
              !matchedBond->hasQuery()) {
            Bond::BondType bondT = GetUnspecifiedBondType(mol, atom1, atom2);
            matchedBond->setBondType(bondT);
          }
          matchedBond->setOwningMol(mol);
          if (matchedBond->getBondType() == Bond::AROMATIC) {
            matchedBond->setIsAromatic(true);
          }

          // add the bond:
          bondIdx = mol->addBond(matchedBond, true);

          // we found a bond, so update the atom's _RingClosures
          // property:
          if (bondIdx > -1) {
            CHECK_INVARIANT(
                atom1->hasProp(common_properties::_RingClosures) &&
                    atom2->hasProp(common_properties::_RingClosures),
                "somehow atom doesn't have _RingClosures property.");
            INT_VECT closures;
            atom1->getProp(common_properties::_RingClosures, closures);
            auto closurePos = std::find(closures.begin(), closures.end(),
                                        -(bookmark.first + 1));
            CHECK_INVARIANT(closurePos != closures.end(),
                            "could not find bookmark in atom _RingClosures");
            *closurePos = bondIdx - 1;
            atom1->setProp(common_properties::_RingClosures, closures);

            atom2->getProp(common_properties::_RingClosures, closures);
            closurePos = std::find(closures.begin(), closures.end(),
                                   -(bookmark.first + 1));
            CHECK_INVARIANT(closurePos != closures.end(),
                            "could not find bookmark in atom _RingClosures");
            *closurePos = bondIdx - 1;
            atom2->setProp(common_properties::_RingClosures, closures);
          }
          bookmarkedAtomsToRemove.push_back(atom1);
          bookmarkedAtomsToRemove.push_back(atom2);
        }
      }
      int mark = bookmark.first;
      ++bookmarkIt;
      for (const auto atom : bookmarkedAtomsToRemove) {
        mol->clearAtomBookmark(mark, atom);
      }
    } else {
      ++bookmarkIt;
    }
  }
};

void CleanupAfterParsing(RWMol *mol) {
  PRECONDITION(mol, "no molecule");
  for (auto atom : mol->atoms()) {
    atom->clearProp(common_properties::_RingClosures);
    atom->clearProp(common_properties::_SmilesStart);
  }
  for (auto bond : mol->bonds()) {
    bond->clearProp(common_properties::_unspecifiedOrder);
    bond->clearProp("_cxsmilesBondIdx");
  }
  for (auto sg : RDKit::getSubstanceGroups(*mol)) {
    sg.clearProp("_cxsmilesindex");
  }
  if (!Chirality::getAllowNontetrahedralChirality()) {
    bool needWarn = false;
    for (auto atom : mol->atoms()) {
      if (atom->hasProp(common_properties::_chiralPermutation)) {
        needWarn = true;
        atom->clearProp(common_properties::_chiralPermutation);
      }
      if (atom->getChiralTag() > Atom::ChiralType::CHI_OTHER) {
        needWarn = true;
        atom->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
      }
    }
    if (needWarn) {
      BOOST_LOG(rdWarningLog)
          << "ignoring non-tetrahedral stereo specification since setAllowNontetrahedralChirality() is false."
          << std::endl;
    }
  }
}

RDKit::QueryBond *getUnspecifiedQueryBond(const RDKit::Atom *a1,
                                          const RDKit::Atom *a2) {
  PRECONDITION(a1, "bad atom pointer");
  QueryBond *newB;
  if (!a1->getIsAromatic() || (a2 && !a2->getIsAromatic())) {
    newB = new QueryBond(Bond::SINGLE);
    newB->setQuery(makeSingleOrAromaticBondQuery());
  } else {
    newB = new QueryBond(Bond::AROMATIC);
    newB->setQuery(makeSingleOrAromaticBondQuery());
  }
  newB->setProp(RDKit::common_properties::_unspecifiedOrder, 1);
  return newB;
}

}  // end of namespace SmilesParseOps
