// $Id$
//
//  Copyright (C) 2001-2010 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Canon.h>
#include <GraphMol/SmilesParse/SmilesParseOps.h>
#include <GraphMol/RDKitQueries.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/hash/hash.hpp>
#include <algorithm>

namespace RDKit {
namespace Canon {
namespace {
bool isUnsaturated(const Atom *atom, const ROMol &mol) {
  ROMol::OEDGE_ITER beg, end;
  boost::tie(beg, end) = mol.getAtomBonds(atom);
  while (beg != end) {
    const Bond *bond = mol[*beg];
    ++beg;
    if (bond->getBondType() != Bond::SINGLE) return true;
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
      std::string descr = (*cIt)->getDescription();
      if (descr == "AtomHCount") {
        if (!(*cIt)->getNegation() &&
            ((ATOM_EQUALS_QUERY *)(*cIt).get())->getVal() == 1) {
          return true;
        }
        return false;
      } else if (descr == "AtomAnd") {
        res = hasSingleHQuery((*cIt).get());
        if (res) return true;
      }
    }
  }
  return res;
}

bool atomHasFourthValence(const Atom *atom) {
  if (atom->getNumExplicitHs() == 1 || atom->getImplicitValence() == 1) {
    return true;
  }
  if (atom->hasQuery()) {
    // the SMARTS [C@@H] produces an atom with a H query, but we also
    // need to treat this like an explicit H for chirality purposes
    // This was Github #1489
    bool res = hasSingleHQuery(atom->getQuery());
    return res;
  }
  return false;
}
}  // end of anonymous namespace

bool chiralAtomNeedsTagInversion(const RDKit::ROMol &mol,
                                 const RDKit::Atom *atom, bool isAtomFirst,
                                 size_t numClosures) {
  PRECONDITION(atom, "bad atom");
  return atom->getDegree() == 3 &&
         ((isAtomFirst && atom->getNumExplicitHs() == 1) ||
          (!atomHasFourthValence(atom) && numClosures == 1 &&
           !isUnsaturated(atom, mol)));
}

struct _possibleCompare
    : public std::binary_function<PossibleType, PossibleType, bool> {
  bool operator()(const PossibleType &arg1, const PossibleType &arg2) const {
    return (arg1.get<0>() < arg2.get<0>());
  }
};

bool checkBondsInSameBranch(MolStack &molStack, Bond *dblBnd, Bond *dirBnd) {
  bool seenDblBond = false;
  int branchCounter = 0;
  for (const auto &item : molStack) {
    switch (item.type) {
      case MOL_STACK_BOND:
        if (item.obj.bond == dirBnd || item.obj.bond == dblBnd) {
          if (seenDblBond) {
            return branchCounter == 0;
          } else {
            seenDblBond = true;
          }
        }
        break;
      case MOL_STACK_BRANCH_OPEN:
        if (seenDblBond) {
          ++branchCounter;
        }
        break;
      case MOL_STACK_BRANCH_CLOSE:
        if (seenDblBond) {
          --branchCounter;
        }
        break;
      default:
        break;
    }
  }
  // We should not ever hit this. But if we do, returning false
  // causes the same behavior as before this patch.
  return false;
}

void switchBondDir(Bond *bond) {
  PRECONDITION(bond, "bad bond");
  PRECONDITION(bond->getBondType() == Bond::SINGLE || bond->getIsAromatic(),
               "bad bond type");
  switch (bond->getBondDir()) {
    case Bond::ENDUPRIGHT:
      bond->setBondDir(Bond::ENDDOWNRIGHT);
      break;
    case Bond::ENDDOWNRIGHT:
      bond->setBondDir(Bond::ENDUPRIGHT);
      break;
    default:
      break;
  }
}

// FIX: this may only be of interest from the SmilesWriter, should we
// move it there?
//
//
void canonicalizeDoubleBond(Bond *dblBond, INT_VECT &bondVisitOrders,
                            INT_VECT &atomVisitOrders, INT_VECT &bondDirCounts,
                            INT_VECT &atomDirCounts, MolStack &molStack) {
  PRECONDITION(dblBond, "bad bond");
  PRECONDITION(dblBond->getBondType() == Bond::DOUBLE, "bad bond order");
  PRECONDITION(dblBond->getStereo() > Bond::STEREOANY, "bad bond stereo");
  PRECONDITION(dblBond->getStereoAtoms().size() >= 2, "bad bond stereo atoms");
  PRECONDITION(atomVisitOrders[dblBond->getBeginAtomIdx()] > 0 ||
                   atomVisitOrders[dblBond->getEndAtomIdx()] > 0,
               "neither end atom traversed");

  // atom1 is the lower numbered atom of the double bond (the one traversed
  // first)
  Atom *atom1, *atom2;
  if (atomVisitOrders[dblBond->getBeginAtomIdx()] <
      atomVisitOrders[dblBond->getEndAtomIdx()]) {
    atom1 = dblBond->getBeginAtom();
    atom2 = dblBond->getEndAtom();
  } else {
    atom1 = dblBond->getEndAtom();
    atom2 = dblBond->getBeginAtom();
  }

  // we only worry about double bonds that begin and end at atoms
  // of degree 2 or 3:
  if ((atom1->getDegree() != 2 && atom1->getDegree() != 3) ||
      (atom2->getDegree() != 2 && atom2->getDegree() != 3)) {
    return;
  }

  Bond *firstFromAtom1 = nullptr, *secondFromAtom1 = nullptr;
  Bond *firstFromAtom2 = nullptr, *secondFromAtom2 = nullptr;

  int firstVisitOrder = 100000;

  ROMol &mol = dblBond->getOwningMol();

  ROMol::OBOND_ITER_PAIR atomBonds;
  // -------------------------------------------------------
  // find the lowest visit order bonds from each end and determine
  // if anything is already constraining our choice of directions:
  bool dir1Set = false, dir2Set = false;
  atomBonds = mol.getAtomBonds(atom1);
  while (atomBonds.first != atomBonds.second) {
    if (mol[*atomBonds.first] != dblBond) {
      int bondIdx = mol[*atomBonds.first]->getIdx();
      if (bondDirCounts[bondIdx] > 0) {
        dir1Set = true;
      }
      if (!firstFromAtom1 || bondVisitOrders[bondIdx] < firstVisitOrder) {
        if (firstFromAtom1) secondFromAtom1 = firstFromAtom1;
        firstFromAtom1 = mol[*atomBonds.first];
        firstVisitOrder = bondVisitOrders[bondIdx];
      } else {
        secondFromAtom1 = mol[*atomBonds.first];
      }
    }
    atomBonds.first++;
  }
  atomBonds = mol.getAtomBonds(atom2);
  firstVisitOrder = 10000;
  while (atomBonds.first != atomBonds.second) {
    if (mol[*atomBonds.first] != dblBond) {
      int bondIdx = mol[*atomBonds.first]->getIdx();
      if (bondDirCounts[bondIdx] > 0) {
        dir2Set = true;
      }
      if (!firstFromAtom2 || bondVisitOrders[bondIdx] < firstVisitOrder) {
        if (firstFromAtom2) secondFromAtom2 = firstFromAtom2;
        firstFromAtom2 = mol[*atomBonds.first];
        firstVisitOrder = bondVisitOrders[bondIdx];
      } else {
        secondFromAtom2 = mol[*atomBonds.first];
      }
    }
    atomBonds.first++;
  }

  // make sure we found everything we need to find:
  CHECK_INVARIANT(firstFromAtom1, "could not find atom1");
  CHECK_INVARIANT(firstFromAtom2, "could not find atom2");
  CHECK_INVARIANT(atom1->getDegree() == 2 || secondFromAtom1,
                  "inconsistency at atom1");
  CHECK_INVARIANT(atom2->getDegree() == 2 || secondFromAtom2,
                  "inconsistency at atom2");

  bool setFromBond1 = true;
  Bond::BondDir atom1Dir = Bond::NONE;
  Bond::BondDir atom2Dir = Bond::NONE;
  Bond *atom1ControllingBond = firstFromAtom1;
  Bond *atom2ControllingBond = firstFromAtom2;
  if (!dir1Set && !dir2Set) {
    // ----------------------------------
    // nothing has touched our bonds so far, so set the
    // directions to "arbitrary" values:

    // the bond we came in on becomes ENDUPRIGHT:
    atom1Dir = Bond::ENDUPRIGHT;
    firstFromAtom1->setBondDir(atom1Dir);
    bondDirCounts[firstFromAtom1->getIdx()] += 1;
    atomDirCounts[atom1->getIdx()] += 1;
  } else if (!dir2Set) {
    // at least one of the bonds on atom1 has its directionality set already:
    if (bondDirCounts[firstFromAtom1->getIdx()] > 0) {
      // The first bond's direction has been set at some earlier point:
      atom1Dir = firstFromAtom1->getBondDir();
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
      // It must be the second bond setting the direction.
      // This happens when the bond dir is set in a branch:
      //        v- this double bond
      //   CC(/C=P/N)=N/O
      //      ^- the second bond sets the direction
      // or when the first bond is a ring closure from an
      // earlier traversed atom:
      //             v- this double bond
      //   NC1=NOC/C1=N\O
      //     ^- this closure ends up being the first bond,
      //        and it does not set the direction.
      //
      // This addresses parts of Issue 185 and sf.net Issue 1842174
      //
      atom1Dir = secondFromAtom1->getBondDir();

      firstFromAtom1->setBondDir(atom1Dir);
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
      atom2Dir = firstFromAtom2->getBondDir();
      bondDirCounts[firstFromAtom2->getIdx()] += 1;
      atomDirCounts[atom2->getIdx()] += 1;
      if (secondFromAtom2) {
        // both bonds have their directionalities set, make sure
        // they are compatible:
        if (firstFromAtom2->getBondDir() == secondFromAtom2->getBondDir() &&
            bondDirCounts[firstFromAtom1->getIdx()]) {
          CHECK_INVARIANT(
              ((firstFromAtom2->getBeginAtomIdx() == atom2->getIdx()) ^
               (secondFromAtom2->getBeginAtomIdx() == atom2->getIdx())),
              "inconsistent state");
        }
      }
    } else {
      // the second bond must be present and setting the direction:
      CHECK_INVARIANT(secondFromAtom2, "inconsistent state");
      CHECK_INVARIANT(bondDirCounts[secondFromAtom2->getIdx()] > 0,
                      "inconsistent state");
      // It must be the second bond setting the direction.
      // This happens when the bond dir is set in a branch:
      //        v- this double bond
      //   CC(/C=P/N)=N/O
      //      ^- the second bond sets the direction
      // or when the first bond is a ring closure from an
      // earlier traversed atom:
      //             v- this double bond
      //   NC1=NOC/C1=N\O
      //     ^- this closure ends up being the first bond,
      //        and it does not set the direction.
      //
      // This addresses parts of Issue 185 and sf.net Issue 1842174
      //
      atom2Dir = secondFromAtom2->getBondDir();

      firstFromAtom2->setBondDir(atom2Dir);
      bondDirCounts[firstFromAtom2->getIdx()] += 1;
      atomDirCounts[atom2->getIdx()] += 2;
      atom2ControllingBond = secondFromAtom2;
    }
    // CHECK_INVARIANT(0,"ring stereochemistry not handled");
  }  // end of the ring stereochemistry if

  // now set the directionality on the other side:
  if (setFromBond1) {
    if (dblBond->getStereo() == Bond::STEREOE ||
        dblBond->getStereo() == Bond::STEREOTRANS) {
      atom2Dir = atom1Dir;
    } else if (dblBond->getStereo() == Bond::STEREOZ ||
               dblBond->getStereo() == Bond::STEREOCIS) {
      atom2Dir = (atom1Dir == Bond::ENDUPRIGHT) ? Bond::ENDDOWNRIGHT
                                                : Bond::ENDUPRIGHT;
    }
    CHECK_INVARIANT(atom2Dir != Bond::NONE, "stereo not set");

    // If we're not looking at the bonds used to determine the
    // stereochemistry, we need to flip the setting on the other bond:
    const INT_VECT &stereoAtoms = dblBond->getStereoAtoms();
    if (atom1->getDegree() == 3 &&
        std::find(stereoAtoms.begin(), stereoAtoms.end(),
                  static_cast<int>(atom1ControllingBond->getOtherAtomIdx(
                      atom1->getIdx()))) == stereoAtoms.end()) {
      atom2Dir = (atom2Dir == Bond::ENDUPRIGHT) ? Bond::ENDDOWNRIGHT
                                                : Bond::ENDUPRIGHT;
    }
    // std::cerr<<" 0 set bond 2: "<<firstFromAtom2->getIdx()<<"
    // "<<atom2Dir<<std::endl;
    if (atom2->getDegree() == 3 &&
        std::find(stereoAtoms.begin(), stereoAtoms.end(),
                  static_cast<int>(firstFromAtom2->getOtherAtomIdx(
                      atom2->getIdx()))) == stereoAtoms.end()) {
      atom2Dir = (atom2Dir == Bond::ENDUPRIGHT) ? Bond::ENDDOWNRIGHT
                                                : Bond::ENDUPRIGHT;
    }
    // std::cerr<<" 1 set bond 2: "<<firstFromAtom2->getIdx()<<"
    // "<<atom2Dir<<std::endl;
    firstFromAtom2->setBondDir(atom2Dir);

    // this block of code is no longer needed
    // if(firstFromAtom2->hasProp(common_properties::_TraversalRingClosureBond)){
    //   // another nice one: we're traversing and come to a ring
    //   // closure bond that has directionality set. This is going to
    //   // have its direction swapped on writing so we need to
    //   // pre-emptively swap it here.
    //   // example situation for this is a non-canonical traversal of
    //   //   C1CCCCN/C=C/1
    //   // starting at atom 0, we hit it on encountering the final bond.
    //   switchBondDir(firstFromAtom2);
    // }

    bondDirCounts[firstFromAtom2->getIdx()] += 1;
    atomDirCounts[atom2->getIdx()] += 1;
  } else {
    // we come before a ring closure:
    if (dblBond->getStereo() == Bond::STEREOZ ||
        dblBond->getStereo() == Bond::STEREOCIS) {
      atom1Dir = atom2Dir;
    } else if (dblBond->getStereo() == Bond::STEREOE ||
               dblBond->getStereo() == Bond::STEREOTRANS) {
      atom1Dir = (atom2Dir == Bond::ENDUPRIGHT) ? Bond::ENDDOWNRIGHT
                                                : Bond::ENDUPRIGHT;
    }
    CHECK_INVARIANT(atom1Dir != Bond::NONE, "stereo not set");
    // If we're not looking at the bonds used to determine the
    // stereochemistry, we need to flip the setting on the other bond:
    const INT_VECT &stereoAtoms = dblBond->getStereoAtoms();
    if (atom2->getDegree() == 3 &&
        std::find(stereoAtoms.begin(), stereoAtoms.end(),
                  static_cast<int>(atom2ControllingBond->getOtherAtomIdx(
                      atom2->getIdx()))) == stereoAtoms.end()) {
      // std::cerr<<"flip 1"<<std::endl;
      atom1Dir = (atom1Dir == Bond::ENDUPRIGHT) ? Bond::ENDDOWNRIGHT
                                                : Bond::ENDUPRIGHT;
    }
    if (atom1->getDegree() == 3 &&
        std::find(stereoAtoms.begin(), stereoAtoms.end(),
                  static_cast<int>(firstFromAtom1->getOtherAtomIdx(
                      atom1->getIdx()))) == stereoAtoms.end()) {
      // std::cerr<<"flip 2"<<std::endl;
      atom1Dir = (atom1Dir == Bond::ENDUPRIGHT) ? Bond::ENDDOWNRIGHT
                                                : Bond::ENDUPRIGHT;
    }

    firstFromAtom1->setBondDir(atom1Dir);
    switchBondDir(firstFromAtom1);
    bondDirCounts[firstFromAtom1->getIdx()] += 1;
    atomDirCounts[atom1->getIdx()] += 1;
  }

  // -----------------------------------
  //
  // Check if there are other bonds from atoms 1 and 2 that need
  // to have their directionalities set:
  ///
  if (atom1->getDegree() == 3) {
    if (!bondDirCounts[secondFromAtom1->getIdx()]) {
      // This bond (the second bond from the starting atom of the double bond)
      // is a special case.  It's going to appear in a branch in the smiles:
      //     X\C(\Y)=C/Z
      //         ^
      //         |- here
      // so it actually needs to go down with the *same* direction as the
      // bond that's already been set (because "pulling the bond out of the
      // branch" reverses its direction).
      // A quick example.  This SMILES:
      //     F/C(\Cl)=C/F
      // is *wrong*. This is the correct form:
      //     F/C(/Cl)=C/F
      // So, since we want this bond to have the opposite direction to the
      // other one, we put it in with the same direction.
      // This was Issue 183

      // UNLESS the bond is not in a branch (in the smiles) (e.g. firstFromAtom1
      // branches off a cycle, and secondFromAtom1 shows up at the end of the
      // cycle). This was Github Issue #2023, see it for an example.
      if (checkBondsInSameBranch(molStack, dblBond, secondFromAtom1)) {
        auto otherDir = (firstFromAtom1->getBondDir() == Bond::ENDUPRIGHT)
                            ? Bond::ENDDOWNRIGHT
                            : Bond::ENDUPRIGHT;
        secondFromAtom1->setBondDir(otherDir);
      } else {
        secondFromAtom1->setBondDir(firstFromAtom1->getBondDir());
      }
    }
    bondDirCounts[secondFromAtom1->getIdx()] += 1;
    atomDirCounts[atom1->getIdx()] += 1;
  }

  if (atom2->getDegree() == 3) {
    if (!bondDirCounts[secondFromAtom2->getIdx()]) {
      // Here we set the bond direction to be opposite the other one (since
      // both come after the atom connected to the double bond).
      Bond::BondDir otherDir;
      if (!secondFromAtom2->hasProp(
              common_properties::_TraversalRingClosureBond)) {
        otherDir = (firstFromAtom2->getBondDir() == Bond::ENDUPRIGHT)
                       ? Bond::ENDDOWNRIGHT
                       : Bond::ENDUPRIGHT;
      } else {
        // another one those irritating little reversal things due to
        // ring closures
        otherDir = firstFromAtom2->getBondDir();
      }
      secondFromAtom2->setBondDir(otherDir);
    }
    bondDirCounts[secondFromAtom2->getIdx()] += 1;
    atomDirCounts[atom2->getIdx()] += 1;
    // std::cerr<<"   other: "<<secondFromAtom2->getIdx()<<"
    // "<<otherDir<<std::endl;
  }

  if (setFromBond1) {
    // This is an odd case... The bonds off the beginning atom are
    // after the start atom in the traversal stack.  These need to
    // have their directions reversed.  An example SMILES (unlikely
    // to actually traverse this way is:
    //   C(=C/O)/F    or C(/F)=C/O
    // That bond is Z, without the reversal, this would come out:
    //   C(=C/O)\F    or C(\F)=C/O
    // which is E.
    //
    // In the case of three-coordinate atoms, we don't need to flip
    // the second bond because the Issue 183 fix (above) already got
    // that one.
    //
    // This was Issue 191 and continued into sf.net issue 1842174
    if (bondVisitOrders[atom1ControllingBond->getIdx()] >
        atomVisitOrders[atom1->getIdx()]) {
      if (bondDirCounts[atom1ControllingBond->getIdx()] == 1) {
        if (!atom1ControllingBond->hasProp(
                common_properties::_TraversalRingClosureBond)) {
          // std::cerr<<"  switcheroo 1"<<std::endl;
          switchBondDir(atom1ControllingBond);
        }
      } else if (bondDirCounts[firstFromAtom2->getIdx()] == 1) {
        // the controlling bond at atom1 is being set by someone else, flip the
        // direction
        // on the atom2 bond instead:
        // std::cerr<<"  switcheroo 2"<<std::endl;
        switchBondDir(firstFromAtom2);
        if (secondFromAtom2 && bondDirCounts[secondFromAtom2->getIdx()] >= 1) {
          switchBondDir(secondFromAtom2);
        }
      }
    }
  }

  // something to watch out for here. For this molecule and traversal order:
  //   0 1 2 3  4 5 6  7 8  <- atom numbers
  //   C/C=C/C(/N=C/C)=C/C
  //        ^  ^
  //        |--|-- these two bonds must match in direction or the SMILES
  //               is inconsistent (according to Daylight, Marvin does ok with
  //               it)
  // That means that the direction of the bond from atom 3->4 needs to be set
  // when the bond from 2->3 is set.
  // Issue2023: But only if 3->4 doesn't have a direction yet?
  //
  // I believe we only need to worry about this for the bonds from atom2.
  const Atom *atom3 = firstFromAtom2->getOtherAtom(atom2);
  if (atom3->getDegree() == 3) {
    Bond *otherAtom3Bond = nullptr;
    bool dblBondPresent = false;
    atomBonds = mol.getAtomBonds(atom3);
    while (atomBonds.first != atomBonds.second) {
      Bond *tbond = mol[*atomBonds.first];
      if (tbond->getBondType() == Bond::DOUBLE &&
          tbond->getStereo() > Bond::STEREOANY) {
        dblBondPresent = true;
      } else if ((tbond->getBondType() == Bond::SINGLE) &&
                 (tbond != firstFromAtom2)) {
        otherAtom3Bond = tbond;
      }
      atomBonds.first++;
    }
    if (dblBondPresent && otherAtom3Bond &&
        otherAtom3Bond->getBondDir() == Bond::NONE) {
      // std::cerr<<"set!"<<std::endl;
      otherAtom3Bond->setBondDir(firstFromAtom2->getBondDir());
      bondDirCounts[otherAtom3Bond->getIdx()] += 1;
      atomDirCounts[atom3->getIdx()] += 1;
    }
  }
}

// finds cycles
void dfsFindCycles(ROMol &mol, int atomIdx, int inBondIdx,
                   std::vector<AtomColors> &colors, const UINT_VECT &ranks,
                   INT_VECT &atomOrders, VECT_INT_VECT &atomRingClosures,
                   const boost::dynamic_bitset<> *bondsInPlay,
                   const std::vector<std::string> *bondSymbols, bool doRandom) {
  Atom *atom = mol.getAtomWithIdx(atomIdx);
  atomOrders.push_back(atomIdx);

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
    if (bondsInPlay && !(*bondsInPlay)[theBond->getIdx()]) continue;
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
        rank = std::rand();
      }
      // std::cerr << "            " << atomIdx << ": " << otherIdx << " " <<
      // rank
      //           << std::endl;
      // std::cerr<<"aIdx: "<< atomIdx <<"   p: "<<otherIdx<<" Rank:
      // "<<ranks[otherIdx] <<" "<<colors[otherIdx]<<"
      // "<<theBond->getBondType()<<" "<<rank<<std::endl;
      possibles.push_back(PossibleType(rank, otherIdx, theBond));
    }
  }

  // ---------------------
  //
  //  Sort on ranks
  //
  // ---------------------
  std::sort(possibles.begin(), possibles.end(), _possibleCompare());
  // if (possibles.size())
  //   std::cerr << " aIdx1: " << atomIdx
  //             << " first: " << possibles.front().get<0>() << " "
  //             << possibles.front().get<1>() << std::endl;
  // // ---------------------
  //
  //  Now work the children
  //
  // ---------------------
  for (auto &possible : possibles) {
    int possibleIdx = possible.get<1>();
    Bond *bond = possible.get<2>();
    switch (colors[possibleIdx]) {
      case WHITE_NODE:
        // -----
        // we haven't seen this node at all before, traverse
        // -----
        dfsFindCycles(mol, possibleIdx, bond->getIdx(), colors, ranks,
                      atomOrders, atomRingClosures, bondsInPlay, bondSymbols,
                      doRandom);
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
                   const UINT_VECT &ranks, INT_VECT &cyclesAvailable,
                   MolStack &molStack, INT_VECT &atomOrders,
                   INT_VECT &bondVisitOrders, VECT_INT_VECT &atomRingClosures,
                   std::vector<INT_LIST> &atomTraversalBondOrder,
                   const boost::dynamic_bitset<> *bondsInPlay,
                   const std::vector<std::string> *bondSymbols, bool doRandom) {
#if 0
    std::cerr<<"traverse from atom: "<<atomIdx<<" via bond "<<inBondIdx<<" num cycles available: "
             <<std::count(cyclesAvailable.begin(),cyclesAvailable.end(),1)<<std::endl;
#endif

  Atom *atom = mol.getAtomWithIdx(atomIdx);
  INT_LIST directTravList, cycleEndList;
  boost::dynamic_bitset<> seenFromHere(mol.getNumAtoms());

  seenFromHere.set(atomIdx);
  molStack.push_back(MolStackElem(atom));
  atomOrders[atom->getIdx()] = rdcast<int>(molStack.size());
  colors[atomIdx] = GREY_NODE;

  INT_LIST travList;
  if (inBondIdx >= 0) travList.push_back(inBondIdx);

  // ---------------------
  //
  //  Add any ring closures
  //
  // ---------------------
  if (atomRingClosures[atomIdx].size()) {
    std::vector<unsigned int> ringsClosed;
    BOOST_FOREACH (int bIdx, atomRingClosures[atomIdx]) {
      travList.push_back(bIdx);
      Bond *bond = mol.getBondWithIdx(bIdx);
      seenFromHere.set(bond->getOtherAtomIdx(atomIdx));
      unsigned int ringIdx;
      if (bond->getPropIfPresent(common_properties::_TraversalRingClosureBond,
                                 ringIdx)) {
        // this is end of the ring closure
        // we can just pull the ring index from the bond itself:
        molStack.push_back(MolStackElem(bond, atomIdx));
        bondVisitOrders[bIdx] = rdcast<int>(molStack.size());
        molStack.push_back(MolStackElem(ringIdx));
        // don't make the ring digit immediately available again: we don't want
        // to have the same
        // ring digit opening and closing rings on an atom.
        ringsClosed.push_back(ringIdx - 1);
      } else {
        // this is the beginning of the ring closure, we need to come up with a
        // ring index:
        INT_VECT::const_iterator cAIt =
            std::find(cyclesAvailable.begin(), cyclesAvailable.end(), 1);
        if (cAIt == cyclesAvailable.end()) {
          throw ValueErrorException(
              "Too many rings open at once. SMILES cannot be generated.");
        }
        unsigned int lowestRingIdx =
            rdcast<unsigned int>(cAIt - cyclesAvailable.begin());
        cyclesAvailable[lowestRingIdx] = 0;
        ++lowestRingIdx;
        bond->setProp(common_properties::_TraversalRingClosureBond,
                      lowestRingIdx);
        molStack.push_back(MolStackElem(lowestRingIdx));
      }
    }
    BOOST_FOREACH (unsigned int ringIdx, ringsClosed) {
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
    if (bondsInPlay && !(*bondsInPlay)[theBond->getIdx()]) continue;
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
        rank = std::rand();
      }

      possibles.push_back(PossibleType(rank, otherIdx, theBond));
    }
  }

  // ---------------------
  //
  //  Sort on ranks
  //
  // ---------------------
  std::sort(possibles.begin(), possibles.end(), _possibleCompare());
  // if (possibles.size())
  //   std::cerr << " aIdx2: " << atomIdx
  //             << " first: " << possibles.front().get<0>() << " "
  //             << possibles.front().get<1>() << std::endl;

  // ---------------------
  //
  //  Now work the children
  //
  // ---------------------
  for (auto possiblesIt = possibles.begin(); possiblesIt != possibles.end();
       possiblesIt++) {
    int possibleIdx = possiblesIt->get<1>();
    if (colors[possibleIdx] != WHITE_NODE) {
      // we're either done or it's a ring-closure, which we already processed...
      // this test isn't strictly required, because we only added WHITE notes to
      // the possibles list, but it seems logical to document it
      continue;
    }
    Bond *bond = possiblesIt->get<2>();
    Atom *otherAtom = mol.getAtomWithIdx(possibleIdx);
    // unsigned int lowestRingIdx;
    INT_VECT::const_iterator cAIt;
    // ww might have some residual data from earlier calls, clean that up:
    if (otherAtom->hasProp(common_properties::_TraversalBondIndexOrder)) {
      otherAtom->clearProp(common_properties::_TraversalBondIndexOrder);
    }
    travList.push_back(bond->getIdx());
    if (possiblesIt + 1 != possibles.end()) {
      // we're branching
      molStack.push_back(
          MolStackElem("(", rdcast<int>(possiblesIt - possibles.begin())));
    }
    molStack.push_back(MolStackElem(bond, atomIdx));
    bondVisitOrders[bond->getIdx()] = rdcast<int>(molStack.size());
    dfsBuildStack(mol, possibleIdx, bond->getIdx(), colors, cycles, ranks,
                  cyclesAvailable, molStack, atomOrders, bondVisitOrders,
                  atomRingClosures, atomTraversalBondOrder, bondsInPlay,
                  bondSymbols, doRandom);
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
                           INT_VECT &cyclesAvailable, MolStack &molStack,
                           INT_VECT &atomOrders, INT_VECT &bondVisitOrders,
                           VECT_INT_VECT &atomRingClosures,
                           std::vector<INT_LIST> &atomTraversalBondOrder,
                           const boost::dynamic_bitset<> *bondsInPlay,
                           const std::vector<std::string> *bondSymbols,
                           bool doRandom) {
  PRECONDITION(colors.size() >= mol.getNumAtoms(), "vector too small");
  PRECONDITION(ranks.size() >= mol.getNumAtoms(), "vector too small");
  PRECONDITION(atomOrders.size() >= mol.getNumAtoms(), "vector too small");
  PRECONDITION(bondVisitOrders.size() >= mol.getNumBonds(), "vector too small");
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
  dfsFindCycles(mol, atomIdx, inBondIdx, tcolors, ranks, atomOrders,
                atomRingClosures, bondsInPlay, bondSymbols, doRandom);
  dfsBuildStack(mol, atomIdx, inBondIdx, colors, cycles, ranks, cyclesAvailable,
                molStack, atomOrders, bondVisitOrders, atomRingClosures,
                atomTraversalBondOrder, bondsInPlay, bondSymbols, doRandom);
}

bool canHaveDirection(const Bond *bond) {
  PRECONDITION(bond, "bad bond");
  Bond::BondType bondType = bond->getBondType();
  return (bondType == Bond::SINGLE || bondType == Bond::AROMATIC);
}

void clearBondDirs(ROMol &mol, Bond *refBond, const Atom *fromAtom,
                   INT_VECT &bondDirCounts, INT_VECT &atomDirCounts,
                   const INT_VECT &bondVisitOrders) {
  RDUNUSED_PARAM(bondVisitOrders);
  PRECONDITION(bondDirCounts.size() >= mol.getNumBonds(), "bad dirCount size");
  PRECONDITION(refBond, "bad bond");
  PRECONDITION(&refBond->getOwningMol() == &mol, "bad bond");
  PRECONDITION(fromAtom, "bad atom");
  PRECONDITION(&fromAtom->getOwningMol() == &mol, "bad bond");

#if 0
    std::copy(bondDirCounts.begin(),bondDirCounts.end(),std::ostream_iterator<int>(std::cerr,", "));
    std::cerr<<"\n";
    std::copy(atomDirCounts.begin(),atomDirCounts.end(),std::ostream_iterator<int>(std::cerr,", "));
    std::cerr<<"\n";
    std::cerr<<"cBD: bond: "<<refBond->getIdx()<<" atom: "<<fromAtom->getIdx()<<": ";
#endif
  ROMol::OEDGE_ITER beg, end;
  boost::tie(beg, end) = mol.getAtomBonds(fromAtom);
  bool nbrPossible = false, adjusted = false;
  while (beg != end) {
    Bond *oBond = mol[*beg];
    // std::cerr<<"  >>"<<oBond->getIdx()<<" "<<canHaveDirection(oBond)<<"
    // "<<bondDirCounts[oBond->getIdx()]<<"-"<<bondDirCounts[refBond->getIdx()]<<"
    // "<<atomDirCounts[oBond->getBeginAtomIdx()]<<"-"<<atomDirCounts[oBond->getEndAtomIdx()]<<std::endl;
    if (oBond != refBond && canHaveDirection(oBond)) {
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
      // std::cerr<<"rb:"<<refBond->getIdx()<<" ";
    }
  }
  // std::cerr<<std::endl;
}

void removeRedundantBondDirSpecs(ROMol &mol, MolStack &molStack,
                                 INT_VECT &bondDirCounts,
                                 INT_VECT &atomDirCounts,
                                 const INT_VECT &bondVisitOrders) {
  PRECONDITION(bondDirCounts.size() >= mol.getNumBonds(), "bad dirCount size");
#if 0
    std::cerr<<"rRBDS: ";
    mol.debugMol(std::cerr);
    std::copy(bondDirCounts.begin(),bondDirCounts.end(),std::ostream_iterator<int>(std::cerr,", "));
    std::cerr<<"\n";
#endif
  // find bonds that have directions indicated that are redundant:
  for (auto &msI : molStack) {
    if (msI.type == MOL_STACK_BOND) {
      Bond *tBond = msI.obj.bond;
      const Atom *canonBeginAtom = mol.getAtomWithIdx(msI.number);
      const Atom *canonEndAtom =
          mol.getAtomWithIdx(tBond->getOtherAtomIdx(msI.number));
      if (canHaveDirection(tBond) && bondDirCounts[tBond->getIdx()] >= 1) {
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

void canonicalizeFragment(ROMol &mol, int atomIdx,
                          std::vector<AtomColors> &colors,
                          const UINT_VECT &ranks, MolStack &molStack,
                          const boost::dynamic_bitset<> *bondsInPlay,
                          const std::vector<std::string> *bondSymbols,
                          bool doIsomericSmiles, bool doRandom) {
  PRECONDITION(colors.size() >= mol.getNumAtoms(), "vector too small");
  PRECONDITION(ranks.size() >= mol.getNumAtoms(), "vector too small");
  PRECONDITION(!bondsInPlay || bondsInPlay->size() >= mol.getNumBonds(),
               "bondsInPlay too small");
  PRECONDITION(!bondSymbols || bondSymbols->size() >= mol.getNumBonds(),
               "bondSymbols too small");
  int nAtoms = mol.getNumAtoms();

  INT_VECT atomVisitOrders(nAtoms, 0);
  INT_VECT bondVisitOrders(mol.getNumBonds(), 0);
  INT_VECT bondDirCounts(mol.getNumBonds(), 0);
  INT_VECT atomDirCounts(nAtoms, 0);
  INT_VECT cyclesAvailable(MAX_CYCLES, 1);

  VECT_INT_VECT cycles(nAtoms);
  for (auto &cycle : cycles) cycle.resize(0);

  boost::dynamic_bitset<> ringStereoChemAdjusted(nAtoms);

  // make sure that we've done the stereo perception:
  if (!mol.hasProp(common_properties::_StereochemDone)) {
    MolOps::assignStereochemistry(mol, false);
  }

  // we need ring information; make sure findSSSR has been called before
  // if not call now
  if (!mol.getRingInfo()->isInitialized()) {
    MolOps::findSSSR(mol);
  }
  mol.getAtomWithIdx(atomIdx)->setProp(common_properties::_TraversalStartPoint,
                                       true);

  VECT_INT_VECT atomRingClosures(nAtoms);
  std::vector<INT_LIST> atomTraversalBondOrder(nAtoms);
  Canon::canonicalDFSTraversal(
      mol, atomIdx, -1, colors, cycles, ranks, cyclesAvailable, molStack,
      atomVisitOrders, bondVisitOrders, atomRingClosures,
      atomTraversalBondOrder, bondsInPlay, bondSymbols, doRandom);

  CHECK_INVARIANT(!molStack.empty(), "Empty stack.");
  CHECK_INVARIANT(molStack.begin()->type == MOL_STACK_ATOM,
                  "Corrupted stack. First element should be an atom.");

  // collect some information about traversal order on chiral atoms
  bool *numSwapsChiralAtoms = (bool *)malloc(nAtoms * sizeof(bool));
  CHECK_INVARIANT(numSwapsChiralAtoms, "failed to allocate memory");
  memset(numSwapsChiralAtoms, 0, nAtoms * sizeof(bool));
  if (doIsomericSmiles) {
    for (ROMol::AtomIterator atomIt = mol.beginAtoms();
         atomIt != mol.endAtoms(); ++atomIt) {
      if ((*atomIt)->getChiralTag() != Atom::CHI_UNSPECIFIED) {
        ROMol::OEDGE_ITER beg, end;
        boost::tie(beg, end) = mol.getAtomBonds(*atomIt);
        while (beg != end) {
          if (bondsInPlay && !(*bondsInPlay)[mol[*beg]->getIdx()]) {
            (*atomIt)->setProp(common_properties::_brokenChirality, true);
            break;
          }
          ++beg;
        }
        if ((*atomIt)->hasProp(common_properties::_brokenChirality)) {
          continue;
        }
        INT_LIST trueOrder = atomTraversalBondOrder[(*atomIt)->getIdx()];
        // Test if the atom is in current fragment
        if (trueOrder.size() > 0) {
          int nSwaps = (*atomIt)->getPerturbationOrder(trueOrder);
          if (chiralAtomNeedsTagInversion(
                  mol, *atomIt,
                  molStack.begin()->obj.atom->getIdx() == (*atomIt)->getIdx(),
                  atomRingClosures[(*atomIt)->getIdx()].size())) {
            // This is a special case. Here's an example:
            //   Our internal representation of a chiral center is equivalent
            //   to:
            //     [C@](F)(O)(C)[H]
            //   we'll be dumping it without the H, which entails a reordering:
            //     [C@@H](F)(O)C
            ++nSwaps;
          }
          if (nSwaps % 2) {
            numSwapsChiralAtoms[(*atomIt)->getIdx()] = 1;
          }
        }
      }
    }
  }

  // remove the current directions on single bonds around double bonds:
  for (ROMol::BondIterator bondIt = mol.beginBonds(); bondIt != mol.endBonds();
       ++bondIt) {
    Bond::BondDir dir = (*bondIt)->getBondDir();
    if (dir == Bond::ENDDOWNRIGHT || dir == Bond::ENDUPRIGHT) {
      (*bondIt)->setBondDir(Bond::NONE);
    }
  }

#if 0
    std::cerr<<"<11111111"<<std::endl;

    std::cerr<<"----------------------------------------->"<<std::endl;
    mol.debugMol(std::cerr);
#endif

  // std::cerr<<"----->\ntraversal stack:"<<std::endl;
  // traverse the stack and canonicalize double bonds and atoms with (ring)
  // stereochemistry
  for (auto &msI : molStack) {
#if 0
      if(msI->type == MOL_STACK_ATOM) std::cerr<<" atom: "<<msI->obj.atom->getIdx()<<std::endl;
      else if(msI->type == MOL_STACK_BOND) std::cerr<<" bond: "<<msI->obj.bond->getIdx()<<" "<<msI->number<<" "<<msI->obj.bond->getBeginAtomIdx()<<"-"<<msI->obj.bond->getEndAtomIdx()<<" order: "<<msI->obj.bond->getBondType()<<std::endl;
      else if(msI->type == MOL_STACK_RING) std::cerr<<" ring: "<<msI->number<<std::endl;
      else if(msI->type == MOL_STACK_BRANCH_OPEN) std::cerr<<" branch open"<<std::endl;
      else if(msI->type == MOL_STACK_BRANCH_CLOSE) std::cerr<<" branch close"<<std::endl;
#endif
    if (msI.type == MOL_STACK_BOND &&
        msI.obj.bond->getBondType() == Bond::DOUBLE &&
        msI.obj.bond->getStereo() > Bond::STEREOANY) {
      if (msI.obj.bond->getStereoAtoms().size() >= 2) {
        Canon::canonicalizeDoubleBond(msI.obj.bond, bondVisitOrders,
                                      atomVisitOrders, bondDirCounts,
                                      atomDirCounts, molStack);
      } else {
        // bad stereo spec:
        msI.obj.bond->setStereo(Bond::STEREONONE);
      }
    }
    if (doIsomericSmiles) {
      if (msI.type == MOL_STACK_ATOM &&
          msI.obj.atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
          !msI.obj.atom->hasProp(common_properties::_brokenChirality)) {
        if (msI.obj.atom->hasProp(common_properties::_ringStereoAtoms)) {
          if (!ringStereoChemAdjusted[msI.obj.atom->getIdx()]) {
            msI.obj.atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
            ringStereoChemAdjusted.set(msI.obj.atom->getIdx());
          }
          const INT_VECT &ringStereoAtoms = msI.obj.atom->getProp<INT_VECT>(
              common_properties::_ringStereoAtoms);
          BOOST_FOREACH (int nbrV, ringStereoAtoms) {
            int nbrIdx = abs(nbrV) - 1;
            // Adjust the chiraliy flag of the ring stereo atoms according to
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
        } else {
          if (msI.obj.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW) {
            if ((numSwapsChiralAtoms[msI.obj.atom->getIdx()])) {
              msI.obj.atom->invertChirality();
            }
          } else if (msI.obj.atom->getChiralTag() ==
                     Atom::CHI_TETRAHEDRAL_CCW) {
            if ((numSwapsChiralAtoms[msI.obj.atom->getIdx()])) {
              msI.obj.atom->invertChirality();
            }
          }
        }
      }
    }
  }
#if 0
    std::cerr<<"<-----"<<std::endl;

    std::cerr<<"----------------------------------------->"<<std::endl;
    mol.debugMol(std::cerr);
#endif
  Canon::removeRedundantBondDirSpecs(mol, molStack, bondDirCounts,
                                     atomDirCounts, bondVisitOrders);
#if 0
    std::cerr<<"----------------------------------------->"<<std::endl;
    mol.debugMol(std::cerr);
    std::cerr<<"----------------------------------------->"<<std::endl;
#endif
  free(numSwapsChiralAtoms);
}
}  // namespace Canon
}  // namespace RDKit
