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
#include <GraphMol/QueryOps.h>
#include <GraphMol/Canon.h>
#include <GraphMol/Rings.h>
#include <GraphMol/SanitException.h>
#include <RDGeneral/RDLog.h>
#include <boost/dynamic_bitset.hpp>
#include <utility>

namespace RDKit {
// Local utility namespace
namespace {

void backTrack(RWMol &mol, INT_INT_DEQ_MAP &, int lastOpt, INT_VECT &done,
               INT_DEQUE &aqueue, boost::dynamic_bitset<> &dBndCands,
               boost::dynamic_bitset<> &dBndAdds) {
  // so we made a wrong turn at the lastOpt
  // remove on done list that comes after the lastOpt including itself

  auto ei = std::find(done.begin(), done.end(), lastOpt);
  INT_VECT tdone;
  tdone.insert(tdone.end(), done.begin(), ei);

  INT_VECT_CRI eri = std::find(done.rbegin(), done.rend(), lastOpt);
  ++eri;
  // and push them back onto the stack
  for (INT_VECT_CRI ri = done.rbegin(); ri != eri; ++ri) {
    aqueue.push_front(*ri);
  }

  // remove any double bonds that were add since we passed through lastOpt
  Bond *bnd;
  unsigned int nbnds = mol.getNumBonds();
  for (unsigned int bi = 0; bi < nbnds; ++bi) {
    if (dBndAdds[bi]) {
      bnd = mol.getBondWithIdx(bi);
      int aid1 = bnd->getBeginAtomIdx();
      int aid2 = bnd->getEndAtomIdx();
      // if one of these atoms has been dealt with before lastOpt
      // we don't have to change the double bond addition
      if ((std::find(tdone.begin(), tdone.end(), aid1) == tdone.end()) &&
          (std::find(tdone.begin(), tdone.end(), aid2) == tdone.end())) {
        // otherwise strip the double bond and set it back to single
        // and add the atoms to candidate for double bonds
        dBndAdds[bi] = 0;
        bnd->setBondType(Bond::SINGLE);
        dBndCands[aid1] = 1;
        dBndCands[aid2] = 1;
      }
    }
  }
  done = tdone;
}

void markDbondCands(RWMol &mol, const INT_VECT &allAtms,
                    boost::dynamic_bitset<> &dBndCands, INT_VECT &questions,
                    INT_VECT &done) {
  // ok this function does more than mark atoms that are candidates for
  // double bonds during kekulization
  // - check that a non-aromatic atom does not have any aromatic bonds
  // - marks all aromatic bonds to single bonds
  // - marks atoms that can take a double bond

  bool hasAromaticOrDummyAtom =
      std::any_of(allAtms.begin(), allAtms.end(), [&mol](int allAtm) {
        return (!mol.getAtomWithIdx(allAtm)->getAtomicNum() ||
                isAromaticAtom(*mol.getAtomWithIdx(allAtm)));
      });
  // if there's not at least one atom in the ring that's
  // marked as being aromatic or a dummy,
  // there's no point in continuing:
  if (!hasAromaticOrDummyAtom) {
    return;
  }
  // mark rings which are not candidates for double bonds
  // i.e. that have at least one atom which is in a single ring
  // and is not aromatic
  boost::dynamic_bitset<> isRingNotCand(mol.getRingInfo()->numRings());
  unsigned int ri = 0;
  for (const auto &aring : mol.getRingInfo()->atomRings()) {
    isRingNotCand.set(ri);
    for (auto ai : aring) {
      const auto at = mol.getAtomWithIdx(ai);
      if (isAromaticAtom(*at) && mol.getRingInfo()->numAtomRings(ai) == 1) {
        isRingNotCand.reset(ri);
        break;
      }
    }
    ++ri;
  }
  std::vector<Bond *> makeSingle;

  boost::dynamic_bitset<> inAllAtms(mol.getNumAtoms());
  for (int allAtm : allAtms) {
    inAllAtms.set(allAtm);
    Atom *at = mol.getAtomWithIdx(allAtm);

    if (at->getAtomicNum() && !isAromaticAtom(*at)) {
      done.push_back(allAtm);
      continue;
    }

    // count the number of neighbors connected with single,
    // double, or aromatic bonds. Along the way, mark
    // bonds that we will later mark as being single:
    int sbo = 0;
    unsigned nToIgnore = 0;
    unsigned int nonArNonDummyNbr = 0;
    for (const auto bond : mol.atomBonds(at)) {
      auto otherAt = bond->getOtherAtom(at);
      if (otherAt->getAtomicNum() && !otherAt->getIsAromatic() &&
          inAllAtms.test(otherAt->getIdx())) {
        ++nonArNonDummyNbr;
      }
      if (bond->getIsAromatic() && (bond->getBondType() == Bond::SINGLE ||
                                    bond->getBondType() == Bond::DOUBLE ||
                                    bond->getBondType() == Bond::AROMATIC)) {
        ++sbo;
        // mark this bond to be marked single later
        // we don't want to do right now because it can screw-up the
        // valence calculation to determine the number of hydrogens below
        makeSingle.push_back(bond);
      } else {
        int bondContrib = std::lround(bond->getValenceContrib(at));
        sbo += bondContrib;
        if (!bondContrib) {
          ++nToIgnore;
        }
      }
    }

    auto numAtomRings = mol.getRingInfo()->numAtomRings(at->getIdx());
    const auto &riVect = mol.getRingInfo()->atomMembers(at->getIdx());
    auto numNonCandRings = std::count_if(
        riVect.begin(), riVect.end(),
        [&isRingNotCand](int ri) { return isRingNotCand.test(ri); });
    if (!at->getAtomicNum() && nonArNonDummyNbr < numAtomRings &&
        numNonCandRings < numAtomRings) {
      // dummies always start as candidates to have a double bond:
      dBndCands[allAtm] = 1;
      // but they don't have to have one, so mark them as questionable:
      questions.push_back(allAtm);
    } else {
      // for non dummies, it's a bit more work to figure out if they
      // can take a double bond:

      sbo += at->getTotalNumHs();
      auto dv =
          PeriodicTable::getTable()->getDefaultValence(at->getAtomicNum());
      auto chrg = at->getFormalCharge();
      if (isEarlyAtom(at->getAtomicNum())) {
        chrg = -chrg;  // fix for GitHub #65
      }
      // special case for carbon - see GitHub #539
      if (at->getAtomicNum() == 6 && chrg > 0) {
        chrg = -chrg;
      }
      dv += chrg;
      int tbo = at->getTotalValence();
      int nRadicals = at->getNumRadicalElectrons();
      int totalDegree = at->getDegree() + at->getImplicitValence() - nToIgnore;

      const auto &valList =
          PeriodicTable::getTable()->getValenceList(at->getAtomicNum());
      unsigned int vi = 1;

      while (tbo > dv && vi < valList.size() && valList[vi] > 0) {
        dv = valList[vi] + chrg;
        ++vi;
      }

      // Kekulize aromatic N-oxides, such as O=n1ccccc1
      // These only reach here if SANITIZE_CLEANUP is disabled.
      if (tbo == 5 && sbo == 4 && dv == 3 && totalDegree == 3 &&
          nRadicals == 0 && chrg == 0 && at->getTotalNumHs() == 0) {
        switch (at->getAtomicNum()) {
          case 7:   // N
          case 15:  // P
          case 33:  // As
            dv = 5;
            break;
        }
      }
      // std::cerr << "  kek: " << at->getIdx() << " tbo:" << tbo << " sbo:" <<
      // sbo
      //           << "  dv : " << dv << " totalDegree : " << totalDegree
      //           << " nRadicals: " << nRadicals << std::endl;
      if (totalDegree + nRadicals >= dv) {
        // if our degree + nRadicals exceeds the default valence,
        // there's no way we can take a double bond, just continue.
        continue;
      }

      // we're a candidate if our total current bond order + nRadicals + 1
      // matches the valence state
      // (including nRadicals here was SF.net issue 3349243)
      if (dv == (sbo + 1 + nRadicals)) {
        dBndCands[allAtm] = 1;
      } else if (!nRadicals && at->getNoImplicit() && dv == (sbo + 2)) {
        // special case: there is currently no radical on the atom, but if
        // if we allow one then this is a candidate:
        dBndCands[allAtm] = 1;
      }
    }
  }  // loop over all atoms in the fused system

  // now turn all the aromatic bond in this fused system to single
  for (auto &bi : makeSingle) {
    bi->setBondType(Bond::SINGLE);
  }
}

bool kekulizeWorker(RWMol &mol, const INT_VECT &allAtms,
                    boost::dynamic_bitset<> dBndCands,
                    boost::dynamic_bitset<> dBndAdds, INT_VECT done,
                    unsigned int maxBackTracks) {
  INT_DEQUE astack;
  INT_INT_DEQ_MAP options;
  int lastOpt = -1;
  boost::dynamic_bitset<> localBondsAdded(mol.getNumBonds());

  // ok the algorithm goes something like this
  // - start with an atom that has been marked aromatic before
  // - check if it can have a double bond
  // - add its neighbors to the stack
  // - check if one of its neighbors can also have a double bond
  // - if yes add a double bond.
  // - if multiple neighbors can have double bonds - add them to a
  //   options stack we may have to retrace out path if we chose the
  //   wrong neighbor to add the double bond
  // - if double bond added update the candidates for double bond
  // - move to the next atom on the stack and repeat the process
  // - if an atom that can have multiple a double bond has no
  //   neighbors that can take double bond - we made a mistake
  //   earlier by picking a wrong candidate for double bond
  // - in this case back track to where we made the mistake

  int curr = -1;
  INT_DEQUE btmoves;
  unsigned int numBT = 0;  // number of back tracks so far
  while ((done.size() < allAtms.size()) || !astack.empty()) {
    // pick a curr atom to work with
    if (astack.size() > 0) {
      curr = astack.front();
      astack.pop_front();
    } else {
      for (int allAtm : allAtms) {
        if (std::find(done.begin(), done.end(), allAtm) == done.end()) {
          curr = allAtm;
          break;
        }
      }
    }
    CHECK_INVARIANT(curr >= 0, "starting point not found");
    done.push_back(curr);

    // loop over the neighbors if we can add double bonds or
    // simply push them onto the stack
    INT_DEQUE opts;
    bool cCand = false;
    if (dBndCands[curr]) {
      cCand = true;
    }
    int ncnd;
    // if we are here because of backtracking
    if (options.find(curr) != options.end()) {
      opts = options[curr];
      CHECK_INVARIANT(opts.size() > 0, "");
    } else {
      INT_DEQUE lstack;
      INT_DEQUE wedgedOpts;
      for (const auto &nbrIdx : boost::make_iterator_range(
               mol.getAtomNeighbors(mol.getAtomWithIdx(curr)))) {
        // ignore if the neighbor has already been dealt with before
        if (std::find(done.begin(), done.end(), nbrIdx) != done.end()) {
          continue;
        }
        // ignore if the neighbor is not part of the fused system
        if (std::find(allAtms.begin(), allAtms.end(), nbrIdx) ==
            allAtms.end()) {
          continue;
        }
        auto nbrBond = mol.getBondBetweenAtoms(curr, nbrIdx);

        // if the neighbor is not on the stack add it
        if (std::find(astack.begin(), astack.end(), nbrIdx) == astack.end()) {
          lstack.push_front(nbrIdx);
        }

        // check if the neighbor is also a candidate for a double bond
        // the refinement that we'll make to the candidate check we've already
        // done is to make sure that the bond is either flagged as aromatic
        // or involves a dummy atom. This was Issue 3525076.
        // This fix is not really 100% of the way there: a situation like
        // that for Issue 3525076 but involving a dummy atom in the cage
        // could lead to the same failure. The full fix would require
        // a fairly detailed analysis of all bonds in the molecule to determine
        // which of them is eligible to be converted.
        if (cCand && dBndCands[nbrIdx] &&
            (nbrBond->getIsAromatic() ||
             mol.getAtomWithIdx(curr)->getAtomicNum() == 0 ||
             mol.getAtomWithIdx(nbrIdx)->getAtomicNum() == 0)) {
          // in order to try and avoid making wedged bonds double, we will add
          // this neighbor at the back of the options after this loop if the
          // bond is wedged. otherwise we append it to the options directly
          if (nbrBond->getBondDir() == Bond::BondDir::BEGINWEDGE ||
              nbrBond->getBondDir() == Bond::BondDir::BEGINDASH) {
            wedgedOpts.push_back(nbrIdx);
          } else {
            opts.push_back(nbrIdx);
          }
        }  // end of curr atoms can have a double bond
      }    // end of looping over neighbors
      // now append to opts the neighbors connected via wedged bonds,
      // if any were found
      opts.insert(opts.end(), wedgedOpts.begin(), wedgedOpts.end());
      astack.insert(astack.end(), lstack.begin(), lstack.end());
    }
    // now add a double bond from current to one of the neighbors if we can
    if (cCand) {
      if (!opts.empty()) {
        ncnd = opts.front();
        opts.pop_front();
        auto bnd = mol.getBondBetweenAtoms(curr, ncnd);
        bnd->setBondType(Bond::DOUBLE);

        // remove current and the neighbor from the dBndCands list
        dBndCands[curr] = 0;
        dBndCands[ncnd] = 0;

        // add them to the list of bonds to which have been made double
        dBndAdds[bnd->getIdx()] = 1;
        localBondsAdded[bnd->getIdx()] = 1;

        // if this is an atom we previously visted and picked we
        // simply tried a different option now, overwrite the options
        // stored for this atoms
        if (options.find(curr) != options.end()) {
          if (opts.size() == 0) {
            options.erase(curr);
            btmoves.pop_back();
            if (btmoves.size() > 0) {
              lastOpt = btmoves.back();
            } else {
              lastOpt = -1;
            }
          } else {
            options[curr] = opts;
          }
        } else {
          // this is new atoms we are trying and have other
          // neighbors as options to add double bond store this to
          // the options stack, we may have made a mistake in
          // which one we chose and have to return here
          if (opts.size() > 0) {
            lastOpt = curr;
            btmoves.push_back(lastOpt);
            options[curr] = opts;
          }
        }

      }  // end of adding a double bond
      else if (mol.getAtomWithIdx(curr)->getAtomicNum()) {
        // we have a non-dummy atom that should be getting a double
        // bond but none of the neighbors can take one. Most likely
        // because of a wrong choice earlier so back track
        if ((lastOpt >= 0) && (numBT < maxBackTracks)) {
          // std::cerr << "PRE BACKTRACK" << std::endl;
          // mol.debugMol(std::cerr);
          backTrack(mol, options, lastOpt, done, astack, dBndCands, dBndAdds);
          // std::cerr << "POST BACKTRACK" << std::endl;
          // mol.debugMol(std::cerr);
          ++numBT;
        } else {
          // undo any remaining changes we made while here
          // this was github #962
          for (unsigned int bidx = 0; bidx < mol.getNumBonds(); ++bidx) {
            if (localBondsAdded[bidx]) {
              mol.getBondWithIdx(bidx)->setBondType(Bond::SINGLE);
            }
          }
          return false;
        }
      }  // end of else try to backtrack
    }    // end of curr atom atom being a cand for double bond
  }      // end of while we are not done with all atoms
  return true;
}

class QuestionEnumerator {
 public:
  QuestionEnumerator(INT_VECT questions)
      : d_questions(std::move(questions)), d_pos(1){};
  INT_VECT next() {
    INT_VECT res;
    if (d_pos >= (0x1u << d_questions.size())) {
      return res;
    }
    for (unsigned int i = 0; i < d_questions.size(); ++i) {
      if (d_pos & (0x1u << i)) {
        res.push_back(d_questions[i]);
      }
    }
    ++d_pos;
    return res;
  };

 private:
  INT_VECT d_questions;
  unsigned int d_pos;
};

bool permuteDummiesAndKekulize(RWMol &mol, const INT_VECT &allAtms,
                               boost::dynamic_bitset<> dBndCands,
                               INT_VECT &questions,
                               unsigned int maxBackTracks) {
  boost::dynamic_bitset<> atomsInPlay(mol.getNumAtoms());
  for (int allAtm : allAtms) {
    atomsInPlay[allAtm] = 1;
  }
  bool kekulized = false;
  QuestionEnumerator qEnum(questions);
  while (!kekulized && questions.size()) {
    boost::dynamic_bitset<> dBndAdds(mol.getNumBonds());
    INT_VECT done;
#if 1
    // reset the state: all aromatic bonds are remarked to single:
    for (const auto bond : mol.bonds()) {
      if (bond->getIsAromatic() && bond->getBondType() != Bond::SINGLE &&
          atomsInPlay[bond->getBeginAtomIdx()] &&
          atomsInPlay[bond->getEndAtomIdx()]) {
        bond->setBondType(Bond::SINGLE);
      }
    }
#endif
    // pick a new permutation of the questionable atoms:
    const auto &switchOff = qEnum.next();
    if (!switchOff.size()) {
      break;
    }
    auto tCands = dBndCands;
    for (int it : switchOff) {
      tCands[it] = 0;
    }
#if 0
        std::cerr<<"permute: ";
        for (boost::dynamic_bitset<>::size_type i = 0; i < tCands.size(); ++i){
          std::cerr << tCands[i];
        }
        std::cerr<<std::endl;
#endif
    // try kekulizing again:
    kekulized =
        kekulizeWorker(mol, allAtms, tCands, dBndAdds, done, maxBackTracks);
  }
  return kekulized;
}

void kekulizeFused(RWMol &mol, const VECT_INT_VECT &arings,
                   unsigned int maxBackTracks) {
  // get all the atoms in the ring system
  INT_VECT allAtms;
  Union(arings, allAtms);
  // get all the atoms that are candidates to receive a double bond
  // also mark atoms in the fused system that are not aromatic to begin with
  // as done. Mark all the bonds that are part of the aromatic system
  // to be single bonds
  INT_VECT done;
  INT_VECT questions;
  auto nats = mol.getNumAtoms();
  auto nbnds = mol.getNumBonds();
  boost::dynamic_bitset<> dBndCands(nats);
  boost::dynamic_bitset<> dBndAdds(nbnds);
  markDbondCands(mol, allAtms, dBndCands, questions, done);
#if 0
      std::cerr << "candidates: ";
      for(int i=0;i<nats;++i) std::cerr << dBndCands[i];
      std::cerr << std::endl;
#endif

  auto kekulized =
      kekulizeWorker(mol, allAtms, dBndCands, dBndAdds, done, maxBackTracks);
  if (!kekulized && questions.size()) {
    // we failed, but there are some dummy atoms we can try permuting.
    kekulized = permuteDummiesAndKekulize(mol, allAtms, dBndCands, questions,
                                          maxBackTracks);
  }
  if (!kekulized) {
    // we exhausted all option (or crossed the allowed
    // number of backTracks) and we still need to backtrack
    // can't kekulize this thing
    std::vector<unsigned int> problemAtoms;
    std::ostringstream errout;
    errout << "Can't kekulize mol.";
    errout << "  Unkekulized atoms:";
    for (unsigned int i = 0; i < nats; ++i) {
      if (dBndCands[i]) {
        errout << " " << i;
        problemAtoms.push_back(i);
      }
    }
    std::string msg = errout.str();
    BOOST_LOG(rdErrorLog) << msg << std::endl;
    throw KekulizeException(msg, problemAtoms);
  }
}
}  // namespace

namespace MolOps {
namespace details {
void KekulizeFragment(RWMol &mol, const boost::dynamic_bitset<> &atomsToUse,
                      boost::dynamic_bitset<> bondsToUse, bool markAtomsBonds,
                      unsigned int maxBackTracks) {
  PRECONDITION(atomsToUse.size() == mol.getNumAtoms(),
               "atomsToUse is wrong size");
  PRECONDITION(bondsToUse.size() == mol.getNumBonds(),
               "bondsToUse is wrong size");
  // if there are no atoms to use we can directly return
  if (atomsToUse.none()) {
    return;
  }

  // there's no point doing kekulization if there are no aromatic bonds
  // without queries:
  bool foundAromatic = false;
  for (const auto bond : mol.bonds()) {
    if (bondsToUse[bond->getIdx()]) {
      if (QueryOps::hasBondTypeQuery(*bond)) {
        // we don't kekulize bonds with bond type queries
        bondsToUse[bond->getIdx()] = 0;
      } else if (bond->getIsAromatic()) {
        foundAromatic = true;
      }
    }
  }

  // before everything do implicit valence calculation and store them
  // we will repeat after kekulization and compare for the sake of error
  // checking
  auto numAtoms = mol.getNumAtoms();
  INT_VECT valences(numAtoms);
  boost::dynamic_bitset<> dummyAts(numAtoms);

  for (auto atom : mol.atoms()) {
    if (!atomsToUse[atom->getIdx()]) {
      continue;
    }
    atom->calcImplicitValence(false);
    valences[atom->getIdx()] = atom->getTotalValence();
    if (isAromaticAtom(*atom)) {
      foundAromatic = true;
    }
    if (!atom->getAtomicNum()) {
      dummyAts[atom->getIdx()] = 1;
    }
  }
  if (!foundAromatic) {
    return;
  }
  // if any bonds to kekulize then give it a try:
  if (bondsToUse.any()) {
    // mark atoms at the beginning of wedged bonds
    boost::dynamic_bitset<> wedgedAtoms(numAtoms);
    for (const auto bond : mol.bonds()) {
      if (bondsToUse[bond->getIdx()] &&
          (bond->getBondDir() == Bond::BEGINWEDGE ||
           bond->getBondDir() == Bond::BEGINDASH)) {
        wedgedAtoms.set(bond->getBeginAtomIdx());
      }
    }

    // A bit on the state of the molecule at this point
    // - aromatic and non aromatic atoms and bonds may be mixed up

    // - for all aromatic bonds it is assumed that that both the following
    //   are true:
    //       - getIsAromatic returns true
    //       - getBondType return aromatic
    // - all aromatic atoms return true for "getIsAromatic"

    // first find all the simple rings in the molecule that are not
    // completely composed of dummy atoms
    VECT_INT_VECT allringsSSSR;
    if (!mol.getRingInfo()->isInitialized()) {
      MolOps::findSSSR(mol, allringsSSSR);
    }
    const VECT_INT_VECT &allrings =
        allringsSSSR.empty() ? mol.getRingInfo()->atomRings() : allringsSSSR;
    std::deque<INT_VECT> tmpRings;
    auto containsNonDummy = [&atomsToUse, &dummyAts](const INT_VECT &ring) {
      bool ringOk = false;
      for (auto ai : ring) {
        if (!atomsToUse[ai]) {
          return false;
        }
        if (!dummyAts[ai]) {
          ringOk = true;
        }
      }
      return ringOk;
    };
    // we can't just copy the rings over: we're going to rearrange them so that
    // we try to favor starting the traversal of any ring from an atom that is
    // at the end of a wedged ring bond. This is part of our attempt to avoid
    // assigning double bonds to bonds with wedging
    for (const auto &ring : allrings) {
      if (containsNonDummy(ring)) {
        unsigned int startPos = 0;
        bool hasWedge = false;
        for (auto ri = 0u; ri < ring.size(); ++ri) {
          if (wedgedAtoms[ring[ri]]) {
            startPos = ri;
            hasWedge = true;
            break;
          }
        }
        INT_VECT nring(ring.size());
        for (auto ri = 0u; ri < ring.size(); ++ri) {
          nring[ri] = ring.at((ri + startPos) % ring.size());
        }
        if (!hasWedge) {
          tmpRings.push_back(nring);
        } else {
          tmpRings.push_front(nring);
        }
      }
    }
    VECT_INT_VECT arings;
    arings.reserve(allrings.size());
    arings.insert(arings.end(), tmpRings.begin(), tmpRings.end());
    VECT_INT_VECT allbrings;
    RingUtils::convertToBonds(arings, allbrings, mol);
    VECT_INT_VECT brings;
    brings.reserve(allbrings.size());
    auto copyBondRingsWithinFragment = [&bondsToUse](const INT_VECT &ring) {
      return std::all_of(ring.begin(), ring.end(), [&bondsToUse](const int bi) {
        return bondsToUse[bi];
      });
    };
    VECT_INT_VECT aringsRemaining;
    aringsRemaining.reserve(arings.size());
    for (unsigned i = 0; i < allbrings.size(); ++i) {
      if (copyBondRingsWithinFragment(allbrings[i])) {
        brings.push_back(allbrings[i]);
        aringsRemaining.push_back(arings[i]);
      }
    }
    arings = std::move(aringsRemaining);

    // make a neighbor map for the rings i.e. a ring is a
    // neighbor to another candidate ring if it shares at least
    // one bond
    // useful to figure out fused systems
    INT_INT_VECT_MAP neighMap;
    RingUtils::makeRingNeighborMap(brings, neighMap);

    int curr = 0;
    int cnrs = rdcast<int>(arings.size());
    boost::dynamic_bitset<> fusDone(cnrs);
    while (curr < cnrs) {
      INT_VECT fused;
      RingUtils::pickFusedRings(curr, neighMap, fused, fusDone);
      VECT_INT_VECT frings(fused.size());
      std::transform(fused.begin(), fused.end(), frings.begin(),
                     [&arings](const int ri) { return arings[ri]; });
      kekulizeFused(mol, frings, maxBackTracks);
      int rix;
      for (rix = 0; rix < cnrs; ++rix) {
        if (!fusDone[rix]) {
          curr = rix;
          break;
        }
      }
      if (rix == cnrs) {
        break;
      }
    }
  }
  if (markAtomsBonds) {
    // if we want the atoms and bonds to be marked non-aromatic do
    // that here.
    if (!mol.getRingInfo()->isInitialized()) {
      MolOps::findSSSR(mol);
    }
    for (auto bond : mol.bonds()) {
      if (bondsToUse[bond->getIdx()]) {
        bond->setIsAromatic(false);
      }
    }
    for (auto atom : mol.atoms()) {
      if (atomsToUse[atom->getIdx()] && atom->getIsAromatic()) {
        // if we're doing the full molecule and there are aromatic atoms not in
        // a ring, throw an exception
        if (atomsToUse.all() && bondsToUse.all() &&
            !mol.getRingInfo()->numAtomRings(atom->getIdx())) {
          std::ostringstream errout;
          errout << "non-ring atom " << atom->getIdx() << " marked aromatic";
          auto msg = errout.str();
          BOOST_LOG(rdErrorLog) << msg << std::endl;
          throw AtomKekulizeException(msg, atom->getIdx());
        }
        atom->setIsAromatic(false);
        // make sure "explicit" Hs on things like pyrroles don't hang around
        // this was Github Issue 141
        if ((atom->getAtomicNum() == 7 || atom->getAtomicNum() == 15) &&
            atom->getFormalCharge() == 0 && atom->getNumExplicitHs() == 1) {
          atom->setNoImplicit(false);
          atom->setNumExplicitHs(0);
          atom->updatePropertyCache(false);
        }
      }
    }
  }

  // ok some error checking here force a implicit valence
  // calculation that should do some error checking by itself. In
  // addition compare them to what they were before kekulizing
  for (auto atom : mol.atoms()) {
    if (!atomsToUse[atom->getIdx()]) {
      continue;
    }
    int val = atom->getTotalValence();
    if (val != valences[atom->getIdx()]) {
      std::ostringstream errout;
      errout << "Kekulization somehow screwed up valence on " << atom->getIdx()
             << ": " << val << "!=" << valences[atom->getIdx()] << std::endl;
      auto msg = errout.str();
      BOOST_LOG(rdErrorLog) << msg << std::endl;
      throw AtomKekulizeException(msg, atom->getIdx());
    }
  }
}
}  // namespace details
void Kekulize(RWMol &mol, bool markAtomsBonds, unsigned int maxBackTracks) {
  boost::dynamic_bitset<> atomsToUse(mol.getNumAtoms());
  atomsToUse.set();
  boost::dynamic_bitset<> bondsToUse(mol.getNumBonds());
  bondsToUse.set();
  details::KekulizeFragment(mol, atomsToUse, bondsToUse, markAtomsBonds,
                            maxBackTracks);
}
bool KekulizeIfPossible(RWMol &mol, bool markAtomsBonds,
                        unsigned int maxBackTracks) {
  boost::dynamic_bitset<> aromaticBonds(mol.getNumBonds());
  for (const auto bond : mol.bonds()) {
    if (bond->getIsAromatic()) {
      aromaticBonds.set(bond->getIdx());
    }
  }
  boost::dynamic_bitset<> aromaticAtoms(mol.getNumAtoms());
  for (const auto atom : mol.atoms()) {
    if (isAromaticAtom(*atom)) {
      aromaticAtoms.set(atom->getIdx());
    }
  }
  bool res = true;
  try {
    Kekulize(mol, markAtomsBonds, maxBackTracks);
  } catch (const MolSanitizeException &) {
    res = false;
    for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {
      if (aromaticBonds[i]) {
        auto bond = mol.getBondWithIdx(i);
        bond->setIsAromatic(true);
        bond->setBondType(Bond::BondType::AROMATIC);
      }
    }
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      if (aromaticAtoms[i]) {
        mol.getAtomWithIdx(i)->setIsAromatic(true);
      }
    }
  }
  return res;
}
}  // namespace MolOps
}  // namespace RDKit
