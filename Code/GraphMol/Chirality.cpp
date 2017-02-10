//
//  Copyright (C) 2004-2017 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/Ranking.h>
#include <GraphMol/new_canon.h>
#include <RDGeneral/types.h>
#include <sstream>
#include <set>
#include <algorithm>
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>

#include <boost/dynamic_bitset.hpp>
#include <Geometry/point.h>

// #define VERBOSE_CANON 1

namespace RDKit {
namespace Chirality {
typedef std::pair<int, int> INT_PAIR;
typedef std::vector<INT_PAIR> INT_PAIR_VECT;
typedef std::vector<INT_PAIR>::iterator INT_PAIR_VECT_I;
typedef std::vector<INT_PAIR>::const_iterator INT_PAIR_VECT_CI;

typedef INT_VECT CIP_ENTRY;
typedef std::vector<CIP_ENTRY> CIP_ENTRY_VECT;

template <typename T>
void debugVect(const std::vector<T> arg) {
  typename std::vector<T>::const_iterator viIt;
  std::stringstream outS;
  for (viIt = arg.begin(); viIt != arg.end(); viIt++) {
    outS << *viIt << " ";
  }
  BOOST_LOG(rdDebugLog) << outS.str() << std::endl;
}

// --------------------------------------------------
//
// Calculates chiral invariants for the atoms of a molecule
//  These are based on Labute's proposal in:
//  "An Efficient Algorithm for the Determination of Topological
//   RS Chirality" Journal of the CCG (1996)
//
// --------------------------------------------------
void buildCIPInvariants(const ROMol &mol, DOUBLE_VECT &res) {
  PRECONDITION(res.size() >= mol.getNumAtoms(), "res vect too small");
  int atsSoFar = 0;
  //
  // NOTE:
  // If you make modifications to this, keep in mind that it is
  // essential that the initial comparison of ranks behave properly.
  // So, though it seems like it would makes sense to include
  // information about the number of Hs (or charge, etc) in the CIP
  // invariants, this will result in bad rankings.  For example, in
  // this molecule: OC[C@H](C)O, including the number of Hs would
  // cause the methyl group (atom 3) to be ranked higher than the CH2
  // connected to O (atom 1).  This is totally wrong.
  //
  // We also don't include any pre-existing stereochemistry information.
  // Though R and S assignments do factor in to the priorities of atoms,
  // we're starting here from scratch and we'll let the R and S stuff
  // be taken into account during the iterations.
  //
  for (ROMol::ConstAtomIterator atIt = mol.beginAtoms(); atIt != mol.endAtoms();
       ++atIt) {
    const unsigned short nMassBits = 10;
    const unsigned short maxMass = 1 << nMassBits;
    Atom const *atom = *atIt;
    unsigned long invariant = 0;
    int num = atom->getAtomicNum() % 128;
    // get an int with the deviation in the mass from the default:
    int mass = 0;
    if (atom->getIsotope()) {
      mass =
          atom->getIsotope() -
          PeriodicTable::getTable()->getMostCommonIsotope(atom->getAtomicNum());
      if (mass >= 0) mass += 1;
    }
    mass += maxMass / 2;
    if (mass < 0)
      mass = 0;
    else
      mass = mass % maxMass;

#if 0
        // NOTE: the inclusion of hybridization in the invariant (as
        // suggested in the original paper), leads to the situation
        // that
        //   C[C@@](O)(C=C)C(C)CC
        // and
        //   C[C@@](O)(C=C)C(C)CO
        // are assigned S chirality even though the rest of the world
        // seems to agree that they ought to be R (atom 3, sp2, is ranked
        // higher than atom 5, sp3, no matter what their environments)
        int hyb=0;
        switch(atom->getHybridization()) {
        case Atom::SP: hyb=6;break;
        case Atom::SP2: hyb=5;break;
        case Atom::SP3: hyb=1;break;
        case Atom::SP3D: hyb=3;break;
        case Atom::SP3D2: hyb=2;break;
        default: break;
        }
#endif

    invariant = num;  // 7 bits here
    invariant = (invariant << nMassBits) | mass;

    int mapnum = -1;
    atom->getPropIfPresent(common_properties::molAtomMapNumber, mapnum);
    mapnum = (mapnum + 1) % 1024;  // increment to allow map numbers of zero
                                   // (though that would be stupid)
    invariant = (invariant << 10) | mapnum;

    res[atsSoFar++] = invariant;
  }
}

void iterateCIPRanks(const ROMol &mol, DOUBLE_VECT &invars, UINT_VECT &ranks,
                     bool seedWithInvars) {
  PRECONDITION(invars.size() == mol.getNumAtoms(), "bad invars size");
  PRECONDITION(ranks.size() >= mol.getNumAtoms(), "bad ranks size");

  unsigned int numAtoms = mol.getNumAtoms();
  CIP_ENTRY_VECT cipEntries(numAtoms);
  INT_LIST allIndices;
  for (unsigned int i = 0; i < numAtoms; ++i) {
    allIndices.push_back(i);
  }
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << "invariants:" << std::endl;
  for (unsigned int i = 0; i < numAtoms; i++) {
    BOOST_LOG(rdDebugLog) << i << ": " << invars[i] << std::endl;
  }
#endif

  // rank those:
  Rankers::rankVect(invars, ranks);
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << "initial ranks:" << std::endl;
  for (unsigned int i = 0; i < numAtoms; ++i) {
    BOOST_LOG(rdDebugLog) << i << ": " << ranks[i] << std::endl;
  }
#endif
  // Start each atom's rank vector with its atomic number:
  //  Note: in general one should avoid the temptation to
  //  use invariants here, those lead to incorrect answers
  for (unsigned int i = 0; i < numAtoms; i++) {
    if (!seedWithInvars) {
      cipEntries[i].push_back(mol[i]->getAtomicNum());
      cipEntries[i].push_back(static_cast<int>(ranks[i]));
    } else {
      cipEntries[i].push_back(static_cast<int>(invars[i]));
    }
  }

  // Loop until either:
  //   1) all classes are uniquified
  //   2) the number of ranks doesn't change from one iteration to
  //      the next
  //   3) we've gone through maxIts times
  //      maxIts is calculated by dividing the number of atoms
  //      by 2. That's a pessimal version of the
  //      maximum number of steps required for two atoms to
  //      "feel" each other (each influences one additional
  //      neighbor shell per iteration).
  unsigned int maxIts = numAtoms / 2 + 1;
  unsigned int numIts = 0;
  int lastNumRanks = -1;
  unsigned int numRanks = *std::max_element(ranks.begin(), ranks.end()) + 1;
  while (numRanks < numAtoms && numIts < maxIts &&
         (lastNumRanks < 0 ||
          static_cast<unsigned int>(lastNumRanks) < numRanks)) {
    unsigned int longestEntry = 0;
    // ----------------------------------------------------
    //
    // for each atom, get a sorted list of its neighbors' ranks:
    //
    for (INT_LIST_I it = allIndices.begin(); it != allIndices.end(); ++it) {
      CIP_ENTRY localEntry;
      localEntry.reserve(16);

      // start by pushing on our neighbors' ranks:
      ROMol::OEDGE_ITER beg, end;
      boost::tie(beg, end) = mol.getAtomBonds(mol[*it].get());
      while (beg != end) {
        const Bond *bond = mol[*beg].get();
        ++beg;
        unsigned int nbrIdx = bond->getOtherAtomIdx(*it);
        const Atom *nbr = mol[nbrIdx].get();

        int rank = ranks[nbrIdx] + 1;
        // put the neighbor in 2N times where N is the bond order as a double.
        // this is to treat aromatic linkages on fair footing. i.e. at least in
        // the
        // first iteration --c(:c):c and --C(=C)-C should look the same.
        // this was part of issue 3009911

        unsigned int count;
        if (bond->getBondType() == Bond::DOUBLE && nbr->getAtomicNum() == 15 &&
            (nbr->getDegree() == 4 || nbr->getDegree() == 3)) {
          // a special case for chiral phophorous compounds
          // (this was leading to incorrect assignment of
          // R/S labels ):
          count = 1;

          // general justification of this is:
          // Paragraph 2.2. in the 1966 article is "Valence-Bond Conventions:
          // Multiple-Bond Unsaturation and Aromaticity". It contains several
          // conventions of which convention (b) is the one applying here:
          // "(b) Contibutions by d orbitals to bonds of quadriligant atoms are
          // neglected."
          // FIX: this applies to more than just P
        } else {
          count = static_cast<unsigned int>(
              floor(2. * bond->getBondTypeAsDouble() + .1));
        }
        CIP_ENTRY::iterator ePos =
            std::lower_bound(localEntry.begin(), localEntry.end(), rank);
        localEntry.insert(ePos, count, rank);
        ++nbr;
      }
      // add a zero for each coordinated H:
      // (as long as we're not a query atom)
      if (!mol[*it]->hasQuery()) {
        localEntry.insert(localEntry.begin(), mol[*it]->getTotalNumHs(), 0);
      }

      // we now have a sorted list of our neighbors' ranks,
      // copy it on in reversed order:
      cipEntries[*it].insert(cipEntries[*it].end(), localEntry.rbegin(),
                             localEntry.rend());
      if (cipEntries[*it].size() > longestEntry) {
        longestEntry = rdcast<unsigned int>(cipEntries[*it].size());
      }
    }
    // ----------------------------------------------------
    //
    // pad the entries so that we compare rounds to themselves:
    //
    for (INT_LIST_I it = allIndices.begin(); it != allIndices.end(); ++it) {
      unsigned int sz = rdcast<unsigned int>(cipEntries[*it].size());
      if (sz < longestEntry) {
        cipEntries[*it].insert(cipEntries[*it].end(), longestEntry - sz, -1);
      }
    }
    // ----------------------------------------------------
    //
    // sort the new ranks and update the list of active indices:
    //
    lastNumRanks = numRanks;

    Rankers::rankVect(cipEntries, ranks);
    numRanks = *std::max_element(ranks.begin(), ranks.end()) + 1;

    // now truncate each vector and stick the rank at the end
    for (unsigned int i = 0; i < numAtoms; ++i) {
      cipEntries[i][numIts + 1] = ranks[i];
      cipEntries[i].erase(cipEntries[i].begin() + numIts + 2,
                          cipEntries[i].end());
    }

    ++numIts;
#ifdef VERBOSE_CANON
    BOOST_LOG(rdDebugLog) << "strings and ranks:" << std::endl;
    for (unsigned int i = 0; i < numAtoms; i++) {
      BOOST_LOG(rdDebugLog) << i << ": " << ranks[i] << " > ";
      debugVect(cipEntries[i]);
    }
#endif
  }
}
// Figure out the CIP ranks for the atoms of a molecule
void assignAtomCIPRanks(const ROMol &mol, UINT_VECT &ranks) {
  PRECONDITION((!ranks.size() || ranks.size() >= mol.getNumAtoms()),
               "bad ranks size");
  if (!ranks.size()) ranks.resize(mol.getNumAtoms());
  unsigned int numAtoms = mol.getNumAtoms();
#ifndef USE_NEW_STEREOCHEMISTRY
  // get the initial invariants:
  DOUBLE_VECT invars(numAtoms, 0);
  buildCIPInvariants(mol, invars);
  iterateCIPRanks(mol, invars, ranks, false);
#else
  Canon::chiralRankMolAtoms(mol, ranks);
#endif

  // copy the ranks onto the atoms:
  for (unsigned int i = 0; i < numAtoms; ++i) {
    mol[i]->setProp(common_properties::_CIPRank, ranks[i], 1);
  }
}

// construct a vector with <atomIdx,direction> pairs for
// neighbors of a given atom.  This list will only be
// non-empty if at least one of the bonds has its direction
// set.
void findAtomNeighborDirHelper(const ROMol &mol, const Atom *atom,
                               const Bond *refBond, UINT_VECT &ranks,
                               INT_PAIR_VECT &neighbors,
                               bool &hasExplicitUnknownStereo) {
  PRECONDITION(atom, "bad atom");
  PRECONDITION(refBond, "bad bond");

  bool seenDir = false;
  ROMol::OEDGE_ITER beg, end;
  boost::tie(beg, end) = mol.getAtomBonds(atom);
  while (beg != end) {
    const BOND_SPTR bond = mol[*beg];
    // check whether this bond is explictly set to have unknown stereo
    if (!hasExplicitUnknownStereo) {
      int explicit_unknown_stereo;
      if (bond->getBondDir() == Bond::UNKNOWN  // there's a squiggle bond
          || (bond->getPropIfPresent<int>(common_properties::_UnknownStereo,
                                          explicit_unknown_stereo) &&
              explicit_unknown_stereo))
        hasExplicitUnknownStereo = true;
    }

    Bond::BondDir dir = bond->getBondDir();
    if (bond->getIdx() != refBond->getIdx()) {
      if (dir == Bond::ENDDOWNRIGHT || dir == Bond::ENDUPRIGHT) {
        seenDir = true;
        // If we're considering the bond "backwards", (i.e. from end
        // to beginning, reverse the effective direction:
        if (atom != bond->getBeginAtom()) {
          if (dir == Bond::ENDDOWNRIGHT)
            dir = Bond::ENDUPRIGHT;
          else
            dir = Bond::ENDDOWNRIGHT;
        }
      }
      Atom *nbrAtom = bond->getOtherAtom(atom);
      neighbors.push_back(std::make_pair(nbrAtom->getIdx(), dir));
    }
    ++beg;
  }
  if (!seenDir) {
    neighbors.clear();
  } else {
    if (neighbors.size() == 2 &&
        ranks[neighbors[0].first] == ranks[neighbors[1].first]) {
      // the two substituents are identical, no stereochemistry here:
      neighbors.clear();
    } else {
      // it's possible that direction was set only one of the bonds, set the
      // other
      // bond's direction to be reversed:
      if (neighbors[0].second != Bond::ENDDOWNRIGHT &&
          neighbors[0].second != Bond::ENDUPRIGHT) {
        CHECK_INVARIANT(neighbors.size() > 1, "too few neighbors");
        neighbors[0].second = neighbors[1].second == Bond::ENDDOWNRIGHT
                                  ? Bond::ENDUPRIGHT
                                  : Bond::ENDDOWNRIGHT;
      } else if (neighbors.size() > 1 &&
                 neighbors[1].second != Bond::ENDDOWNRIGHT &&
                 neighbors[1].second != Bond::ENDUPRIGHT) {
        neighbors[1].second = neighbors[0].second == Bond::ENDDOWNRIGHT
                                  ? Bond::ENDUPRIGHT
                                  : Bond::ENDDOWNRIGHT;
      }
    }
  }
}

// find the neighbors for an atoms that are not connected by single bond that is
// not refBond
// if checkDir is true only neighbor atoms with bonds marked with a direction
// will be returned
void findAtomNeighborsHelper(const ROMol &mol, const Atom *atom,
                             const Bond *refBond, UINT_VECT &neighbors,
                             bool checkDir = false) {
  PRECONDITION(atom, "bad atom");
  PRECONDITION(refBond, "bad bond");
  neighbors.clear();
  ROMol::OEDGE_ITER beg, end;
  boost::tie(beg, end) = mol.getAtomBonds(atom);
  while (beg != end) {
    const BOND_SPTR bond = mol[*beg];
    Bond::BondDir dir = bond->getBondDir();
    if (bond->getBondType() == Bond::SINGLE &&
        bond->getIdx() != refBond->getIdx()) {
      if (checkDir) {
        if ((dir != Bond::ENDDOWNRIGHT) && (dir != Bond::ENDUPRIGHT)) {
          ++beg;
          continue;
        }
      }
      Atom *nbrAtom = bond->getOtherAtom(atom);
      neighbors.push_back(nbrAtom->getIdx());
    }
    ++beg;
  }
}

// conditions for an atom to be a candidate for ring stereochem:
//   1) two non-ring neighbors that have different ranks
//   2) one non-ring neighbor and two ring neighbors (the ring neighbors will
//      have the same rank)
//   3) four ring neighbors with three different ranks
//   4) three ring neighbors with two different ranks
//     example for this last one: C[C@H]1CC2CCCC3CCCC(C1)[C@@H]23
bool atomIsCandidateForRingStereochem(const ROMol &mol, const Atom *atom) {
  PRECONDITION(atom, "bad atom");
  bool res = false;
  std::set<unsigned int> nbrRanks;
  if (!atom->getPropIfPresent(common_properties::_ringStereochemCand, res)) {
    const RingInfo *ringInfo = mol.getRingInfo();
    if (ringInfo->isInitialized() && ringInfo->numAtomRings(atom->getIdx())) {
      ROMol::OEDGE_ITER beg, end;
      boost::tie(beg, end) = mol.getAtomBonds(atom);
      std::vector<const Atom *> nonRingNbrs;
      std::vector<const Atom *> ringNbrs;
      while (beg != end) {
        const BOND_SPTR bond = mol[*beg];
        if (!ringInfo->numBondRings(bond->getIdx())) {
          nonRingNbrs.push_back(bond->getOtherAtom(atom));
        } else {
          const Atom *nbr = bond->getOtherAtom(atom);
          ringNbrs.push_back(nbr);
          unsigned int rnk = 0;
          nbr->getPropIfPresent(common_properties::_CIPRank, rnk);
          nbrRanks.insert(rnk);
        }
        ++beg;
      }

      unsigned int rank1 = 0, rank2 = 0;
      switch (nonRingNbrs.size()) {
        case 2:
          if (nonRingNbrs[0]->getPropIfPresent(common_properties::_CIPRank,
                                               rank1) &&
              nonRingNbrs[1]->getPropIfPresent(common_properties::_CIPRank,
                                               rank2)) {
            if (rank1 == rank2) {
              res = false;
            } else {
              res = true;
            }
          }
          break;
        case 1:
          if (ringNbrs.size() == 2) res = true;
          break;
        case 0:
          if (ringNbrs.size() == 4 && nbrRanks.size() == 3) {
            res = true;
          } else if (ringNbrs.size() == 3 && nbrRanks.size() == 2) {
            res = true;
          } else {
            res = false;
          }
          break;
        default:
          res = false;
      }
    }
    atom->setProp(common_properties::_ringStereochemCand, res, 1);
  }
  return res;
}

// finds all possible chiral special cases.
// at the moment this is just candidates for ring stereochemistry
void findChiralAtomSpecialCases(ROMol &mol,
                                boost::dynamic_bitset<> &possibleSpecialCases) {
  PRECONDITION(possibleSpecialCases.size() >= mol.getNumAtoms(),
               "bit vector too small");
  possibleSpecialCases.reset();
  if (!mol.getRingInfo()->isInitialized()) {
    VECT_INT_VECT sssrs;
    MolOps::symmetrizeSSSR(mol, sssrs);
  }
  boost::dynamic_bitset<> atomsSeen(mol.getNumAtoms());
  boost::dynamic_bitset<> atomsUsed(mol.getNumAtoms());
  boost::dynamic_bitset<> bondsSeen(mol.getNumBonds());

  for (ROMol::AtomIterator ait = mol.beginAtoms(); ait != mol.endAtoms();
       ++ait) {
    const Atom *atom = *ait;
    if (atomsSeen[atom->getIdx()]) continue;
    if (atom->getChiralTag() == Atom::CHI_UNSPECIFIED ||
        atom->hasProp(common_properties::_CIPCode) ||
        !mol.getRingInfo()->numAtomRings(atom->getIdx()) ||
        !atomIsCandidateForRingStereochem(mol, atom)) {
      continue;
    }
    // do a BFS from this ring atom along ring bonds and find other
    // stereochemistry candidates.
    std::list<const Atom *> nextAtoms;
    // start with finding viable neighbors
    ROMol::OEDGE_ITER beg, end;
    boost::tie(beg, end) = mol.getAtomBonds(atom);
    while (beg != end) {
      unsigned int bidx = mol[*beg]->getIdx();
      if (!bondsSeen[bidx]) {
        bondsSeen.set(bidx);
        if (mol.getRingInfo()->numBondRings(bidx)) {
          const Atom *oatom = mol[*beg]->getOtherAtom(atom);
          if (!atomsSeen[oatom->getIdx()]) {
            nextAtoms.push_back(oatom);
            atomsUsed.set(oatom->getIdx());
          }
        }
      }
      ++beg;
    }
    INT_VECT ringStereoAtoms(0);
    if (!nextAtoms.empty()) {
      atom->getPropIfPresent(common_properties::_ringStereoAtoms,
                             ringStereoAtoms);
    }

    while (!nextAtoms.empty()) {
      const Atom *ratom = nextAtoms.front();
      nextAtoms.pop_front();
      atomsSeen.set(ratom->getIdx());
      if (ratom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
          !ratom->hasProp(common_properties::_CIPCode) &&
          atomIsCandidateForRingStereochem(mol, ratom)) {
        int same = (ratom->getChiralTag() == atom->getChiralTag()) ? 1 : -1;
        ringStereoAtoms.push_back(same * (ratom->getIdx() + 1));
        INT_VECT oringatoms(0);
        ratom->getPropIfPresent(common_properties::_ringStereoAtoms,
                                oringatoms);
        oringatoms.push_back(same * (atom->getIdx() + 1));
        ratom->setProp(common_properties::_ringStereoAtoms, oringatoms, true);
        possibleSpecialCases.set(ratom->getIdx());
        possibleSpecialCases.set(atom->getIdx());
      }
      // now push this atom's neighbors
      boost::tie(beg, end) = mol.getAtomBonds(ratom);
      while (beg != end) {
        unsigned int bidx = mol[*beg]->getIdx();
        if (!bondsSeen[bidx]) {
          bondsSeen.set(bidx);
          if (mol.getRingInfo()->numBondRings(bidx)) {
            const Atom *oatom = mol[*beg]->getOtherAtom(ratom);
            if (!atomsSeen[oatom->getIdx()] && !atomsUsed[oatom->getIdx()]) {
              nextAtoms.push_back(oatom);
              atomsUsed.set(oatom->getIdx());
            }
          }
        }
        ++beg;
      }
    }  // end of BFS
    if (ringStereoAtoms.size() != 0) {
      atom->setProp(common_properties::_ringStereoAtoms, ringStereoAtoms, true);
      // because we're only going to hit each ring atom once, the first atom we
      // encounter in a ring is going to end up with all the other atoms set as
      // stereoAtoms, but each of them will only have the first atom present. We
      // need to fix that. because the traverse from the first atom only
      // followed ring bonds, these things are all by definition in one ring
      // system. (Q: is this true if there's a spiro center in there?)
      INT_VECT same(mol.getNumAtoms(), 0);
      BOOST_FOREACH (int ringAtomEntry, ringStereoAtoms) {
        int ringAtomIdx =
            ringAtomEntry < 0 ? -ringAtomEntry - 1 : ringAtomEntry - 1;
        same[ringAtomIdx] = ringAtomEntry;
      }
      for (INT_VECT_CI rae = ringStereoAtoms.begin();
           rae != ringStereoAtoms.end(); ++rae) {
        int ringAtomEntry = *rae;
        int ringAtomIdx =
            ringAtomEntry < 0 ? -ringAtomEntry - 1 : ringAtomEntry - 1;
        INT_VECT lringatoms(0);
        mol.getAtomWithIdx(ringAtomIdx)
            ->getPropIfPresent(common_properties::_ringStereoAtoms, lringatoms);
        CHECK_INVARIANT(lringatoms.size() > 0, "no other ring atoms found.");
        for (INT_VECT_CI orae = rae + 1; orae != ringStereoAtoms.end();
             ++orae) {
          int oringAtomEntry = *orae;
          int oringAtomIdx =
              oringAtomEntry < 0 ? -oringAtomEntry - 1 : oringAtomEntry - 1;
          int theseDifferent = (ringAtomEntry < 0) ^ (oringAtomEntry < 0);
          lringatoms.push_back(theseDifferent ? -(oringAtomIdx + 1)
                                              : (oringAtomIdx + 1));
          INT_VECT olringatoms(0);
          mol.getAtomWithIdx(oringAtomIdx)
              ->getPropIfPresent(common_properties::_ringStereoAtoms,
                                 olringatoms);
          CHECK_INVARIANT(olringatoms.size() > 0, "no other ring atoms found.");
          olringatoms.push_back(theseDifferent ? -(ringAtomIdx + 1)
                                               : (ringAtomIdx + 1));
          mol.getAtomWithIdx(oringAtomIdx)
              ->setProp(common_properties::_ringStereoAtoms, olringatoms);
        }
        mol.getAtomWithIdx(ringAtomIdx)
            ->setProp(common_properties::_ringStereoAtoms, lringatoms);
      }

    } else {
      possibleSpecialCases.reset(atom->getIdx());
    }
    atomsSeen.set(atom->getIdx());
  }
}

std::pair<bool, bool> isAtomPotentialChiralCenter(
    const Atom *atom, const ROMol &mol, const UINT_VECT &ranks,
    Chirality::INT_PAIR_VECT &nbrs) {
  // loop over all neighbors and form a decorated list of their
  // ranks:
  bool legalCenter = true;
  bool hasDupes = false;

  if (atom->getTotalDegree() > 4) {
    // we only know tetrahedral chirality
    legalCenter = false;
  } else {
    boost::dynamic_bitset<> codesSeen(mol.getNumAtoms());
    ROMol::OEDGE_ITER beg, end;
    boost::tie(beg, end) = mol.getAtomBonds(atom);
    while (beg != end) {
      unsigned int otherIdx = mol[*beg]->getOtherAtom(atom)->getIdx();
      CHECK_INVARIANT(ranks[otherIdx] < mol.getNumAtoms(),
                      "CIP rank higher than the number of atoms.");
      // watch for neighbors with duplicate ranks, which would mean
      // that we cannot be chiral:
      if (codesSeen[ranks[otherIdx]]) {
        // we've already seen this code, it's a dupe
        hasDupes = true;
        break;
      }
      codesSeen[ranks[otherIdx]] = 1;
      nbrs.push_back(std::make_pair(ranks[otherIdx], mol[*beg]->getIdx()));
      ++beg;
    }

    // figure out if this is a legal chiral center or not:
    if (!hasDupes) {
      if (nbrs.size() < 3) {
        // less than three neighbors is never stereogenic
        legalCenter = false;
      } else if (nbrs.size() == 3) {
        // three-coordinate with a single H we'll accept automatically:
        if (atom->getTotalNumHs() != 1) {
          // otherwise we default to not being a legal center
          legalCenter = false;
          // but there are a few special cases we'll accept
          // sulfur or selenium with either a positive charge or a double
          // bond:
          if ((atom->getAtomicNum() == 16 || atom->getAtomicNum() == 34) &&
              (atom->getExplicitValence() == 4 ||
               (atom->getExplicitValence() == 3 &&
                atom->getFormalCharge() == 1))) {
            legalCenter = true;
          }
        }
      }
    }
  }
  return std::make_pair(legalCenter, hasDupes);
}

// returns a pair:
//   1) are there unassigned stereoatoms
//   2) did we assign any?
std::pair<bool, bool> assignAtomChiralCodes(ROMol &mol, UINT_VECT &ranks,
                                            bool flagPossibleStereoCenters) {
  PRECONDITION((!ranks.size() || ranks.size() == mol.getNumAtoms()),
               "bad rank vector size");
  bool atomChanged = false;
  unsigned int unassignedAtoms = 0;

  // ------------------
  // now loop over each atom and, if it's marked as chiral,
  //  figure out the appropriate CIP label:
  for (ROMol::AtomIterator atIt = mol.beginAtoms(); atIt != mol.endAtoms();
       ++atIt) {
    Atom *atom = *atIt;
    Atom::ChiralType tag = atom->getChiralTag();

    // only worry about this atom if it has a marked chirality
    // we understand:
    if (flagPossibleStereoCenters ||
        (tag != Atom::CHI_UNSPECIFIED && tag != Atom::CHI_OTHER)) {
      if (atom->hasProp(common_properties::_CIPCode)) {
        continue;
      }

      if (!ranks.size()) {
        //  if we need to, get the "CIP" ranking of each atom:
        assignAtomCIPRanks(mol, ranks);
      }
      Chirality::INT_PAIR_VECT nbrs;
      bool legalCenter, hasDupes;
      boost::tie(legalCenter, hasDupes) =
          isAtomPotentialChiralCenter(atom, mol, ranks, nbrs);
      if (legalCenter) {
        ++unassignedAtoms;
      }
      if (legalCenter && !hasDupes && flagPossibleStereoCenters) {
        atom->setProp(common_properties::_ChiralityPossible, 1);
      }

      if (legalCenter && !hasDupes && tag != Atom::CHI_UNSPECIFIED &&
          tag != Atom::CHI_OTHER) {
        // stereochem is possible and we have no duplicate neighbors, assign
        // a CIP code:
        atomChanged = true;
        --unassignedAtoms;

        // sort the list of neighbors by their CIP ranks:
        std::sort(nbrs.begin(), nbrs.end(), Rankers::pairLess<int, int>());

        // collect the list of neighbor indices:
        std::list<int> nbrIndices;
        for (Chirality::INT_PAIR_VECT_CI nbrIt = nbrs.begin();
             nbrIt != nbrs.end(); ++nbrIt) {
          nbrIndices.push_back((*nbrIt).second);
        }
        // ask the atom how many swaps we have to make:
        int nSwaps = atom->getPerturbationOrder(nbrIndices);

        // if the atom has 3 neighbors and a hydrogen, add a swap:
        if (nbrIndices.size() == 3 && atom->getTotalNumHs() == 1) {
          ++nSwaps;
        }

        // if that number is odd, we'll change our chirality:
        if (nSwaps % 2) {
          if (tag == Atom::CHI_TETRAHEDRAL_CCW)
            tag = Atom::CHI_TETRAHEDRAL_CW;
          else
            tag = Atom::CHI_TETRAHEDRAL_CCW;
        }
        // now assign the CIP code:
        std::string cipCode;
        if (tag == Atom::CHI_TETRAHEDRAL_CCW)
          cipCode = "S";
        else
          cipCode = "R";
        atom->setProp(common_properties::_CIPCode, cipCode);
      }
    }
  }
  return std::make_pair((unassignedAtoms > 0), atomChanged);
}

// returns a pair:
//   1) are there unassigned stereo bonds?
//   2) did we assign any?
std::pair<bool, bool> assignBondStereoCodes(ROMol &mol, UINT_VECT &ranks) {
  PRECONDITION((!ranks.size() || ranks.size() == mol.getNumAtoms()),
               "bad rank vector size");
  bool assignedABond = false;
  unsigned int unassignedBonds = 0;

  // find the double bonds:
  for (ROMol::BondIterator bondIt = mol.beginBonds(); bondIt != mol.endBonds();
       ++bondIt) {
    if ((*bondIt)->getBondType() == Bond::DOUBLE) {
      Bond *dblBond = *bondIt;
      if (dblBond->getStereo() != Bond::STEREONONE) {
        continue;
      }
      if (!ranks.size()) {
        assignAtomCIPRanks(mol, ranks);
      }
      dblBond->getStereoAtoms().clear();

      // at the moment we are ignoring stereochem on ring bonds with less than
      // 8
      // members.
      if (!mol.getRingInfo()->numBondRings(dblBond->getIdx()) ||
          mol.getRingInfo()->minBondRingSize(dblBond->getIdx()) > 7) {
        const Atom *begAtom = dblBond->getBeginAtom();
        const Atom *endAtom = dblBond->getEndAtom();
        // we're only going to handle 2 or three coordinate atoms:
        if ((begAtom->getDegree() == 2 || begAtom->getDegree() == 3) &&
            (endAtom->getDegree() == 2 || endAtom->getDegree() == 3)) {
          ++unassignedBonds;

          // look around each atom and see if it has at least one bond with
          // direction marked:

          // the pairs here are: atomrank,bonddir
          Chirality::INT_PAIR_VECT begAtomNeighbors, endAtomNeighbors;
          bool hasExplicitUnknownStereo = false;
          int bgn_stereo = false, end_stereo = false;
          if ((dblBond->getBeginAtom()->getPropIfPresent(
                   common_properties::_UnknownStereo, bgn_stereo) &&
               bgn_stereo) ||
              (dblBond->getEndAtom()->getPropIfPresent(
                   common_properties::_UnknownStereo, end_stereo) &&
               end_stereo)) {
            hasExplicitUnknownStereo = true;
          }
          Chirality::findAtomNeighborDirHelper(mol, begAtom, dblBond, ranks,
                                               begAtomNeighbors,
                                               hasExplicitUnknownStereo);
          Chirality::findAtomNeighborDirHelper(mol, endAtom, dblBond, ranks,
                                               endAtomNeighbors,
                                               hasExplicitUnknownStereo);

          if (begAtomNeighbors.size() && endAtomNeighbors.size()) {
            // Each atom has at least one neighboring bond with marked
            // directionality.  Find the highest-ranked directionality
            // on each side:

            int begDir, endDir, endNbrAid, begNbrAid;
            if (begAtomNeighbors.size() == 1 ||
                ranks[begAtomNeighbors[0].first] >
                    ranks[begAtomNeighbors[1].first]) {
              begDir = begAtomNeighbors[0].second;
              begNbrAid = begAtomNeighbors[0].first;
            } else {
              begDir = begAtomNeighbors[1].second;
              begNbrAid = begAtomNeighbors[1].first;
            }
            if (endAtomNeighbors.size() == 1 ||
                ranks[endAtomNeighbors[0].first] >
                    ranks[endAtomNeighbors[1].first]) {
              endDir = endAtomNeighbors[0].second;
              endNbrAid = endAtomNeighbors[0].first;
            } else {
              endDir = endAtomNeighbors[1].second;
              endNbrAid = endAtomNeighbors[1].first;
            }
            dblBond->getStereoAtoms().push_back(begNbrAid);
            dblBond->getStereoAtoms().push_back(endNbrAid);
            if (hasExplicitUnknownStereo) {
              dblBond->setStereo(Bond::STEREOANY);
              assignedABond = true;
            } else if (begDir == endDir) {
              // In findAtomNeighborDirHelper, we've set up the
              // bond directions here so that they correspond to
              // having both single bonds START at the double bond.
              // This means that if the single bonds point in the same
              // direction, the bond is cis, "Z"
              dblBond->setStereo(Bond::STEREOZ);
              assignedABond = true;
            } else {
              dblBond->setStereo(Bond::STEREOE);
              assignedABond = true;
            }
            --unassignedBonds;
          }
        }
      }
    }
  }
  return std::make_pair(unassignedBonds > 0, assignedABond);
}

// reassign atom ranks by supplementing the current ranks
// with information about known chirality
void rerankAtoms(const ROMol &mol, UINT_VECT &ranks) {
  PRECONDITION(ranks.size() == mol.getNumAtoms(), "bad rank vector size");
  unsigned int factor = 100;
  while (factor < mol.getNumAtoms()) factor *= 10;

#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << "rerank PRE: " << std::endl;
  for (int i = 0; i < mol.getNumAtoms(); i++) {
    BOOST_LOG(rdDebugLog) << "  " << i << ": " << ranks[i] << std::endl;
  }
#endif

  DOUBLE_VECT invars(mol.getNumAtoms());
  // and now supplement them:
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    invars[i] = ranks[i] * factor;
    const Atom *atom = mol.getAtomWithIdx(i);
    // Priority order: R > S > nothing
    std::string cipCode;
    if (atom->getPropIfPresent(common_properties::_CIPCode, cipCode)) {
      if (cipCode == "S") {
        invars[i] += 10;
      } else if (cipCode == "R") {
        invars[i] += 20;
      }
    }
    ROMol::OEDGE_ITER beg, end;
    boost::tie(beg, end) = mol.getAtomBonds(atom);
    while (beg != end) {
      const BOND_SPTR oBond = mol[*beg];
      if (oBond->getBondType() == Bond::DOUBLE) {
        if (oBond->getStereo() == Bond::STEREOE) {
          invars[i] += 1;
        } else if (oBond->getStereo() == Bond::STEREOZ) {
          invars[i] += 2;
        }
      }
      ++beg;
    }
  }
  iterateCIPRanks(mol, invars, ranks, true);
  // copy the ranks onto the atoms:
  for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    mol.getAtomWithIdx(i)->setProp(common_properties::_CIPRank, ranks[i]);
  }

#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << "   post: " << std::endl;
  for (int i = 0; i < mol.getNumAtoms(); i++) {
    BOOST_LOG(rdDebugLog) << "  " << i << ": " << ranks[i] << std::endl;
  }
#endif
}
}  // end of chirality namespace

namespace MolOps {

/*
    We're going to do this iteratively:
      1) assign atom stereochemistry
      2) assign bond stereochemistry
      3) if there are still unresolved atoms or bonds
         repeat the above steps as necessary
 */
void assignStereochemistry(ROMol &mol, bool cleanIt, bool force,
                           bool flagPossibleStereoCenters) {
  if (!force && mol.hasProp(common_properties::_StereochemDone)) {
    return;
  }

  // later we're going to need ring information, get it now if we don't
  // have it already:
  if (!mol.getRingInfo()->isInitialized()) {
    MolOps::fastFindRings(mol);
  }

#if 0
  std::cerr << ">>>>>>>>>>>>>\n";
  std::cerr << "assign stereochem\n";
  mol.debugMol(std::cerr);
#endif

  // as part of the preparation, we'll loop over the atoms and
  // bonds to see if anything has stereochemistry
  // indicated. There's no point in doing the work here if there
  // are neither stereocenters nor bonds that we need to consider.
  // The exception to this is when flagPossibleStereoCenters is
  // true; then we always need to do the work
  bool hasStereoAtoms = flagPossibleStereoCenters;
  for (ROMol::AtomIterator atIt = mol.beginAtoms(); atIt != mol.endAtoms();
       ++atIt) {
    if (cleanIt) {
      if ((*atIt)->hasProp(common_properties::_CIPCode)) {
        (*atIt)->clearProp(common_properties::_CIPCode);
      }
      if ((*atIt)->hasProp(common_properties::_ChiralityPossible)) {
        (*atIt)->clearProp(common_properties::_ChiralityPossible);
      }
    }
    if (!hasStereoAtoms && (*atIt)->getChiralTag() != Atom::CHI_UNSPECIFIED &&
        (*atIt)->getChiralTag() != Atom::CHI_OTHER) {
      hasStereoAtoms = true;
    }
  }
  bool hasStereoBonds = false;
  for (ROMol::BondIterator bondIt = mol.beginBonds(); bondIt != mol.endBonds();
       ++bondIt) {
    if (cleanIt) {
      if ((*bondIt)->getBondType() == Bond::DOUBLE) {
        if ((*bondIt)->getBondDir() == Bond::EITHERDOUBLE) {
          (*bondIt)->setStereo(Bond::STEREOANY);
        } else if ((*bondIt)->getStereo() != Bond::STEREOANY) {
          (*bondIt)->setStereo(Bond::STEREONONE);
          (*bondIt)->getStereoAtoms().clear();
        }
      }
    }
    if (!hasStereoBonds && (*bondIt)->getBondType() == Bond::DOUBLE) {
      ROMol::OEDGE_ITER beg, end;
      boost::tie(beg, end) = mol.getAtomBonds((*bondIt)->getBeginAtom());
      while (!hasStereoBonds && beg != end) {
        const BOND_SPTR nbond = mol[*beg];
        ++beg;
        if (nbond->getBondDir() == Bond::ENDDOWNRIGHT ||
            nbond->getBondDir() == Bond::ENDUPRIGHT) {
          hasStereoBonds = true;
        }
      }
      boost::tie(beg, end) = mol.getAtomBonds((*bondIt)->getEndAtom());
      while (!hasStereoBonds && beg != end) {
        const BOND_SPTR nbond = mol[*beg];
        ++beg;
        if (nbond->getBondDir() == Bond::ENDDOWNRIGHT ||
            nbond->getBondDir() == Bond::ENDUPRIGHT) {
          hasStereoBonds = true;
        }
      }
    }
  }
  UINT_VECT atomRanks;
  bool keepGoing = hasStereoAtoms | hasStereoBonds;
  bool changedStereoAtoms, changedStereoBonds;
  while (keepGoing) {
    if (hasStereoAtoms) {
      boost::tie(hasStereoAtoms, changedStereoAtoms) =
          Chirality::assignAtomChiralCodes(mol, atomRanks,
                                           flagPossibleStereoCenters);
    } else {
      changedStereoAtoms = false;
    }
    if (hasStereoBonds) {
      boost::tie(hasStereoBonds, changedStereoBonds) =
          Chirality::assignBondStereoCodes(mol, atomRanks);
    } else {
      changedStereoBonds = false;
    }
    keepGoing = (hasStereoAtoms || hasStereoBonds) &&
                (changedStereoAtoms || changedStereoBonds);

    if (keepGoing) {
      // update the atom ranks based on the new information we have:
      Chirality::rerankAtoms(mol, atomRanks);
    }
#if 0
    std::cout << "*************** done iteration " << keepGoing
              << " ***********" << std::endl;
    mol.debugMol(std::cout);
    std::cout << "*************** done iteration " << keepGoing
              << " ***********" << std::endl;
#endif
  }

  if (cleanIt) {
    for (ROMol::AtomIterator atIt = mol.beginAtoms(); atIt != mol.endAtoms();
         ++atIt) {
      if ((*atIt)->hasProp(common_properties::_ringStereochemCand))
        (*atIt)->clearProp(common_properties::_ringStereochemCand);
      if ((*atIt)->hasProp(common_properties::_ringStereoAtoms))
        (*atIt)->clearProp(common_properties::_ringStereoAtoms);
    }
    boost::dynamic_bitset<> possibleSpecialCases(mol.getNumAtoms());
    Chirality::findChiralAtomSpecialCases(mol, possibleSpecialCases);

    for (ROMol::AtomIterator atIt = mol.beginAtoms(); atIt != mol.endAtoms();
         ++atIt) {
      Atom *atom = *atIt;
      if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
          !atom->hasProp(common_properties::_CIPCode) &&
          (possibleSpecialCases[atom->getIdx()] ||
           atom->hasProp(common_properties::_ringStereoAtoms))) {
      }

      if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
          !atom->hasProp(common_properties::_CIPCode) &&
          (!possibleSpecialCases[atom->getIdx()] ||
           !atom->hasProp(common_properties::_ringStereoAtoms))) {
        atom->setChiralTag(Atom::CHI_UNSPECIFIED);

        // If the atom has an explicit hydrogen and no charge, that H
        // was probably put there solely because of the chirality.
        // So we'll go ahead and remove it.
        // This was Issue 194
        if (atom->getNumExplicitHs() == 1 && atom->getFormalCharge() == 0 &&
            !atom->getIsAromatic()) {
          atom->setNumExplicitHs(0);
          atom->setNoImplicit(false);
          atom->calcExplicitValence(false);
          atom->calcImplicitValence(false);
        }
      }
    }
    for (ROMol::BondIterator bondIt = mol.beginBonds();
         bondIt != mol.endBonds(); ++bondIt) {
      // wedged bonds to atoms that have no stereochem
      // should be removed. (github issue 87)
      if (((*bondIt)->getBondDir() == Bond::BEGINWEDGE ||
           (*bondIt)->getBondDir() == Bond::BEGINDASH) &&
          (*bondIt)->getBeginAtom()->getChiralTag() == Atom::CHI_UNSPECIFIED &&
          (*bondIt)->getEndAtom()->getChiralTag() == Atom::CHI_UNSPECIFIED) {
        (*bondIt)->setBondDir(Bond::NONE);
      }
    }
  }
  mol.setProp(common_properties::_StereochemDone, 1, true);

#if 0
      std::cerr<<"---\n";
      mol.debugMol(std::cerr);
      std::cerr<<"<<<<<<<<<<<<<<<<\n";
#endif
}

// Find bonds than can be cis/trans in a molecule and mark them as
// Bond::STEREOANY.
void findPotentialStereoBonds(ROMol &mol, bool cleanIt) {
  // FIX: The earlier thought was to provide an optional argument to ignore or
  // consider
  //  double bonds in a ring. But I am removing this optional argument and
  //  ignoring ring bonds
  //  completely for now. This is because finding a potential stereo bond in a
  //  ring involves
  //  more than just checking the CIPranks for the neighbors - SP 05/04/04

  // make this function callable multiple times
  if ((mol.hasProp(common_properties::_BondsPotentialStereo)) && (!cleanIt)) {
    return;
  } else {
    UINT_VECT ranks;
    ranks.resize(mol.getNumAtoms());
    bool cipDone = false;

    ROMol::BondIterator bondIt;
    for (bondIt = mol.beginBonds(); bondIt != mol.endBonds(); ++bondIt) {
      if ((*bondIt)->getBondType() == Bond::DOUBLE &&
          !(mol.getRingInfo()->numBondRings((*bondIt)->getIdx()))) {
        // we are ignoring ring bonds here - read the FIX above
        Bond *dblBond = *bondIt;
        // if the bond is flagged as EITHERDOUBLE, we ignore it:
        if (dblBond->getBondDir() == Bond::EITHERDOUBLE ||
            dblBond->getStereo() == Bond::STEREOANY) {
          break;
        }
        // proceed only if we either want to clean the stereocode on this bond
        // or if none is set on it yet
        if (cleanIt || dblBond->getStereo() == Bond::STEREONONE) {
          dblBond->setStereo(Bond::STEREONONE);
          const Atom *begAtom = dblBond->getBeginAtom(),
                     *endAtom = dblBond->getEndAtom();
          // we're only going to handle 2 or three coordinate atoms:
          if ((begAtom->getDegree() == 2 || begAtom->getDegree() == 3) &&
              (endAtom->getDegree() == 2 || endAtom->getDegree() == 3)) {
            // ------------------
            // get the CIP ranking of each atom if we need it:
            if (!cipDone) {
              if (!begAtom->hasProp(common_properties::_CIPRank)) {
                Chirality::assignAtomCIPRanks(mol, ranks);
              } else {
                // no need to recompute if we don't need to recompute. :-)
                for (unsigned int ai = 0; ai < mol.getNumAtoms(); ++ai) {
                  ranks[ai] = mol.getAtomWithIdx(ai)->getProp<unsigned int>(
                      common_properties::_CIPRank);
                }
              }
              cipDone = true;
            }
            // find the neighbors for the begin atom and the endAtom
            UINT_VECT begAtomNeighbors, endAtomNeighbors;
            Chirality::findAtomNeighborsHelper(mol, begAtom, dblBond,
                                               begAtomNeighbors);
            Chirality::findAtomNeighborsHelper(mol, endAtom, dblBond,
                                               endAtomNeighbors);
            if (begAtomNeighbors.size() > 0 && endAtomNeighbors.size() > 0) {
              if ((begAtomNeighbors.size() == 2) &&
                  (endAtomNeighbors.size() == 2)) {
// if both of the atoms have 2 neighbors (other than the one
// connected
// by the double bond) and ....
#if 0
                std::cerr << "Bond: " << dblBond->getIdx() << " "
                          << begAtom->getIdx() << "=" << endAtom->getIdx()
                          << std::endl;
                std::cerr << "   " << begAtomNeighbors[0] << "="
                          << ranks[begAtomNeighbors[0]] << ":";
                std::cerr << "   " << begAtomNeighbors[1] << "="
                          << ranks[begAtomNeighbors[1]] << std::endl;
                std::cerr << "   " << endAtomNeighbors[0] << "="
                          << ranks[endAtomNeighbors[0]] << ":";
                std::cerr << "   " << endAtomNeighbors[1] << "="
                          << ranks[endAtomNeighbors[1]] << std::endl;
#endif
                if ((ranks[begAtomNeighbors[0]] !=
                     ranks[begAtomNeighbors[1]]) &&
                    (ranks[endAtomNeighbors[0]] !=
                     ranks[endAtomNeighbors[1]])) {
                  // the neighbors ranks are different at both the ends,
                  // this bond can be part of a cis/trans system
                  if (ranks[begAtomNeighbors[0]] > ranks[begAtomNeighbors[1]]) {
                    dblBond->getStereoAtoms().push_back(begAtomNeighbors[0]);
                  } else {
                    dblBond->getStereoAtoms().push_back(begAtomNeighbors[1]);
                  }
                  if (ranks[endAtomNeighbors[0]] > ranks[endAtomNeighbors[1]]) {
                    dblBond->getStereoAtoms().push_back(endAtomNeighbors[0]);
                  } else {
                    dblBond->getStereoAtoms().push_back(endAtomNeighbors[1]);
                  }
                }
              } else if (begAtomNeighbors.size() == 2) {
                // if the begAtom has two neighbors and ....
                if (ranks[begAtomNeighbors[0]] != ranks[begAtomNeighbors[1]]) {
                  // their ranks are different
                  if (ranks[begAtomNeighbors[0]] > ranks[begAtomNeighbors[1]]) {
                    dblBond->getStereoAtoms().push_back(begAtomNeighbors[0]);
                  } else {
                    dblBond->getStereoAtoms().push_back(begAtomNeighbors[1]);
                  }
                  dblBond->getStereoAtoms().push_back(endAtomNeighbors[0]);
                }
              } else if (endAtomNeighbors.size() == 2) {
                // if the endAtom has two neighbors and ...
                if (ranks[endAtomNeighbors[0]] != ranks[endAtomNeighbors[1]]) {
                  // their ranks are different
                  dblBond->getStereoAtoms().push_back(begAtomNeighbors[0]);
                  if (ranks[endAtomNeighbors[0]] > ranks[endAtomNeighbors[1]]) {
                    dblBond->getStereoAtoms().push_back(endAtomNeighbors[0]);
                  } else {
                    dblBond->getStereoAtoms().push_back(endAtomNeighbors[1]);
                  }
                }
              } else {
                // end and beg atoms has only one neighbor each, it doesn't
                // matter what the ranks are:
                dblBond->getStereoAtoms().push_back(begAtomNeighbors[0]);
                dblBond->getStereoAtoms().push_back(endAtomNeighbors[0]);
              }  // end of different number of neighbors on beg and end atoms

              // mark this double bond as a potential stereo bond
              if (!dblBond->getStereoAtoms().empty()) {
                dblBond->setStereo(Bond::STEREOANY);
              }
            }  // end of check that beg and end atoms have at least 1
               // neighbor:
          }    // end of 2 and 3 coordinated atoms only
        }      // end of we want it or CIP code is not set
      }        // end of double bond
    }          // end of for loop over all bonds
    mol.setProp(common_properties::_BondsPotentialStereo, 1, true);
  }
}

// removes chirality markers from sp and sp2 hybridized centers:
void cleanupChirality(RWMol &mol) {
  for (ROMol::AtomIterator atomIt = mol.beginAtoms(); atomIt != mol.endAtoms();
       ++atomIt) {
    if ((*atomIt)->getChiralTag() != Atom::CHI_UNSPECIFIED &&
        (*atomIt)->getHybridization() < Atom::SP3) {
      (*atomIt)->setChiralTag(Atom::CHI_UNSPECIFIED);
    }
  }
}

void assignChiralTypesFrom3D(ROMol &mol, int confId, bool replaceExistingTags) {
  const double ZERO_VOLUME_TOL = 0.1;
  if (!mol.getNumConformers()) return;
  const Conformer &conf = mol.getConformer(confId);
  if (!conf.is3D()) return;

  // if the molecule already has stereochemistry
  // perceived, remove the flags that indicate
  // this... what we're about to do will require
  // that we go again.
  if (mol.hasProp(common_properties::_StereochemDone)) {
    mol.clearProp(common_properties::_StereochemDone);
  }

  for (ROMol::AtomIterator atomIt = mol.beginAtoms(); atomIt != mol.endAtoms();
       ++atomIt) {
    Atom *atom = *atomIt;
    // if we aren't replacing existing tags and the atom is already tagged,
    // punt:
    if (!replaceExistingTags && atom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
      continue;
    }
    atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // additional reasons to skip the atom:
    if (atom->getDegree() < 3 || atom->getTotalDegree() > 4) {
      // not enough explicit neighbors or too many total neighbors
      continue;
    } else {
      int anum = atom->getAtomicNum();
      if (anum != 16 && anum != 34 &&  // S or Se are special
                                       // (just using the InChI list for now)
          (atom->getTotalDegree() != 4 ||  // not enough total neighbors
           atom->getTotalNumHs(true) > 1)) {
        continue;
      }
    }
    const RDGeom::Point3D &p0 = conf.getAtomPos(atom->getIdx());
    ROMol::ADJ_ITER nbrIdx, endNbrs;
    boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
    const RDGeom::Point3D &p1 = conf.getAtomPos(*nbrIdx);
    ++nbrIdx;
    const RDGeom::Point3D &p2 = conf.getAtomPos(*nbrIdx);
    ++nbrIdx;
    const RDGeom::Point3D &p3 = conf.getAtomPos(*nbrIdx);

    RDGeom::Point3D v1 = p1 - p0;
    RDGeom::Point3D v2 = p2 - p0;
    RDGeom::Point3D v3 = p3 - p0;

    double chiralVol = v1.dotProduct(v2.crossProduct(v3));
    if (chiralVol < -ZERO_VOLUME_TOL) {
      atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
    } else if (chiralVol > ZERO_VOLUME_TOL) {
      atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
    } else {
      atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    }
  }
}

void removeStereochemistry(ROMol &mol) {
  if (mol.hasProp(common_properties::_StereochemDone)) {
    mol.clearProp(common_properties::_StereochemDone);
  }
  for (ROMol::AtomIterator atIt = mol.beginAtoms(); atIt != mol.endAtoms();
       ++atIt) {
    (*atIt)->setChiralTag(Atom::CHI_UNSPECIFIED);
    if ((*atIt)->hasProp(common_properties::_CIPCode)) {
      (*atIt)->clearProp(common_properties::_CIPCode);
    }
    if ((*atIt)->hasProp(common_properties::_CIPRank)) {
      (*atIt)->clearProp(common_properties::_CIPRank);
    }
  }
  for (ROMol::BondIterator bondIt = mol.beginBonds(); bondIt != mol.endBonds();
       ++bondIt) {
    if ((*bondIt)->getBondType() == Bond::DOUBLE) {
      (*bondIt)->setStereo(Bond::STEREONONE);
      (*bondIt)->getStereoAtoms().clear();
    } else if ((*bondIt)->getBondType() == Bond::SINGLE) {
      (*bondIt)->setBondDir(Bond::NONE);
    }
  }
}
}  // end of namespace MolOps
}  // end of namespace RDKit
