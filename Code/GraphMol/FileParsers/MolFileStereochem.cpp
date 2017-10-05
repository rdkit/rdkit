//
//  Copyright (C) 2004-2017 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
#include <list>
#include <RDGeneral/RDLog.h>
#include "MolFileStereochem.h"
#include <Geometry/point.h>
#include <boost/dynamic_bitset.hpp>
#include <algorithm>
#include "MolFileStereochem.h"
#include <RDGeneral/Ranking.h>

namespace RDKit {
typedef std::list<double> DOUBLE_LIST;

// ----------------------------------- -----------------------------------
// This algorithm is identical to that used in the CombiCode Mol file
//  parser (also developed by RD).
//
//
// SUMMARY:
//   Derive a chiral code for an atom that has a wedged (or dashed) bond
//   drawn to it.
//
// RETURNS:
//   The chiral type
//
// CAVEATS:
//   This is careful to ensure that the central atom has 4 neighbors and
//   only single bonds to it, but that's about it.
//
// NOTE: this isn't careful at all about checking to make sure that
// things actually *should* be chiral. e.g. if the file has a
// 3-coordinate N with a wedged bond, it will make some erroneous
// assumptions about the chirality.
//
// ----------------------------------- -----------------------------------

Atom::ChiralType FindAtomStereochemistry(const RWMol &mol, const Bond *bond,
                                         const Conformer *conf) {
  PRECONDITION(bond, "no bond");
  PRECONDITION(conf, "no conformer");
  Bond::BondDir bondDir = bond->getBondDir();
  PRECONDITION(bondDir == Bond::BEGINWEDGE || bondDir == Bond::BEGINDASH,
               "bad bond direction");

  // NOTE that according to the CT file spec, wedging assigns chirality
  // to the atom at the point of the wedge, (atom 1 in the bond).
  const Atom *atom = bond->getBeginAtom();
  PRECONDITION(atom, "no atom");

  // we can't do anything with atoms that have more than 4 neighbors:
  if (atom->getDegree() > 4) {
    return Atom::CHI_UNSPECIFIED;
  }
  const Atom *bondAtom = bond->getEndAtom();

  Atom::ChiralType res = Atom::CHI_UNSPECIFIED;

  INT_LIST neighborBondIndices;
  RDGeom::Point3D centerLoc, tmpPt;
  centerLoc = conf->getAtomPos(atom->getIdx());
  tmpPt = conf->getAtomPos(bondAtom->getIdx());
  centerLoc.z = 0.0;
  tmpPt.z = 0.0;

  RDGeom::Point3D refVect = centerLoc.directionVector(tmpPt);

  //----------------------------------------------------------
  //
  //  start by ensuring that all the bonds to neighboring atoms
  //  are single bonds and collecting a list of neighbor indices:
  //
  //----------------------------------------------------------
  bool hSeen = false;

  neighborBondIndices.push_back(bond->getIdx());
  if (bondAtom->getAtomicNum() == 1 && bondAtom->getIsotope() == 0)
    hSeen = true;

  bool allSingle = true;
  ROMol::OEDGE_ITER beg, end;
  boost::tie(beg, end) = mol.getAtomBonds(atom);
  while (beg != end) {
    Bond *nbrBond = mol[*beg].get();
    if (nbrBond->getBondType() != Bond::SINGLE) {
      allSingle = false;
      // break;
    }
    if (nbrBond != bond) {
      if ((nbrBond->getOtherAtom(atom)->getAtomicNum() == 1 &&
           nbrBond->getOtherAtom(atom)->getIsotope() == 0))
        hSeen = true;
      neighborBondIndices.push_back(nbrBond->getIdx());
    }
    ++beg;
  }
  size_t nNbrs = neighborBondIndices.size();

  //----------------------------------------------------------
  //
  //  Return now if there aren't at least 3 non-H bonds to the atom.
  //  (we can implicitly add a single H to 3 coordinate atoms, but
  //  we're horked otherwise).
  //
  //----------------------------------------------------------
  if (nNbrs < 3 || (hSeen && nNbrs < 4)) {
    return Atom::CHI_UNSPECIFIED;
  }

  //----------------------------------------------------------
  //
  //  Continue if there are all single bonds or if we're considering
  //  4-coordinate P or S
  //
  //----------------------------------------------------------
  if (allSingle || atom->getAtomicNum() == 15 || atom->getAtomicNum() == 16) {
    //------------------------------------------------------------
    //
    //  Here we need to figure out the rotation direction between
    //  the neighbor bonds and the wedged bond:
    //
    //------------------------------------------------------------
    bool isCCW = true;
    double angle0, angle1, angle2;
    const Bond *bond1, *bond2, *bond3;
    RDGeom::Point3D atomVect0, atomVect1, atomVect2;
    INT_LIST::const_iterator bondIter = neighborBondIndices.begin();
    ++bondIter;
    bond1 = mol.getBondWithIdx(*bondIter);
    int oaid = bond1->getOtherAtom(atom)->getIdx();
    tmpPt = conf->getAtomPos(oaid);
    tmpPt.z = 0;
    atomVect0 = centerLoc.directionVector(tmpPt);
    angle0 = refVect.signedAngleTo(atomVect0);
    if (angle0 < 0) angle0 += 2. * M_PI;

    ++bondIter;
    bond2 = mol.getBondWithIdx(*bondIter);
    oaid = bond2->getOtherAtom(atom)->getIdx();
    tmpPt = conf->getAtomPos(oaid);
    tmpPt.z = 0;
    atomVect1 = centerLoc.directionVector(tmpPt);
    angle1 = refVect.signedAngleTo(atomVect1);
    if (angle1 < 0) angle1 += 2. * M_PI;

    // We proceed differently for 3 and 4 coordinate atoms:
    double firstAngle, secondAngle;
    if (nNbrs == 4) {
      bool flipIt = false;
      // grab the angle to the last neighbor:
      ++bondIter;
      bond3 = mol.getBondWithIdx(*bondIter);
      oaid = bond3->getOtherAtom(atom)->getIdx();
      tmpPt = conf->getAtomPos(oaid);
      tmpPt.z = 0;
      atomVect2 = centerLoc.directionVector(tmpPt);
      angle2 = refVect.signedAngleTo(atomVect2);
      if (angle2 < 0) angle2 += 2. * M_PI;

      // find the lowest and second-lowest angle and keep track of
      // whether or not we have to do a non-cyclic permutation to
      // get there:
      if (angle0 < angle1) {
        if (angle1 < angle2) {
          // order is angle0 -> angle1 -> angle2
          firstAngle = angle0;
          secondAngle = angle1;
        } else if (angle0 < angle2) {
          // order is angle0 -> angle2 -> angle1
          firstAngle = angle0;
          secondAngle = angle2;
          flipIt = true;
        } else {
          // order is angle2 -> angle0 -> angle1
          firstAngle = angle2;
          secondAngle = angle0;
        }
      } else if (angle0 < angle2) {
        // order is angle1 -> angle0 -> angle2
        firstAngle = angle1;
        secondAngle = angle0;
        flipIt = true;
      } else {
        if (angle1 < angle2) {
          // order is angle1 -> angle2 -> angle0
          firstAngle = angle1;
          secondAngle = angle2;
        } else {
          // order is angle2 -> angle1 -> angle0
          firstAngle = angle2;
          secondAngle = angle1;
          flipIt = true;
        }
      }
      if (flipIt) {
        isCCW = !isCCW;
      }
    } else {
      // it's three coordinate.  Things are a bit different here
      // because we have to at least kind of figure out where the
      // hydrogen might be.

      // before getting started with that, use some of the inchi rules
      // for contradictory stereochemistry
      // (Table 10 in the InChi v1 technical manual)

      angle2 = atomVect0.signedAngleTo(atomVect1);
      if (angle2 < 0) angle2 += 2. * M_PI;

      //  this one is never allowed:
      //     0   2
      //      \ /
      //       C
      //       *
      //       1
      if (angle0 < (M_PI - 1e-3) && angle1 < (M_PI - 1e-3) &&
          angle2 < (M_PI - 1e-3)) {
        if ((bond1->getBondDir() != Bond::NONE &&
             bond1->getBeginAtomIdx() == bond->getBeginAtomIdx() &&
             (bond1->getBondDir() != bond->getBondDir() ||
              (bond2->getBondDir() != Bond::NONE &&
               bond2->getBeginAtomIdx() == bond->getBeginAtomIdx() &&
               bond2->getBondDir() != bond1->getBondDir()))) ||
            (bond2->getBondDir() != Bond::NONE &&
             bond2->getBeginAtomIdx() == bond->getBeginAtomIdx() &&
             bond2->getBondDir() != bond->getBondDir())) {
          BOOST_LOG(rdWarningLog)
              << "Warning: conflicting stereochemistry at atom "
              << bond->getBeginAtomIdx() << " ignored."
              << std::endl;  // by rule 1." << std::endl;
          return Atom::CHI_UNSPECIFIED;
        }
      }
      if (bond1->getBondDir() != Bond::NONE &&
          bond1->getBeginAtomIdx() == bond->getBeginAtomIdx()) {
        if (!(bond2->getBondDir() != Bond::NONE &&
              bond2->getBeginAtomIdx() == bond->getBeginAtomIdx())) {
          BOOST_LOG(rdWarningLog)
              << "Warning: conflicting stereochemistry at atom "
              << bond->getBeginAtomIdx() << " ignored."
              << std::endl;  // by rule 2a." << std::endl;
        }
        if (bond1->getBondDir() != bond->getBondDir()) {
          // bond1 has a spec and does not match the bond0 spec.
          // the only cases this is allowed are:
          //      1        0 1 2
          //      *         \*/
          //  0 - C - 2      C
          //    and
          //      1        2 1 0
          //      *         \*/
          //  2 - C - 0      C
          //
          if ((angle0 > M_PI && angle0 < angle1) ||
              (angle0 < M_PI && angle0 > angle1)) {
            BOOST_LOG(rdWarningLog)
                << "Warning: conflicting stereochemistry at atom "
                << bond->getBeginAtomIdx() << " ignored."
                << std::endl;  // by rule 2b." << std::endl;
            return Atom::CHI_UNSPECIFIED;
          }
        } else {
          // bond1 matches, what about bond2 ?
          if (bond2->getBondDir() != bond->getBondDir()) {
            // the only cases this is allowed are:
            //      2        0 2 1
            //      *         \*/
            //  0 - C - 1      C
            //    and
            //      2        1 2 0
            //      *         \*/
            //  1 - C - 0      C
            //
            if ((angle1 > M_PI && angle1 < angle0) ||
                (angle1 < M_PI && angle1 > angle0)) {
              BOOST_LOG(rdWarningLog)
                  << "Warning: conflicting stereochemistry at atom "
                  << bond->getBeginAtomIdx() << " ignored."
                  << std::endl;  // by rule 2c." << std::endl;
              return Atom::CHI_UNSPECIFIED;
            }
          }
        }
      } else if (bond2->getBondDir() != Bond::NONE &&
                 bond2->getBeginAtomIdx() == bond->getBeginAtomIdx() &&
                 bond2->getBondDir() != bond->getBondDir()) {
        // bond2 has a spec and does not match the bond0 spec, but bond1
        // is not set: this is never allowed.
        BOOST_LOG(rdWarningLog)
            << "Warning: conflicting stereochemistry at atom "
            << bond->getBeginAtomIdx() << " ignored."
            << std::endl;  // by rule 3." << std::endl;
        return Atom::CHI_UNSPECIFIED;
      }

      if (angle0 < angle1) {
        firstAngle = angle0;
        secondAngle = angle1;
        isCCW = true;
      } else {
        firstAngle = angle1;
        secondAngle = angle0;
        isCCW = false;
      }
      if (secondAngle - firstAngle >= (M_PI - 1e-4)) {
        // it's a situation like one of these:
        //
        //      0        1 0 2
        //      *         \*/
        //  1 - C - 2      C
        //
        // In each of these cases, the implicit H is between atoms 1
        // and 2, so we need to flip the rotation direction (go
        // around the back).
        isCCW = !isCCW;
      }
    }
    // reverse the rotation direction if the reference is wedged down:
    if (bondDir == Bond::BEGINDASH) {
      isCCW = !isCCW;
    }

    // ----------------
    //
    // We now have the rotation direction using mol-file order.
    // We need to convert that into the appropriate label for the
    // central atom
    //
    // ----------------
    int nSwaps = atom->getPerturbationOrder(neighborBondIndices);
    if (nSwaps % 2) isCCW = !isCCW;
    if (isCCW)
      res = Atom::CHI_TETRAHEDRAL_CCW;
    else
      res = Atom::CHI_TETRAHEDRAL_CW;
  }

  return res;
}

void WedgeMolBonds(ROMol &mol, const Conformer *conf) {
  PRECONDITION(conf, "no conformer");
  INT_MAP_INT wedgeBonds = pickBondsToWedge(mol);
  for (ROMol::BondIterator bondIt = mol.beginBonds(); bondIt != mol.endBonds();
       ++bondIt) {
    Bond *bond = *bondIt;
    if (bond->getBondType() == Bond::SINGLE) {
      Bond::BondDir dir = DetermineBondWedgeState(bond, wedgeBonds, conf);
      if (dir == Bond::BEGINWEDGE || dir == Bond::BEGINDASH) {
        bond->setBondDir(dir);
      }
    }
  }
}

INT_MAP_INT pickBondsToWedge(const ROMol &mol) {
  // we need ring information; make sure findSSSR has been called before
  // if not call now
  if (!mol.getRingInfo()->isInitialized()) {
    MolOps::findSSSR(mol);
  }

  static int noNbrs = 100;
  INT_VECT nChiralNbrs(mol.getNumAtoms(), noNbrs);

  // start by looking for bonds that are already wedged
  for (ROMol::ConstBondIterator cbi = mol.beginBonds(); cbi != mol.endBonds();
       ++cbi) {
    const Bond *bond = *cbi;
    if (bond->getBondDir() == Bond::BEGINWEDGE ||
        bond->getBondDir() == Bond::BEGINDASH ||
        bond->getBondDir() == Bond::UNKNOWN) {
      if (bond->getBeginAtom()->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
          bond->getBeginAtom()->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW)
        nChiralNbrs[bond->getBeginAtomIdx()] = noNbrs + 1;
      else if (bond->getEndAtom()->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
               bond->getEndAtom()->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW)
        nChiralNbrs[bond->getEndAtomIdx()] = noNbrs + 1;
    }
  }

  // now rank atoms by the number of chiral neighbors or Hs they have:
  bool chiNbrs = false;
  for (ROMol::ConstAtomIterator cai = mol.beginAtoms(); cai != mol.endAtoms();
       ++cai) {
    const Atom *at = *cai;
    if (nChiralNbrs[at->getIdx()] > noNbrs) {
      // std::cerr << " SKIPPING1: " << at->getIdx() << std::endl;
      continue;
    }
    Atom::ChiralType type = at->getChiralTag();
    if (type != Atom::CHI_TETRAHEDRAL_CW && type != Atom::CHI_TETRAHEDRAL_CCW)
      continue;
    nChiralNbrs[at->getIdx()] = 0;
    chiNbrs = true;
    ROMol::ADJ_ITER nbrIdx, endNbrs;
    boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(at);
    while (nbrIdx != endNbrs) {
      const ATOM_SPTR nat = mol[*nbrIdx];
      ++nbrIdx;
      if (nat->getAtomicNum() == 1) {
        // special case: it's an H... we weight these especially high:
        nChiralNbrs[at->getIdx()] -= 10;
        continue;
      }
      type = nat->getChiralTag();
      if (type != Atom::CHI_TETRAHEDRAL_CW && type != Atom::CHI_TETRAHEDRAL_CCW)
        continue;
      nChiralNbrs[at->getIdx()] -= 1;
    }
  }
  std::vector<unsigned int> indices(mol.getNumAtoms());
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) indices[i] = i;
  if (chiNbrs) {
    std::sort(indices.begin(), indices.end(),
              Rankers::argless<INT_VECT>(nChiralNbrs));
  }
#if 0
  std::cerr << "  nbrs: ";
  std::copy(nChiralNbrs.begin(), nChiralNbrs.end(),
            std::ostream_iterator<int>(std::cerr, " "));
  std::cerr << std::endl;
  std::cerr << "  order: ";
  std::copy(indices.begin(), indices.end(),
            std::ostream_iterator<int>(std::cerr, " "));
  std::cerr << std::endl;
#endif
  // picks a bond for each atom that we will wedge when we write the mol file
  // here is what we are going to do
  // - at each chiral center look for a bond that is begins at the atom and
  //   is not yet picked to be wedged for a different chiral center, preferring
  //   bonds to Hs
  // - if we do not find a bond that begins at the chiral center - we will take
  //   the first bond that is not yet picked by any other chiral centers
  // we use the orders calculated above to determine which order to do the
  // wedging
  INT_MAP_INT res;
  BOOST_FOREACH (unsigned int idx, indices) {
    if (nChiralNbrs[idx] > noNbrs) {
      // std::cerr << " SKIPPING2: " << idx << std::endl;
      continue;  // already have a wedged bond here
    }
    const Atom *atom = mol.getAtomWithIdx(idx);
    Atom::ChiralType type = atom->getChiralTag();
    // the indices are ordered such that all chiral atoms come first. If
    // this has no chiral flag, we can stop the whole loop:
    if (type != Atom::CHI_TETRAHEDRAL_CW && type != Atom::CHI_TETRAHEDRAL_CCW)
      break;
    RDKit::ROMol::OBOND_ITER_PAIR atomBonds = mol.getAtomBonds(atom);
    std::vector<std::pair<int, int>> nbrScores;
    while (atomBonds.first != atomBonds.second) {
      const Bond *bond = mol[*atomBonds.first].get();
      ++atomBonds.first;

      // can only wedge single bonds:
      if (bond->getBondType() != Bond::SINGLE) continue;

      int bid = bond->getIdx();
      if (res.find(bid) == res.end()) {
        // very strong preference for Hs:
        if (bond->getOtherAtom(atom)->getAtomicNum() == 1) {
          nbrScores.push_back(std::make_pair(
              -1000000, bid));  // lower than anything else can be
          continue;
        }
        // prefer lower atomic numbers with lower degrees and no specified
        // chirality:
        const Atom *oatom = bond->getOtherAtom(atom);
        int nbrScore = oatom->getAtomicNum() + 10 * oatom->getDegree() +
                       100 * ((oatom->getChiralTag() != Atom::CHI_UNSPECIFIED));
        // prefer neighbors that are nonchiral or have as few chiral neighbors
        // as possible:
        int oIdx = oatom->getIdx();
        if (nChiralNbrs[oIdx] < noNbrs) {
          // the counts are negative, so we have to subtract them off
          nbrScore -= 10000 * nChiralNbrs[oIdx];
        }
        // prefer bonds to non-ring atoms:
        nbrScore += 1000 * mol.getRingInfo()->numAtomRings(oIdx);
        // prefer non-ring bonds;
        nbrScore += 1000 * mol.getRingInfo()->numBondRings(bid);
        // std::cerr << "    nrbScore: " << idx << " - " << oIdx << " : "
        //           << nbrScore << " nChiralNbrs: " << nChiralNbrs[oIdx]
        //           << std::endl;
        nbrScores.push_back(std::make_pair(nbrScore, bid));
      }
    }
    // There's still one situation where this whole thing can fail: an unlucky
    // situation where all neighbors of all neighbors of an atom are chiral and
    // that atom ends up being the last one picked for stereochem assignment.
    //
    // We'll catch that as an error here and hope that it's as unlikely to occur
    // as it seems like it is. (I'm going into this knowing that it's bound to
    // happen; I'll kick myself and do the hard solution at that point.)
    CHECK_INVARIANT(nbrScores.size(),
                    "no eligible neighbors for chiral center");
    std::sort(nbrScores.begin(), nbrScores.end(),
              Rankers::pairLess<int, int>());
    res[nbrScores[0].second] = idx;
  }
  return res;
}

//
// Determine bond wedge state
///
Bond::BondDir DetermineBondWedgeState(const Bond *bond,
                                      const INT_MAP_INT &wedgeBonds,
                                      const Conformer *conf) {
  PRECONDITION(bond, "no bond");
  PRECONDITION(bond->getBondType() == Bond::SINGLE,
               "bad bond order for wedging");
  const ROMol *mol = &(bond->getOwningMol());
  PRECONDITION(mol, "no mol");

  Bond::BondDir res = bond->getBondDir();
  if (!conf) {
    return res;
  }

  int bid = bond->getIdx();
  auto wbi = wedgeBonds.find(bid);
  if (wbi == wedgeBonds.end()) {
    return res;
  }

  unsigned int waid = wbi->second;

  Atom *atom, *bondAtom;  // = bond->getBeginAtom();
  if (bond->getBeginAtom()->getIdx() == waid) {
    atom = bond->getBeginAtom();
    bondAtom = bond->getEndAtom();
  } else {
    atom = bond->getEndAtom();
    bondAtom = bond->getBeginAtom();
  }

  Atom::ChiralType chiralType = atom->getChiralTag();
  CHECK_INVARIANT(chiralType == Atom::CHI_TETRAHEDRAL_CW ||
                      chiralType == Atom::CHI_TETRAHEDRAL_CCW,
                  "");

  // if we got this far, we really need to think about it:
  INT_LIST neighborBondIndices;
  DOUBLE_LIST neighborBondAngles;
  RDGeom::Point3D centerLoc, tmpPt;
  centerLoc = conf->getAtomPos(atom->getIdx());
  tmpPt = conf->getAtomPos(bondAtom->getIdx());
  centerLoc.z = 0.0;
  tmpPt.z = 0.0;
  RDGeom::Point3D refVect = centerLoc.directionVector(tmpPt);

  neighborBondIndices.push_back(bond->getIdx());
  neighborBondAngles.push_back(0.0);

  ROMol::OEDGE_ITER beg, end;
  boost::tie(beg, end) = mol->getAtomBonds(atom);
  while (beg != end) {
    Bond *nbrBond = (*mol)[*beg].get();
    Atom *otherAtom = nbrBond->getOtherAtom(atom);
    if (nbrBond != bond) {
      tmpPt = conf->getAtomPos(otherAtom->getIdx());
      tmpPt.z = 0.0;
      RDGeom::Point3D tmpVect = centerLoc.directionVector(tmpPt);
      double angle = refVect.signedAngleTo(tmpVect);
      if (angle < 0.0) angle += 2. * M_PI;
      auto nbrIt = neighborBondIndices.begin();
      auto angleIt = neighborBondAngles.begin();
      // find the location of this neighbor in our angle-sorted list
      // of neighbors:
      while (angleIt != neighborBondAngles.end() && angle > (*angleIt)) {
        ++angleIt;
        ++nbrIt;
      }
      neighborBondAngles.insert(angleIt, angle);
      neighborBondIndices.insert(nbrIt, nbrBond->getIdx());
    }
    ++beg;
  }

  // at this point, neighborBondIndices contains a list of bond
  // indices from the central atom.  They are arranged starting
  // at the reference bond in CCW order (based on the current
  // depiction).
  int nSwaps = atom->getPerturbationOrder(neighborBondIndices);

  // in the case of three-coordinated atoms we may have to worry about
  // the location of the implicit hydrogen - Issue 209
  // Check if we have one of these situation
  //
  //      0        1 0 2
  //      *         \*/
  //  1 - C - 2      C
  //
  // here the hydrogen will be between 1 and 2 and we need to add an additional
  // swap
  if (neighborBondAngles.size() == 3) {
    // three coordinated
    auto angleIt = neighborBondAngles.begin();
    ++angleIt;  // the first is the 0 (or reference bond - we will ignoire that
    double angle1 = (*angleIt);
    ++angleIt;
    double angle2 = (*angleIt);
    if (angle2 - angle1 >= (M_PI - 1e-4)) {
      // we have the above situation
      nSwaps++;
    }
  }

#ifdef VERBOSE_STEREOCHEM
  BOOST_LOG(rdDebugLog) << "--------- " << nSwaps << std::endl;
  std::copy(neighborBondIndices.begin(), neighborBondIndices.end(),
            std::ostream_iterator<int>(BOOST_LOG(rdDebugLog), " "));
  BOOST_LOG(rdDebugLog) << std::endl;
  std::copy(neighborBondAngles.begin(), neighborBondAngles.end(),
            std::ostream_iterator<double>(BOOST_LOG(rdDebugLog), " "));
  BOOST_LOG(rdDebugLog) << std::endl;
#endif
  if (chiralType == Atom::CHI_TETRAHEDRAL_CCW) {
    if (nSwaps % 2 == 1) {  // ^ reverse) {
      res = Bond::BEGINDASH;
    } else {
      res = Bond::BEGINWEDGE;
    }
  } else {
    if (nSwaps % 2 == 1) {  // ^ reverse) {
      res = Bond::BEGINWEDGE;
    } else {
      res = Bond::BEGINDASH;
    }
  }

  return res;
}

// handles stereochem markers set by the Mol file parser and
// converts them to the RD standard:
void DetectAtomStereoChemistry(RWMol &mol, const Conformer *conf) {
  PRECONDITION(conf, "no conformer");

  for (RWMol::BondIterator bondIt = mol.beginBonds(); bondIt != mol.endBonds();
       ++bondIt) {
    Bond *bond = *bondIt;
    if (bond->getBondDir() != Bond::UNKNOWN) {
      Bond::BondDir dir = bond->getBondDir();
      // the bond is marked as chiral:
      if (dir == Bond::BEGINWEDGE || dir == Bond::BEGINDASH) {
        Atom *atom = bond->getBeginAtom();
        if (atom->getImplicitValence() == -1) {
          atom->calcExplicitValence(false);
          atom->calcImplicitValence(false);
        }
        Atom::ChiralType code = FindAtomStereochemistry(mol, bond, conf);
        atom->setChiralTag(code);
        // within the RD representation, if a three-coordinate atom
        // is chiral and has an implicit H, that H needs to be made explicit:
        if (atom->getDegree() == 3 && !atom->getNumExplicitHs() &&
            atom->getNumImplicitHs() == 1) {
          atom->setNumExplicitHs(1);
          // recalculated number of implicit Hs:
          atom->updatePropertyCache();
        }
      }
    }
  }
}

void ClearSingleBondDirFlags(ROMol &mol) {
  for (RWMol::BondIterator bondIt = mol.beginBonds(); bondIt != mol.endBonds();
       ++bondIt) {
    if ((*bondIt)->getBondType() == Bond::SINGLE) {
      if ((*bondIt)->getBondDir() == Bond::UNKNOWN)
        (*bondIt)->setProp(common_properties::_UnknownStereo, 1);
      (*bondIt)->setBondDir(Bond::NONE);
    }
  }
}

void DetectBondStereoChemistry(ROMol &mol, const Conformer *conf) {
  PRECONDITION(conf, "no conformer");
  PRECONDITION(&(conf->getOwningMol()) == &mol,
               "conformer does not belong to molecule");
  MolOps::detectBondStereochemistry(mol, conf->getId());
}
}
