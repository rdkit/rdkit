//
//  Copyright (C) 2004-2018 Greg Landrum and Rational Discovery LLC
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
#include "Chirality.h"

// #define VERBOSE_CANON 1

namespace RDKit {

namespace {
bool shouldDetectDoubleBondStereo(const Bond *bond) {
  const RingInfo *ri = bond->getOwningMol().getRingInfo();
  return (!ri->numBondRings(bond->getIdx()) ||
          ri->minBondRingSize(bond->getIdx()) > 7);
}

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

Atom::ChiralType atomChiralTypeFromBondDir(const ROMol &mol, const Bond *bond,
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
  if (bondAtom->getAtomicNum() == 1 && bondAtom->getIsotope() == 0) {
    hSeen = true;
  }

  bool allSingle = true;
  for(auto *nbrBond : atom->bonds()) {
    if (nbrBond->getBondType() != Bond::SINGLE) {
      allSingle = false;
      // break;
    }
    if (nbrBond != bond) {
      if ((nbrBond->getOtherAtom(atom)->getAtomicNum() == 1 &&
           nbrBond->getOtherAtom(atom)->getIsotope() == 0)) {
        hSeen = true;
      }
      neighborBondIndices.push_back(nbrBond->getIdx());
    }
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
    INT_LIST::const_iterator bondIter = neighborBondIndices.begin();
    ++bondIter;
    auto bond1 = mol.getBondWithIdx(*bondIter);
    int oaid = bond1->getOtherAtom(atom)->getIdx();
    tmpPt = conf->getAtomPos(oaid);
    tmpPt.z = 0;
    auto atomVect0 = centerLoc.directionVector(tmpPt);
    auto angle01 = refVect.signedAngleTo(atomVect0);

    ++bondIter;
    auto bond2 = mol.getBondWithIdx(*bondIter);
    oaid = bond2->getOtherAtom(atom)->getIdx();
    tmpPt = conf->getAtomPos(oaid);
    tmpPt.z = 0;
    auto atomVect1 = centerLoc.directionVector(tmpPt);
    auto angle02 = refVect.signedAngleTo(atomVect1);

    // order everything so that looking from the top in a CCW direction we
    // have 0, 1, 2
    bool swappedIt = false;
    if (angle01 > angle02) {
      // std::cerr << " swap because " << angle01 << " " << angle02 << " "
      //           << bond1->getIdx() << "->" << bond2->getIdx() << std::endl;
      std::swap(angle01, angle02);
      std::swap(bond1, bond2);
      std::swap(atomVect0, atomVect1);
      swappedIt = true;
    }

    double firstAngle, secondAngle;
    // We proceed differently for 3 and 4 coordinate atoms:
    double angle12;
    if (nNbrs == 4) {
      bool flipIt = false;
      // grab the angle to the last neighbor:
      ++bondIter;
      auto bond3 = mol.getBondWithIdx(*bondIter);
      oaid = bond3->getOtherAtom(atom)->getIdx();
      tmpPt = conf->getAtomPos(oaid);
      tmpPt.z = 0;
      auto atomVect2 = centerLoc.directionVector(tmpPt);
      angle12 = refVect.signedAngleTo(atomVect2);

      // find the lowest and second-lowest angle and keep track of
      // whether or not we have to do a non-cyclic permutation to
      // get there:
      if (angle01 < angle02) {
        if (angle02 < angle12) {
          // order is angle01 -> angle02 -> angle12
          firstAngle = angle01;
          secondAngle = angle02;
        } else if (angle01 < angle12) {
          // order is angle01 -> angle12 -> angle02
          firstAngle = angle01;
          secondAngle = angle12;
          flipIt = true;
        } else {
          // order is angle12 -> angle01 -> angle02
          firstAngle = angle12;
          secondAngle = angle01;
        }
      } else if (angle01 < angle12) {
        // order is angle02 -> angle01 -> angle12
        firstAngle = angle02;
        secondAngle = angle01;
        flipIt = true;
      } else {
        if (angle02 < angle12) {
          // order is angle02 -> angle12 -> angle01
          firstAngle = angle02;
          secondAngle = angle12;
        } else {
          // order is angle12 -> angle02 -> angle01
          firstAngle = angle12;
          secondAngle = angle02;
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

      angle12 = atomVect0.signedAngleTo(atomVect1);
      double angle20 = atomVect1.signedAngleTo(refVect);

      // to simplify the code below, pick out the directions of the bonds
      // if and only if they start at our atom:
      auto dir0 = bondDir;
      auto dir1 = (bond1->getBeginAtomIdx() == bond->getBeginAtomIdx())
                      ? bond1->getBondDir()
                      : Bond::NONE;
      auto dir2 = (bond2->getBeginAtomIdx() == bond->getBeginAtomIdx())
                      ? bond2->getBondDir()
                      : Bond::NONE;

      // we know bond 0 has the direction set

      //  this one is never allowed with different directions:
      //     0   2
      //      \ /
      //       C
      //       *
      //       1
      if (angle01 < (M_PI - 1e-3) && angle12 < (M_PI - 1e-3) &&
          angle20 < (M_PI - 1e-3)) {
        if ((dir1 != Bond::NONE && dir1 != dir0) ||
            (dir2 != Bond::NONE && dir2 != dir0)) {
          BOOST_LOG(rdWarningLog)
              << "Warning: conflicting stereochemistry at atom "
              << bond->getBeginAtomIdx() << " ignored."
              << " by rule 1a." << std::endl;
          return Atom::CHI_UNSPECIFIED;
        }
      } else {
        // otherwise they cannot all be the same
        if (dir1 == dir0 && dir1 == dir2) {
          BOOST_LOG(rdWarningLog)
              << "Warning: conflicting stereochemistry at atom "
              << bond->getBeginAtomIdx() << " ignored."
              << " by rule 1b." << std::endl;
          return Atom::CHI_UNSPECIFIED;
        }

        // the remaining cases where stereo is allowed are:
        //      0        1 0 2
        //      *         \*/
        //  1 - C - 2      C      for these two bond1 and bond2 must match
        //    and
        //      1        2 1 0
        //      *         \*/
        //  2 - C - 0      C      for these two bond0 and bond2 must match
        //    and
        //      2        0 2 1
        //      *         \*/
        //  0 - C - 1      C      for these two bond0 and bond1 must match
        //

        if (dir1 != Bond::NONE && dir1 != dir0) {
          if ((angle01 >= M_PI) ||  // last two examples
              (angle02 > M_PI && dir2 != Bond::NONE &&
               dir1 != dir2) ||  // top two examples
              (angle02 <= M_PI && dir2 != Bond::NONE &&
               dir0 != dir2)) {  // middle two examples
            BOOST_LOG(rdWarningLog)
                << "Warning: conflicting stereochemistry at atom "
                << bond->getBeginAtomIdx() << " ignored."
                << " by rule 2a." << std::endl;
            return Atom::CHI_UNSPECIFIED;
          }
        } else if (dir0 == dir2) {
          // all bonds are the same and we already removed the "all
          // angles less than 180" case above
          BOOST_LOG(rdWarningLog)
              << "Warning: conflicting stereochemistry at atom "
              << bond->getBeginAtomIdx() << " ignored."
              << " by rule 2b." << std::endl;
          return Atom::CHI_UNSPECIFIED;
        }

        if (angle01 < angle02) {
          firstAngle = angle01;
          secondAngle = angle02;
          isCCW = true;
        } else {
          firstAngle = angle02;
          secondAngle = angle01;
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
    }
    // reverse the rotation direction if the reference is wedged down:
    if (bondDir == Bond::BEGINDASH) {
      isCCW = !isCCW;
    }
    if (swappedIt) {
      // we swapped the order of the bonds to simplify the analysis above:
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
    if (nSwaps % 2) {
      isCCW = !isCCW;
    }
    if (isCCW) {
      res = Atom::CHI_TETRAHEDRAL_CCW;
    } else {
      res = Atom::CHI_TETRAHEDRAL_CW;
    }
  }

  return res;
}  // namespace RDKit

Bond::BondDir getOppositeBondDir(Bond::BondDir dir) {
  PRECONDITION(dir == Bond::ENDDOWNRIGHT || dir == Bond::ENDUPRIGHT,
               "bad bond direction");
  switch (dir) {
    case Bond::ENDDOWNRIGHT:
      return Bond::ENDUPRIGHT;
    case Bond::ENDUPRIGHT:
      return Bond::ENDDOWNRIGHT;
    default:
      return Bond::NONE;
  }
}

void setBondDirRelativeToAtom(Bond *bond, Atom *atom, Bond::BondDir dir,
                              bool reverse, boost::dynamic_bitset<> &needsDir) {
  PRECONDITION(bond, "bad bond");
  PRECONDITION(atom, "bad atom");
  PRECONDITION(dir == Bond::ENDUPRIGHT || dir == Bond::ENDDOWNRIGHT, "bad dir");
  PRECONDITION(atom == bond->getBeginAtom() || atom == bond->getEndAtom(),
               "atom doesn't belong to bond");
  RDUNUSED_PARAM(needsDir);

  if (bond->getBeginAtom() != atom) {
    reverse = !reverse;
  }

  if (reverse) {
    dir = (dir == Bond::ENDUPRIGHT ? Bond::ENDDOWNRIGHT : Bond::ENDUPRIGHT);
  }
  // to ensure maximum compatibility, even when a bond has unknown stereo (set
  // explicitly and recorded in _UnknownStereo property), I will still let a
  // direction to be computed. You must check the _UnknownStereo property to
  // make sure whether this bond is explicitly set to have no direction info.
  // This makes sense because the direction info are all derived from
  // coordinates, the _UnknownStereo property is like extra metadata to be
  // used with the direction info.
  bond->setBondDir(dir);
}

bool isLinearArrangement(const RDGeom::Point3D &v1, const RDGeom::Point3D &v2,
                         double tol = 0.035) {  // tolerance of 2 degrees
  return fabs(v2.angleTo(v1) - M_PI) < tol;
}
void updateDoubleBondNeighbors(ROMol &mol, Bond *dblBond, const Conformer *conf,
                               boost::dynamic_bitset<> &needsDir,
                               std::vector<unsigned int> &singleBondCounts,
                               const VECT_INT_VECT &singleBondNbrs) {
  // we want to deal only with double bonds:
  PRECONDITION(dblBond, "bad bond");
  PRECONDITION(dblBond->getBondType() == Bond::DOUBLE, "not a double bond");
  if (!needsDir[dblBond->getIdx()]) {
    return;
  }
  needsDir.set(dblBond->getIdx(), 0);
#if 0
  std::cerr << "**********************\n";
  std::cerr << "**********************\n";
  std::cerr << "**********************\n";
  std::cerr << "UDBN: " << dblBond->getIdx() << " "
            << dblBond->getBeginAtomIdx() << "=" << dblBond->getEndAtomIdx()
            << "\n";
#endif

  std::vector<Bond *> followupBonds;

  Bond *bond1 = nullptr, *obond1 = nullptr;
  bool squiggleBondSeen = false;
  bool doubleBondSeen = false;
    for(auto *tBond : dblBond->getBeginAtom()->bonds()) {
    if (tBond == dblBond) {
      continue;
    }
    if (tBond->getBondType() == Bond::SINGLE ||
        tBond->getBondType() == Bond::AROMATIC) {
      // prefer bonds that already have their directionality set
      // or that are adjacent to more double bonds:
      if (!bond1) {
        bond1 = tBond;
      } else if (needsDir[tBond->getIdx()]) {
        if (singleBondCounts[tBond->getIdx()] >
            singleBondCounts[bond1->getIdx()]) {
          obond1 = bond1;
          bond1 = tBond;
        } else {
          obond1 = tBond;
        }
      } else {
        obond1 = bond1;
        bond1 = tBond;
      }
    } else if (tBond->getBondType() == Bond::DOUBLE) {
      doubleBondSeen = true;
    }
    int explicit_unknown_stereo;
    if ((tBond->getBondType() == Bond::SINGLE ||
         tBond->getBondType() == Bond::AROMATIC) &&
        (tBond->getBondDir() == Bond::UNKNOWN ||
         ((tBond->getPropIfPresent<int>(common_properties::_UnknownStereo,
                                        explicit_unknown_stereo) &&
           explicit_unknown_stereo)))) {
      squiggleBondSeen = true;
      break;
    }
  }
  // Don't do any direction setting if we've seen a squiggle bond, but do mark
  // the double bond as a crossed bond and return
  if (!bond1 || squiggleBondSeen || doubleBondSeen) {
    if (!doubleBondSeen) {
      // FIX: This is the fix for #2649, but it will need to be modified once we
      // decide to properly handle allenese
      dblBond->setBondDir(Bond::EITHERDOUBLE);
    }
    return;
  }
  Bond *bond2 = nullptr, *obond2 = nullptr;
  for(auto *tBond : dblBond->getEndAtom()->bonds()) {
    if (tBond == dblBond) {
      continue;
    }
    if (tBond->getBondType() == Bond::SINGLE ||
        tBond->getBondType() == Bond::AROMATIC) {
      if (!bond2) {
        bond2 = tBond;
      } else if (needsDir[tBond->getIdx()]) {
        if (singleBondCounts[tBond->getIdx()] >
            singleBondCounts[bond2->getIdx()]) {
          obond2 = bond2;
          bond2 = tBond;
        } else {
          obond2 = tBond;
        }
      } else {
        // we already had a bond2 and we don't need to set the direction
        // on the new one, so swap.
        obond2 = bond2;
        bond2 = tBond;
      }
    } else if (tBond->getBondType() == Bond::DOUBLE) {
      doubleBondSeen = true;
    }
    int explicit_unknown_stereo;
    if (tBond->getBondType() == Bond::SINGLE &&
        (tBond->getBondDir() == Bond::UNKNOWN ||
         ((tBond->getPropIfPresent<int>(common_properties::_UnknownStereo,
                                        explicit_unknown_stereo) &&
           explicit_unknown_stereo)))) {
      squiggleBondSeen = true;
      break;
    }
  }
  // Don't do any direction setting if we've seen a squiggle bond, but do mark
  // the double bond as a crossed bond and return
  if (!bond2 || squiggleBondSeen || doubleBondSeen) {
    if (!doubleBondSeen) {
      // FIX: This is the fix for #2649, but it will need to be modified once we
      // decide to properly handle allenese
      dblBond->setBondDir(Bond::EITHERDOUBLE);
    }
    return;
  }

  CHECK_INVARIANT(bond1 && bond2, "no bonds found");

  bool sameTorsionDir = false;
  if (conf) {
    RDGeom::Point3D beginP = conf->getAtomPos(dblBond->getBeginAtomIdx());
    RDGeom::Point3D endP = conf->getAtomPos(dblBond->getEndAtomIdx());
    RDGeom::Point3D bond1P =
        conf->getAtomPos(bond1->getOtherAtomIdx(dblBond->getBeginAtomIdx()));
    RDGeom::Point3D bond2P =
        conf->getAtomPos(bond2->getOtherAtomIdx(dblBond->getEndAtomIdx()));
    // check for a linear arrangement of atoms on either end:
    bool linear = false;
    RDGeom::Point3D p1;
    RDGeom::Point3D p2;
    p1 = bond1P - beginP;
    p2 = endP - beginP;
    if (isLinearArrangement(p1, p2)) {
      if (!obond1) {
        linear = true;
      } else {
        // one of the bonds was linear; what about the other one?
        Bond *tBond = bond1;
        bond1 = obond1;
        obond1 = tBond;
        bond1P = conf->getAtomPos(
            bond1->getOtherAtomIdx(dblBond->getBeginAtomIdx()));
        p1 = bond1P - beginP;
        if (isLinearArrangement(p1, p2)) {
          linear = true;
        }
      }
    }
    if (!linear) {
      p1 = bond2P - endP;
      p2 = beginP - endP;
      if (isLinearArrangement(p1, p2)) {
        if (!obond2) {
          linear = true;
        } else {
          Bond *tBond = bond2;
          bond2 = obond2;
          obond2 = tBond;
          bond2P = conf->getAtomPos(
              bond2->getOtherAtomIdx(dblBond->getEndAtomIdx()));
          p1 = bond2P - beginP;
          if (isLinearArrangement(p1, p2)) {
            linear = true;
          }
        }
      }
    }
    if (linear) {
      dblBond->setBondDir(Bond::EITHERDOUBLE);
      return;
    }

    double ang = RDGeom::computeDihedralAngle(bond1P, beginP, endP, bond2P);
    if (ang < M_PI / 2) {
      sameTorsionDir = false;
    } else {
      sameTorsionDir = true;
    }
    // std::cerr << "   angle: " << ang << " sameTorsionDir: " << sameTorsionDir
    // << "\n";
  } else {
    if (dblBond->getStereo() == Bond::STEREOCIS ||
        dblBond->getStereo() == Bond::STEREOZ) {
      sameTorsionDir = false;
    } else if (dblBond->getStereo() == Bond::STEREOTRANS ||
               dblBond->getStereo() == Bond::STEREOE) {
      sameTorsionDir = true;
    } else {
      return;
    }
    // if bond1 or bond2 are not to the stereo-controlling atoms, flip
    // our expections of the torsion dir
    int bond1AtomIdx = bond1->getOtherAtomIdx(dblBond->getBeginAtomIdx());
    if (bond1AtomIdx != dblBond->getStereoAtoms()[0] &&
        bond1AtomIdx != dblBond->getStereoAtoms()[1]) {
      sameTorsionDir = !sameTorsionDir;
    }
    int bond2AtomIdx = bond2->getOtherAtomIdx(dblBond->getEndAtomIdx());
    if (bond2AtomIdx != dblBond->getStereoAtoms()[0] &&
        bond2AtomIdx != dblBond->getStereoAtoms()[1]) {
      sameTorsionDir = !sameTorsionDir;
    }
  }

  /*
     Time for some clarificatory text, because this gets really
     confusing really fast.

     The dihedral angle analysis above is based on viewing things
     with an atom order as follows:

     1
      \
       2 = 3
            \
             4

     so dihedrals > 90 correspond to sameDir=true

     however, the stereochemistry representation is
     based on something more like this:

     2
      \
       1 = 3
            \
             4
     (i.e. we consider the direction-setting single bonds to be
      starting at the double-bonded atom)

  */
  bool reverseBondDir = sameTorsionDir;

  Atom *atom1 = dblBond->getBeginAtom(), *atom2 = dblBond->getEndAtom();
  if (needsDir[bond1->getIdx()]) {
    BOOST_FOREACH (int bidx, singleBondNbrs[bond1->getIdx()]) {
      // std::cerr << "       neighbor from: " << bond1->getIdx() << " " << bidx
      //           << ": " << needsDir[bidx] << std::endl;
      if (needsDir[bidx]) {
        followupBonds.push_back(mol.getBondWithIdx(bidx));
      }
    }
  }
  if (needsDir[bond2->getIdx()]) {
    BOOST_FOREACH (int bidx, singleBondNbrs[bond2->getIdx()]) {
      // std::cerr << "       neighbor from: " << bond2->getIdx() << " " << bidx
      //           << ": " << needsDir[bidx] << std::endl;
      if (needsDir[bidx]) {
        followupBonds.push_back(mol.getBondWithIdx(bidx));
      }
    }
  }
  if (!needsDir[bond1->getIdx()]) {
    if (!needsDir[bond2->getIdx()]) {
      // check that we agree
    } else {
      if (bond1->getBeginAtom() != atom1) {
        reverseBondDir = !reverseBondDir;
      }
      setBondDirRelativeToAtom(bond2, atom2, bond1->getBondDir(),
                               reverseBondDir, needsDir);
    }
  } else if (!needsDir[bond2->getIdx()]) {
    if (bond2->getBeginAtom() != atom2) {
      reverseBondDir = !reverseBondDir;
    }
    setBondDirRelativeToAtom(bond1, atom1, bond2->getBondDir(), reverseBondDir,
                             needsDir);
  } else {
    setBondDirRelativeToAtom(bond1, atom1, Bond::ENDDOWNRIGHT, false, needsDir);
    setBondDirRelativeToAtom(bond2, atom2, Bond::ENDDOWNRIGHT, reverseBondDir,
                             needsDir);
  }
  needsDir[bond1->getIdx()] = 0;
  needsDir[bond2->getIdx()] = 0;
  if (obond1 && needsDir[obond1->getIdx()]) {
    setBondDirRelativeToAtom(obond1, atom1, bond1->getBondDir(),
                             bond1->getBeginAtom() == atom1, needsDir);
    needsDir[obond1->getIdx()] = 0;
  }
  if (obond2 && needsDir[obond2->getIdx()]) {
    setBondDirRelativeToAtom(obond2, atom2, bond2->getBondDir(),
                             bond2->getBeginAtom() == atom2, needsDir);
    needsDir[obond2->getIdx()] = 0;
  }
#if 0
  std::cerr << "  1:" << bond1->getIdx() << " ";
  if (obond1)
    std::cerr << obond1->getIdx() << std::endl;
  else
    std::cerr << "N/A" << std::endl;
  std::cerr << "  2:" << bond2->getIdx() << " ";
  if (obond2)
    std::cerr << obond2->getIdx() << std::endl;
  else
    std::cerr << "N/A" << std::endl;
  std::cerr << "**********************\n";
  std::cerr << "**********************\n";
  std::cerr << "**********************\n";
#endif
  BOOST_FOREACH (Bond *oDblBond, followupBonds) {
    // std::cerr << "FOLLOWUP: " << oDblBond->getIdx() << " "
    //           << needsDir[oDblBond->getIdx()] << std::endl;
    updateDoubleBondNeighbors(mol, oDblBond, conf, needsDir, singleBondCounts,
                              singleBondNbrs);
  }
}

bool isBondCandidateForStereo(const Bond *bond) {
  PRECONDITION(bond, "no bond");
  if (bond->getBondType() == Bond::DOUBLE &&
      bond->getStereo() != Bond::STEREOANY &&
      bond->getBondDir() != Bond::EITHERDOUBLE &&
      bond->getBeginAtom()->getDegree() > 1 &&
      bond->getEndAtom()->getDegree() > 1 &&
      shouldDetectDoubleBondStereo(bond)) {
    return true;
  }
  return false;
}

const Atom *findHighestCIPNeighbor(const Atom *atom, const Atom *skipAtom) {
  PRECONDITION(atom, "bad atom");

  unsigned bestCipRank = 0;
  const Atom *bestCipRankedAtom = nullptr;
  const auto &mol = atom->getOwningMol();
  
  for(auto *neighbor : atom->nbrs()) {
    if (neighbor == skipAtom) {
      continue;
    }
    unsigned cip = 0;
    if (!neighbor->getPropIfPresent(common_properties::_CIPRank, cip)) {
      // If at least one of the atoms doesn't have a CIP rank, the highest rank
      // does not make sense, so return a nullptr.
      return nullptr;
    } else if (cip > bestCipRank || bestCipRankedAtom == nullptr) {
      bestCipRank = cip;
      bestCipRankedAtom = neighbor;
    } else if (cip == bestCipRank) {
      // This also doesn't make sense if there is a tie (if that's possible).
      // We still keep the best CIP rank in case something better comes around
      // (also not sure if that's possible).
      BOOST_LOG(rdWarningLog)
          << "Warning: duplicate CIP ranks found in findHighestCIPNeighbor()"
          << std::endl;
      bestCipRankedAtom = nullptr;
    }
  }
  return bestCipRankedAtom;
}

}  // namespace

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
      if (mass >= 0) {
        mass += 1;
      }
    }
    mass += maxMass / 2;
    if (mass < 0) {
      mass = 0;
    } else {
      mass = mass % maxMass;
    }

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
      cipEntries[i].push_back(mol.atoms()[i]->getAtomicNum());
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
    for (int &index : allIndices) {
      CIP_ENTRY localEntry;
      localEntry.reserve(16);

      // start by pushing on our neighbors' ranks:
      for(auto *bond : mol.atoms()[index]->bonds()) {
        unsigned int nbrIdx = bond->getOtherAtomIdx(index);
        const Atom *nbr = mol.getAtomWithIdx(nbrIdx);
        int rank = ranks[nbrIdx] + 1;
        // put the neighbor in 2N times where N is the bond order as a double.
        // this is to treat aromatic linkages on fair footing. i.e. at least in
        // the
        // first iteration --c(:c):c and --C(=C)-C should look the same.
        // this was part of issue 3009911

        unsigned int count;
        if (bond->getBondType() == Bond::DOUBLE && nbr->getAtomicNum() == 15 &&
            (nbr->getDegree() == 4 || nbr->getDegree() == 3)) {
          // a special case for chiral phosphorous compounds
          // (this was leading to incorrect assignment of
          // R/S labels ):
          count = 1;

          // general justification of this is:
          // Paragraph 2.2. in the 1966 article is "Valence-Bond Conventions:
          // Multiple-Bond Unsaturation and Aromaticity". It contains several
          // conventions of which convention (b) is the one applying here:
          // "(b) Contributions by d orbitals to bonds of quadriligant atoms are
          // neglected."
          // FIX: this applies to more than just P
        } else {
          count = static_cast<unsigned int>(
              floor(2. * bond->getBondTypeAsDouble() + .1));
        }
        auto ePos =
            std::lower_bound(localEntry.begin(), localEntry.end(), rank);
        localEntry.insert(ePos, count, rank);
        ++nbr;
      }
      // add a zero for each coordinated H:
      // (as long as we're not a query atom)
      if (!mol.getAtomWithIdx(index)->hasQuery()) {
        localEntry.insert(localEntry.begin(), mol.getAtomWithIdx(index)->getTotalNumHs(), 0);
      }

      // we now have a sorted list of our neighbors' ranks,
      // copy it on in reversed order:
      cipEntries[index].insert(cipEntries[index].end(), localEntry.rbegin(),
                               localEntry.rend());
      if (cipEntries[index].size() > longestEntry) {
        longestEntry = rdcast<unsigned int>(cipEntries[index].size());
      }
    }
    // ----------------------------------------------------
    //
    // pad the entries so that we compare rounds to themselves:
    //
    for (int &index : allIndices) {
      auto sz = rdcast<unsigned int>(cipEntries[index].size());
      if (sz < longestEntry) {
        cipEntries[index].insert(cipEntries[index].end(), longestEntry - sz,
                                 -1);
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
  if (!ranks.size()) {
    ranks.resize(mol.getNumAtoms());
  }
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
    mol.atoms()[i]->setProp(common_properties::_CIPRank, ranks[i], 1);
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
  for(auto *bond : atom->bonds()) {
    // check whether this bond is explicitly set to have unknown stereo
    if (!hasExplicitUnknownStereo) {
      int explicit_unknown_stereo;
      if (bond->getBondDir() == Bond::UNKNOWN  // there's a squiggle bond
          || (bond->getPropIfPresent<int>(common_properties::_UnknownStereo,
                                          explicit_unknown_stereo) &&
              explicit_unknown_stereo)) {
        hasExplicitUnknownStereo = true;
      }
    }

    Bond::BondDir dir = bond->getBondDir();
    if (bond->getIdx() != refBond->getIdx()) {
      if (dir == Bond::ENDDOWNRIGHT || dir == Bond::ENDUPRIGHT) {
        seenDir = true;
        // If we're considering the bond "backwards", (i.e. from end
        // to beginning, reverse the effective direction:
        if (atom != bond->getBeginAtom()) {
          if (dir == Bond::ENDDOWNRIGHT) {
            dir = Bond::ENDUPRIGHT;
          } else {
            dir = Bond::ENDDOWNRIGHT;
          }
        }
      }
      Atom *nbrAtom = bond->getOtherAtom(atom);
      neighbors.push_back(std::make_pair(nbrAtom->getIdx(), dir));
    }
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
                             bool checkDir = false,
                             bool includeAromatic = false) {
  PRECONDITION(atom, "bad atom");
  PRECONDITION(refBond, "bad bond");
  neighbors.clear();

  for(auto *bond : atom->bonds()) {
    Bond::BondDir dir = bond->getBondDir();
    if ((bond->getBondType() == Bond::SINGLE ||
         (includeAromatic && bond->getBondType() == Bond::AROMATIC)) &&
        bond->getIdx() != refBond->getIdx()) {
      if (checkDir) {
        if ((dir != Bond::ENDDOWNRIGHT) && (dir != Bond::ENDUPRIGHT)) {
          continue;
        }
      }
      Atom *nbrAtom = bond->getOtherAtom(atom);
      neighbors.push_back(nbrAtom->getIdx());
    }
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
      std::vector<const Atom *> nonRingNbrs;
      std::vector<const Atom *> ringNbrs;
      for(auto *bond : atom->bonds()) {
        if (!ringInfo->numBondRings(bond->getIdx())) {
          nonRingNbrs.push_back(bond->getOtherAtom(atom));
        } else {
          const Atom *nbr = bond->getOtherAtom(atom);
          ringNbrs.push_back(nbr);
          unsigned int rnk = 0;
          nbr->getPropIfPresent(common_properties::_CIPRank, rnk);
          nbrRanks.insert(rnk);
        }
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
          if (ringNbrs.size() >= 2) {
            res = true;
          }
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
    if (atomsSeen[atom->getIdx()]) {
      continue;
    }
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
    for(auto *bond : atom->bonds()) {
      unsigned int bidx = bond->getIdx();
      if (!bondsSeen[bidx]) {
        bondsSeen.set(bidx);
        if (mol.getRingInfo()->numBondRings(bidx)) {
          const Atom *oatom = bond->getOtherAtom(atom);
          if (!atomsSeen[oatom->getIdx()]) {
            nextAtoms.push_back(oatom);
            atomsUsed.set(oatom->getIdx());
          }
        }
      }
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
      for(auto *bond : ratom->bonds()) {
        unsigned int bidx = bond->getIdx();
        if (!bondsSeen[bidx]) {
          bondsSeen.set(bidx);
          if (mol.getRingInfo()->numBondRings(bidx)) {
            const Atom *oatom = bond->getOtherAtom(ratom);
            if (!atomsSeen[oatom->getIdx()] && !atomsUsed[oatom->getIdx()]) {
              nextAtoms.push_back(oatom);
              atomsUsed.set(oatom->getIdx());
            }
          }
        }
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
        for (auto orae = rae + 1; orae != ringStereoAtoms.end(); ++orae) {
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
    for(auto *bond : atom->bonds()) {
      unsigned int otherIdx = bond->getOtherAtom(atom)->getIdx();
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
      nbrs.push_back(std::make_pair(ranks[otherIdx], bond->getIdx()));

    }

    // figure out if this is a legal chiral center or not:
    if (!hasDupes) {
      if (nbrs.size() < 3 &&
          (atom->getAtomicNum() != 15 && atom->getAtomicNum() != 33)) {
        // less than three neighbors is never stereogenic
        // unless it is a phosphine/arsine with implicit H
        legalCenter = false;
      } else if (atom->getAtomicNum() == 15 || atom->getAtomicNum() == 33) {
        // from logical flow: nbrs.size is 3 or 4, or 2 (implicit H)
        // Since InChI Software v. 1.02-standard (2009), phosphines and arsines
        // are always treated as stereogenic even with H atom neighbors.
        // Accept automatically.
        legalCenter = true;
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
          } else if (atom->getAtomicNum() == 7 &&
                     mol.getRingInfo()->isAtomInRingOfSize(atom->getIdx(), 3)) {
            // N in a three-membered ring is another one of the InChI special
            // cases
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
          if (tag == Atom::CHI_TETRAHEDRAL_CCW) {
            tag = Atom::CHI_TETRAHEDRAL_CW;
          } else {
            tag = Atom::CHI_TETRAHEDRAL_CCW;
          }
        }
        // now assign the CIP code:
        std::string cipCode;
        if (tag == Atom::CHI_TETRAHEDRAL_CCW) {
          cipCode = "S";
        } else {
          cipCode = "R";
        }
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
  boost::dynamic_bitset<> bondsToClear(mol.getNumBonds());
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
      if (shouldDetectDoubleBondStereo(dblBond)) {
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

            bool conflictingBegin =
                (begAtomNeighbors.size() == 2 &&
                 begAtomNeighbors[0].second == begAtomNeighbors[1].second);
            bool conflictingEnd =
                (endAtomNeighbors.size() == 2 &&
                 endAtomNeighbors[0].second == endAtomNeighbors[1].second);
            if (conflictingBegin || conflictingEnd) {
              dblBond->setStereo(Bond::STEREONONE);
              BOOST_LOG(rdWarningLog) << "Conflicting single bond directions "
                                         "around double bond at index "
                                      << dblBond->getIdx() << "." << std::endl;
              BOOST_LOG(rdWarningLog) << "  BondStereo set to STEREONONE and "
                                         "single bond directions set to NONE."
                                      << std::endl;
              assignedABond = true;
              if (conflictingBegin) {
                bondsToClear[mol.getBondBetweenAtoms(begAtomNeighbors[0].first,
                                                     begAtom->getIdx())
                                 ->getIdx()] = 1;
                bondsToClear[mol.getBondBetweenAtoms(begAtomNeighbors[1].first,
                                                     begAtom->getIdx())
                                 ->getIdx()] = 1;
              }
              if (conflictingEnd) {
                bondsToClear[mol.getBondBetweenAtoms(endAtomNeighbors[0].first,
                                                     endAtom->getIdx())
                                 ->getIdx()] = 1;
                bondsToClear[mol.getBondBetweenAtoms(endAtomNeighbors[1].first,
                                                     endAtom->getIdx())
                                 ->getIdx()] = 1;
              }
            } else {
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
            }
            --unassignedBonds;
          }
        }
      }
    }
  }

  for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {
    if (bondsToClear[i]) {
      mol.getBondWithIdx(i)->setBondDir(Bond::NONE);
    }
  }

  return std::make_pair(unassignedBonds > 0, assignedABond);
}

// reassign atom ranks by supplementing the current ranks
// with information about known chirality
void rerankAtoms(const ROMol &mol, UINT_VECT &ranks) {
  PRECONDITION(ranks.size() == mol.getNumAtoms(), "bad rank vector size");
  unsigned int factor = 100;
  while (factor < mol.getNumAtoms()) {
    factor *= 10;
  }

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
    for(auto *oBond : atom->bonds()) {
      if (oBond->getBondType() == Bond::DOUBLE) {
        if (oBond->getStereo() == Bond::STEREOE) {
          invars[i] += 1;
        } else if (oBond->getStereo() == Bond::STEREOZ) {
          invars[i] += 2;
        }
      }
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

bool hasStereoBondDir(const Bond *bond) {
  PRECONDITION(bond, "no bond");
  return bond->getBondDir() == Bond::BondDir::ENDDOWNRIGHT ||
         bond->getBondDir() == Bond::BondDir::ENDUPRIGHT;
}

const Bond *getNeighboringDirectedBond(const ROMol &mol, const Atom *atom) {
  PRECONDITION(atom, "no atom");
  for(auto *bond : atom->bonds()) {
    if (bond->getBondType() != Bond::BondType::DOUBLE &&
        hasStereoBondDir(bond)) {
      return bond;
    }
  }
  return nullptr;
}

Bond::BondStereo translateEZLabelToCisTrans(Bond::BondStereo label) {
  switch (label) {
    case Bond::STEREOE:
      return Bond::STEREOTRANS;
    case Bond::STEREOZ:
      return Bond::STEREOCIS;
    default:
      return label;
  }
}

INT_VECT findStereoAtoms(const Bond *bond) {
  PRECONDITION(bond, "bad bond");
  PRECONDITION(bond->hasOwningMol(), "no mol");
  PRECONDITION(bond->getBondType() == Bond::DOUBLE, "not double bond");
  PRECONDITION(bond->getStereo() > Bond::BondStereo::STEREOANY,
               "no defined stereo");

  if (!bond->getStereoAtoms().empty()) {
    return bond->getStereoAtoms();
  }
  if (bond->getStereo() == Bond::BondStereo::STEREOE ||
      bond->getStereo() == Bond::BondStereo::STEREOZ) {
    const Atom *startStereoAtom =
        findHighestCIPNeighbor(bond->getBeginAtom(), bond->getEndAtom());
    const Atom *endStereoAtom =
        findHighestCIPNeighbor(bond->getEndAtom(), bond->getBeginAtom());

    if (startStereoAtom == nullptr || endStereoAtom == nullptr) {
      return {};
    }

    int startStereoAtomIdx = static_cast<int>(startStereoAtom->getIdx());
    int endStereoAtomIdx = static_cast<int>(endStereoAtom->getIdx());

    return {startStereoAtomIdx, endStereoAtomIdx};
  } else {
    BOOST_LOG(rdWarningLog) << "Unable to assign stereo atoms for bond "
                            << bond->getIdx() << std::endl;
    return {};
  }
}

}  // namespace Chirality

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
      for(auto *nbond : (*bondIt)->getBeginAtom()->bonds()) {
        if (nbond->getBondDir() == Bond::ENDDOWNRIGHT ||
            nbond->getBondDir() == Bond::ENDUPRIGHT) {
          hasStereoBonds = true;
          break;
        }
      }
      for(auto *nbond : (*bondIt)->getEndAtom()->bonds()) {
        if (nbond->getBondDir() == Bond::ENDDOWNRIGHT ||
            nbond->getBondDir() == Bond::ENDUPRIGHT) {
          hasStereoBonds = true;
          break;
        }
      }
    }
    if (!cleanIt && hasStereoBonds) {
      break;  // no reason to keep iterating if we've already
              // determined there are stereo bonds to consider
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
    // if the ranks are needed again, this will force them to be
    // re-calculated based on the stereo calculated above.
    // atomRanks.clear();

    for (ROMol::AtomIterator atIt = mol.beginAtoms(); atIt != mol.endAtoms();
         ++atIt) {
      if ((*atIt)->hasProp(common_properties::_ringStereochemCand)) {
        (*atIt)->clearProp(common_properties::_ringStereochemCand);
      }
      if ((*atIt)->hasProp(common_properties::_ringStereoAtoms)) {
        (*atIt)->clearProp(common_properties::_ringStereoAtoms);
      }
    }
    boost::dynamic_bitset<> possibleSpecialCases(mol.getNumAtoms());
    Chirality::findChiralAtomSpecialCases(mol, possibleSpecialCases);

    for (auto atom : mol.atoms()) {
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
    for (auto bond : mol.bonds()) {
      // wedged bonds to atoms that have no stereochem
      // should be removed. (github issue 87)
      if ((bond->getBondDir() == Bond::BEGINWEDGE ||
           bond->getBondDir() == Bond::BEGINDASH) &&
          bond->getBeginAtom()->getChiralTag() == Atom::CHI_UNSPECIFIED &&
          bond->getEndAtom()->getChiralTag() == Atom::CHI_UNSPECIFIED) {
        bond->setBondDir(Bond::NONE);
      }

      // check for directionality on single bonds around
      // double bonds without stereo. This was github #2422
      if (bond->getBondType() == Bond::DOUBLE &&
          (bond->getStereo() == Bond::STEREOANY ||
           bond->getStereo() == Bond::STEREONONE)) {
        std::vector<Atom *> batoms = {bond->getBeginAtom(), bond->getEndAtom()};
        for (auto batom : batoms) {
            for (const auto &nbrBndI : batom->bonds()) {
            if ((nbrBndI->getBondDir() == Bond::ENDDOWNRIGHT ||
                 nbrBndI->getBondDir() == Bond::ENDUPRIGHT) &&
                (nbrBndI->getBondType() == Bond::SINGLE ||
                 nbrBndI->getBondType() == Bond::AROMATIC)) {
              // direction is set, and we know it's not because of our
              // bond. What about other neighbors?
              bool okToClear = true;
              for (const auto &nbrBndJ : nbrBndI->getOtherAtom(batom)->bonds()) {
                if (nbrBndJ->getBondType() == Bond::DOUBLE &&
                    nbrBndJ->getStereo() != Bond::STEREOANY &&
                    nbrBndJ->getStereo() != Bond::STEREONONE) {
                  okToClear = false;
                  break;
                }
              }
              if (okToClear) {
                nbrBndI->setBondDir(Bond::NONE);
              }
            }
          }
        }
      }
#if 0
      // make sure CIS/TRANS assignments are actually stereo bonds
      if ((*bondIt)->getBondType() == Bond::DOUBLE) {
        if ((*bondIt)->getStereo() == Bond::STEREOCIS ||
            (*bondIt)->getStereo() == Bond::STEREOTRANS) {
          if (!atomRanks.size()) {
            Chirality::assignAtomCIPRanks(mol, atomRanks);
          }

          const Atom *begAtom = (*bondIt)->getBeginAtom(),
                     *endAtom = (*bondIt)->getEndAtom();
          UINT_VECT begAtomNeighbors, endAtomNeighbors;
          Chirality::findAtomNeighborsHelper(mol, begAtom, *bondIt,
                                             begAtomNeighbors);
          Chirality::findAtomNeighborsHelper(mol, endAtom, *bondIt,
                                             endAtomNeighbors);

          // Note, this relies on this being a hydrogen-suppressed
          // graph as the 'Note' in the doc string of this function
          // indicates is a pre-condition.
          if ((begAtomNeighbors.size() == 2 &&
               atomRanks[begAtomNeighbors[0]] ==
                   atomRanks[begAtomNeighbors[1]]) ||
              (endAtomNeighbors.size() == 2 &&
               atomRanks[endAtomNeighbors[0]] ==
                   atomRanks[endAtomNeighbors[1]])) {
            (*bondIt)->setStereo(Bond::STEREONONE);
            (*bondIt)->getStereoAtoms().clear();
          }
        }
      }
#endif
    }
  }
  mol.setProp(common_properties::_StereochemDone, 1, true);

#if 0
  std::cerr << "---\n";
  mol.debugMol(std::cerr);
  std::cerr << "<<<<<<<<<<<<<<<<\n";
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
        // We ignore bonds flagged as EITHERDOUBLE or STEREOANY which have
        // stereo atoms set.
        if (dblBond->getBondDir() == Bond::EITHERDOUBLE ||
            (dblBond->getStereo() == Bond::STEREOANY &&
             dblBond->getStereoAtoms().size() == 2)) {
          continue;
        }
        // proceed only if we either want to clean the stereocode on this bond,
        // if none is set on it yet, or it is STEREOANY and we need to find
        // stereoatoms
        if (cleanIt || dblBond->getStereo() == Bond::STEREONONE ||
            dblBond->getStereo() == Bond::STEREOANY) {
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
            bool checkDir = false;
            bool includeAromatic = true;
            Chirality::findAtomNeighborsHelper(mol, begAtom, dblBond,
                                               begAtomNeighbors, checkDir,
                                               includeAromatic);
            Chirality::findAtomNeighborsHelper(mol, endAtom, dblBond,
                                               endAtomNeighbors, checkDir,
                                               includeAromatic);
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
  if (!mol.getNumConformers()) {
    return;
  }
  const Conformer &conf = mol.getConformer(confId);
  if (!conf.is3D()) {
    return;
  }

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
    const RDGeom::Point3D &p1 = conf.getAtomPos((*nbrIdx)->getIdx());
    ++nbrIdx;
    const RDGeom::Point3D &p2 = conf.getAtomPos((*nbrIdx)->getIdx());
    ++nbrIdx;
    const RDGeom::Point3D &p3 = conf.getAtomPos((*nbrIdx)->getIdx());

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

void assignChiralTypesFromMolParity(ROMol &mol, bool replaceExistingTags) {
  static const std::vector<Atom::ChiralType> chiralTypeVect{
      Atom::CHI_TETRAHEDRAL_CW, Atom::CHI_TETRAHEDRAL_CCW};
  // if the molecule already has stereochemistry
  // perceived, remove the flags that indicate
  // this... what we're about to do will require
  // that we go again.
  if (mol.hasProp(common_properties::_StereochemDone)) {
    mol.clearProp(common_properties::_StereochemDone);
  }
  // Atom-based parity
  // Number the atoms surrounding the stereo center with 1, 2, 3, and 4
  // in order of increasing atom number (position in the atom block)
  // (an implicit hydrogen should be considered the highest numbered atom).
  // View the center from a position such that the bond connecting the
  // highest-numbered atom (4) projects behind the plane formed by
  // atoms 1, 2, and 3.
  //
  // Parity 1 (CW)        Parity 2 (CCW)
  //     3   1                3   2
  //      \ /                  \ /
  //       |                    |
  //       2                    1
  //
  for (auto atom : mol.atoms()) {
    // if we aren't replacing existing tags and the atom is already tagged,
    // punt:
    if (!replaceExistingTags && atom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
      continue;
    }
    int parity = 0;
    atom->getPropIfPresent(common_properties::molParity, parity);
    if (parity <= 0 || parity > 2 || atom->getDegree() < 3) {
      atom->setChiralTag(Atom::CHI_UNSPECIFIED);
      continue;
    }
    // if we are here, parity was 1 (CW) or 2 (CCW)
    // now we set parity 0 to be CW and 1 to be CCW
    --parity;
    RDKit::ROMol::OBOND_ITER_PAIR nbrBonds = mol.getAtomBonds(atom);
    INT_LIST nbrBondIdxList;
    std::transform(
        nbrBonds.first, nbrBonds.second, std::back_inserter(nbrBondIdxList),
        [mol](const ROMol::edge_descriptor &e) { return e->getIdx(); });
    unsigned int atomIdx = atom->getIdx();
    nbrBondIdxList.sort([mol, atomIdx](const int ai, const int bi) {
      return (mol.getBondWithIdx(ai)->getOtherAtomIdx(atomIdx) <
              mol.getBondWithIdx(bi)->getOtherAtomIdx(atomIdx));
    });
    int nSwaps = atom->getPerturbationOrder(nbrBondIdxList);
    if (nSwaps % 2) {
      parity = 1 - parity;
    }
    atom->setChiralTag(chiralTypeVect[parity]);
    if (atom->getImplicitValence() == -1) {
      atom->calcExplicitValence(false);
      atom->calcImplicitValence(false);
    }
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

void setDoubleBondNeighborDirections(ROMol &mol, const Conformer *conf) {
  // used to store the number of single bonds a given
  // single bond is adjacent to
  std::vector<unsigned int> singleBondCounts(mol.getNumBonds(), 0);
  std::vector<Bond *> bondsInPlay;
  // keeps track of which single bonds are adjacent to each double bond:
  VECT_INT_VECT dblBondNbrs(mol.getNumBonds());
  // keeps track of which double bonds are adjacent to each single bond:
  VECT_INT_VECT singleBondNbrs(mol.getNumBonds());
  // keeps track of which single bonds need a dir set and which double bonds
  // need to have their neighbors' dirs set
  boost::dynamic_bitset<> needsDir(mol.getNumBonds());

  // find double bonds that should be considered for
  // stereochemistry
  // NOTE that we are explicitly excluding double bonds in rings
  // with this test.
  bool resetRings = false;
  if (!mol.getRingInfo()->isInitialized()) {
    resetRings = true;
    MolOps::fastFindRings(mol);
  }

  for (RWMol::BondIterator bondIt = mol.beginBonds(); bondIt != mol.endBonds();
       ++bondIt) {
    if (isBondCandidateForStereo(*bondIt)) {
      const Atom *a1 = (*bondIt)->getBeginAtom();
      const Atom *a2 = (*bondIt)->getEndAtom();

      for(auto *nbrBond : a1->bonds())
        if (nbrBond->getBondType() == Bond::SINGLE ||
            nbrBond->getBondType() == Bond::AROMATIC) {
          singleBondCounts[nbrBond->getIdx()] += 1;
          auto nbrDir = nbrBond->getBondDir();
          if (nbrDir == Bond::BondDir::NONE ||
              nbrDir == Bond::BondDir::ENDDOWNRIGHT ||
              nbrDir == Bond::BondDir::ENDUPRIGHT) {
            needsDir[nbrBond->getIdx()] = 1;
          }
          needsDir[(*bondIt)->getIdx()] = 1;
          dblBondNbrs[(*bondIt)->getIdx()].push_back(nbrBond->getIdx());
          // the search may seem inefficient, but these vectors are going to
          // be at most 2 long (with very few exceptions). It's just not worth
          // using a different data structure
          if (std::find(singleBondNbrs[nbrBond->getIdx()].begin(),
                        singleBondNbrs[nbrBond->getIdx()].end(),
                        (*bondIt)->getIdx()) ==
              singleBondNbrs[nbrBond->getIdx()].end()) {
            singleBondNbrs[nbrBond->getIdx()].push_back((*bondIt)->getIdx());
          }
        }
      
      //boost::tie(beg, end) = mol.getAtomBonds(a2);
      //while (beg != end) {
      //  const Bond *nbrBond = mol[*beg];
      for(auto *nbrBond : a2->bonds()) {
        if (nbrBond->getBondType() == Bond::SINGLE ||
            nbrBond->getBondType() == Bond::AROMATIC) {
          singleBondCounts[nbrBond->getIdx()] += 1;
          auto nbrDir = nbrBond->getBondDir();
          if (nbrDir == Bond::BondDir::NONE ||
              nbrDir == Bond::BondDir::ENDDOWNRIGHT ||
              nbrDir == Bond::BondDir::ENDUPRIGHT) {
            needsDir[nbrBond->getIdx()] = 1;
          }
          needsDir[(*bondIt)->getIdx()] = 1;
          dblBondNbrs[(*bondIt)->getIdx()].push_back(nbrBond->getIdx());

          // the search may seem inefficient, but these vectors are going to
          // be at most 2 long (with very few exceptions). It's just not worth
          // using a different data structure
          if (std::find(singleBondNbrs[nbrBond->getIdx()].begin(),
                        singleBondNbrs[nbrBond->getIdx()].end(),
                        (*bondIt)->getIdx()) ==
              singleBondNbrs[nbrBond->getIdx()].end()) {
            singleBondNbrs[nbrBond->getIdx()].push_back((*bondIt)->getIdx());
          }
        }
      }
      bondsInPlay.push_back(*bondIt);
    }
  }

  if (!bondsInPlay.size()) {
    if (resetRings) {
      mol.getRingInfo()->reset();
    }
    return;
  }

  // order the double bonds based on the singleBondCounts of their neighbors:
  std::vector<std::pair<unsigned int, Bond *>> orderedBondsInPlay;
  for (auto dblBond : bondsInPlay) {
    unsigned int countHere =
        std::accumulate(dblBondNbrs[dblBond->getIdx()].begin(),
                        dblBondNbrs[dblBond->getIdx()].end(), 0);
    // and favor double bonds that are *not* in rings. The combination of
    // using the sum above (instead of the max) and this ring-membershipt test
    // seem to fix sf.net issue 3009836
    if (!(mol.getRingInfo()->numBondRings(dblBond->getIdx()))) {
      countHere *= 10;
    }
    orderedBondsInPlay.push_back(std::make_pair(countHere, dblBond));
  }
  std::sort(orderedBondsInPlay.begin(), orderedBondsInPlay.end());

  // oof, now loop over the double bonds in that order and
  // update their neighbor directionalities:
  std::vector<std::pair<unsigned int, Bond *>>::reverse_iterator pairIter;
  for (pairIter = orderedBondsInPlay.rbegin();
       pairIter != orderedBondsInPlay.rend(); ++pairIter) {
    // std::cerr << "RESET?: " << pairIter->second->getIdx() << " "
    //           << pairIter->second->getStereo() << std::endl;
    updateDoubleBondNeighbors(mol, pairIter->second, conf, needsDir,
                              singleBondCounts, singleBondNbrs);
  }
  if (resetRings) {
    mol.getRingInfo()->reset();
  }
}

void detectBondStereochemistry(ROMol &mol, int confId) {
  if (!mol.getNumConformers()) {
    return;
  }
  const Conformer &conf = mol.getConformer(confId);
  setDoubleBondNeighborDirections(mol, &conf);
}

void setBondStereoFromDirections(ROMol &mol) {
  for (Bond *bond : mol.bonds()) {
    if (bond->getBondType() == Bond::DOUBLE) {
      const Atom *stereoBondBeginAtom = bond->getBeginAtom();
      const Atom *stereoBondEndAtom = bond->getEndAtom();

      const Bond *directedBondAtBegin =
          Chirality::getNeighboringDirectedBond(mol, stereoBondBeginAtom);
      const Bond *directedBondAtEnd =
          Chirality::getNeighboringDirectedBond(mol, stereoBondEndAtom);

      if (directedBondAtBegin != nullptr && directedBondAtEnd != nullptr) {
        unsigned beginSideStereoAtom =
            directedBondAtBegin->getOtherAtomIdx(stereoBondBeginAtom->getIdx());
        unsigned endSideStereoAtom =
            directedBondAtEnd->getOtherAtomIdx(stereoBondEndAtom->getIdx());

        bond->setStereoAtoms(beginSideStereoAtom, endSideStereoAtom);

        auto beginSideBondDirection = directedBondAtBegin->getBondDir();
        if (directedBondAtBegin->getBeginAtom() == stereoBondBeginAtom) {
          beginSideBondDirection = getOppositeBondDir(beginSideBondDirection);
        }

        auto endSideBondDirection = directedBondAtEnd->getBondDir();
        if (directedBondAtEnd->getEndAtom() == stereoBondEndAtom) {
          endSideBondDirection = getOppositeBondDir(endSideBondDirection);
        }

        if (beginSideBondDirection == endSideBondDirection) {
          bond->setStereo(Bond::STEREOTRANS);
        } else {
          bond->setStereo(Bond::STEREOCIS);
        }
      }
    }
  }
}

void assignStereochemistryFrom3D(ROMol &mol, int confId,
                                 bool replaceExistingTags) {
  if (!mol.getNumConformers() || !mol.getConformer(confId).is3D()) {
    return;
  }

  detectBondStereochemistry(mol, confId);
  assignChiralTypesFrom3D(mol, confId, replaceExistingTags);
  bool force = true;
  bool flagPossibleStereoCenters = true;
  assignStereochemistry(mol, replaceExistingTags, force,
                        flagPossibleStereoCenters);
}

void assignChiralTypesFromBondDirs(ROMol &mol, const int confId,
                                   const bool replaceExistingTags) {
  if (!mol.getNumConformers()) {
    return;
  }
  auto conf = mol.getConformer(confId);
  boost::dynamic_bitset<> atomsSet(mol.getNumAtoms(), 0);
  for (auto &bond : mol.bonds()) {
    const Bond::BondDir dir = bond->getBondDir();
    if (dir != Bond::UNKNOWN) {
      // the bond is marked as chiral:
      if (dir == Bond::BEGINWEDGE || dir == Bond::BEGINDASH) {
        Atom *atom = bond->getBeginAtom();
        if (atomsSet[atom->getIdx()] ||
            (!replaceExistingTags &&
             atom->getChiralTag() != Atom::CHI_UNSPECIFIED)) {
          continue;
        }
        if (atom->getImplicitValence() == -1) {
          atom->calcExplicitValence(false);
          atom->calcImplicitValence(false);
        }
        Atom::ChiralType code = atomChiralTypeFromBondDir(mol, bond, &conf);
        if (code != Atom::ChiralType::CHI_UNSPECIFIED) {
          atomsSet.set(atom->getIdx());
          //   std::cerr << "atom " << atom->getIdx() << " code " << code
          //             << " from bond " << bond->getIdx() << std::endl;
        }
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
