//
//  Copyright (C) 2004-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Chirality.h"

#include <Geometry/point.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/Ranking.h>
#include <GraphMol/new_canon.h>
#include <GraphMol/Atropisomers.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Ranking.h>
#include <RDGeneral/types.h>
#include <RDGeneral/utils.h>

#include <boost/dynamic_bitset.hpp>
#include <boost/algorithm/string.hpp>

#include <algorithm>
#include <cstdlib>
#include <optional>
#include <set>
#include <sstream>
#include <utility>

// #define VERBOSE_CANON 1

namespace RDKit {

namespace {
bool shouldDetectDoubleBondStereo(const Bond *bond) {
  const RingInfo *ri = bond->getOwningMol().getRingInfo();
  return (!ri->numBondRings(bond->getIdx()) ||
          ri->minBondRingSize(bond->getIdx()) >=
              Chirality::minRingSizeForDoubleBondStereo);
}

bool getValFromEnvironment(const char *var, bool defVal) {
  auto evar = std::getenv(var);
  if (evar != nullptr) {
    if (!strcmp(evar, "0")) {
      return false;
    } else {
      return true;
    }
  }
  return defVal;
}

bool is_regular_h(const Atom &atom) {
  return atom.getAtomicNum() == 1 && atom.getIsotope() == 0;
}

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
                              bool reverse, boost::dynamic_bitset<> &) {
  PRECONDITION(bond, "bad bond");
  PRECONDITION(atom, "bad atom");
  PRECONDITION(dir == Bond::ENDUPRIGHT || dir == Bond::ENDDOWNRIGHT, "bad dir");
  PRECONDITION(atom == bond->getBeginAtom() || atom == bond->getEndAtom(),
               "atom doesn't belong to bond");

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

bool isLinearArrangement(const RDGeom::Point3D &v1, const RDGeom::Point3D &v2) {
  double lsq = v1.lengthSq() * v2.lengthSq();

  // treat zero length vectors as linear
  if (lsq < 1.0e-6) {
    return true;
  }

  double dotProd = v1.dotProduct(v2);

  double cos178 =
      -0.999388;  // == cos(M_PI-0.035), corresponds to a tolerance of 2 degrees
  return dotProd < cos178 * sqrt(lsq);
}

void controllingBondFromAtom(const ROMol &mol,
                             const boost::dynamic_bitset<> &needsDir,
                             const std::vector<unsigned int> &singleBondCounts,
                             const Bond *dblBond, const Atom *atom, Bond *&bond,
                             Bond *&obond, bool &squiggleBondSeen,
                             bool &doubleBondSeen) {
  bond = nullptr;
  obond = nullptr;
  for (const auto tBond : mol.atomBonds(atom)) {
    if (tBond == dblBond) {
      continue;
    }
    if ((tBond->getBondType() == Bond::SINGLE ||
         tBond->getBondType() == Bond::AROMATIC) &&
        (tBond->getBondDir() == Bond::BondDir::NONE ||
         tBond->getBondDir() == Bond::BondDir::ENDDOWNRIGHT ||
         tBond->getBondDir() == Bond::BondDir::ENDUPRIGHT)) {
      // prefer bonds that already have their directionality set
      // or that are adjacent to more double bonds:
      if (!bond) {
        bond = tBond;
      } else if (needsDir[tBond->getIdx()]) {
        if (singleBondCounts[tBond->getIdx()] >
            singleBondCounts[bond->getIdx()]) {
          obond = bond;
          bond = tBond;
        } else {
          obond = tBond;
        }
      } else {
        obond = bond;
        bond = tBond;
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

  controllingBondFromAtom(mol, needsDir, singleBondCounts, dblBond,
                          dblBond->getBeginAtom(), bond1, obond1,
                          squiggleBondSeen, doubleBondSeen);

  // Don't do any direction setting if we've seen a squiggle bond, but do mark
  // the double bond as a crossed bond and return
  if (squiggleBondSeen) {
    Chirality::detail::setStereoForBond(mol, dblBond, Bond::STEREOANY);
    return;
  }
  if (!bond1) {
    return;
  }

  Bond *bond2 = nullptr, *obond2 = nullptr;
  controllingBondFromAtom(mol, needsDir, singleBondCounts, dblBond,
                          dblBond->getEndAtom(), bond2, obond2,
                          squiggleBondSeen, doubleBondSeen);

  // Don't do any direction setting if we've seen a squiggle bond, but do mark
  // the double bond as a crossed bond and return
  if (squiggleBondSeen) {
    Chirality::detail::setStereoForBond(mol, dblBond, Bond::STEREOANY);
    return;
  }
  if (!bond2) {
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
      Chirality::detail::setStereoForBond(mol, dblBond, Bond::STEREOANY);
      return;
    }

    double ang = RDGeom::computeDihedralAngle(bond1P, beginP, endP, bond2P);
    sameTorsionDir = ang >= M_PI / 2;
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
    for (auto bidx : singleBondNbrs[bond1->getIdx()]) {
      // std::cerr << "       neighbor from: " << bond1->getIdx() << " " << bidx
      //           << ": " << needsDir[bidx] << std::endl;
      if (needsDir[bidx]) {
        followupBonds.push_back(mol.getBondWithIdx(bidx));
      }
    }
  }
  if (needsDir[bond2->getIdx()]) {
    for (auto bidx : singleBondNbrs[bond2->getIdx()]) {
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
  for (Bond *oDblBond : followupBonds) {
    // std::cerr << "FOLLOWUP: " << oDblBond->getIdx() << " "
    //           << needsDir[oDblBond->getIdx()] << std::endl;
    updateDoubleBondNeighbors(mol, oDblBond, conf, needsDir, singleBondCounts,
                              singleBondNbrs);
  }
}

bool isBondCandidateForStereo(const Bond *bond) {
  PRECONDITION(bond, "no bond");
  return bond->getBondType() == Bond::DOUBLE &&
         bond->getStereo() != Bond::STEREOANY &&
         bond->getBondDir() != Bond::EITHERDOUBLE &&
         bond->getBeginAtom()->getDegree() > 1u &&
         bond->getEndAtom()->getDegree() > 1u &&
         shouldDetectDoubleBondStereo(bond);
}

const Atom *findHighestCIPNeighbor(const Atom *atom, const Atom *skipAtom) {
  PRECONDITION(atom, "bad atom");

  unsigned bestCipRank = 0;
  const Atom *bestCipRankedAtom = nullptr;
  const auto &mol = atom->getOwningMol();

  for (const auto neighbor : mol.atomNeighbors(atom)) {
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

std::optional<Atom::ChiralType> atomChiralTypeFromBondDirPseudo3D(
    const ROMol &mol, const Bond *bond, const Conformer *conf,
    double pseudo3DOffset = 0.1, double volumeTolerance = 0.01) {
  PRECONDITION(bond, "no bond");
  PRECONDITION(conf, "no conformer");
  auto bondDir = bond->getBondDir();
  PRECONDITION(bondDir == Bond::BEGINWEDGE || bondDir == Bond::BEGINDASH,
               "bad bond direction");
  constexpr double coordZeroTol = 1e-4;
  constexpr double zeroTol = 1e-3;

  // NOTE that according to the CT file spec, wedging assigns chirality
  // to the atom at the point of the wedge, (atom 1 in the bond).
  const auto atom = bond->getBeginAtom();
  PRECONDITION(atom, "no atom");

  // we can't do anything with atoms that have more than 4 neighbors:
  if (atom->getDegree() > 4) {
    return Atom::CHI_UNSPECIFIED;
  }
  const auto bondAtom = bond->getEndAtom();

  Atom::ChiralType res = Atom::CHI_UNSPECIFIED;

  auto centerLoc = conf->getAtomPos(atom->getIdx());
  centerLoc.z = 0.0;
  auto refPt = conf->getAtomPos(bondAtom->getIdx());
  refPt.z =
      bondDir == Bond::BondDir::BEGINWEDGE ? pseudo3DOffset : -pseudo3DOffset;

  //----------------------------------------------------------
  //
  //  collect indices and bond vectors of neighbors and track whether or
  //  not there's an H neighbor and if all bonds are single
  //
  //  at the end of this process bond 0 is the input wedged bond
  //
  //----------------------------------------------------------
  bool hSeen = false;

  INT_VECT neighborBondIndices;
  if (is_regular_h(*bondAtom)) {
    hSeen = true;
  }

  unsigned int refIdx = mol.getNumBonds() + 1;
  std::vector<RDGeom::Point3D> bondVects;
  bool allSingle = true;
  unsigned int nbrIdx = 0;
  for (const auto nbrBond : mol.atomBonds(atom)) {
    const auto oAtom = nbrBond->getOtherAtom(atom);
    auto tmpPt = conf->getAtomPos(oAtom->getIdx());
    if (nbrBond == bond) {
      refIdx = nbrIdx;
      tmpPt = refPt;
    } else {
      // theoretically we could confirm that this is a single bond,
      // but it's not impossible that at some point in the future we
      // could allow wedged multiple bonds for things like atropisomers
      if (nbrBond->getBeginAtomIdx() == atom->getIdx() &&
          (nbrBond->getBondDir() == Bond::BondDir::BEGINWEDGE ||
           nbrBond->getBondDir() == Bond::BondDir::BEGINDASH)) {
        tmpPt.z = nbrBond->getBondDir() == Bond::BondDir::BEGINWEDGE
                      ? pseudo3DOffset
                      : -pseudo3DOffset;
      } else {
        tmpPt.z = 0;
      }
      // check for overly short bonds. Note that we're doing this check *after*
      // adjusting the z coordinate.
      //    We want to allow atoms to overlap in x-y space if they are connected
      //    via a wedged bond.
      if ((centerLoc - tmpPt).lengthSq() < zeroTol) {
        BOOST_LOG(rdWarningLog)
            << "Warning: ambiguous stereochemistry - zero-length (or near zero-length) bond - at atom "
            << atom->getIdx() << " ignored." << std::endl;
        return std::nullopt;
      }
    }
    ++nbrIdx;
    if (nbrBond->getBondType() != Bond::SINGLE) {
      allSingle = false;
    }
    bondVects.push_back(centerLoc.directionVector(tmpPt));
    if (is_regular_h(*oAtom)) {
      hSeen = true;
    }
    neighborBondIndices.push_back(nbrBond->getIdx());
  }
  CHECK_INVARIANT(refIdx < mol.getNumBonds(),
                  "could not find reference bond in neighbors");

  auto nNbrs = bondVects.size();

  //----------------------------------------------------------
  //
  //  Return now if there aren't at least 3 non-H bonds to the atom.
  //  (we can implicitly add a single H to 3 coordinate atoms, but
  //  we're horked otherwise).
  //
  //----------------------------------------------------------
  if (nNbrs < 3 || nNbrs > 4 || (hSeen && nNbrs < 4)) {
    return std::nullopt;
  }

  //----------------------------------------------------------
  //  Check for neighbor atoms which overlap
  //----------------------------------------------------------
  for (auto i = 0u; i < nNbrs; ++i) {
    for (auto j = 0u; j < i; ++j) {
      if ((bondVects[i] - bondVects[j]).lengthSq() < zeroTol) {
        BOOST_LOG(rdWarningLog)
            << "Warning: ambiguous stereochemistry - overlapping neighbors  - at atom "
            << atom->getIdx() << " ignored" << std::endl;
        return std::nullopt;
      }
    }
  }

  //----------------------------------------------------------
  //
  //  Continue if there are all single bonds or if we're considering
  //  4-coordinate P or S
  //
  //----------------------------------------------------------
  if (allSingle || atom->getAtomicNum() == 15 || atom->getAtomicNum() == 16) {
    double vol;
    unsigned int order[4] = {0, 1, 2, 3};
    double prefactor = 1;
    if (refIdx != 0) {
      // bring the wedged bond to the front so that we always consider it
      std::swap(order[0], order[refIdx]);
      prefactor *= -1;
    }

    // check for the case that bonds 1 and 2 are co-linear but 1 and 0 are
    // not:
    if (nNbrs > 3 &&
        bondVects[order[1]].crossProduct(bondVects[order[2]]).lengthSq() <
            10 * zeroTol &&
        bondVects[order[1]].crossProduct(bondVects[order[0]]).lengthSq() >
            10 * zeroTol) {
      bondVects[order[1]].z = bondVects[order[0]].z * -1;
      // that bondVect is no longer normalized, but this hopefully won't break
      // anything
    }

    //----------------------------------------------------------
    //
    // order the bonds so that the rotation order is:
    //   0 - 1 - 2        for three coordinate
    // or
    //   0 - 1 - 2 - 3    for four coordinate
    //
    // this makes the rest of the code a lot simpler
    //
    //----------------------------------------------------------

    // checks to see if the vectors 1 and 2 need to have their order
    //    relative to vector 0 swapped.
    // we don't actually pass the vectors in, but use their cross products
    // and dot products to vector 0 to figure out if they need to be swapped
#if defined(__clang__)
// Clang apparently doesn't need to capture the constexpr zeroTol, and complains
// about it being specified, but MSVC does need it, and removing it will break
// the build
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-lambda-capture"
#endif
    auto needsSwap = [&zeroTol](const RDGeom::Point3D &cp01,
                                const RDGeom::Point3D &cp02, double dp01,
                                double dp02) -> bool {
      if (fabs(dp01) - 1 > -zeroTol) {
        if (cp02.z < 0) {
          return true;
        }
        return false;
      }
      if (fabs(dp02) - 1 > -zeroTol) {
        if (cp01.z < 0) {
          return true;
        }
      }

      if ((cp01.z * cp02.z) < -zeroTol) {
        if (cp01.z < cp02.z) {
          return true;
        }
        return false;
      }
      if (dp01 * dp02 < -zeroTol) {
        if (dp01 < dp02) {
          return true;
        }
        return false;
      }
      return fabs(dp01) > fabs(dp02);
    };
#if defined(__clang__)
#pragma GCC diagnostic pop
#endif

    if (nNbrs == 3) {
      // this case is simple, we either need to swap vectors 1 and 2 or we
      // don't:
      auto cp01 = bondVects[order[0]].crossProduct(bondVects[order[1]]);
      auto cp02 = bondVects[order[0]].crossProduct(bondVects[order[2]]);
      auto dp01 = bondVects[order[0]].dotProduct(bondVects[order[1]]);
      auto dp02 = bondVects[order[0]].dotProduct(bondVects[order[2]]);
      if (needsSwap(cp01, cp02, dp01, dp02)) {
        std::swap(order[1], order[2]);
        prefactor *= -1;
      }
    } else if (nNbrs > 3) {
      // here there are more permutations. Rather than hand-coding all of them
      // we'll just sort bonds 1, 2, and 3 based on their cross- and dot-
      // products to bond 0
      std::vector<std::tuple<double, double, unsigned>> orderedBonds(3);
      for (auto i = 1u; i < 4; ++i) {
        auto cp0i = bondVects[order[0]].crossProduct(bondVects[order[i]]);
        auto sgn = cp0i.z < -zeroTol ? -1 : 1;
        auto dp0i = bondVects[order[0]].dotProduct(bondVects[order[i]]);
        orderedBonds[i - 1] = std::make_tuple(sgn, sgn * dp0i, order[i]);
      }
      std::sort(orderedBonds.rbegin(), orderedBonds.rend());

      // update the order array and figure out whether or not we've done a
      // cyclic permutation
      auto nChanged = 0;
      for (auto i = 1u; i < 4; ++i) {
        auto ni = std::get<2>(orderedBonds[i - 1]);
        if (order[i] != ni) {
          order[i] = ni;
          ++nChanged;
        }
      }
      if (nChanged == 2) {
        // this is always an acyclic permutation
        prefactor *= -1;
      }
    }

    // std::cerr<<"ORDER "<<neighborBondIndices[order[0]]<<"
    // "<<neighborBondIndices[order[1]]<<" "<<neighborBondIndices[order[2]]<<"
    // "<<neighborBondIndices[order[3]]<<std::endl;

    // check for opposing bonds with opposite wedging
    for (auto i = 0u; i < nNbrs; ++i) {
      for (auto j = i + 1; j < nNbrs; ++j) {
        if (bondVects[order[i]].z * bondVects[order[j]].z < -zeroTol) {
          auto cp =
              bondVects[order[i]].crossProduct(bondVects[order[j]]).lengthSq();
          if (cp < 0.01) {
            // exception to our rejection of these structures: in some horrible
            // pseudo-3D drawings of things like sugars the ring substituents
            // are drawn 180 degrees apart and with opposite wedging. Let that
            // one pass.
            if (nNbrs == 4 &&
                fabs(bondVects[order[i]].dotProduct(bondVects[order[j]]) + 1) <
                    zeroTol) {
              // this is allowed for neighboring bonds
              if (j - i == 1 || (i == 0 && j == 3)) {
                // std::cerr << " skip it " << std::endl;
                bondVects[order[j]].z = 0.0;
                continue;
              }
            }
            BOOST_LOG(rdWarningLog)
                << "Warning: ambiguous stereochemistry - opposing bonds have opposite wedging - at atom "
                << atom->getIdx() << " ignored." << std::endl;
            return std::nullopt;
          }
        }
      }
    }

    // three-coordinate special cases where chirality cannot be determined
    //
    //  Case 1:
    //  this one is never allowed with different directions for the bonds to 1
    //  and 2
    //     0   2
    //      \ /
    //       C
    //       *
    //       1
    //   This is ST-1.2.10 in the IUPAC guidelines
    //
    //  Case 2: all bonds are wedged in the same direction
    if (nNbrs == 3) {
      bool conflict = false;
      if (bondVects[order[1]].z * bondVects[order[0]].z < -coordZeroTol &&
          fabs(bondVects[order[2]].z) < coordZeroTol) {
        conflict = bondVects[order[2]].crossProduct(bondVects[order[0]]).z *
                       bondVects[order[2]].crossProduct(bondVects[order[1]]).z <
                   -1e-4;
      } else if (bondVects[order[2]].z * bondVects[order[0]].z <
                     -coordZeroTol &&
                 fabs(bondVects[order[1]].z) < coordZeroTol) {
        conflict = bondVects[order[1]].crossProduct(bondVects[order[0]]).z *
                       bondVects[order[1]].crossProduct(bondVects[order[2]]).z <
                   -coordZeroTol;
      }
      if (conflict) {
        BOOST_LOG(rdWarningLog)
            << "Warning: conflicting stereochemistry - bond wedging contradiction - at atom "
            << atom->getIdx() << " ignored" << std::endl;
        return std::nullopt;
      }
    }
    // for the purposes of the cross products we ignore any pseudo-3D
    // coordinates
    auto bv1 = bondVects[order[1]];
    bv1.z = 0;
    auto bv2 = bondVects[order[2]];
    bv2.z = 0;
    auto crossp1 = bv1.crossProduct(bv2);
    // catch linear arrangements
    if (nNbrs == 3) {
      if (crossp1.lengthSq() < 5 * zeroTol) {
        // in a linear relationship with three neighbors we assume that the
        // two perpendicular bonds are wedged in the other direction from the
        // one that was provided.
        // that's this situation:
        //
        //              0
        //              |   <- wedged up
        //           1--C--2
        //
        //  here we assume that bonds C-1 and C-2 are wedged down
        //
        // ST-1.2.12 of the IUPAC guidelines says that this form is wrong since
        // it's for a "T-shaped" configuration instead of a tetrahedron, but it
        // shows up fairly frequently, particularly with fused ring systems
        bv1.z = -bondVects[order[0]].z;
        bv2.z = -bondVects[order[0]].z;
        crossp1 = bv1.crossProduct(bv2);
      }
    } else if (crossp1.lengthSq() < 10 * zeroTol) {
      // if the other bond is flat:
      if (fabs(bondVects[order[3]].z) < coordZeroTol) {
        // By construction this is a neighboring bond, so make it the opposite
        // wedging from us.
        bondVects[order[3]].z = -1 * bondVects[order[0]].z;
        // that bondVect is no longer normalized, but this hopefully won't break
        // anything
      }
    }
    vol = crossp1.dotProduct(bondVects[order[0]]);
    if (nNbrs == 4) {
      const auto dotp1 = bondVects[order[1]].dotProduct(bondVects[order[2]]);
      // for the purposes of the cross products we ignore any pseudo-3D
      // coordinates
      auto bv3 = bondVects[order[3]];
      bv3.z = 0;
      const auto crossp2 = bv1.crossProduct(bv3);
      const auto dotp2 = bondVects[order[1]].dotProduct(bondVects[order[3]]);
      auto vol2 = crossp2.dotProduct(bondVects[order[0]]);
#if 0
      std::cerr << neighborBondIndices[order[0]] << " " <<
      bondVects[order[0]]
                << std::endl;
      std::cerr << neighborBondIndices[order[1]] << " " <<
      bondVects[order[1]]
                << std::endl;
      std::cerr << neighborBondIndices[order[2]] << " " <<
      bondVects[order[2]]
                << std::endl;
      std::cerr << neighborBondIndices[order[3]] << " " <<
      bondVects[order[3]]
                << std::endl;
      std::cerr << "------------" << std::endl;
      std::cerr << crossp1 << " l2=" << crossp1.lengthSq()<<" " << std::endl;
      std::cerr << crossp2 << " l2=" << crossp2.lengthSq()<< " " << dotp2 << std::endl;
      std::cerr << " !!! " << vol << " " << vol2 << std::endl;
#endif

      // detect the case where there's no chiral volume for the default
      // evaluation
      if (fabs(vol) < zeroTol) {
        // and check the other evaluation:
        if (fabs(vol2) < zeroTol) {
          BOOST_LOG(rdWarningLog)
              << "Warning: ambiguous stereochemistry - no chiral volume - at atom "
              << atom->getIdx() << " ignored" << std::endl;
          return std::nullopt;
        }
        vol = vol2;
        prefactor *= -1;
      } else if (vol * vol2 > 0 && fabs(vol2) > volumeTolerance &&
                 dotp1 < dotp2) {
        // both volumes give the same answer, but in the second case the cross
        // product is between two bonds with a better dot product
        vol = vol2;
        prefactor *= -1;
      } else if (fabs(vol) < volumeTolerance && fabs(vol2) > volumeTolerance) {
        // if the first volume is too small, but the second isn't, take the
        // second
        if (vol * vol2 < 0) {
          prefactor *= -1;
        }
        vol = vol2;
      }
    }
    vol *= prefactor;
    // std::cerr << " final " << vol << std::endl;

    // at this point we can assign our atomic stereo based on the sign of the
    // chiral volume
    if (vol > volumeTolerance) {
      res = Atom::ChiralType::CHI_TETRAHEDRAL_CCW;
    } else if (vol < -volumeTolerance) {
      res = Atom::ChiralType::CHI_TETRAHEDRAL_CW;
    } else {
      BOOST_LOG(rdWarningLog)
          << "Warning: ambiguous stereochemistry - zero final chiral volume - at atom "
          << atom->getIdx() << " ignored" << std::endl;
      return std::nullopt;
    }
  }

  return res;
}

#ifdef _WIN32
int setenv(const char *name, const char *value, int) {
  return _putenv_s(name, value);
}
#endif

void setAllowNontetrahedralChirality(bool val) {
  if (val) {
    setenv(nonTetrahedralStereoEnvVar, "1", 1);
  } else {
    setenv(nonTetrahedralStereoEnvVar, "0", 1);
  }
}
bool getAllowNontetrahedralChirality() {
  return getValFromEnvironment(nonTetrahedralStereoEnvVar,
                               nonTetrahedralStereoDefaultVal);
}

void setUseLegacyStereoPerception(bool val) {
  if (val) {
    setenv(useLegacyStereoEnvVar, "1", 1);
  } else {
    setenv(useLegacyStereoEnvVar, "0", 1);
  }
}
bool getUseLegacyStereoPerception() {
  return getValFromEnvironment(useLegacyStereoEnvVar,
                               useLegacyStereoDefaultVal);
}

namespace detail {
bool bondAffectsAtomChirality(const Bond *bond, const Atom *atom) {
  // FIX consider how to handle organometallics
  PRECONDITION(bond, "bad bond pointer");
  PRECONDITION(atom, "bad atom pointer");
  if (bond->getBondType() == Bond::BondType::UNSPECIFIED ||
      bond->getBondType() == Bond::BondType::ZERO ||
      (bond->getBondType() == Bond::BondType::DATIVE &&
       bond->getBeginAtomIdx() == atom->getIdx())) {
    return false;
  }
  return true;
}
unsigned int getAtomNonzeroDegree(const Atom *atom) {
  PRECONDITION(atom, "bad pointer");
  PRECONDITION(atom->hasOwningMol(), "no owning molecule");
  unsigned int res = 0;
  for (auto bond : atom->getOwningMol().atomBonds(atom)) {
    if (!bondAffectsAtomChirality(bond, atom)) {
      continue;
    }
    ++res;
  }
  return res;
}

bool has_protium_neighbor(const ROMol &mol, const Atom *atom) {
  for (const auto nbr : mol.atomNeighbors(atom)) {
    if (is_regular_h(*nbr)) {
      return true;
    }
  }
  return false;
}

void setStereoForBond(ROMol &mol, Bond *bond, Bond::BondStereo stereo) {
  // NOTE:  moved from parse_doublebond_stereo CXSmilesOps
  // the cis/trans/unknown marker is relative to the lowest numbered atom
  // connected to the lowest numbered double bond atom and the
  // highest-numbered atom connected to the highest-numbered double bond
  // atom find those
  auto begAtom = bond->getBeginAtom();
  auto endAtom = bond->getEndAtom();
  if (begAtom->getIdx() > endAtom->getIdx()) {
    std::swap(begAtom, endAtom);
  }
  if (begAtom->getDegree() > 1 && endAtom->getDegree() > 1) {
    unsigned int begControl = mol.getNumAtoms();
    for (auto nbr : mol.atomNeighbors(begAtom)) {
      if (nbr == endAtom) {
        continue;
      }
      begControl = std::min(nbr->getIdx(), begControl);
    }
    unsigned int endControl = 0;
    for (auto nbr : mol.atomNeighbors(endAtom)) {
      if (nbr == begAtom) {
        continue;
      }
      endControl = std::max(nbr->getIdx(), endControl);
    }
    if (begAtom != bond->getBeginAtom()) {
      std::swap(begControl, endControl);
    }
    bond->setStereoAtoms(begControl, endControl);
    bond->setStereo(stereo);
    mol.setProp("_needsDetectBondStereo", 1);
  }
}
}  // namespace detail

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
  for (const auto atom : mol.atoms()) {
    const unsigned short nMassBits = 10;
    const unsigned short maxMass = 1 << nMassBits;
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

void iterateCIPRanks(const ROMol &mol, const DOUBLE_VECT &invars,
                     UINT_VECT &ranks, bool seedWithInvars) {
  PRECONDITION(invars.size() == mol.getNumAtoms(), "bad invars size");
  PRECONDITION(ranks.size() >= mol.getNumAtoms(), "bad ranks size");

  unsigned int numAtoms = mol.getNumAtoms();
  CIP_ENTRY_VECT cipEntries(numAtoms);
  for (auto &vec : cipEntries) {
    vec.reserve(16);
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
    if (seedWithInvars) {
      cipEntries[i].push_back(static_cast<int>(invars[i]));
    } else {
      cipEntries[i].push_back(mol[i]->getAtomicNum());
      cipEntries[i].push_back(static_cast<int>(ranks[i]));
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
  std::vector<unsigned int> counts(ranks.size());
  std::vector<unsigned int> updatedNbrIdxs;
  updatedNbrIdxs.reserve(8);
  while (numRanks < numAtoms && numIts < maxIts &&
         (lastNumRanks < 0 ||
          static_cast<unsigned int>(lastNumRanks) < numRanks)) {
    unsigned int longestEntry = 0;
    // ----------------------------------------------------
    //
    // for each atom, get a sorted list of its neighbors' ranks:
    //
    for (unsigned int index = 0; index < numAtoms; ++index) {
      // Note: counts is cleaned up when we drain into cipEntries.
      updatedNbrIdxs.clear();

      // start by pushing on our neighbors' ranks:
      for (const auto bond : mol.atomBonds(mol[index])) {
        unsigned int nbrIdx = bond->getOtherAtomIdx(index);
        updatedNbrIdxs.push_back(nbrIdx);

        // put the neighbor in 2N times where N is the bond order as a double.
        // this is to treat aromatic linkages on fair footing. i.e. at least in
        // the first iteration --c(:c):c and --C(=C)-C should look the same.
        // this was part of issue 3009911

        // a special case for chiral phosphorus compounds
        // (this was leading to incorrect assignment of R/S labels ):
        bool isChiralPhosphorusSpecialCase = false;
        if (bond->getBondType() == Bond::DOUBLE) {
          const Atom *nbr = mol[nbrIdx];
          if (nbr->getAtomicNum() == 15) {
            unsigned int nbrDeg = nbr->getDegree();
            isChiralPhosphorusSpecialCase = nbrDeg == 3 || nbrDeg == 4;
          }
        };

        // general justification of this is:
        // Paragraph 2.2. in the 1966 article is "Valence-Bond Conventions:
        // Multiple-Bond Unsaturation and Aromaticity". It contains several
        // conventions of which convention (b) is the one applying here:
        // "(b) Contributions by d orbitals to bonds of quadriligant atoms are
        // neglected."
        // FIX: this applies to more than just P
        if (isChiralPhosphorusSpecialCase) {
          counts[nbrIdx] += 1;
        } else {
          counts[nbrIdx] += getTwiceBondType(*bond);
        }
      }

      // For each of our neighbors' ranks weighted by bond type, copy it N times
      // to our cipEntry in reverse rank order, where N is the weight.
      if (updatedNbrIdxs.size() > 1) {  // compare vs 1 for performance.
        std::sort(std::begin(updatedNbrIdxs), std::end(updatedNbrIdxs),
                  [&ranks](unsigned int idx1, unsigned int idx2) {
                    return ranks[idx1] > ranks[idx2];
                  });
      }
      auto &cipEntry = cipEntries[index];
      for (auto nbrIdx : updatedNbrIdxs) {
        unsigned int count = counts[nbrIdx];
        cipEntry.insert(cipEntry.end(), count, ranks[nbrIdx] + 1);
        counts[nbrIdx] = 0;
      }
      // add a zero for each coordinated H as long as we're not a query atom
      if (!mol[index]->hasQuery()) {
        cipEntry.insert(cipEntry.end(), mol[index]->getTotalNumHs(), 0);
      }

      if (cipEntry.size() > longestEntry) {
        longestEntry = rdcast<unsigned int>(cipEntry.size());
      }
    }
    // ----------------------------------------------------
    //
    // pad the entries so that we compare rounds to themselves:
    //
    for (unsigned int index = 0; index < numAtoms; ++index) {
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
    if (static_cast<unsigned int>(lastNumRanks) != numRanks) {
      for (unsigned int i = 0; i < numAtoms; ++i) {
        cipEntries[i][numIts + 1] = ranks[i];
        cipEntries[i].erase(cipEntries[i].begin() + numIts + 2,
                            cipEntries[i].end());
      }
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
  for (const auto bond : mol.atomBonds(atom)) {
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
  for (const auto bond : mol.atomBonds(atom)) {
    if (bond == refBond) {
      continue;
    }
    Bond::BondDir dir = bond->getBondDir();
    if (bond->getBondType() == Bond::SINGLE ||
        (includeAromatic && bond->getBondType() == Bond::AROMATIC)) {
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
// Note that N atoms are only candidates if they are in a 3-ring
bool atomIsCandidateForRingStereochem(const ROMol &mol, const Atom *atom) {
  PRECONDITION(atom, "bad atom");
  bool res = false;
  std::set<unsigned int> nbrRanks;
  if (!atom->getPropIfPresent(common_properties::_ringStereochemCand, res)) {
    const RingInfo *ringInfo = mol.getRingInfo();
    if (ringInfo->isInitialized() && ringInfo->numAtomRings(atom->getIdx())) {
      // three-coordinate N additional requirements:
      //   in a ring of size 3  (from InChI)
      // OR
      //   a bridgehead (RDKit extension)
      if (atom->getAtomicNum() == 7 && atom->getTotalDegree() == 3 &&
          !ringInfo->isAtomInRingOfSize(atom->getIdx(), 3) &&
          !queryIsAtomBridgehead(atom)) {
        return false;
      }
      std::vector<const Atom *> nonRingNbrs;
      std::vector<const Atom *> ringNbrs;
      for (const auto bond : mol.atomBonds(atom)) {
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
            res = rank1 != rank2;
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

  for (const auto atom : mol.atoms()) {
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
    for (const auto bond : mol.atomBonds(atom)) {
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
      for (const auto bond : mol.atomBonds(ratom)) {
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
      for (auto ringAtomEntry : ringStereoAtoms) {
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

  auto nzDegree = Chirality::detail::getAtomNonzeroDegree(atom);
  auto tnzDegree = nzDegree + atom->getTotalNumHs();
  if (tnzDegree > 4) {
    // we only know tetrahedral chirality
    legalCenter = false;
  } else {
    // cases we can exclude immediately without having to look at neighbors
    // ranks:
    if (tnzDegree < 3) {
      legalCenter = false;
    } else if (nzDegree < 3 &&
               (atom->getAtomicNum() != 15 && atom->getAtomicNum() != 33)) {
      // less than three neighbors is never stereogenic
      // unless it is a phosphine/arsine with implicit H (this is from InChI)
      legalCenter = false;
    } else if (nzDegree == 3) {
      if (atom->getTotalNumHs() == 1) {
        // three-coordinate with more than one H is never stereogenic
        if (detail::has_protium_neighbor(mol, atom)) {
          legalCenter = false;
        }
      } else {
        // assume something that's really three coordinate isn't potentially
        // chiral, then look for exceptions
        legalCenter = false;
        if (atom->getAtomicNum() == 7) {
          // three-coordinate N additional requirements:
          //   in a ring of size 3  (from InChI)
          // OR
          /// is a bridgehead atom (RDKit extension)
          if (mol.getRingInfo()->isAtomInRingOfSize(atom->getIdx(), 3) ||
              queryIsAtomBridgehead(atom)) {
            legalCenter = true;
          }
        } else if (atom->getAtomicNum() == 15 || atom->getAtomicNum() == 33) {
          // three-coordinate phosphines and arsines
          // are always treated as stereogenic even with H atom neighbors.
          // (this is from InChI)
          legalCenter = true;
        } else if (atom->getAtomicNum() == 16 || atom->getAtomicNum() == 34) {
          if (atom->getExplicitValence() == 4 ||
              (atom->getExplicitValence() == 3 &&
               atom->getFormalCharge() == 1)) {
            // we also accept sulfur or selenium with either a positive charge
            // or a double bond:
            legalCenter = true;
          }
        }
      }
    }

    if (legalCenter) {
      boost::dynamic_bitset<> codesSeen(mol.getNumAtoms());
      for (const auto bond : mol.atomBonds(atom)) {
        unsigned int otherIdx = bond->getOtherAtom(atom)->getIdx();
        nbrs.push_back(std::make_pair(ranks[otherIdx], bond->getIdx()));
        if (!Chirality::detail::bondAffectsAtomChirality(bond, atom)) {
          continue;
        }
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
  for (auto atom : mol.atoms()) {
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
      // note that hasDupes is only evaluated if legalCenter==true
      auto [legalCenter, hasDupes] =
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
        std::sort(nbrs.begin(), nbrs.end(), Rankers::pairLess);

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
  for (auto dblBond : mol.bonds()) {
    if (dblBond->getBondType() == Bond::BondType::DOUBLE) {
      if (dblBond->getStereo() != Bond::BondStereo::STEREONONE) {
        continue;
      }
      if (!ranks.size()) {
        assignAtomCIPRanks(mol, ranks);
      }
      dblBond->getStereoAtoms().clear();

      // at the moment we are ignoring stereochem on ring bonds with less than
      // 8 members.
      if (shouldDetectDoubleBondStereo(dblBond)) {
        const Atom *begAtom = dblBond->getBeginAtom();
        const Atom *endAtom = dblBond->getEndAtom();
        // we're only going to handle 2 or three coordinate atoms:
        if ((begAtom->getDegree() == 2 || begAtom->getDegree() == 3) &&
            (endAtom->getDegree() == 2 || endAtom->getDegree() == 3)) {
          ++unassignedBonds;

          // look around each atom and see if it has at least one bond with
          // direction marked:

          // the pairs here are: atomIdx,bonddir
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

void assignLegacyCIPLabels(ROMol &mol, bool flagPossibleStereoCenters) {
  std::vector<unsigned int> atomRanks;
  assignAtomChiralCodes(mol, atomRanks, flagPossibleStereoCenters);

  // reset any already-specfied double bonds:
  for (auto bond : mol.bonds()) {
    if (bond->getBondType() == Bond::BondType::DOUBLE &&
        bond->getStereo() > Bond::BondStereo::STEREOANY) {
      bond->setStereo(Bond::BondStereo::STEREONONE);
    }
  }
  assignBondStereoCodes(mol, atomRanks);
}

void assignBondCisTrans(ROMol &mol, const StereoInfo &sinfo) {
  if (sinfo.type != StereoType::Bond_Double ||
      sinfo.specified != StereoSpecified::Unspecified ||
      sinfo.controllingAtoms.size() != 4 ||
      ((sinfo.controllingAtoms[0] == StereoInfo::NOATOM &&
        sinfo.controllingAtoms[1] == StereoInfo::NOATOM) ||
       (sinfo.controllingAtoms[2] == StereoInfo::NOATOM &&
        sinfo.controllingAtoms[3] == StereoInfo::NOATOM))) {
    return;
  }

  auto dblBond = mol.getBondWithIdx(sinfo.centeredOn);

  bool begFirstNeighbor = true;
  auto begBond = mol.getBondBetweenAtoms(dblBond->getBeginAtomIdx(),
                                         sinfo.controllingAtoms[0]);
  CHECK_INVARIANT(begBond, "no initial bond found");
  auto begDir = begBond->getBondDir();
  if (begDir != Bond::BondDir::ENDDOWNRIGHT &&
      begDir != Bond::BondDir::ENDUPRIGHT) {
    begFirstNeighbor = false;
    if (sinfo.controllingAtoms[1] != StereoInfo::NOATOM) {
      begBond = mol.getBondBetweenAtoms(dblBond->getBeginAtomIdx(),
                                        sinfo.controllingAtoms[1]);
      CHECK_INVARIANT(begBond, "no initial bond found");
      begDir = begBond->getBondDir();
    }
  }
  // no direction found at beginning
  if (begDir != Bond::BondDir::ENDDOWNRIGHT &&
      begDir != Bond::BondDir::ENDUPRIGHT) {
    return;
  }
  if (begBond->getBeginAtomIdx() != dblBond->getBeginAtomIdx()) {
    begDir = begDir == Bond::BondDir::ENDDOWNRIGHT
                 ? Bond::BondDir::ENDUPRIGHT
                 : Bond::BondDir::ENDDOWNRIGHT;
  }

  bool endFirstNeighbor = true;
  auto endBond = mol.getBondBetweenAtoms(dblBond->getEndAtomIdx(),
                                         sinfo.controllingAtoms[2]);
  CHECK_INVARIANT(endBond, "no final bond found");
  auto endDir = endBond->getBondDir();
  if (endDir != Bond::BondDir::ENDDOWNRIGHT &&
      endDir != Bond::BondDir::ENDUPRIGHT) {
    endFirstNeighbor = false;
    if (sinfo.controllingAtoms[3] != StereoInfo::NOATOM) {
      endBond = mol.getBondBetweenAtoms(dblBond->getEndAtomIdx(),
                                        sinfo.controllingAtoms[3]);
      CHECK_INVARIANT(endBond, "no final bond found");
      endDir = endBond->getBondDir();
    }
  }
  // no direction found at end
  if (endDir != Bond::BondDir::ENDDOWNRIGHT &&
      endDir != Bond::BondDir::ENDUPRIGHT) {
    return;
  }
  if (endBond->getBeginAtomIdx() != dblBond->getEndAtomIdx()) {
    endDir = endDir == Bond::BondDir::ENDDOWNRIGHT
                 ? Bond::BondDir::ENDUPRIGHT
                 : Bond::BondDir::ENDDOWNRIGHT;
  }

  // we've set up the bond directions here so that they correspond to having
  // both single bonds START at the double bond. This means that if the single
  // bonds point in the same direction, the bond is cis
  bool sameDir = begDir == endDir;

  // if either the direction bond at the beginning or the direction bond at the
  // end wasn't to the first neighbor on that side (but not both), then we need
  // to swap
  if (begFirstNeighbor ^ endFirstNeighbor) {
    sameDir = !sameDir;
  }

  dblBond->setStereoAtoms(sinfo.controllingAtoms[0], sinfo.controllingAtoms[2]);
  if (sameDir) {
    dblBond->setStereo(Bond::BondStereo::STEREOCIS);
  } else {
    dblBond->setStereo(Bond::BondStereo::STEREOTRANS);
  }
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
    for (const auto oBond : mol.atomBonds(atom)) {
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
  for (const auto &bondIdx :
       boost::make_iterator_range(mol.getAtomBonds(atom))) {
    const Bond *bond = mol[bondIdx];

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

void cleanupStereoGroups(ROMol &mol) {
  std::vector<StereoGroup> newsgs;
  for (auto sg : mol.getStereoGroups()) {
    std::vector<Atom *> okatoms;
    std::vector<Bond *> okbonds;
    bool keep = true;
    for (const auto atom : sg.getAtoms()) {
      if (atom->getChiralTag() == Atom::ChiralType::CHI_UNSPECIFIED) {
        keep = false;
      } else {
        okatoms.push_back(atom);
      }
    }
    for (const auto bond : sg.getBonds()) {
      if (bond->getStereo() != Bond::BondStereo::STEREOATROPCCW &&
          bond->getStereo() != Bond::BondStereo::STEREOATROPCW) {
        keep = false;
      } else {
        okbonds.push_back(bond);
      }
    }

    if (keep) {
      newsgs.push_back(sg);
    } else if (!okatoms.empty()) {
      newsgs.emplace_back(sg.getGroupType(), std::move(okatoms),
                          std::move(okbonds), sg.getReadId());
    }
  }
  mol.setStereoGroups(std::move(newsgs));
}

// ****************************************************************************
std::ostream &operator<<(std::ostream &oss, const StereoType &s) {
  switch (s) {
    case StereoType::Unspecified:
      oss << "Unspecified";
      break;
    case StereoType::Atom_Tetrahedral:
      oss << "Atom_Tetrahedral";
      break;
    case StereoType::Atom_SquarePlanar:
      oss << "Atom_SquarePlanar";
      break;
    case StereoType::Atom_TrigonalBipyramidal:
      oss << "Atom_TrigonalBipyramidal";
      break;
    case StereoType::Atom_Octahedral:
      oss << "Atom_Octahedral";
      break;
    case StereoType::Bond_Double:
      oss << "Bond_Double";
      break;
    case StereoType::Bond_Cumulene_Even:
      oss << "Bond_Cumulene_Even";
      break;
    case StereoType::Bond_Atropisomer:
      oss << "Bond_Atropisomer";
      break;
  }
  return oss;
}

// ****************************************************************************
std::ostream &operator<<(std::ostream &oss, const StereoSpecified &s) {
  switch (s) {
    case StereoSpecified::Unspecified:
      oss << "Unspecified";
      break;
    case StereoSpecified::Specified:
      oss << "Specified";
      break;
    case StereoSpecified::Unknown:
      oss << "Unknown";
      break;
  }
  return oss;
}

/*
    We're going to do this iteratively:
      1) assign atom stereochemistry
      2) assign bond stereochemistry
      3) if there are still unresolved atoms or bonds
         repeat the above steps as necessary
 */
void legacyStereoPerception(ROMol &mol, bool cleanIt,
                            bool flagPossibleStereoCenters) {
  mol.clearProp("_needsDetectBondStereo");

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
  for (auto atom : mol.atoms()) {
    if (cleanIt) {
      if (atom->hasProp(common_properties::_CIPCode)) {
        atom->clearProp(common_properties::_CIPCode);
      }
      if (atom->hasProp(common_properties::_ChiralityPossible)) {
        atom->clearProp(common_properties::_ChiralityPossible);
      }
    }
    if (!hasStereoAtoms && atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
        atom->getChiralTag() != Atom::CHI_OTHER) {
      hasStereoAtoms = true;
    }
  }
  bool hasStereoBonds = false;
  for (auto bond : mol.bonds()) {
    if (cleanIt) {
      // enforce no stereo on small rings
      if ((bond->getBondType() == Bond::DOUBLE ||
           bond->getBondType() == Bond::AROMATIC) &&
          !shouldDetectDoubleBondStereo(bond)) {
        if (bond->getBondDir() == Bond::EITHERDOUBLE) {
          bond->setBondDir(Bond::NONE);
        }
        if (bond->getStereo() != Bond::STEREONONE) {
          bond->setStereo(Bond::STEREONONE);
          bond->getStereoAtoms().clear();
        }
        continue;
      } else if (bond->getBondType() == Bond::DOUBLE) {
        if (bond->getBondDir() == Bond::EITHERDOUBLE) {
          bond->setStereo(Bond::STEREOANY);
          bond->getStereoAtoms().clear();
          bond->setBondDir(Bond::NONE);
        } else if (bond->getStereo() != Bond::STEREOANY) {
          bond->setStereo(Bond::STEREONONE);
          bond->getStereoAtoms().clear();
        }
      }
    }
    if (!hasStereoBonds && bond->getBondType() == Bond::DOUBLE) {
      for (auto nbond : mol.atomBonds(bond->getBeginAtom())) {
        if (nbond->getBondDir() == Bond::ENDDOWNRIGHT ||
            nbond->getBondDir() == Bond::ENDUPRIGHT) {
          hasStereoBonds = true;
          break;
        }
      }
      if (!hasStereoBonds) {
        for (auto nbond : mol.atomBonds(bond->getEndAtom())) {
          if (nbond->getBondDir() == Bond::ENDDOWNRIGHT ||
              nbond->getBondDir() == Bond::ENDUPRIGHT) {
            hasStereoBonds = true;
            break;
          }
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
      std::tie(hasStereoAtoms, changedStereoAtoms) =
          Chirality::assignAtomChiralCodes(mol, atomRanks,
                                           flagPossibleStereoCenters);
    } else {
      changedStereoAtoms = false;
    }
    if (hasStereoBonds) {
      std::tie(hasStereoBonds, changedStereoBonds) =
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

    for (auto atom : mol.atoms()) {
      if (atom->hasProp(common_properties::_ringStereochemCand)) {
        atom->clearProp(common_properties::_ringStereochemCand);
      }
      if (atom->hasProp(common_properties::_ringStereoAtoms)) {
        atom->clearProp(common_properties::_ringStereoAtoms);
      }
    }
    boost::dynamic_bitset<> possibleSpecialCases(mol.getNumAtoms());
    Chirality::findChiralAtomSpecialCases(mol, possibleSpecialCases);

    for (auto atom : mol.atoms()) {
      if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
          !Chirality::hasNonTetrahedralStereo(atom) &&
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
          for (const auto nbrBndI : mol.atomBonds(batom)) {
            if (nbrBndI == bond) {
              continue;
            }
            if ((nbrBndI->getBondDir() == Bond::ENDDOWNRIGHT ||
                 nbrBndI->getBondDir() == Bond::ENDUPRIGHT) &&
                (nbrBndI->getBondType() == Bond::SINGLE ||
                 nbrBndI->getBondType() == Bond::AROMATIC)) {
              // direction is set, and we know it's not because of our
              // bond. What about other neighbors?
              bool okToClear = true;
              for (const auto nbrBndJ :
                   mol.atomBonds(nbrBndI->getOtherAtom(batom))) {
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
    }
    Chirality::cleanupStereoGroups(mol);
  }
}

void updateDoubleBondStereo(ROMol &mol, const std::vector<StereoInfo> &sinfo,
                            bool cleanIt) {
  boost::dynamic_bitset<> bondsTouched(mol.getNumBonds(), 0);
  for (const auto &si : sinfo) {
    if (si.type == Chirality::StereoType::Bond_Double) {
      auto bond = mol.getBondWithIdx(si.centeredOn);
      bondsTouched.set(bond->getIdx());
      bond->setStereo(Bond::BondStereo::STEREONONE);
      if (si.specified == Chirality::StereoSpecified::Specified) {
        TEST_ASSERT(si.controllingAtoms.size() == 4);
        bond->setStereoAtoms(si.controllingAtoms[0], si.controllingAtoms[2]);
        switch (si.descriptor) {
          case Chirality::StereoDescriptor::Bond_Cis:
            bond->setStereo(Bond::BondStereo::STEREOCIS);
            break;
          case Chirality::StereoDescriptor::Bond_Trans:
            bond->setStereo(Bond::BondStereo::STEREOTRANS);
            break;
          default:
            BOOST_LOG(rdWarningLog)
                << "unrecognized bond stereo type" << std::endl;
        }
      } else if (si.specified == Chirality::StereoSpecified::Unknown) {
        bond->setStereo(Bond::BondStereo::STEREOANY);
        bond->getStereoAtoms().clear();
        bond->setBondDir(Bond::BondDir::NONE);
      } else if (si.specified == Chirality::StereoSpecified::Unspecified) {
        assignBondCisTrans(mol, si);
      }
    }
  }
  if (cleanIt) {
    for (auto bond : mol.bonds()) {
      if (bondsTouched[bond->getIdx()] ||
          bond->getBondType() != Bond::BondType::DOUBLE) {
        continue;
      }
      // we didn't see it above, so it can't have stereo:
      bond->setStereo(Bond::BondStereo::STEREONONE);
      bond->setBondDir(Bond::BondDir::NONE);
      bond->getStereoAtoms().clear();
    }
  }
}
void stereoPerception(ROMol &mol, bool cleanIt,
                      bool flagPossibleStereoCenters) {
  if (cleanIt) {
    for (auto atom : mol.atoms()) {
      atom->clearProp(common_properties::_CIPCode);
      atom->clearProp(common_properties::_ChiralityPossible);
    }
    for (auto bond : mol.bonds()) {
      if (bond->getBondDir() == Bond::BondDir::EITHERDOUBLE) {
        bond->setStereo(Bond::BondStereo::STEREOANY);
        bond->getStereoAtoms().clear();
        bond->setBondDir(Bond::BondDir::NONE);
      }
    }
  }
  // we need cis/trans markers on the double bonds... set those now:
  MolOps::setBondStereoFromDirections(mol);

  // do the actual perception
  auto sinfo =
      Chirality::findPotentialStereo(mol, cleanIt, flagPossibleStereoCenters);

  if (flagPossibleStereoCenters) {
    for (const auto &si : sinfo) {
      if (si.type == Chirality::StereoType::Atom_Tetrahedral ||
          si.type == Chirality::StereoType::Atom_SquarePlanar ||
          si.type == Chirality::StereoType::Atom_TrigonalBipyramidal ||
          si.type == Chirality::StereoType::Atom_Octahedral) {
        mol.getAtomWithIdx(si.centeredOn)
            ->setProp(common_properties::_ChiralityPossible, 1);
      }
    }
  }
  // populate double bond stereo info:
  updateDoubleBondStereo(mol, sinfo, cleanIt);
  if (cleanIt) {
    Chirality::cleanupStereoGroups(mol);
  }
}

bool canBeStereoBond(const Bond *bond) {
  PRECONDITION(bond, "no bond");
  if (bond->getBondType() != Bond::BondType::DOUBLE &&
      bond->getBondType() != Bond::BondType::AROMATIC) {
    return false;
  }
  auto beginAtom = bond->getBeginAtom();
  auto endAtom = bond->getEndAtom();
  for (const auto atom : {beginAtom, endAtom}) {
    std::vector<int> nbrRanks;
    for (auto nbrBond : bond->getOwningMol().atomBonds(atom)) {
      if (nbrBond == bond) {
        continue;  // a bond is NOT its own neighbor
      }

      if (nbrBond->getBondType() == Bond::SINGLE) {
        // if a neighbor has a wedge or hash bond, do NOT mark it as double
        // crossed
        if (nbrBond->getBondDir() == Bond::ENDUPRIGHT ||
            nbrBond->getBondDir() == Bond::ENDDOWNRIGHT) {
          return false;
        }

        // if a neighbor has a wiggle bond, do NOT mark it as crossed (although
        // it is unknown
        if (nbrBond->getBondDir() == Bond::BondDir::UNKNOWN &&
            nbrBond->getBeginAtom() == atom) {
          return false;
        }

        // if two neighbors havr the same CIP ranking, this is not stereo
        const auto otherAtom = nbrBond->getOtherAtom(atom);
        int rank;
        if (RDKit::Chirality::getUseLegacyStereoPerception()) {
          if (!otherAtom->getPropIfPresent(common_properties::_CIPRank, rank)) {
            rank = -1;
          }
        } else {  // NOT legacy stereo
          if (!otherAtom->getPropIfPresent(common_properties::_ChiralAtomRank,
                                           rank)) {
            rank = -1;
          }
        }
        if (rank >= 0) {
          if (std::find(nbrRanks.begin(), nbrRanks.end(), rank) !=
              nbrRanks.end()) {
            return false;
          } else {
            nbrRanks.push_back(rank);
          }
        }
      }
    }
  }
  return true;
}

bool shouldBeACrossedBond(const Bond *bond) {
  PRECONDITION(bond, "");

  // double bond stereochemistry -
  // if the bond isn't specified, then it should go in the mol block
  // as "any", this was sf.net issue 2963522.
  // two caveats to this:
  // 1) if it's a ring bond, we'll only put the "any"
  //    in the mol block if the user specifically asked for it.
  //    Constantly seeing crossed bonds in rings, though maybe
  //    technically correct, is irritating.
  // 2) if it's a terminal bond (where there's no chance of
  //    stereochemistry anyway), we also skip the any.
  //    this was sf.net issue 3009756

  if (bond->getStereo() == Bond::STEREOANY) {
    // see if any of the neighbors have a wiggle bond - if so, do NOT make this
    // one a cross bond
    for (auto nbrBond : bond->getOwningMol().atomBonds(bond->getBeginAtom())) {
      if (nbrBond->getBondDir() == Bond::UNKNOWN &&
          nbrBond->getBeginAtom()->getIdx() == bond->getBeginAtom()->getIdx()) {
        return false;
      }
    }
    for (auto nbrBond : bond->getOwningMol().atomBonds(bond->getEndAtom())) {
      if (nbrBond->getBondDir() == Bond::UNKNOWN &&
          nbrBond->getBeginAtom()->getIdx() == bond->getEndAtom()->getIdx()) {
        return false;
      }
    }

    return true;  // crossed double bond
  }
  if (bond->getStereo() != Bond::BondStereo::STEREONONE) {
    return false;
  }

  // if it is in a ring it is not makred as stereo.
  // If either end is terminal, it is not stereo

  if (!Chirality::detail::isBondPotentialStereoBond(bond)) {
    return false;
  }
  // we don't know that it's explicitly unspecified (covered above with
  // the ==STEREOANY check)

  if (bond->getBondDir() == Bond::EITHERDOUBLE) {
    return true;  // crossed double bond
  }

  const auto beginAtom = bond->getBeginAtom();
  const auto endAtom = bond->getEndAtom();
  if (beginAtom->getDegree() > 1 && endAtom->getDegree() > 1 &&
      (beginAtom->getTotalValence() - beginAtom->getTotalDegree()) == 1 &&
      (endAtom->getTotalValence() - endAtom->getTotalDegree()) == 1) {
    // we only do this if each atom only has one unsaturation
    // FIX: this is the fix for github #2649, but we will need to
    // change it once we start handling allenes properly

    if (canBeStereoBond(bond)) {
      return true;  // crossed double bond
    }
  }

  return false;  // NOT crossed double bond
}

// only valid for single or aromatic  bonds
int BondGetDirCode(const Bond::BondDir dir) {
  int res = 0;
  switch (dir) {
    case Bond::NONE:
      res = 0;
      break;
    case Bond::BEGINWEDGE:
      res = 1;
      break;
    case Bond::BEGINDASH:
      res = 6;
      break;
    case Bond::UNKNOWN:
      res = 4;
      break;
    case Bond::BondDir::EITHERDOUBLE:
      res = 3;
      break;
    default:
      break;
  }
  return res;
}

void GetMolFileBondStereoInfo(
    const Bond *bond,
    const std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
        &wedgeBonds,
    const Conformer *conf, Bond::BondDir &dir, bool &reverse) {
  PRECONDITION(bond, "");
  reverse = false;
  dir = Bond::NONE;
  if (canHaveDirection(*bond)) {
    // single bond stereo chemistry

    dir = Chirality::detail::determineBondWedgeState(bond, wedgeBonds, conf);

    // if this bond needs to be wedged it is possible that this
    // wedging was determined by a chiral atom at the end of the
    // bond (instead of at the beginning). In this case we need to
    // reverse the begin and end atoms for the bond when we write
    // the mol file
    if ((dir == Bond::BEGINDASH) ||
        (dir == Bond::BEGINWEDGE || dir == Bond::UNKNOWN)) {
      auto wbi = wedgeBonds.find(bond->getIdx());
      if (wbi != wedgeBonds.end() &&
          wbi->second->getType() ==
              Chirality::WedgeInfoType::WedgeInfoTypeChiral &&
          static_cast<unsigned int>(wbi->second->getIdx()) !=
              bond->getBeginAtomIdx()) {
        reverse = true;
      }
    }
  } else if (bond->getBondType() == Bond::DOUBLE) {
    if (Chirality::shouldBeACrossedBond(bond)) {
      dir = Bond::BondDir::EITHERDOUBLE;
    }
  }
}

void GetMolFileBondStereoInfo(
    const Bond *bond,
    const std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
        &wedgeBonds,
    const Conformer *conf, int &dirCode, bool &reverse) {
  Bond::BondDir dir;
  GetMolFileBondStereoInfo(bond, wedgeBonds, conf, dir, reverse);
  dirCode = BondGetDirCode(dir);
}

void removeNonExplicit3DChirality(ROMol &mol) {
  for (auto atom : mol.atoms()) {
    if (atom->hasProp(common_properties::_NonExplicit3DChirality)) {
      atom->clearProp(common_properties::_NonExplicit3DChirality);
      atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    }
  }
}

void addStereoAnnotations(ROMol &mol, std::string absLabel, std::string orLabel,
                          std::string andLabel, std::string cipLabel,
                          std::string bondLabel) {
  auto sgs = mol.getStereoGroups();
  assignStereoGroupIds(sgs);
  boost::dynamic_bitset<> doneAts(mol.getNumAtoms());
  for (const auto &sg : sgs) {
    std::string gid = std::to_string(sg.getWriteId());
    for (const auto atom : sg.getAtoms()) {
      if (doneAts[atom->getIdx()]) {
        BOOST_LOG(rdWarningLog) << "Warning: atom " << atom->getIdx()
                                << " is in more than one stereogroup. Only the "
                                   "label from the first group will be used."
                                << std::endl;
        continue;
      }
      std::string cip;
      atom->getPropIfPresent(common_properties::_CIPCode, cip);

      std::string lab;
      switch (sg.getGroupType()) {
        case StereoGroupType::STEREO_ABSOLUTE:
          lab = absLabel;
          doneAts.set(atom->getIdx());
          break;
        case StereoGroupType::STEREO_OR:
          lab = orLabel;
          doneAts.set(atom->getIdx());
          break;
        case StereoGroupType::STEREO_AND:
          lab = andLabel;
          doneAts.set(atom->getIdx());
          break;
        default:
          break;
      }

      if (!lab.empty()) {
        boost::algorithm::replace_all(lab, "{id}", gid);
        if (!cip.empty()) {
          boost::algorithm::replace_all(lab, "{cip}", cip);
        }
        atom->setProp(common_properties::atomNote, lab);
      }
    }
  }
  if (!cipLabel.empty()) {
    for (auto atom : mol.atoms()) {
      std::string cip;
      if (!doneAts[atom->getIdx()] &&
          atom->getPropIfPresent(common_properties::_CIPCode, cip)) {
        std::string lab = cipLabel;
        boost::algorithm::replace_all(lab, "{cip}", cip);
        atom->setProp(common_properties::atomNote, lab);
      }
    }
  }
  if (!bondLabel.empty()) {
    for (auto bond : mol.bonds()) {
      std::string cip;
      if (!bond->getPropIfPresent(common_properties::_CIPCode, cip)) {
        if (bond->getStereo() == Bond::STEREOE) {
          cip = "E";
        } else if (bond->getStereo() == Bond::STEREOZ) {
          cip = "Z";
        }
      }
      if (!cip.empty()) {
        std::string lab = bondLabel;
        boost::algorithm::replace_all(lab, "{cip}", cip);
        bond->setProp(common_properties::bondNote, lab);
      }
    }
  }
}

}  // namespace Chirality

namespace MolOps {

void assignStereochemistry(ROMol &mol, bool cleanIt, bool force,
                           bool flagPossibleStereoCenters) {
  if (!force && mol.hasProp(common_properties::_StereochemDone)) {
    return;
  }

  if (!Chirality::getUseLegacyStereoPerception()) {
    Chirality::stereoPerception(mol, cleanIt, flagPossibleStereoCenters);
  } else {
    Chirality::legacyStereoPerception(mol, cleanIt, flagPossibleStereoCenters);
  }
  mol.setProp(common_properties::_StereochemDone, 1, true);
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
        // proceed only if we either want to clean the stereocode on this bond,
        // if none is set on it yet, or it is STEREOANY and we need to find
        // stereoatoms
        if (cleanIt || dblBond->getStereo() == Bond::STEREONONE ||
            (dblBond->getStereo() == Bond::STEREOANY &&
             dblBond->getStereoAtoms().size() != 2)) {
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
  unsigned int degree, perm;
  for (auto atom : mol.atoms()) {
    switch (atom->getChiralTag()) {
      case Atom::CHI_TETRAHEDRAL_CW:
      case Atom::CHI_TETRAHEDRAL_CCW:
        if (atom->getHybridization() != Atom::SP3) {
          atom->setChiralTag(Atom::CHI_UNSPECIFIED);
        }
        break;

      case Atom::CHI_TETRAHEDRAL:
        if (atom->getHybridization() != Atom::SP3) {
          atom->setChiralTag(Atom::CHI_UNSPECIFIED);
        } else {
          perm = 0;
          atom->getPropIfPresent(common_properties::_chiralPermutation, perm);
          if (perm > 2) {
            perm = 0;
            atom->setProp(common_properties::_chiralPermutation, perm);
          }
        }
        break;

      case Atom::CHI_SQUAREPLANAR:
        degree = atom->getTotalDegree();
        if (degree < 2 || degree > 4) {
          atom->setChiralTag(Atom::CHI_UNSPECIFIED);
        } else {
          perm = 0;
          atom->getPropIfPresent(common_properties::_chiralPermutation, perm);
          if (perm > 3) {
            perm = 0;
            atom->setProp(common_properties::_chiralPermutation, perm);
          }
        }
        break;

      case Atom::CHI_TRIGONALBIPYRAMIDAL:
        degree = atom->getTotalDegree();
        if (degree < 2 || degree > 5) {
          atom->setChiralTag(Atom::CHI_UNSPECIFIED);
        } else {
          perm = 0;
          atom->getPropIfPresent(common_properties::_chiralPermutation, perm);
          if (perm > 20) {
            perm = 0;
            atom->setProp(common_properties::_chiralPermutation, perm);
          }
        }
        break;

      case Atom::CHI_OCTAHEDRAL:
        degree = atom->getTotalDegree();
        if (degree < 2 || degree > 6) {
          atom->setChiralTag(Atom::CHI_UNSPECIFIED);
        } else {
          perm = 0;
          atom->getPropIfPresent(common_properties::_chiralPermutation, perm);
          if (perm > 30) {
            perm = 0;
            atom->setProp(common_properties::_chiralPermutation, perm);
          }
        }
        break;

      default:
        /* ??? Handle other types in future.  */
        break;
    }
  }
}

#define VOLTEST(X, Y, Z) (v[X].dotProduct(v[Y].crossProduct(v[Z])) >= 0.0)

static unsigned int OctahedralPermFrom3D(unsigned char *pair,
                                         const RDGeom::Point3D *v) {
  switch (pair[0]) {
    case 2:  // a-b
      switch (pair[2]) {
        case 4:
          return VOLTEST(0, 3, 4) ? 28 : 27;
        case 5:
          return VOLTEST(0, 2, 3) ? 25 : 30;
        default:  // 0 or 6
          return VOLTEST(0, 2, 3) ? 26 : 29;
      }
      break;
    case 3:  // a-c
      switch (pair[1]) {
        case 4:
          return VOLTEST(0, 3, 4) ? 22 : 21;
        case 5:
          return VOLTEST(0, 1, 3) ? 19 : 24;
        default:  // 0 or 6
          return VOLTEST(0, 1, 3) ? 20 : 23;
      }
      break;
    case 4:  // a-d
      switch (pair[1]) {
        case 3:
          return VOLTEST(0, 2, 4) ? 13 : 12;
        case 5:
          return VOLTEST(0, 1, 2) ? 6 : 18;
        default:  // 0 or 6
          return VOLTEST(0, 1, 2) ? 7 : 17;
      }
      break;
    case 5:  // a-e
      switch (pair[1]) {
        case 3:
          return VOLTEST(0, 2, 3) ? 11 : 9;
        case 4:
          return VOLTEST(0, 1, 2) ? 3 : 16;
        default:  // 0 or 6
          return VOLTEST(0, 1, 2) ? 5 : 15;
      }
      break;
    default:  // 0 or 6  a-f
      switch (pair[1]) {
        case 3:
          return VOLTEST(0, 2, 3) ? 10 : 8;
        case 4:
          return VOLTEST(0, 1, 2) ? 1 : 2;
        default:  // 5
          return VOLTEST(0, 1, 2) ? 4 : 14;
      }
  }
  // unreachable
  return 0;
}

bool isWigglyBond(const Bond *bond, const Atom *atom) {
  int hasWigglyBond = 0;
  if (bond->getBeginAtomIdx() == atom->getIdx() &&
      bond->getBondType() == Bond::BondType::SINGLE &&
      (bond->getBondDir() == Bond::BondDir::UNKNOWN ||
       (bond->getPropIfPresent<int>(common_properties::_UnknownStereo,
                                    hasWigglyBond) &&
        hasWigglyBond))) {
    return true;
  }
  return false;
}
// The tolerance here is pretty high in order to accomodate things coming from
// the dgeom code As we get more experience with real-world structures and/or
// improve the dgeom code, we can think about lowering this.
static bool assignNontetrahedralChiralTypeFrom3D(ROMol &mol,
                                                 const Conformer &conf,
                                                 Atom *atom,
                                                 double tolerance = 0.1) {
  // FIX: add tests for dative and zero order bonds
  // Fail fast check for non-tetrahedral elements
  if (atom->getAtomicNum() < 15) {
    return false;
  }

  // check for wiggly bonds
  for (const auto bond : mol.atomBonds(atom)) {
    if (isWigglyBond(bond, atom)) {
      return false;
    }
  }
  RDGeom::Point3D cen = conf.getAtomPos(atom->getIdx());
  RDGeom::Point3D v[6];
  unsigned int count = 0;

  ROMol::ADJ_ITER nbrIdx, endNbrs;
  boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
  while (nbrIdx != endNbrs) {
    if (count == 6) {
      return false;
    }
    RDGeom::Point3D p = conf.getAtomPos(*nbrIdx);
    v[count] = cen.directionVector(p);
    ++count;
    ++nbrIdx;
  }

  if (count < 3) {
    return false;
  }

  unsigned char pair[6];
  memset(pair, 0, 6);

  unsigned int pairs = 0;
  for (unsigned int i = 0; i < count; i++) {
    for (unsigned int j = i + 1; j < count; j++) {
      if (v[i].dotProduct(v[j]) < -(1 - tolerance)) {
        if (pair[i] || pair[j]) {
          return false;
        }
        pair[i] = j + 1;
        pair[j] = i + 1;
        pairs++;
      }
    }
  }

#if 0
  printf("count=%u pairs=%u [%u,%u,%u,%u,%u,%u]\n", count, pairs,
         pair[0], pair[1], pair[2], pair[3], pair[4], pair[5]);
#endif

  Atom::ChiralType tag;
  unsigned int perm;
  bool res = false;
  switch (pairs) {
    case 0:
      break;
    case 1:
      switch (count) {
        case 3: /* T-shape */
          atom->setChiralTag(Atom::ChiralType::CHI_SQUAREPLANAR);
          res = true;
          if (pair[0] == 0) {
            perm = 3;  // Z
          } else if (pair[0] == 2) {
            perm = 2;  // 4
          } else /* pair[0] == 3 */ {
            perm = 1;  // U
          }
          atom->setProp(common_properties::_chiralPermutation, perm);
          break;
        case 4:                /* See-saw */
          if (pair[0] == 2) {  // a b
            if (v[2].angleTo(v[3]) < 100 * M_PI / 180.0) {
              tag = Atom::ChiralType::CHI_OCTAHEDRAL;
              perm = VOLTEST(0, 2, 3) ? 25 : 29;
            } else {
              tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
              perm = VOLTEST(0, 2, 3) ? 7 : 8;
            }
          } else if (pair[0] == 3) {  // a c
            if (v[1].angleTo(v[3]) < 100 * M_PI / 180.0) {
              tag = Atom::ChiralType::CHI_OCTAHEDRAL;
              perm = VOLTEST(0, 1, 3) ? 19 : 23;
            } else {
              tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
              perm = VOLTEST(0, 1, 3) ? 5 : 6;
            }
          } else if (pair[0] == 4) {  // a d
            if (v[1].angleTo(v[2]) < 100 * M_PI / 180.0) {
              tag = Atom::ChiralType::CHI_OCTAHEDRAL;
              perm = VOLTEST(0, 1, 2) ? 6 : 17;
            } else {
              tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
              perm = VOLTEST(0, 1, 2) ? 3 : 4;
            }
          } else if (pair[1] == 3) {  // b c
            if (v[0].angleTo(v[3]) < 100 * M_PI / 180.0) {
              tag = Atom::ChiralType::CHI_OCTAHEDRAL;
              perm = VOLTEST(0, 1, 3) ? 10 : 8;
            } else {
              tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
              perm = VOLTEST(1, 0, 3) ? 13 : 14;
            }
          } else if (pair[1] == 4) {  // b d
            if (v[0].angleTo(v[2]) < 100 * M_PI / 180.0) {
              tag = Atom::ChiralType::CHI_OCTAHEDRAL;
              perm = VOLTEST(0, 1, 3) ? 1 : 2;
            } else {
              tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
              perm = VOLTEST(1, 0, 2) ? 10 : 12;
            }
          } else /* pair[2] == 4 */ {  // c d
            if (v[0].angleTo(v[1]) < 100 * M_PI / 180.0) {
              tag = Atom::ChiralType::CHI_OCTAHEDRAL;
              perm = VOLTEST(0, 1, 3) ? 4 : 14;
            } else {
              tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
              perm = VOLTEST(3, 0, 1) ? 16 : 19;
            }
          }
          atom->setChiralTag(tag);
          res = true;
          atom->setProp(common_properties::_chiralPermutation, perm);
          break;
        case 5: /* Trigonal bipyramidal */
          atom->setChiralTag(Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL);
          res = true;
          if (pair[0] == 2) {
            perm = VOLTEST(0, 2, 3) ? 7 : 8;  // a b
          } else if (pair[0] == 3) {
            perm = VOLTEST(0, 1, 3) ? 5 : 6;  // a c
          } else if (pair[0] == 4) {
            perm = VOLTEST(0, 1, 2) ? 3 : 4;  // a d
          } else if (pair[0] == 5) {
            perm = VOLTEST(0, 1, 2) ? 1 : 2;  // a e
          } else if (pair[1] == 3) {
            perm = VOLTEST(1, 0, 3) ? 13 : 14;  // b c
          } else if (pair[1] == 4) {
            perm = VOLTEST(1, 0, 2) ? 10 : 12;  // b d
          } else if (pair[1] == 5) {
            perm = VOLTEST(1, 0, 2) ? 9 : 11;  // b e
          } else if (pair[2] == 4) {
            perm = VOLTEST(2, 0, 1) ? 16 : 19;  // c d
          } else if (pair[2] == 5) {
            perm = VOLTEST(2, 0, 1) ? 15 : 20;  // c e
          } else /* pair[2] == 4 */ {
            perm = VOLTEST(3, 0, 1) ? 17 : 18;  // d e
          }
          atom->setProp(common_properties::_chiralPermutation, perm);
          break;
      }
      break;
    case 2:
      if (count == 4) {
        /* Square planar */
        atom->setChiralTag(Atom::ChiralType::CHI_SQUAREPLANAR);
        res = true;
        if (pair[0] == 2) {
          perm = 2;  // 4
        } else if (pair[0] == 3) {
          perm = 1;  // U
        } else /* pair[1] == 4 */ {
          perm = 3;  // Z
        }
        atom->setProp(common_properties::_chiralPermutation, perm);
      } else if (count == 5) {
        /* Square pyramidal */
        atom->setChiralTag(Atom::ChiralType::CHI_OCTAHEDRAL);
        res = true;
        perm = OctahedralPermFrom3D(pair, v);
        atom->setProp(common_properties::_chiralPermutation, perm);
      }
      break;
    case 3:
      if (count == 6) {
        /* Octahedral */
        atom->setChiralTag(Atom::ChiralType::CHI_OCTAHEDRAL);
        res = true;
        perm = OctahedralPermFrom3D(pair, v);
        atom->setProp(common_properties::_chiralPermutation, perm);
      }
      break;
  }
  return res;
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

  auto allowNontetrahedralStereo = Chirality::getAllowNontetrahedralChirality();

  boost::dynamic_bitset<> explicitAtoms;
  explicitAtoms.resize(mol.getNumAtoms(), 0);
  for (auto bond : mol.bonds()) {
    auto bondDir = bond->getBondDir();
    if (bondDir == Bond::BondDir::BEGINWEDGE ||
        bondDir == Bond::BondDir::BEGINDASH) {
      explicitAtoms[bond->getBeginAtom()->getIdx()] = 1;
    }
  }

  for (auto atom : mol.atoms()) {
    if (atom->getChiralTag() != Atom::ChiralType::CHI_UNSPECIFIED) {
      explicitAtoms[atom->getIdx()] = 1;
    }
  }

  for (auto atom : mol.atoms()) {
    // see if only the explicitly wedged atoms are to be used

    // if we aren't replacing existing tags and the atom is already tagged,
    // punt:
    if (!replaceExistingTags && atom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
      continue;
    }
    atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // additional reasons to skip the atom:
    auto nzDegree = Chirality::detail::getAtomNonzeroDegree(atom);
    auto tnzDegree = nzDegree + atom->getTotalNumHs();
    if (nzDegree < 3 || tnzDegree > 6) {
      // not enough explicit neighbors or too many total neighbors
      continue;
    }
    if (allowNontetrahedralStereo &&
        assignNontetrahedralChiralTypeFrom3D(mol, conf, atom)) {
      if (explicitAtoms[atom->getIdx()] == 0) {
        atom->setProp(common_properties::_NonExplicit3DChirality, 1);
      }
      continue;
    }
    /* We're only doing tetrahedral cases here */
    if (tnzDegree > 4) {
      continue;
    }
    int anum = atom->getAtomicNum();
    if (anum != 16 && anum != 34 &&  // S or Se are special
                                     // (just using the InChI list for now)
        tnzDegree != 4               // not enough total neighbors
    ) {
      continue;
    }

    const auto &p0 = conf.getAtomPos(atom->getIdx());
    const RDGeom::Point3D *nbrs[4];
    unsigned int nbrIdx = 0;
    int hasWigglyBond = 0;
    for (const auto bond : mol.atomBonds(atom)) {
      hasWigglyBond = isWigglyBond(bond, atom);
      if (hasWigglyBond) {
        break;
      }
      if (!Chirality::detail::bondAffectsAtomChirality(bond, atom)) {
        continue;
      }
      nbrs[nbrIdx++] = &conf.getAtomPos(bond->getOtherAtomIdx(atom->getIdx()));
    }
    if (hasWigglyBond) {
      continue;
    }
    auto v1 = *nbrs[0] - p0;
    auto v2 = *nbrs[1] - p0;
    auto v3 = *nbrs[2] - p0;

    double chiralVol = v1.dotProduct(v2.crossProduct(v3));
    bool chiralitySet = false;
    if (chiralVol < -ZERO_VOLUME_TOL) {
      atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
      chiralitySet = true;
    } else if (chiralVol > ZERO_VOLUME_TOL) {
      atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
      chiralitySet = true;
    } else if (nbrIdx == 4) {
      // The first three neighbors are on the same plane as the chiral atom (or
      // very close to it). If a 4th neighbor is present, let's see if this one
      // determines a chiral volume

      auto v4 = *nbrs[3] - p0;
      // v4 would be in the opposite direction to v3
      chiralVol = -v1.dotProduct(v2.crossProduct(v4));
      if (chiralVol < -ZERO_VOLUME_TOL) {
        atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
        chiralitySet = true;
      } else if (chiralVol > ZERO_VOLUME_TOL) {
        atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
        chiralitySet = true;
      } else {
        atom->setChiralTag(Atom::CHI_UNSPECIFIED);
      }
    } else {
      atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    }

    if (chiralitySet && explicitAtoms[atom->getIdx()] == 0) {
      atom->setProp<int>(common_properties::_NonExplicit3DChirality, 1);
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
        [&mol](const ROMol::edge_descriptor &e) { return mol[e]->getIdx(); });
    unsigned int atomIdx = atom->getIdx();
    nbrBondIdxList.sort([&mol, atomIdx](const int ai, const int bi) {
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
  if (!mol.getRingInfo()->isSymmSssr()) {
    RDKit::MolOps::symmetrizeSSSR(mol);
  }

  for (auto bond : mol.bonds()) {
    if (isBondCandidateForStereo(bond)) {
      bool isCandidate = true;
      for (const auto bondAtom : {bond->getBeginAtom(), bond->getEndAtom()}) {
        for (const auto nbrBond : mol.atomBonds(bondAtom)) {
          if (nbrBond->getBondType() == Bond::SINGLE ||
              nbrBond->getBondType() == Bond::AROMATIC) {
            singleBondCounts[nbrBond->getIdx()] += 1;
            auto nbrDir = nbrBond->getBondDir();
            int hasUnknownStereo = 0;
            if (nbrBond->getBeginAtom() == bondAtom &&
                nbrDir == Bond::BondDir::UNKNOWN &&
                nbrBond->getPropIfPresent(common_properties::_UnknownStereo,
                                          hasUnknownStereo) &&
                hasUnknownStereo) {
              // if there's a wiggly bond starting here, then we're not a
              // candidate for stereo
              isCandidate = false;
            } else {
              needsDir[bond->getIdx()] = 1;
              if (nbrDir == Bond::BondDir::NONE ||
                  nbrDir == Bond::BondDir::ENDDOWNRIGHT ||
                  nbrDir == Bond::BondDir::ENDUPRIGHT) {
                needsDir[nbrBond->getIdx()] = 1;
                dblBondNbrs[bond->getIdx()].push_back(nbrBond->getIdx());
                // the search may seem inefficient, but these vectors are
                // going to be at most 2 long (with very few exceptions). It's
                // just not worth using a different data structure
                if (std::find(singleBondNbrs[nbrBond->getIdx()].begin(),
                              singleBondNbrs[nbrBond->getIdx()].end(),
                              bond->getIdx()) ==
                    singleBondNbrs[nbrBond->getIdx()].end()) {
                  singleBondNbrs[nbrBond->getIdx()].push_back(bond->getIdx());
                }
              }
            }
          }
          if (!isCandidate) {
            break;
          }
        }
        if (!isCandidate) {
          break;
        }
      }
      if (isCandidate) {
        bondsInPlay.push_back(bond);
      }
    }
  }

  if (!bondsInPlay.size()) {
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
}

void detectBondStereochemistry(ROMol &mol, int confId) {
  if (!mol.getNumConformers()) {
    return;
  }
  const Conformer &conf = mol.getConformer(confId);
  setDoubleBondNeighborDirections(mol, &conf);
}

void clearSingleBondDirFlags(ROMol &mol, bool onlyWedgeFlags) {
  for (auto bond : mol.bonds()) {
    if (bond->getBondType() == Bond::SINGLE) {
      if (bond->getBondDir() == Bond::UNKNOWN) {
        bond->setProp(common_properties::_UnknownStereo, 1);
      }

      if (!onlyWedgeFlags ||
          (bond->getBondDir() != Bond::BondDir::ENDDOWNRIGHT &&
           bond->getBondDir() != Bond::BondDir::ENDUPRIGHT)) {
        bond->setBondDir(Bond::NONE);
      }
    }
  }
}

void clearDirFlags(ROMol &mol, bool onlyWedgeTypeBondDirs) {
  for (auto bond : mol.bonds()) {
    if (bond->getBondDir() == Bond::UNKNOWN ||
        bond->getBondDir() == Bond::BondDir::EITHERDOUBLE) {
      bond->setProp(common_properties::_UnknownStereo, 1);
    }

    if (onlyWedgeTypeBondDirs == false ||
        (bond->getBondDir() != Bond::BondDir::ENDDOWNRIGHT &&
         bond->getBondDir() != Bond::BondDir::ENDUPRIGHT)) {
      bond->setBondDir(Bond::NONE);
    }
  }
}

void clearAllBondDirFlags(ROMol &mol) { clearDirFlags(mol, false); }

void setBondStereoFromDirections(ROMol &mol) {
  mol.clearProp("_needsDetectBondStereo");
  for (Bond *bond : mol.bonds()) {
    if (bond->getBondType() == Bond::DOUBLE &&
        bond->getStereo() != Bond::STEREOANY) {
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
  if (mol.needsUpdatePropertyCache()) {
    mol.updatePropertyCache(false);
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
    Atom *atom = bond->getBeginAtom();
    if (dir == Bond::UNKNOWN) {
      if (atomsSet[atom->getIdx()] || replaceExistingTags) {
        atom->setChiralTag(Atom::CHI_UNSPECIFIED);
        atomsSet.set(atom->getIdx());
      }
    } else {
      // the bond is marked as chiral:
      if (dir == Bond::BEGINWEDGE || dir == Bond::BEGINDASH) {
        if (atomsSet[atom->getIdx()] ||
            (!replaceExistingTags &&
             atom->getChiralTag() != Atom::CHI_UNSPECIFIED)) {
          continue;
        }
        if (atom->getImplicitValence() == -1) {
          atom->calcExplicitValence(false);
          atom->calcImplicitValence(false);
        }
        Atom::ChiralType code =
            Chirality::atomChiralTypeFromBondDirPseudo3D(mol, bond, &conf)
                .value_or(Atom::ChiralType::CHI_UNSPECIFIED);
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
  for (auto atom : mol.atoms()) {
    atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    if (atom->hasProp(common_properties::_CIPCode)) {
      atom->clearProp(common_properties::_CIPCode);
    }
    if (atom->hasProp(common_properties::_CIPRank)) {
      atom->clearProp(common_properties::_CIPRank);
    }
  }
  for (auto bond : mol.bonds()) {
    if (bond->getBondType() == Bond::DOUBLE) {
      bond->setStereo(Bond::BondStereo::STEREONONE);
      bond->getStereoAtoms().clear();
      bond->setBondDir(Bond::BondDir::NONE);
    } else if (bond->getBondType() == Bond::SINGLE) {
      bond->setBondDir(Bond::BondDir::NONE);
    }
  }
  std::vector<StereoGroup> sgs;
  static_cast<RWMol &>(mol).setStereoGroups(std::move(sgs));
}

}  // namespace MolOps
}  // namespace RDKit
