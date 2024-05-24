//
//  Copyright (C) 2004-2021 Tad hurst/CDD  and other RDKit contributors
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
#include <GraphMol/Chirality.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/Atropisomers.h>
#include <Geometry/point.h>
#include <boost/dynamic_bitset.hpp>
#include <algorithm>
#include <RDGeneral/Ranking.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <RDGeneral/BoostEndInclude.h>

constexpr double REALLY_SMALL_BOND_LEN = 0.0000001;

namespace RDKit {
namespace Atropisomers {

bool getAtropisomerAtomsAndBonds(const Bond *bond,
                                 AtropAtomAndBondVec atomsAndBondVects[2],
                                 const ROMol &mol) {
  PRECONDITION(bond, "no bond");
  atomsAndBondVects[0].first = bond->getBeginAtom();
  atomsAndBondVects[1].first = bond->getEndAtom();

  // get the one or two bonds on each end

  for (int bondAtomIndex = 0; bondAtomIndex < 2; ++bondAtomIndex) {
    for (const auto nbrBond :
         mol.atomBonds(atomsAndBondVects[bondAtomIndex].first)) {
      if (nbrBond == bond) {
        continue;  // a bond is NOT its own neighbor
      }
      atomsAndBondVects[bondAtomIndex].second.push_back(nbrBond);
    }
    if (atomsAndBondVects[bondAtomIndex].second.size() == 0) {
      return false;  // no neighbor bonds found
    }

    // make sure the bond with this lowest atom is is first

    if (atomsAndBondVects[bondAtomIndex].second.size() == 2 &&
        atomsAndBondVects[bondAtomIndex]
                .second[1]
                ->getOtherAtom(atomsAndBondVects[bondAtomIndex].first)
                ->getIdx() <
            atomsAndBondVects[bondAtomIndex]
                .second[0]
                ->getOtherAtom(atomsAndBondVects[bondAtomIndex].first)
                ->getIdx()) {
      std::swap(atomsAndBondVects[bondAtomIndex].second[0],
                atomsAndBondVects[bondAtomIndex].second[1]);
    }
  }

  return true;
}

bool getBondFrameOfReference(const Bond *bond, const Conformer *conf,
                             RDGeom::Point3D &xAxis, RDGeom::Point3D &yAxis,
                             RDGeom::Point3D &zAxis) {
  // create a frame of reference that has its X-axis along the atrop bond
  // for 2D confs, the yAxis is in the 2D plane and the zAxis is perpendicular
  // to that plane) for 3D confs  the yAxis and the zAxis are arbitrary.

  PRECONDITION(bond, "bad bond");

  xAxis = conf->getAtomPos(bond->getEndAtom()->getIdx()) -
          conf->getAtomPos(bond->getBeginAtom()->getIdx());
  if (xAxis.length() < REALLY_SMALL_BOND_LEN) {
    return false;  // bond len is xero
  }
  xAxis.normalize();
  if (!conf->is3D()) {
    yAxis = RDGeom::Point3D(-xAxis.y, xAxis.x, 0);
    yAxis.normalize();
    zAxis = RDGeom::Point3D(0.0, 0.0, 1.0);
    return true;
  }

  // here for 3D conf

  if (xAxis.x > REALLY_SMALL_BOND_LEN || xAxis.y > REALLY_SMALL_BOND_LEN) {
    zAxis = RDGeom::Point3D(-xAxis.y, xAxis.x,
                            0);  // temp z axis - used to find yAxis
  } else {
    zAxis = RDGeom::Point3D(xAxis.z, xAxis.z,
                            0);  // temp z axis - used to find yAxis
  }

  yAxis = zAxis.crossProduct(xAxis);
  zAxis = xAxis.crossProduct(yAxis);
  yAxis.normalize();
  zAxis.normalize();

  return true;
}

Bond::BondDir getBondDirForAtropisomerNoConf(Bond::BondStereo bondStereo,
                                             unsigned int whichEnd,
                                             unsigned int whichBond) {
  // the convention is that in the absence of coords, the coordiates are choosen
  // with the lowest numbered atom of the atrop bond down, and the other atom
  // straight up.
  // On each end, the lowest numbered connecting atom is on the left
  //
  //              a      b
  //               \   /
  //                 c
  //                 |
  //                 d
  //               /   \     aaa
  //              e      f
  //
  // where  c > d
  //        a < b
  //        e < f

  PRECONDITION(whichEnd <= 1, "whichEnd must be 0 or 1");
  PRECONDITION(whichBond <= 1, "whichBond must be 0 or 1");
  PRECONDITION(bondStereo == Bond::BondStereo::STEREOATROPCW ||
                   bondStereo == Bond::BondStereo::STEREOATROPCCW,
               "bondStereo must be BondAtropisomerCW or BondAtropisomerCCW");

  int flips = 0;
  if (bondStereo == Bond::BondStereo::STEREOATROPCW) {
    ++flips;
  }
  if (whichBond == 1) {
    ++flips;
  }
  if (whichEnd == 1) {
    ++flips;
  }

  return flips % 2 ? Bond::BEGINDASH : Bond::BEGINWEDGE;
}

Bond::BondDir getBondDirForAtropisomer2d(RDGeom::Point3D bondVecs[2],
                                         Bond::BondStereo bondStereo,
                                         unsigned int whichEnd,
                                         unsigned int whichBond) {
  PRECONDITION(whichEnd <= 1, "whichEnd must be 0 or 1");
  PRECONDITION(whichBond <= 1, "whichBond must be 0 or 1");
  PRECONDITION(bondStereo == Bond::BondStereo::STEREOATROPCW ||
                   bondStereo == Bond::BondStereo::STEREOATROPCCW,
               "bondStereo must be BondAtropisomerCW or BondAtropisomerCCW");

  int flips = 0;
  if (bondStereo == Bond::BondStereo::STEREOATROPCCW) {
    ++flips;
  }
  if (whichBond == 1) {
    ++flips;
  }
  if (whichEnd == 1) {
    ++flips;
  }
  if (bondVecs[1 - whichEnd].y < 0) {
    ++flips;  // if the OTHER end is negative for the low index bond vec, it
              // is a flip
  }

  return flips % 2 ? Bond::BEGINWEDGE : Bond::BEGINDASH;
}

Bond::BondDir getBondDirForAtropisomer3d(Bond *whichBond,
                                         const Conformer *conf) {
  // for 3D we mark it as wedge or hash depending on the z-value of the bond
  // vector
  //  IT really doesn't matter since we ignore these except as MARKERS for
  //  which bonds are atropisomer bonds
  if ((conf->getAtomPos(whichBond->getEndAtom()->getIdx()).z -
       conf->getAtomPos(whichBond->getBeginAtom()->getIdx()).z) >
      REALLY_SMALL_BOND_LEN) {
    return Bond::BondDir::BEGINWEDGE;
  } else {
    return Bond::BondDir::BEGINDASH;
  }
}

bool getAtropIsomerEndVect(const AtropAtomAndBondVec &atomAndBondVec,
                           const RDGeom::Point3D &yAxis,
                           const RDGeom::Point3D &zAxis, const Conformer *conf,
                           RDGeom::Point3D &bondVec) {
  PRECONDITION(
      atomAndBondVec.second.size() > 0 && atomAndBondVec.second.size() < 3,
      "bad bond size");
  PRECONDITION(atomAndBondVec.second[0], "bad first bond");
  PRECONDITION(atomAndBondVec.second.size() == 1 || atomAndBondVec.second[1],
               "bad second bond");

  bondVec = conf->getAtomPos(atomAndBondVec.second[0]
                                 ->getOtherAtom(atomAndBondVec.first)
                                 ->getIdx()) -
            conf->getAtomPos(
                atomAndBondVec.first->getIdx());  // in old frame of reference

  bondVec = RDGeom::Point3D(0.0, bondVec.dotProduct(yAxis),
                            bondVec.dotProduct(zAxis));  // in new frame

  // make sure the other atom is on the other side

  if (atomAndBondVec.second.size() == 2) {
    RDGeom::Point3D otherVec =
        conf->getAtomPos(atomAndBondVec.second[1]
                             ->getOtherAtom(atomAndBondVec.first)
                             ->getIdx()) -
        conf->getAtomPos(
            atomAndBondVec.first->getIdx());  // in old frame of reference
    otherVec = RDGeom::Point3D(0.0, otherVec.dotProduct(yAxis),
                               otherVec.dotProduct(zAxis));  // in new frame

    if (bondVec.length() < REALLY_SMALL_BOND_LEN) {
      bondVec = -otherVec;  // put it on the other side of otherVec
    } else if (bondVec.dotProduct(otherVec) > REALLY_SMALL_BOND_LEN) {
      // the product of dotproducts (y-values) should be
      // negative (or at least zero)
      BOOST_LOG(rdWarningLog)
          << "Both bonds on one end of an atropisomer are on the same side - atoms is : "
          << atomAndBondVec.first->getIdx() << std::endl;
      return false;
    }
  }
  if (bondVec.length() < REALLY_SMALL_BOND_LEN) {
    BOOST_LOG(rdWarningLog)
        << "Could not find a bond on one end of an atropisomer that is not co-linear - atoms are : "
        << atomAndBondVec.first->getIdx() << std::endl;
    return false;
  }

  bondVec.normalize();
  return true;
}

std::pair<bool, Bond::BondDir> getBondDir(
    const Bond *bond, const AtropAtomAndBondVec &atomAndBondVec) {
  // get the wedge dir for this end of the bond
  // if the first bond 1 has a bondDir, use it
  // if the second bond has a bond dir use the opposite of if
  // if both bonds have a dir, make sure they are different

  auto bond1Dir = atomAndBondVec.second[0]->getBondDir();
  if (bond1Dir != Bond::BEGINWEDGE && bond1Dir != Bond::BEGINDASH) {
    bond1Dir = Bond::NONE;  //  we dont care if it any thing else
  }
  auto bond2Dir = atomAndBondVec.second.size() == 2
                      ? atomAndBondVec.second[1]->getBondDir()
                      : Bond::NONE;
  if (bond2Dir != Bond::BEGINWEDGE && bond2Dir != Bond::BEGINDASH) {
    bond2Dir = Bond::NONE;
  }

  // if both are set to a direction, they must NOT be the same - one
  // must be a dash and the other a hash

  if (bond1Dir != Bond::NONE && bond2Dir != Bond::NONE &&
      bond1Dir == bond2Dir) {
    BOOST_LOG(rdWarningLog)
        << "The bonds on one end of an atropisomer are both UP or both DOWN - atoms are: "
        << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    return {false, Bond::BondDir::NONE};
  }

  if (bond1Dir == Bond::BEGINWEDGE || bond2Dir == Bond::BEGINDASH) {
    return {true, Bond::BondDir::BEGINWEDGE};
  }
  if (bond1Dir == Bond::BEGINDASH || bond2Dir == Bond::BEGINWEDGE) {
    return {true, Bond::BondDir::BEGINDASH};
  }
  return {true, Bond::BondDir::NONE};
}

bool DetectAtropisomerChiralityOneBond(Bond *bond, ROMol &mol,
                                       const Conformer *conf) {
  // the approach is this:
  // we will view the system along the line from the potential atropisomer
  // bond, from atom1 to atom 2 and we do a coordinate transformation to
  // the plane of reference where that vector, from a1 to a2, is the x-AXIS.
  // For 2D, the y axis is in the 2D plane, and the zaxis is perpendicaul to
  // the 2D plane For 3D, the Y and Z axes are taken arbitrarily to form
  // a right-handed system with the X-axis.
  //  atoms 1 and 2 each have one or two bonds out from the main potential
  //  atrop bond. for each end of the main bond, we find a vector to reprent
  //  the neighbor atom with the smallest index as its projection onto the
  //  x=0 plane.
  // (In 2d, this projection is on the y-AXIS for the end that does NOT have
  // a wedge/hash bond, and  on the z axis - out of the plane - for the end
  // that does have a wedge/hash). The chirality is recorded as the
  // direction we rotate from, atom 1's projection to atom2's proejection -
  // either clockwise or counter clockwise

  PRECONDITION(bond, "bad bond");

  // one vector for each end - each one - should end up with 1 or 2 entries
  AtropAtomAndBondVec atomAndBondVecs[2];
  if (!getAtropisomerAtomsAndBonds(bond, atomAndBondVecs, mol)) {
    return false;  // not an atropisomer
  }

  // make sure we do not have wiggle bonds

  for (auto atomAndBondVec : atomAndBondVecs) {
    for (auto endBond : atomAndBondVec.second) {
      if (endBond->getBondDir() == Bond::UNKNOWN) {
        return false;  // not an atropisomer
      }
    }
  }

  // the convention is that in the absence of coords, the coordiates are choosen
  // with the lowest numbered atom of the atrop bond down, and the other atom
  // straight up.
  // On each end, the lowest numbered connecting atom is on the left
  //
  //              a      b
  //               \   /
  //                 c
  //                 |
  //                 d
  //               /   \     aaa
  //              e      f
  //
  // where  c > d
  //        a < b
  //        e < f

  if (conf == nullptr) {
    std::pair<bool, Bond::BondDir> bond1DirResult;
    bond1DirResult = getBondDir(bond, atomAndBondVecs[0]);
    if (!bond1DirResult.first) {
      return false;
    }
    std::pair<bool, Bond::BondDir> bond2DirResult;
    bond2DirResult = getBondDir(bond, atomAndBondVecs[1]);
    if (!bond2DirResult.first) {
      return false;
    }
    if (bond1DirResult.second == bond2DirResult.second) {
      BOOST_LOG(rdWarningLog)
          << "inconsistent bond wedging for an atropisomer.  Atoms are: "
          << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
          << std::endl;
      return false;
    }
    if (bond1DirResult.second == Bond::BEGINWEDGE ||
        bond2DirResult.second == Bond::BEGINDASH) {
      bond->setStereo(Bond::BondStereo::STEREOATROPCCW);
    } else if (bond1DirResult.second == Bond::BEGINDASH ||
               bond2DirResult.second == Bond::BEGINWEDGE) {
      bond->setStereo(Bond::BondStereo::STEREOATROPCW);
    }

    return true;
  }

  // create a frame of reference that has its X-axis along the atrop bond

  RDGeom::Point3D xAxis, yAxis, zAxis;
  if (!getBondFrameOfReference(bond, conf, xAxis, yAxis, zAxis)) {
    // connot percieve atroisomer
    BOOST_LOG(rdWarningLog)
        << "Failed to get a frame of reference along an atropisomer bond - atoms are: "
        << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    return false;
  }
  RDGeom::Point3D bondVecs[2];  // one bond vector from each end of the
                                // potential atropisomer bond

  for (int bondAtomIndex = 0; bondAtomIndex < 2; ++bondAtomIndex) {
    // if the conf is 2D, we use the wedge bonds to set the coords for the
    // projected vector onto the xAxis perpendicular plane (looking down
    // the atrop bond )

    if (!conf->is3D()) {
      // get the wedge dir for this end of the bond
      // if the first bond 1 has a bondDir, use it
      // if the second bond has a bond dir use the opposite of if
      // if both bonds have a dir, make sure they are different

      std::pair<bool, Bond::BondDir> bondDirResult;

      bondDirResult = getBondDir(bond, atomAndBondVecs[bondAtomIndex]);
      if (!bondDirResult.first) {
        return false;
      }

      if (!getAtropIsomerEndVect(atomAndBondVecs[bondAtomIndex], yAxis, zAxis,
                                 conf, bondVecs[bondAtomIndex])) {
        return false;
      }

      if (bondDirResult.second == Bond::BEGINWEDGE) {
        bondVecs[bondAtomIndex].y *= 0.707;
        bondVecs[bondAtomIndex].z = fabs(bondVecs[bondAtomIndex].y);
      } else if (bondDirResult.second == Bond::BEGINDASH) {
        bondVecs[bondAtomIndex].y *= 0.707;
        bondVecs[bondAtomIndex].z = -fabs(bondVecs[bondAtomIndex].y);
      }
    } else {  // the conf is 3D
      // to be considered, one or more neighbor bonds must have a wedge or
      // hash

      // find the projection of the bond(s) on this end in the frame of
      // reference's  x=0  plane
      RDGeom::Point3D tempBondVec =
          conf->getAtomPos(
              atomAndBondVecs[bondAtomIndex]
                  .second[0]
                  ->getOtherAtom(atomAndBondVecs[bondAtomIndex].first)
                  ->getIdx()) -
          conf->getAtomPos(atomAndBondVecs[bondAtomIndex].first->getIdx());
      bondVecs[bondAtomIndex] = RDGeom::Point3D(
          0.0, tempBondVec.dotProduct(yAxis), tempBondVec.dotProduct(zAxis));

      if (atomAndBondVecs[bondAtomIndex].second.size() == 2) {
        tempBondVec =
            conf->getAtomPos(
                atomAndBondVecs[bondAtomIndex]
                    .second[1]
                    ->getOtherAtom(atomAndBondVecs[bondAtomIndex].first)
                    ->getIdx()) -
            conf->getAtomPos(atomAndBondVecs[bondAtomIndex].first->getIdx());

        // get the projection of the 2nd bond on the x=0 plane

        RDGeom::Point3D otherBondVec = RDGeom::Point3D(
            0.0, tempBondVec.dotProduct(yAxis), tempBondVec.dotProduct(zAxis));

        // if the first atom is co-linear with the main atrop bond, use
        // the opposite of the 2nd atom

        if (bondVecs[bondAtomIndex].length() < REALLY_SMALL_BOND_LEN) {
          bondVecs[bondAtomIndex] =
              -otherBondVec;  // note - it might still be co-linear-
                              // this is checked below
        } else if (bondVecs[bondAtomIndex].dotProduct(otherBondVec) >
                   REALLY_SMALL_BOND_LEN) {
          BOOST_LOG(rdWarningLog)
              << "Both bonds on one end of an atropisomer are on the same side - atoms are: "
              << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
              << std::endl;
          return false;
        }
      }

      if (bondVecs[bondAtomIndex].length() < REALLY_SMALL_BOND_LEN) {
        BOOST_LOG(rdWarningLog)
            << "Failed to find a bond on one end of an atropisomer that is NOT co-linear - atoms are: "
            << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
            << std::endl;
        return false;
      }
    }
  }

  auto crossProduct = bondVecs[1].crossProduct(bondVecs[0]);

  if (crossProduct.x > REALLY_SMALL_BOND_LEN) {
    bond->setStereo(Bond::BondStereo::STEREOATROPCCW);
  } else if (crossProduct.x < -REALLY_SMALL_BOND_LEN) {
    bond->setStereo(Bond::BondStereo::STEREOATROPCW);
  } else {
    BOOST_LOG(rdWarningLog)
        << "The 2 defining bonds for an atropisomer are co-planar - atoms are: "
        << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    return false;
  }

  return true;
}

void cleanupAtropisomerStereoGroups(ROMol &mol) {
  std::vector<StereoGroup> newsgs;
  for (auto sg : mol.getStereoGroups()) {
    std::vector<Atom *> okatoms;
    std::vector<Bond *> okbonds;

    for (auto atom : sg.getAtoms()) {
      bool foundAtrop = false;
      for (auto bndI : boost::make_iterator_range(mol.getAtomBonds(atom))) {
        auto bond = (mol)[bndI];
        if (bond->getStereo() == Bond::BondStereo::STEREOATROPCCW ||
            bond->getStereo() == Bond::BondStereo::STEREOATROPCW) {
          foundAtrop = true;
          if (std::find(okbonds.begin(), okbonds.end(), bond) ==
              okbonds.end()) {
            okbonds.push_back(bond);
          }
        }
      }

      if (!foundAtrop) {
        okatoms.push_back(atom);
      }
    }

    if (okbonds.empty()) {
      newsgs.push_back(sg);
    } else {
      newsgs.emplace_back(sg.getGroupType(), std::move(okatoms),
                          std::move(okbonds));
    }
  }
  mol.setStereoGroups(std::move(newsgs));
}

void detectAtropisomerChirality(ROMol &mol, const Conformer *conf) {
  PRECONDITION(conf == nullptr || &(conf->getOwningMol()) == &mol,
               "conformer does not belong to molecule");

  std::set<Bond *> bondsToTry;

  for (auto bond : mol.bonds()) {
    if (canHaveDirection(*bond) &&
        (bond->getBondDir() == Bond::BondDir::BEGINDASH ||
         bond->getBondDir() == Bond::BondDir::BEGINWEDGE)) {
      for (const auto &nbrBond : mol.atomBonds(bond->getBeginAtom())) {
        if (nbrBond == bond) {
          continue;  // a bond is NOT its own neighbor
        }
        bondsToTry.insert(nbrBond);
      }
    }
  }

  bool foundAtrop = false;
  for (auto bondToTry : bondsToTry) {
    if (bondToTry->getBeginAtom()->getImplicitValence() == -1) {
      bondToTry->getBeginAtom()->calcExplicitValence(false);
      bondToTry->getBeginAtom()->calcImplicitValence(false);
    }
    if (bondToTry->getEndAtom()->getImplicitValence() == -1) {
      bondToTry->getEndAtom()->calcExplicitValence(false);
      bondToTry->getEndAtom()->calcImplicitValence(false);
    }
    if (bondToTry->getBondType() != Bond::SINGLE ||
        bondToTry->getStereo() == Bond::BondStereo::STEREOANY ||
        bondToTry->getBeginAtom()->getTotalDegree() < 2 ||
        bondToTry->getEndAtom()->getTotalDegree() < 2 ||
        bondToTry->getBeginAtom()->getTotalDegree() > 3 ||
        bondToTry->getEndAtom()->getTotalDegree() > 3) {
      continue;
    }

    if (DetectAtropisomerChiralityOneBond(bondToTry, mol, conf)) {
      foundAtrop = true;
    }
  }

  if (foundAtrop) {
    cleanupAtropisomerStereoGroups(mol);
  }
}
void getAllAtomIdsForStereoGroup(
    const ROMol &mol, const StereoGroup &group,
    std::vector<unsigned int> &atomIds,
    const std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
        &wedgeBonds) {
  atomIds.clear();
  for (auto &&atom : group.getAtoms()) {
    atomIds.push_back(atom->getIdx());
  }

  for (auto &&bond : group.getBonds()) {
    // figure out which atoms of the bond get wedge/hash indications
    // mark the atom with the wedge/hash

    for (auto atom : {bond->getBeginAtom(), bond->getEndAtom()}) {
      for (const auto atomBond : mol.atomBonds(atom)) {
        if (atomBond->getIdx() == bond->getIdx()) {
          continue;
        }

        if (atomBond->getBondDir() == Bond::BEGINWEDGE ||
            atomBond->getBondDir() == Bond::BEGINDASH ||
            (wedgeBonds.find(atomBond->getIdx()) != wedgeBonds.end() &&
             (wedgeBonds.at(atomBond->getIdx())->getType()) ==
                 Chirality::WedgeInfoType::WedgeInfoTypeAtropisomer)) {
          if (std::find(atomIds.begin(), atomIds.end(), atom->getIdx()) ==
              atomIds.end()) {
            atomIds.push_back(atom->getIdx());
          }
        }
      }
    }
  }
}

bool WedgeBondFromAtropisomerOneBondNoConf(
    Bond *bond, const ROMol &mol,
    std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
        &wedgeBonds) {
  PRECONDITION(bond, "no bond");

  AtropAtomAndBondVec atomAndBondVecs[2];
  if (!getAtropisomerAtomsAndBonds(bond, atomAndBondVecs, mol)) {
    return false;  // not an atropisomer
  }

  //  make sure we do not have wiggle bonds

  for (auto atomAndBondVec : atomAndBondVecs) {
    for (auto endBond : atomAndBondVec.second) {
      if (endBond->getBondDir() == Bond::UNKNOWN) {
        return false;  // not an atropisomer)
      }
    }
  }

  // first see if any candidate bond is already set to a wedge or hash
  // if so, we will use that bond as a wedge or hash

  std::vector<int> useBondsAtEnd[2];
  bool foundBondDir = false;

  for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    for (unsigned int whichBond = 0;
         whichBond < atomAndBondVecs[whichEnd].second.size(); ++whichBond) {
      auto bondDir = atomAndBondVecs[whichEnd].second[whichBond]->getBondDir();

      // see if it is a wedge or hash and its origin is the atom in the
      // main bond

      if ((bondDir == Bond::BEGINWEDGE || bondDir == Bond::BEGINDASH) &&
          atomAndBondVecs[whichEnd].second[whichBond]->getBeginAtom() ==
              atomAndBondVecs[whichEnd].first &&
          canHaveDirection(*bond)) {
        useBondsAtEnd[whichEnd].push_back(whichBond);
        foundBondDir = true;
      }
    }
  }

  if (foundBondDir) {
    for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
      for (unsigned int whichBondIndex = 0;
           whichBondIndex < useBondsAtEnd[whichEnd].size(); ++whichBondIndex) {
        atomAndBondVecs[whichEnd]
            .second[useBondsAtEnd[whichEnd][whichBondIndex]]
            ->setBondDir(getBondDirForAtropisomerNoConf(
                bond->getStereo(), whichEnd,
                useBondsAtEnd[whichEnd][whichBondIndex]));
      }
    }

    return true;
  }

  // did not find a good bond dir - pick one to use
  // we would like to have one that is not in a ring, and will be a wedge

  const RingInfo *ri = bond->getOwningMol().getRingInfo();

  int bestBondEnd = -1, bestBondNumber = -1;
  bool bestBondIsSingle = false;
  unsigned int bestRingCount = INT_MAX;
  Bond::BondDir bestBondDir = Bond::BondDir::NONE;
  for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    for (unsigned int whichBond = 0;
         whichBond < atomAndBondVecs[whichEnd].second.size(); ++whichBond) {
      auto bondToTry = atomAndBondVecs[whichEnd].second[whichBond];

      if (!canHaveDirection(*bondToTry) ||
          wedgeBonds.find(bondToTry->getIdx()) != wedgeBonds.end()) {
        continue;  // must be a single OR aromatic bond and not already
                   // spoken for by a chiral center
      }

      if (bondToTry->getBondDir() != Bond::BondDir::NONE) {
        if (bondToTry->getBeginAtom()->getIdx() ==
            atomAndBondVecs[whichEnd].first->getIdx()) {
          BOOST_LOG(rdWarningLog)
              << "Wedge or hash bond found on atropisomer where not expected - atoms are: "
              << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
              << std::endl;
          return false;
        } else {
          continue;  // wedge or hash bond affecting the OTHER atom
                     // = perhaps a chiral center
        }
      }
      auto ringCount = ri->numBondRings(bondToTry->getIdx());
      if (ringCount > bestRingCount) {
        continue;
      }

      else if (ringCount < bestRingCount) {
        bestBondEnd = whichEnd;
        bestBondNumber = whichBond;
        bestRingCount = ringCount;
        bestBondIsSingle = (bondToTry->getBondType() == Bond::BondType::SINGLE);
        bestBondDir = getBondDirForAtropisomerNoConf(bond->getStereo(),
                                                     whichEnd, whichBond);
      } else if (bestBondIsSingle &&
                 bondToTry->getBondType() != Bond::BondType::SINGLE) {
        continue;

      } else if (!bestBondIsSingle &&
                 bondToTry->getBondType() == Bond::BondType::SINGLE) {
        bestBondEnd = whichEnd;
        bestBondNumber = whichBond;
        bestRingCount = ringCount;
        bestBondIsSingle = true;
        bestBondDir = getBondDirForAtropisomerNoConf(bond->getStereo(),
                                                     whichEnd, whichBond);

      } else {
        auto bondDir = getBondDirForAtropisomerNoConf(bond->getStereo(),
                                                      whichEnd, whichBond);
        if (bestBondDir == Bond::BondDir::NONE ||
            (bestBondDir == Bond::BondDir::BEGINDASH &&
             bondDir == Bond::BondDir::BEGINWEDGE)) {
          bestBondEnd = whichEnd;
          bestBondNumber = whichBond;
          bestRingCount = ringCount;
          bestBondIsSingle =
              (bondToTry->getBondType() == Bond::BondType::SINGLE);
          bestBondDir = bondDir;
        }
      }
    }
  }

  if (bestBondEnd >= 0)  // we found a good one
  {
    // make sure the atoms on the bond are in the right order for the
    // wedge/hash the atom on the end of the main bond must be listed
    // first for the wedge/has bond

    auto bestBond = atomAndBondVecs[bestBondEnd].second[bestBondNumber];
    if (bestBond->getBeginAtom() != atomAndBondVecs[bestBondEnd].first) {
      bestBond->setEndAtom(bestBond->getBeginAtom());
      bestBond->setBeginAtom(atomAndBondVecs[bestBondEnd].first);
    }

    bestBond->setBondDir(bestBondDir);

    auto newWedgeInfo = std::unique_ptr<RDKit::Chirality::WedgeInfoBase>(
        new RDKit::Chirality::WedgeInfoAtropisomer(bond->getIdx(),
                                                   bestBondDir));
    wedgeBonds[bestBond->getIdx()] = std::move(newWedgeInfo);
  } else {
    BOOST_LOG(rdWarningLog)
        << "Failed to find a good bond to set as UP or DOWN for an atropisomer - atoms are: "
        << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    return false;
  }

  return true;
}

bool WedgeBondFromAtropisomerOneBond2d(
    Bond *bond, const ROMol &mol, const Conformer *conf,
    std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
        &wedgeBonds) {
  PRECONDITION(bond, "no bond");

  AtropAtomAndBondVec atomAndBondVecs[2];
  if (!getAtropisomerAtomsAndBonds(bond, atomAndBondVecs, mol)) {
    return false;  // not an atropisomer
  }

  //  make sure we do not have wiggle bonds

  for (auto atomAndBondVec : atomAndBondVecs) {
    for (auto endBond : atomAndBondVec.second) {
      if (endBond->getBondDir() == Bond::UNKNOWN) {
        return false;  // not an atropisomer)
      }
    }
  }

  // create a frame of reference that has its X-axis along the atrop bond

  RDGeom::Point3D xAxis, yAxis, zAxis;

  if (!getBondFrameOfReference(bond, conf, xAxis, yAxis, zAxis)) {
    // connot percieve atroisomer bond

    BOOST_LOG(rdWarningLog)
        << "Cound not get a frame of reference for an atropisomer bond - atoms are: "
        << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    return false;
  }

  RDGeom::Point3D bondVecs[2];  // one bond vector from each end of the
                                // potential atropisome bond

  for (int bondAtomIndex = 0; bondAtomIndex < 2; ++bondAtomIndex) {
    // find a vector to represent the lowest numbered atom on each end
    // this vector is NOT the bond vector, but is y-value in the bond
    // frame or reference

    if (!getAtropIsomerEndVect(atomAndBondVecs[bondAtomIndex], yAxis, zAxis,
                               conf, bondVecs[bondAtomIndex])) {
      return false;
    }

    if (bondVecs[bondAtomIndex].length() < REALLY_SMALL_BOND_LEN) {
      // did not find a non-colinear bond

      BOOST_LOG(rdWarningLog)
          << "Failed to get a representative vector for the defining bond of an atropisomer - atoms are: "
          << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
          << std::endl;
      return false;
    }
  }

  // first see if any candidate bond is already set to a wedge or hash
  // if so, we will use that bond as a wedge or hash

  std::vector<int> useBondsAtEnd[2];
  bool foundBondDir = false;

  for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    for (unsigned int whichBond = 0;
         whichBond < atomAndBondVecs[whichEnd].second.size(); ++whichBond) {
      auto bondDir = atomAndBondVecs[whichEnd].second[whichBond]->getBondDir();

      // see if it is a wedge or hash and its origin is the atom in the
      // main bond

      if ((bondDir == Bond::BEGINWEDGE || bondDir == Bond::BEGINDASH) &&
          atomAndBondVecs[whichEnd].second[whichBond]->getBeginAtom() ==
              atomAndBondVecs[whichEnd].first &&
          canHaveDirection(*bond)) {
        useBondsAtEnd[whichEnd].push_back(whichBond);
        foundBondDir = true;
      }
    }
  }

  if (foundBondDir) {
    for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
      for (unsigned int whichBondIndex = 0;
           whichBondIndex < useBondsAtEnd[whichEnd].size(); ++whichBondIndex) {
        atomAndBondVecs[whichEnd]
            .second[useBondsAtEnd[whichEnd][whichBondIndex]]
            ->setBondDir(getBondDirForAtropisomer2d(
                bondVecs, bond->getStereo(), whichEnd,
                useBondsAtEnd[whichEnd][whichBondIndex]));
      }
    }

    return true;
  }

  // did not find a good bond dir - pick one to use
  // we would like to have one that is in a ring, and will favor it being a
  // wedge

  // We favor rings here because wedging non-ring bonds makes it too likely that
  // we'll end up accidentally creating new atropisomeric bonds. This was github
  // issue 7371

  const RingInfo *ri = bond->getOwningMol().getRingInfo();

  int bestBondEnd = -1, bestBondNumber = -1;
  bool bestBondIsSingle = false;
  unsigned int bestRingCount = INT_MAX;
  unsigned int largestRingSize = 0;
  Bond::BondDir bestBondDir = Bond::BondDir::NONE;
  for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    for (unsigned int whichBond = 0;
         whichBond < atomAndBondVecs[whichEnd].second.size(); ++whichBond) {
      auto bondToTry = atomAndBondVecs[whichEnd].second[whichBond];

      if (!canHaveDirection(*bondToTry) ||
          wedgeBonds.find(bondToTry->getIdx()) != wedgeBonds.end()) {
        continue;  // must be a single OR aromatic bond and not already
                   // spoken for by a chiral center
      }

      if (bondToTry->getBondDir() != Bond::BondDir::NONE) {
        if (bondToTry->getBeginAtom()->getIdx() ==
            atomAndBondVecs[whichEnd].first->getIdx()) {
          if (bondToTry->getBondDir() == Bond::BEGINWEDGE ||
              bondToTry->getBondDir() == Bond::BEGINDASH) {
            BOOST_LOG(rdWarningLog)
                << "Wedge or hash bond found on atropisomer where not expected - atoms are: "
                << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
                << std::endl;
            return false;
          } else {
            continue;  // probably a slash up or down for a double bond
          }
        } else {
          continue;  // wedge or hash bond affecting the OTHER atom
                     // = perhaps a chiral center
        }
      }
      auto ringCount = ri->numBondRings(bondToTry->getIdx());
      unsigned int ringSize = 0;
      if (!ringCount) {
        ringCount = 10;
      } else {
        // we're going to prefer to put wedges in larger rings, but don't want
        // to end up wedging macrocyles if it's avoidable.
        ringSize = ri->minBondRingSize(bondToTry->getIdx());
        if (ringSize > 8) {
          ringSize = 0;
        }
      }
      if (ringCount > bestRingCount) {
        continue;
      } else if (ringCount < bestRingCount || ringSize > largestRingSize) {
        bestBondEnd = whichEnd;
        bestBondNumber = whichBond;
        bestRingCount = ringCount;
        largestRingSize = ringSize;
        bestBondIsSingle = (bondToTry->getBondType() == Bond::BondType::SINGLE);
        bestBondDir = getBondDirForAtropisomer2d(bondVecs, bond->getStereo(),
                                                 whichEnd, whichBond);
      } else if (bestBondIsSingle &&
                 bondToTry->getBondType() != Bond::BondType::SINGLE) {
        continue;

      } else if (!bestBondIsSingle &&
                 bondToTry->getBondType() == Bond::BondType::SINGLE) {
        bestBondEnd = whichEnd;
        bestBondNumber = whichBond;
        bestRingCount = ringCount;
        bestBondIsSingle = true;
        bestBondDir = getBondDirForAtropisomer2d(bondVecs, bond->getStereo(),
                                                 whichEnd, whichBond);

      } else {
        auto bondDir = getBondDirForAtropisomer2d(bondVecs, bond->getStereo(),
                                                  whichEnd, whichBond);
        if (bestBondDir == Bond::BondDir::NONE ||
            (bestBondDir == Bond::BondDir::BEGINDASH &&
             bondDir == Bond::BondDir::BEGINWEDGE)) {
          bestBondEnd = whichEnd;
          bestBondNumber = whichBond;
          bestRingCount = ringCount;
          bestBondIsSingle =
              (bondToTry->getBondType() == Bond::BondType::SINGLE);
          bestBondDir = bondDir;
        }
      }
    }
  }

  if (bestBondEnd >= 0) {
    // we found a good one
    // make sure the atoms on the bond are in the right order for the
    // wedge/hash the atom on the end of the main bond must be listed
    // first for the wedge/has bond

    auto bestBond = atomAndBondVecs[bestBondEnd].second[bestBondNumber];
    if (bestBond->getBeginAtom() != atomAndBondVecs[bestBondEnd].first) {
      bestBond->setEndAtom(bestBond->getBeginAtom());
      bestBond->setBeginAtom(atomAndBondVecs[bestBondEnd].first);
    }
    bestBond->setBondDir(bestBondDir);

    auto newWedgeInfo = std::unique_ptr<RDKit::Chirality::WedgeInfoBase>(
        new RDKit::Chirality::WedgeInfoAtropisomer(bond->getIdx(),
                                                   bestBondDir));
    wedgeBonds[bestBond->getIdx()] = std::move(newWedgeInfo);

  } else {
    BOOST_LOG(rdWarningLog)
        << "Failed to find a good bond to set as UP or DOWN for an atropisomer - atoms are: "
        << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    return false;
  }

  return true;
}

bool WedgeBondFromAtropisomerOneBond3d(
    Bond *bond, const ROMol &mol, const Conformer *conf,
    std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
        &wedgeBonds) {
  PRECONDITION(bond, "bad bond");

  AtropAtomAndBondVec atomAndBondVecs[2];
  if (!getAtropisomerAtomsAndBonds(bond, atomAndBondVecs, mol)) {
    return false;  // not an atropisomer
  }

  //  make sure we do not have wiggle bonds

  for (auto atomAndBondVecs : atomAndBondVecs) {
    for (auto endBond : atomAndBondVecs.second) {
      if (endBond->getBondDir() == Bond::UNKNOWN) {
        return false;  // not an atropisomer)
      }
    }
  }

  // first see if any candidate bond is already set to a wedge or hash
  // if so, we will use that bond as a wedge or hash

  std::vector<Bond *> useBonds;

  for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    for (unsigned int whichBond = 0;
         whichBond < atomAndBondVecs[whichEnd].second.size(); ++whichBond) {
      auto bond = atomAndBondVecs[whichEnd].second[whichBond];
      auto bondDir = bond->getBondDir();

      // see if it is a wedge or hash and its origin is the atom in the
      // main bond

      if ((bondDir == Bond::BEGINWEDGE || bondDir == Bond::BEGINDASH) &&
          bond->getBeginAtom() == atomAndBondVecs[whichEnd].first &&
          canHaveDirection(*bond)) {
        useBonds.push_back(bond);
      }
    }
  }

  // the following may seem redundant, since we just found the useBonds
  // based on their bond dir PRESENCE, but this endures that the values are
  // correct.

  if (useBonds.size() > 0) {
    for (auto useBond : useBonds) {
      useBond->setBondDir(getBondDirForAtropisomer3d(useBond, conf));
    }

    return true;
  }

  // did not find a used bond dir - pick one to use
  // we would like to have one that is not in a ring, and will be a dash

  const RingInfo *ri = bond->getOwningMol().getRingInfo();

  Bond *bestBond = nullptr;
  int bestBondEnd = -1;
  unsigned int bestRingCount = UINT_MAX;
  unsigned int largestRingSize = 0;
  Bond::BondDir bestBondDir = Bond::BondDir::NONE;
  bool bestBondIsSingle = false;
  for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    for (unsigned int whichBond = 0;
         whichBond < atomAndBondVecs[whichEnd].second.size(); ++whichBond) {
      auto bondToTry = atomAndBondVecs[whichEnd].second[whichBond];

      // cannot use a bond that is not single, nor if it is already slated
      // to be used for a chiral center

      if (!canHaveDirection(*bondToTry) ||
          wedgeBonds.find(bond->getIdx()) != wedgeBonds.end()) {
        continue;  // must be a single bond and not already spoken
                   // for by a chiral center
      }

      // make sure the atoms on the bond are in the right order for the
      // wedge/hash the atom on the end of the main bond must be listed
      // first

      if (bondToTry->getBeginAtom() != atomAndBondVecs[whichEnd].first) {
        bondToTry->setEndAtom(bondToTry->getBeginAtom());
        bondToTry->setBeginAtom(atomAndBondVecs[whichEnd].first);
      }

      if (bondToTry->getBondDir() != Bond::BondDir::NONE) {
        if (bondToTry->getBeginAtom()->getIdx() ==
            atomAndBondVecs[whichEnd].first->getIdx()) {
          BOOST_LOG(rdWarningLog)
              << "Wedge or hash bond found on atropisomer where not expected - atoms are: "
              << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
              << std::endl;
          return false;
        } else {
          continue;  // wedge or hash bond affecting the OTHER atom
                     // = perhaps a chiral center
        }
      }
      auto ringCount = ri->numBondRings(bondToTry->getIdx());
      unsigned int ringSize = 0;
      if (!ringCount) {
        ringCount = 10;
      } else {
        // we're going to prefer to put wedges in larger rings, but don't want
        // to end up wedging macrocyles if it's avoidable.
        ringSize = ri->minBondRingSize(bondToTry->getIdx());
        if (ringSize > 8) {
          ringSize = 0;
        }
      }
      if (ringCount > bestRingCount) {
        continue;
      } else if (ringCount < bestRingCount || ringSize > largestRingSize) {
        bestBond = bondToTry;
        bestBondEnd = whichEnd;
        bestRingCount = ringCount;
        largestRingSize = ringSize;
        bestBondIsSingle = (bondToTry->getBondType() == Bond::BondType::SINGLE);
        bestBondDir = getBondDirForAtropisomer3d(bondToTry, conf);
      } else if (bestBondIsSingle &&
                 bondToTry->getBondType() != Bond::BondType::SINGLE) {
        continue;
      } else if (!bestBondIsSingle &&
                 bondToTry->getBondType() == Bond::BondType::SINGLE) {
        bestBondEnd = whichEnd;
        bestBond = bondToTry;
        bestRingCount = ringCount;
        bestBondIsSingle = true;
        bestBondDir = getBondDirForAtropisomer3d(bondToTry, conf);
      } else {
        auto bondDir = getBondDirForAtropisomer3d(bondToTry, conf);
        if (bestBondDir == Bond::BondDir::NONE ||
            (bestBondDir == Bond::BondDir::BEGINDASH &&
             bondDir == Bond::BondDir::BEGINWEDGE)) {
          bestBond = bondToTry;
          bestBondEnd = whichEnd;
          bestRingCount = ringCount;
          bestBondIsSingle =
              (bondToTry->getBondType() == Bond::BondType::SINGLE);

          bestBondDir = bondDir;
        }
      }
    }
  }

  if (bestBond != nullptr) {
    // we found a good one

    // make sure the atoms on the bond are in the right order for the
    // wedge/hash the atom on the end of the main bond must be listed
    // first for the wedge/has bond

    if (bestBond->getBeginAtom() != atomAndBondVecs[bestBondEnd].first) {
      bestBond->setEndAtom(bestBond->getBeginAtom());
      bestBond->setBeginAtom(atomAndBondVecs[bestBondEnd].first);
    }
    bestBond->setBondDir(bestBondDir);
    auto newWedgeInfo = std::unique_ptr<RDKit::Chirality::WedgeInfoBase>(
        new RDKit::Chirality::WedgeInfoAtropisomer(bond->getIdx(),
                                                   bestBondDir));

    wedgeBonds[bestBond->getIdx()] = std::move(newWedgeInfo);
  } else {
    BOOST_LOG(rdWarningLog)
        << "Failed to find a good bond to set as UP or DOWN for an atropisomer - atoms are: "
        << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    return false;
  }

  return true;
}

void wedgeBondsFromAtropisomers(
    const ROMol &mol, const Conformer *conf,
    std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
        &wedgeBonds) {
  PRECONDITION(conf == nullptr || &(conf->getOwningMol()) == &mol,
               "conformer does not belong to molecule");

  // WedgeBondFromAtropisomerOneBond 2d/3d requires ring bond counts
  if (!mol.getRingInfo()->isSssrOrBetter()) {
    RDKit::MolOps::findSSSR(mol);
  }

  for (auto bond : mol.bonds()) {
    auto bondStereo = bond->getStereo();

    if (bond->getBondType() != Bond::BondType::SINGLE ||
        (bondStereo != Bond::BondStereo::STEREOATROPCW &&
         bondStereo != Bond::BondStereo::STEREOATROPCCW) ||
        bond->getBeginAtom()->getTotalDegree() < 2 ||
        bond->getEndAtom()->getTotalDegree() < 2 ||
        bond->getBeginAtom()->getTotalDegree() > 3 ||
        bond->getEndAtom()->getTotalDegree() > 3) {
      continue;
    }

    if (conf) {
      if (conf->is3D()) {
        WedgeBondFromAtropisomerOneBond3d(bond, mol, conf, wedgeBonds);
      } else {
        WedgeBondFromAtropisomerOneBond2d(bond, mol, conf, wedgeBonds);
      }
    } else {  // no conformer
      WedgeBondFromAtropisomerOneBondNoConf(bond, mol, wedgeBonds);
    }
  }
}

bool doesMolHaveAtropisomers(const ROMol &mol) {
  for (auto bond : mol.bonds()) {
    auto bondStereo = bond->getStereo();

    if (bondStereo == Bond::BondStereo::STEREOATROPCW ||
        bondStereo == Bond::BondStereo::STEREOATROPCCW) {
      return true;
    }
  }
  return false;
}
}  // namespace Atropisomers
}  // namespace RDKit
