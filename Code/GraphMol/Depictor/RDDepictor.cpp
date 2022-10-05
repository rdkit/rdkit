//
//  Copyright (C) 2003-2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "RDDepictor.h"
#include "EmbeddedFrag.h"

#ifdef RDK_BUILD_COORDGEN_SUPPORT
#include <CoordGen/CoordGen.h>
#endif

#include <RDGeneral/types.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/Chirality.h>
#include <cmath>
#include <GraphMol/MolOps.h>
#include <GraphMol/Rings.h>
#include <Geometry/point.h>
#include <Geometry/Transform2D.h>
#include <Geometry/Transform3D.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/Substruct/SubstructUtils.h>
#include "EmbeddedFrag.h"
#include "DepictUtils.h"
#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include <algorithm>

namespace RDDepict {

bool preferCoordGen = false;

namespace DepictorLocal {

constexpr auto ISQRT2 = 0.707107;
constexpr auto SQRT3_2 = 0.866025;

std::vector<const RDKit::Atom *> getRankedAtomNeighbors(
    const RDKit::ROMol &mol, const RDKit::Atom *atom,
    const std::vector<int> &atomRanks) {
  std::vector<const RDKit::Atom *> nbrs;
  for (auto nbr : mol.atomNeighbors(atom)) {
    nbrs.push_back(nbr);
  }
  std::sort(nbrs.begin(), nbrs.end(),
            [&atomRanks](const auto e1, const auto e2) {
              return atomRanks[e1->getIdx()] < atomRanks[e2->getIdx()];
            });
  return nbrs;
}

void embedSquarePlanar(const RDKit::ROMol &mol, const RDKit::Atom *atom,
                       std::list<EmbeddedFrag> &efrags,
                       const std::vector<int> &atomRanks) {
  static const RDGeom::Point2D idealPoints[] = {
      RDGeom::Point2D(ISQRT2 * BOND_LEN, ISQRT2 * BOND_LEN),
      RDGeom::Point2D(ISQRT2 * BOND_LEN, -ISQRT2 * BOND_LEN),
      RDGeom::Point2D(-ISQRT2 * BOND_LEN, -ISQRT2 * BOND_LEN),
      RDGeom::Point2D(-ISQRT2 * BOND_LEN, ISQRT2 * BOND_LEN),
  };
  PRECONDITION(atom, "bad atom");
  if (atom->getChiralTag() != RDKit::Atom::ChiralType::CHI_SQUAREPLANAR) {
    return;
  }
  auto nbrs = getRankedAtomNeighbors(mol, atom, atomRanks);
  RDGeom::INT_POINT2D_MAP coordMap;
  coordMap[atom->getIdx()] = RDGeom::Point2D(0., 0.);
  coordMap[nbrs[0]->getIdx()] = idealPoints[0];
  bool q2Full = false;
  for (const auto nbr : nbrs) {
    if (nbr == nbrs.front()) {
      continue;
    }
    auto angle =
        RDKit::Chirality::getIdealAngleBetweenLigands(atom, nbrs.front(), nbr);
    if (fabs(angle - 180) < 0.1) {
      coordMap[nbr->getIdx()] = idealPoints[2];
    } else {
      if (!q2Full) {
        coordMap[nbr->getIdx()] = idealPoints[1];
        q2Full = true;
      } else {
        coordMap[nbr->getIdx()] = idealPoints[3];
      }
    }
  }
  efrags.emplace_back(&mol, coordMap);
}

void embedTBP(const RDKit::ROMol &mol, const RDKit::Atom *atom,
              std::list<EmbeddedFrag> &efrags,
              const std::vector<int> &atomRanks) {
  static const RDGeom::Point2D idealPoints[] = {
      RDGeom::Point2D(0, BOND_LEN),                        // axial
      RDGeom::Point2D(0, -BOND_LEN),                       // axial
      RDGeom::Point2D(-SQRT3_2 * BOND_LEN, BOND_LEN / 2),  // equatorial
      RDGeom::Point2D(-SQRT3_2 * BOND_LEN,
                      -BOND_LEN / 2),  // equatorial
      RDGeom::Point2D(BOND_LEN, 0),    // equatorial
  };
  PRECONDITION(atom, "bad atom");
  if (atom->getChiralTag() !=
      RDKit::Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL) {
    return;
  }
  auto nbrs = getRankedAtomNeighbors(mol, atom, atomRanks);
  RDGeom::INT_POINT2D_MAP coordMap;
  coordMap[atom->getIdx()] = RDGeom::Point2D(0., 0.);
  const RDKit::Atom *axial1 =
      RDKit::Chirality::getTrigonalBipyramidalAxialAtom(atom);
  const RDKit::Atom *axial2 =
      RDKit::Chirality::getTrigonalBipyramidalAxialAtom(atom, -1);
  if (axial1) {
    coordMap[axial1->getIdx()] = idealPoints[0];
  }
  if (axial2) {
    coordMap[axial2->getIdx()] = idealPoints[1];
  }
  unsigned whichEq = 2;
  for (const auto nbr : nbrs) {
    if (nbr != axial1 && nbr != axial2) {
      coordMap[nbr->getIdx()] = idealPoints[whichEq++];
    }
  }
  efrags.emplace_back(&mol, coordMap);
}

void embedOctahedral(const RDKit::ROMol &mol, const RDKit::Atom *atom,
                     std::list<EmbeddedFrag> &efrags,
                     const std::vector<int> &atomRanks) {
  static const RDGeom::Point2D idealPoints[] = {
      RDGeom::Point2D(0, BOND_LEN),                         // axial
      RDGeom::Point2D(0, -BOND_LEN),                        // axial
      RDGeom::Point2D(SQRT3_2 * BOND_LEN, BOND_LEN / 2),    // equatorial
      RDGeom::Point2D(SQRT3_2 * BOND_LEN, -BOND_LEN / 2),   // equatorial
      RDGeom::Point2D(-SQRT3_2 * BOND_LEN, -BOND_LEN / 2),  // equatorial
      RDGeom::Point2D(-SQRT3_2 * BOND_LEN, BOND_LEN / 2),   // equatorial
  };
  PRECONDITION(atom, "bad atom");
  if (atom->getChiralTag() != RDKit::Atom::ChiralType::CHI_OCTAHEDRAL) {
    return;
  }
  auto nbrs = getRankedAtomNeighbors(mol, atom, atomRanks);
  RDGeom::INT_POINT2D_MAP coordMap;
  coordMap[atom->getIdx()] = RDGeom::Point2D(0., 0.);
  const RDKit::Atom *axial1 = nullptr;
  const RDKit::Atom *axial2 = nullptr;
  for (auto i = 0u; i < nbrs.size(); ++i) {
    bool all90 = true;
    for (auto j = i + 1; j < nbrs.size(); ++j) {
      if (fabs(RDKit::Chirality::getIdealAngleBetweenLigands(atom, nbrs[i],
                                                             nbrs[j]) -
               180) < 0.1) {
        axial1 = nbrs[i];
        axial2 = nbrs[j];
        all90 = false;
        break;
      } else if (fabs(RDKit::Chirality::getIdealAngleBetweenLigands(
                          atom, nbrs[i], nbrs[j]) -
                      90) > 0.1) {
        all90 = false;
      }
    }
    if (all90) {
      axial1 = nbrs[i];
    }
    if (axial1) {
      break;
    }
  }
  if (axial1) {
    coordMap[axial1->getIdx()] = idealPoints[0];
  }
  if (axial2) {
    coordMap[axial2->getIdx()] = idealPoints[1];
  }
  const RDKit::Atom *refEqAtom1 = nullptr;
  const RDKit::Atom *refEqAtom2 = nullptr;
  for (const auto nbr : nbrs) {
    if (nbr != axial1 && nbr != axial2) {
      if (!refEqAtom1) {
        refEqAtom1 = nbr;
        coordMap[nbr->getIdx()] = idealPoints[2];
        refEqAtom2 = RDKit::Chirality::getChiralAcrossAtom(atom, nbr);
        if (refEqAtom2) {
          coordMap[refEqAtom2->getIdx()] = idealPoints[4];
        }
      } else {
        if (nbr == refEqAtom2 || nbr == refEqAtom1) {
          continue;
        }
        coordMap[nbr->getIdx()] = idealPoints[3];
        const auto acrossAtom2 =
            RDKit::Chirality::getChiralAcrossAtom(atom, nbr);
        if (acrossAtom2) {
          coordMap[acrossAtom2->getIdx()] = idealPoints[5];
        }
        break;
      }
    }
  }
  efrags.emplace_back(&mol, coordMap);
}

void embedNontetrahedralStereo(const RDKit::ROMol &mol,
                               std::list<EmbeddedFrag> &efrags,
                               const std::vector<int> &atomRanks) {
  boost::dynamic_bitset<> consider(mol.getNumAtoms());
  for (const auto atm : mol.atoms()) {
    if (RDKit::Chirality::hasNonTetrahedralStereo(atm)) {
      consider[atm->getIdx()] = 1;
    }
  }
  if (consider.empty()) {
    return;
  }
  for (const auto atm : mol.atoms()) {
    if (!consider[atm->getIdx()]) {
      continue;
    }
    switch (atm->getChiralTag()) {
      case RDKit::Atom::ChiralType::CHI_SQUAREPLANAR:
        embedSquarePlanar(mol, atm, efrags, atomRanks);
        break;
      case RDKit::Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL:
        embedTBP(mol, atm, efrags, atomRanks);
        break;
      case RDKit::Atom::ChiralType::CHI_OCTAHEDRAL:
        embedOctahedral(mol, atm, efrags, atomRanks);
        break;
      default:
        break;
    }
  }
}

// arings: indices of atoms in rings
void embedFusedSystems(const RDKit::ROMol &mol,
                       const RDKit::VECT_INT_VECT &arings,
                       std::list<EmbeddedFrag> &efrags) {
  RDKit::INT_INT_VECT_MAP neighMap;
  RingUtils::makeRingNeighborMap(arings, neighMap);

  auto cnrs = arings.size();
  boost::dynamic_bitset<> fusDone(cnrs);

  auto curr = 0u;
  while (curr < cnrs) {
    // embed all ring and fused ring systems
    RDKit::INT_VECT fused;
    RingUtils::pickFusedRings(curr, neighMap, fused, fusDone);
    RDKit::VECT_INT_VECT frings;
    frings.reserve(fused.size());
    for (auto rid : fused) {
      frings.push_back(arings.at(rid));
    }
    EmbeddedFrag efrag(&mol, frings);
    efrag.setupNewNeighs();
    efrags.push_back(efrag);
    size_t rix;
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

void embedCisTransSystems(const RDKit::ROMol &mol,
                          std::list<EmbeddedFrag> &efrags) {
  for (auto bond : mol.bonds()) {
    // check if this bond is in a cis/trans double bond
    // and it is not a ring bond
    if ((bond->getBondType() == RDKit::Bond::DOUBLE)  // this is a double bond
        && (bond->getStereo() >
            RDKit::Bond::STEREOANY)  // and has stereo chemistry specified
        && (!bond->getOwningMol().getRingInfo()->numBondRings(
               bond->getIdx()))) {  // not in a ring
      if (bond->getStereoAtoms().size() != 2) {
        BOOST_LOG(rdWarningLog)
            << "WARNING: bond found with stereo spec but no stereo atoms"
            << std::endl;
        continue;
      }
      EmbeddedFrag efrag(bond);
      efrag.setupNewNeighs();
      efrags.push_back(efrag);
    }
  }
}

RDKit::INT_LIST getNonEmbeddedAtoms(const RDKit::ROMol &mol,
                                    const std::list<EmbeddedFrag> &efrags) {
  RDKit::INT_LIST res;
  boost::dynamic_bitset<> done(mol.getNumAtoms());
  for (const auto &efrag : efrags) {
    const auto &oatoms = efrag.GetEmbeddedAtoms();
    for (const auto &oatom : oatoms) {
      done[oatom.first] = 1;
    }
  }
  for (auto aid = 0u; aid < mol.getNumAtoms(); ++aid) {
    if (!done[aid]) {
      res.push_back(aid);
    }
  }
  return res;
}

// find the largest fragments that is not done yet (
//  i.e. merged with the master fragments)
// if do not find anything we return efrags.end()
std::list<EmbeddedFrag>::iterator _findLargestFrag(
    std::list<EmbeddedFrag> &efrags) {
  std::list<EmbeddedFrag>::iterator mfri;
  int msiz = 0;
  for (auto efri = efrags.begin(); efri != efrags.end(); ++efri) {
    if ((!efri->isDone()) && (efri->Size() > msiz)) {
      msiz = efri->Size();
      mfri = efri;
    }
  }
  if (msiz == 0) {
    mfri = efrags.end();
  }
  return mfri;
}

void _shiftCoords(std::list<EmbeddedFrag> &efrags) {
  // shift the coordinates if there are multiple fragments
  // so that the fragments do not overlap each other
  if (efrags.empty()) {
    return;
  }
  for (auto &efrag : efrags) {
    efrag.computeBox();
  }
  auto eri = efrags.begin();
  auto xmax = eri->getBoxPx();
  auto xmin = eri->getBoxNx();
  auto ymax = eri->getBoxPy();
  auto ymin = eri->getBoxNy();

  ++eri;
  while (eri != efrags.end()) {
    bool xshift = true;

    if (xmax + xmin > ymax + ymin) {
      xshift = false;
    }
    auto xn = eri->getBoxNx();
    auto xp = eri->getBoxPx();
    auto yn = eri->getBoxNy();
    auto yp = eri->getBoxPy();
    RDGeom::Point2D shift(0.0, 0.0);
    if (xshift) {
      shift.x = xmax + xn + 1.0;
      shift.y = 0.0;
      xmax += xp + xn + 1.0;
    } else {
      shift.x = 0.0;
      shift.y = ymax + yn + 1.0;
      ymax += yp + yn + 1.0;
    }
    eri->Translate(shift);

    ++eri;
  }
}

// we do not use std::copysign as we need a tolerance
double copySign(double to, double from, double tol) {
  return (from < -tol ? -fabs(to) : fabs(to));
}

struct ThetaBin {
  double d_thetaAvg = 0.0;
  std::vector<double> thetaValues;
};
}  // namespace DepictorLocal

void computeInitialCoords(RDKit::ROMol &mol,
                          const RDGeom::INT_POINT2D_MAP *coordMap,
                          std::list<EmbeddedFrag> &efrags) {
  std::vector<int> atomRanks;
  atomRanks.resize(mol.getNumAtoms());
  for (auto i = 0u; i < mol.getNumAtoms(); ++i) {
    atomRanks[i] = getAtomDepictRank(mol.getAtomWithIdx(i));
  }
  RDKit::VECT_INT_VECT arings;

  // first find all the rings
  RDKit::MolOps::symmetrizeSSSR(mol, arings);

  // do stereochemistry
  RDKit::MolOps::assignStereochemistry(mol, false);

  efrags.clear();

  // user-specified coordinates exist
  bool preSpec = false;
  // first embed any atoms for which the coordinates have been specified.
  if ((coordMap) && (coordMap->size() > 1)) {
    EmbeddedFrag efrag(&mol, *coordMap);
    // add this to the list of embedded fragments
    efrags.push_back(efrag);
    preSpec = true;
  }

  if (arings.size() > 0) {
    // first deal with the fused rings
    DepictorLocal::embedFusedSystems(mol, arings, efrags);
  }

  // do non-tetrahedral stereo
  DepictorLocal::embedNontetrahedralStereo(mol, efrags, atomRanks);

  // deal with any cis/trans systems
  DepictorLocal::embedCisTransSystems(mol, efrags);
  // now get the atoms that are not yet embedded in either a cis/trans system
  // or a ring system (or simply the first atom)
  auto nratms = DepictorLocal::getNonEmbeddedAtoms(mol, efrags);
  std::list<EmbeddedFrag>::iterator mri;
  if (preSpec) {
    // if the user specified coordinates on some of the atoms use that as
    // as the starting fragment and it should be at the beginning of the vector
    mri = efrags.begin();
  } else {
    // otherwise - find the largest fragment that was embedded
    mri = DepictorLocal::_findLargestFrag(efrags);
  }

  while ((mri != efrags.end()) || (nratms.size() > 0)) {
    if (mri == efrags.end()) {
      // we are out of embedded fragments, if there are any
      // non embedded atoms use them to start a fragment
      auto mrank = static_cast<int>(RDKit::MAX_INT);
      RDKit::INT_LIST_I mnri;
      for (auto nri = nratms.begin(); nri != nratms.end(); ++nri) {
        auto rank = atomRanks.at(*nri);
        rank *= mol.getNumAtoms();
        // use the atom index as well so that we at least
        // get reproducible depictions in cases where things
        // have identical ranks.
        rank += *nri;
        if (rank < mrank) {
          mrank = rank;
          mnri = nri;
        }
      }
      EmbeddedFrag efrag((*mnri), &mol);
      nratms.erase(mnri);
      efrags.push_back(efrag);
      mri = efrags.end();
      --mri;
    }
    mri->markDone();
    mri->expandEfrag(nratms, efrags);
    mri = DepictorLocal::_findLargestFrag(efrags);
  }
  // at this point any remaining efrags should belong individual fragments in
  // the molecule
}

unsigned int copyCoordinate(RDKit::ROMol &mol, std::list<EmbeddedFrag> &efrags,
                            bool clearConfs) {
  // create a conformation to store the coordinates and add it to the molecule
  auto *conf = new RDKit::Conformer(mol.getNumAtoms());
  conf->set3D(false);
  std::list<EmbeddedFrag>::iterator eri;
  for (const auto &efrag : efrags) {
    for (const auto &eai : efrag.GetEmbeddedAtoms()) {
      const auto &cr = eai.second.loc;
      RDGeom::Point3D fcr(cr.x, cr.y, 0.0);
      conf->setAtomPos(eai.first, fcr);
    }
  }
  unsigned int confId = 0;
  if (clearConfs) {
    // clear all the conformation on the molecules and assign conf ID 0 to this
    // conformation
    mol.clearConformers();
    conf->setId(confId);
    // conf ID has already been set in this case to 0 - not other
    // confs on the molecule at this point
    mol.addConformer(conf);
  } else {
    // let add conf assign a conformation ID for the conformation
    confId = mol.addConformer(conf, true);
  }
  return confId;
}
//
//
// 50,000 foot algorithm:
//   1) Find rings
//   2) Find fused systems
//   3) embed largest fused system
//   4) for each unfinished atom:
//      1) find neighbors
//      2) if neighbor is non-ring atom, embed it; otherwise merge the
//         ring system
//      3) add all atoms just merged/embedded to unfinished atom list
//
//
unsigned int compute2DCoords(RDKit::ROMol &mol,
                             const RDGeom::INT_POINT2D_MAP *coordMap,
                             bool canonOrient, bool clearConfs,
                             unsigned int nFlipsPerSample,
                             unsigned int nSamples, int sampleSeed,
                             bool permuteDeg4Nodes, bool forceRDKit) {
  if (mol.needsUpdatePropertyCache()) {
    mol.updatePropertyCache(false);
  }
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  // default to use CoordGen if we have it installed
  if (!forceRDKit && preferCoordGen) {
    RDKit::CoordGen::CoordGenParams params;
    if (coordMap) {
      params.coordMap = *coordMap;
    }
    auto cid = RDKit::CoordGen::addCoords(mol, &params);
    return cid;
  };
#endif
  // storage for pieces of a molecule/s that are embedded in 2D
  std::list<EmbeddedFrag> efrags;
  computeInitialCoords(mol, coordMap, efrags);

#if 1
  // perform random sampling here to improve the density
  for (auto &eri : efrags) {
    // either sample the 2D space by randomly flipping rotatable
    // bonds in the structure or flip only bonds along the shortest
    // path between colliding atoms - don't do both
    if ((nSamples > 0) && (nFlipsPerSample > 0)) {
      eri.randomSampleFlipsAndPermutations(nFlipsPerSample, nSamples,
                                           sampleSeed, nullptr, 0.0,
                                           permuteDeg4Nodes);
    } else {
      eri.removeCollisionsBondFlip();
    }
  }
  for (auto &eri : efrags) {
    // if there are any remaining collisions
    eri.removeCollisionsOpenAngles();
    eri.removeCollisionsShortenBonds();
  }
  if (!coordMap || !coordMap->size()) {
    if (canonOrient && efrags.size()) {
      // if we do not have any prespecified coordinates - canonicalize
      // the orientation of the fragment so that the longest axes fall
      // along the x-axis etc.
      for (auto &eri : efrags) {
        eri.canonicalizeOrientation();
      }
    }
  }
  DepictorLocal::_shiftCoords(efrags);
#endif
  // create a conformation on the molecule and copy the coordinates
  auto cid = copyCoordinate(mol, efrags, clearConfs);

  // special case for a single-atom coordMap template
  if ((coordMap) && (coordMap->size() == 1)) {
    auto &conf = mol.getConformer(cid);
    auto cRef = coordMap->begin();
    const auto &confPos = conf.getAtomPos(cRef->first);
    auto refPos = cRef->second;
    refPos.x -= confPos.x;
    refPos.y -= confPos.y;
    for (auto i = 0u; i < conf.getNumAtoms(); ++i) {
      auto confPos = conf.getAtomPos(i);
      confPos.x += refPos.x;
      confPos.y += refPos.y;
      conf.setAtomPos(i, confPos);
    }
  }

  return cid;
}

//! \brief Compute the 2D coordinates such that the interatom distances
//!        mimic those in a distance matrix
/*!
  This function generates 2D coordinates such that the inter atom
  distance mimic those specified via dmat. This is done by randomly
  sampling(flipping) the rotatable bonds in the molecule and
  evaluating a cost function which contains two components. The
  first component is the sum of inverse of the squared inter-atom
  distances, this helps in spreading the atoms far from each
  other. The second component is the sum of squares of the
  difference in distance between those in dmat and the generated
  structure.  The user can adjust the relative importance of the two
  components via an adjustable parameter (see below)

  ARGUMENTS:
  \param mol - molecule involved in the fragment

  \param dmat - the distance matrix we want to mimic, this is
                symmetric N by N matrix when N is the number of
                atoms in mol. All negative entries in dmat are
                ignored.

  \param canonOrient - canonicalize the orientation after the 2D
                       embedding is done

  \param clearConfs - clear any previously existing conformations on
                      mol before adding a conformation

  \param weightDistMat - A value between 0.0 and 1.0, this
                         determines the importance of mimicking the
                         inter atoms distances in dmat. (1.0 -
                         weightDistMat) is the weight associated to
                         spreading out the structure (density) in
                         the cost function

  \param nFlipsPerSample - the number of rotatable bonds that are
                           randomly flipped for each sample

  \param nSample - the number of samples

  \param sampleSeed - seed for the random sampling process
*/
unsigned int compute2DCoordsMimicDistMat(
    RDKit::ROMol &mol, const DOUBLE_SMART_PTR *dmat, bool canonOrient,
    bool clearConfs, double weightDistMat, unsigned int nFlipsPerSample,
    unsigned int nSamples, int sampleSeed, bool permuteDeg4Nodes, bool) {
  // storage for pieces of a molecule/s that are embedded in 2D
  std::list<EmbeddedFrag> efrags;
  computeInitialCoords(mol, nullptr, efrags);

  // now perform random flips of rotatable bonds so that we can sample the space
  // and try to mimic the distances in dmat
  std::list<EmbeddedFrag>::iterator eri;
  for (auto &eri : efrags) {
    eri.randomSampleFlipsAndPermutations(nFlipsPerSample, nSamples, sampleSeed,
                                         dmat, weightDistMat, permuteDeg4Nodes);
  }
  if (canonOrient && efrags.size()) {
    // canonicalize the orientation of the fragment so that the
    // longest axes fall along the x-axis etc.
    for (auto &eri : efrags) {
      eri.canonicalizeOrientation();
    }
  }

  DepictorLocal::_shiftCoords(efrags);
  // create a conformation on the molecule and copy the coordinates
  return copyCoordinate(mol, efrags, clearConfs);
}

//! \brief Compute 2D coordinates where a piece of the molecule is
//   constrained to have the same coordinates as a reference;
//   correspondences between reference and molecule atom indices
//   are determined by refMatchVect
void generateDepictionMatching2DStructure(
    RDKit::ROMol &mol, const RDKit::ROMol &reference,
    const RDKit::MatchVectType &refMatchVect, int confId, bool forceRDKit) {
  if (refMatchVect.size() > reference.getNumAtoms()) {
    throw RDDepict::DepictException(
        "When a refMatchVect is provided, it must have size "
        "<= number of atoms in the reference");
  }
  RDGeom::INT_POINT2D_MAP coordMap;
  const RDKit::Conformer &conf = reference.getConformer(confId);
  for (const auto &mv : refMatchVect) {
    if (mv.first > static_cast<int>(reference.getNumAtoms())) {
      throw RDDepict::DepictException(
          "Reference atom index in refMatchVect out of range");
    }
    if (mv.second > static_cast<int>(mol.getNumAtoms())) {
      throw RDDepict::DepictException(
          "Molecule atom index in refMatchVect out of range");
    }
    const auto &pt3 = conf.getAtomPos(mv.first);
    RDGeom::Point2D pt2(pt3.x, pt3.y);
    coordMap[mv.second] = pt2;
  }
  RDDepict::compute2DCoords(mol, &coordMap, false /* canonOrient */,
                            true /* clearConfs */, 0, 0, 0, false, forceRDKit);
}

//! \brief Compute 2D coordinates where a piece of the molecule is
//   constrained to have the same coordinates as a reference.
RDKit::MatchVectType generateDepictionMatching2DStructure(
    RDKit::ROMol &mol, const RDKit::ROMol &reference, int confId,
    const RDKit::ROMol *referencePattern, bool acceptFailure, bool forceRDKit,
    bool allowOptionalAttachments) {
  std::unique_ptr<RDKit::ROMol> referenceHs;
  std::vector<int> refMatch;
  RDKit::MatchVectType matchVect;
  std::vector<RDKit::MatchVectType> multiRefMatchVect;
  RDKit::MatchVectType singleRefMatchVect;
  auto &refMatchVectRef = singleRefMatchVect;
  const RDKit::ROMol &query =
      (referencePattern ? *referencePattern : reference);
  if (allowOptionalAttachments) {
    // we do not need the allowOptionalAttachments logic if there are no
    // terminal dummy atoms
    allowOptionalAttachments = false;
    for (const auto queryAtom : query.atoms()) {
      if (queryAtom->getAtomicNum() == 0 && queryAtom->getDegree() == 1) {
        allowOptionalAttachments = true;
        break;
      }
    }
  }
  if (referencePattern) {
    if (allowOptionalAttachments &&
        referencePattern->getNumAtoms() > reference.getNumAtoms()) {
      referenceHs.reset(RDKit::MolOps::addHs(reference));
      CHECK_INVARIANT(referenceHs, "addHs returned a nullptr");
      multiRefMatchVect =
          RDKit::SubstructMatch(*referenceHs, *referencePattern);
      if (!multiRefMatchVect.empty()) {
        refMatchVectRef = RDKit::getMostSubstitutedCoreMatch(
            *referenceHs, *referencePattern, multiRefMatchVect);
      }
    } else if (referencePattern->getNumAtoms() <= reference.getNumAtoms()) {
      RDKit::SubstructMatch(reference, *referencePattern, singleRefMatchVect);
    }
    if (refMatchVectRef.empty()) {
      throw RDDepict::DepictException(
          "Reference pattern does not map to reference.");
    }
    refMatch.resize(query.getNumAtoms(), -1);
    for (auto &i : refMatchVectRef) {
      // skip indices corresponding to added Hs
      if (allowOptionalAttachments &&
          referenceHs->getAtomWithIdx(i.second)->getAtomicNum() == 1) {
        continue;
      }
      refMatch[i.first] = i.second;
    }
  } else {
    refMatch.resize(reference.getNumAtoms());
    std::iota(refMatch.begin(), refMatch.end(), 0);
  }
  if (allowOptionalAttachments) {
    std::unique_ptr<RDKit::ROMol> molHs(RDKit::MolOps::addHs(mol));
    CHECK_INVARIANT(molHs, "addHs returned a nullptr");
    auto matches = SubstructMatch(*molHs, query);
    if (matches.empty()) {
      allowOptionalAttachments = false;
    } else {
      for (const auto &pair :
           getMostSubstitutedCoreMatch(*molHs, query, matches)) {
        if (molHs->getAtomWithIdx(pair.second)->getAtomicNum() != 1 &&
            refMatch.at(pair.first) >= 0) {
          matchVect.push_back(pair);
        }
      }
    }
  }
  if (!allowOptionalAttachments) {
    RDKit::SubstructMatch(mol, query, matchVect);
  }
  if (matchVect.empty() && !acceptFailure) {
    throw RDDepict::DepictException(
        "Substructure match with reference not found.");
  }
  for (auto &pair : matchVect) {
    pair.first = refMatch.at(pair.first);
  }
  generateDepictionMatching2DStructure(mol, reference, matchVect, confId,
                                       forceRDKit);
  return matchVect;
}

//! \brief Generate a 2D depiction for a molecule where all or part of
//   it mimics the coordinates of a 3D reference structure.
void generateDepictionMatching3DStructure(RDKit::ROMol &mol,
                                          const RDKit::ROMol &reference,
                                          int confId,
                                          RDKit::ROMol *referencePattern,
                                          bool acceptFailure, bool forceRDKit) {
  auto num_ats = mol.getNumAtoms();
  if (!referencePattern && reference.getNumAtoms() < num_ats) {
    if (acceptFailure) {
      RDDepict::compute2DCoords(mol);
      return;
    } else {
      throw RDDepict::DepictException(
          "Reference molecule not compatible with target molecule.");
    }
  }

  std::vector<int> mol_to_ref(num_ats, -1);
  if (referencePattern && referencePattern->getNumAtoms()) {
    RDKit::MatchVectType molMatchVect, refMatchVect;
    RDKit::SubstructMatch(mol, *referencePattern, molMatchVect);
    RDKit::SubstructMatch(reference, *referencePattern, refMatchVect);
    if (molMatchVect.empty() || refMatchVect.empty()) {
      if (acceptFailure) {
        RDDepict::compute2DCoords(mol);
        return;
      } else {
        throw RDDepict::DepictException(
            "Reference pattern didn't match molecule or reference.");
      }
    }
    for (size_t i = 0; i < molMatchVect.size(); ++i) {
      mol_to_ref[molMatchVect[i].second] = refMatchVect[i].second;
    }

  } else {
    for (unsigned int i = 0; i < num_ats; ++i) {
      mol_to_ref[i] = i;
    }
  }

  const RDKit::Conformer &conf = reference.getConformer(confId);
  // the distance matrix is a triangular representation
  RDDepict::DOUBLE_SMART_PTR dmat(new double[num_ats * (num_ats - 1) / 2]);
  // negative distances are ignored, so initialise to -1.0 so subset by
  // referencePattern works.
  std::fill(dmat.get(), dmat.get() + num_ats * (num_ats - 1) / 2, -1.0);
  for (unsigned int i = 0; i < num_ats; ++i) {
    if (-1 == mol_to_ref[i]) {
      continue;
    }
    RDGeom::Point3D cds_i = conf.getAtomPos(i);
    for (unsigned int j = i + 1; j < num_ats; ++j) {
      if (-1 == mol_to_ref[j]) {
        continue;
      }
      RDGeom::Point3D cds_j = conf.getAtomPos(mol_to_ref[j]);
      dmat[(j * (j - 1) / 2) + i] = (cds_i - cds_j).length();
    }
  }

  RDDepict::compute2DCoordsMimicDistMat(mol, &dmat, false, true, 0.5, 3, 100,
                                        25, true, forceRDKit);
}

void straightenDepiction(RDKit::ROMol &mol, int confId, bool minimizeRotation) {
  if (!mol.getNumBonds()) {
    return;
  }
  constexpr double RAD2DEG = 180. / M_PI;
  constexpr double DEG2RAD = M_PI / 180.;
  constexpr double ALMOST_ZERO = 1.e-5;
  constexpr double INCR_DEG = 30.;
  constexpr double HALF_INCR_DEG = 0.5 * INCR_DEG;
  constexpr double QUARTER_INCR_DEG = 0.25 * INCR_DEG;
  auto &conf = mol.getConformer(confId);
  auto &pos = conf.getPositions();
  std::unordered_map<int, DepictorLocal::ThetaBin> thetaBins;
  for (const auto b : mol.bonds()) {
    auto bi = b->getBeginAtomIdx();
    auto ei = b->getEndAtomIdx();
    auto bv = pos.at(bi) - pos.at(ei);
    bv.x = (bv.x < 0.) ? std::min(-ALMOST_ZERO, bv.x)
                       : std::max(ALMOST_ZERO, bv.x);
    auto theta = RAD2DEG * atan(bv.y / bv.x);
    auto d_theta = fmod(-theta, INCR_DEG);
    if (fabs(d_theta) > HALF_INCR_DEG) {
      d_theta -= DepictorLocal::copySign(INCR_DEG, d_theta, ALMOST_ZERO);
    }
    int thetaKey = static_cast<int>(
        d_theta + DepictorLocal::copySign(0.5, d_theta, ALMOST_ZERO));
    auto &thetaBin = thetaBins[thetaKey];
    thetaBin.d_thetaAvg += d_theta;
    thetaBin.thetaValues.push_back(theta);
  }
  CHECK_INVARIANT(!thetaBins.empty(), "");
  double d_thetaSmallest = std::numeric_limits<double>::max();
  for (auto &it : thetaBins) {
    auto &thetaBin = it.second;
    thetaBin.d_thetaAvg /= static_cast<double>(thetaBin.thetaValues.size());
    if (fabs(thetaBin.d_thetaAvg) < fabs(d_thetaSmallest)) {
      d_thetaSmallest = thetaBin.d_thetaAvg;
    }
  }
  const auto &minRotationBin =
      std::max_element(
          thetaBins.begin(), thetaBins.end(),
          [](const auto &a, const auto &b) {
            const auto &aBin = a.second;
            const auto &bBin = b.second;
            return (aBin.thetaValues.size() < bBin.thetaValues.size() ||
                    (aBin.thetaValues.size() == bBin.thetaValues.size() &&
                     fabs(aBin.d_thetaAvg) > fabs(bBin.d_thetaAvg)));
          })
          ->second;
  double d_thetaMin = minRotationBin.d_thetaAvg;
  // unless we want to preserve as much as possible the initial orientation,
  // we try to orient the molecule such that the majority of bonds have
  // an angle of 30 or 90 degrees with the X axis
  if (!minimizeRotation) {
    unsigned int count60vs30[2] = {0, 0};
    for (auto theta : minRotationBin.thetaValues) {
      theta += d_thetaMin;
      auto idx = static_cast<unsigned int>((fabs(theta) + 0.5) / INCR_DEG) % 2;
      CHECK_INVARIANT(idx < 2, "");
      ++count60vs30[idx];
    }
    if (count60vs30[0] > count60vs30[1]) {
      d_thetaMin -= DepictorLocal::copySign(INCR_DEG, d_thetaMin, ALMOST_ZERO);
    }
  } else if (fabs(d_thetaSmallest) < ALMOST_ZERO ||
             (fabs(d_thetaSmallest) < fabs(d_thetaMin) &&
              fabs(d_thetaMin) > QUARTER_INCR_DEG)) {
    d_thetaMin = d_thetaSmallest;
  }
  if (fabs(d_thetaMin) > ALMOST_ZERO) {
    d_thetaMin *= DEG2RAD;
    RDGeom::Transform3D trans;
    trans.SetRotation(d_thetaMin, RDGeom::Z_Axis);
    MolTransforms::transformConformer(conf, trans);
  }
}

double normalizeDepiction(RDKit::ROMol &mol, int confId, int canonicalize,
                          double scaleFactor) {
  constexpr double SCALE_FACTOR_THRESHOLD = 1.e-5;
  if (!mol.getNumBonds()) {
    return -1.;
  }
  auto &conf = mol.getConformer(confId);
  if (scaleFactor < 0.0) {
    constexpr double RDKIT_BOND_LEN = 1.5;
    int mostCommonBondLengthInt = -1;
    unsigned int maxCount = 0;
    std::unordered_map<int, unsigned int> binnedBondLengths;
    for (const auto b : mol.bonds()) {
      int bondLength =
          static_cast<int>(MolTransforms::getBondLength(
                               conf, b->getBeginAtomIdx(), b->getEndAtomIdx()) *
                               10.0 +
                           0.5);
      auto it = binnedBondLengths.find(bondLength);
      if (it == binnedBondLengths.end()) {
        it = binnedBondLengths.emplace(bondLength, 0U).first;
      }
      ++it->second;
      if (it->second > maxCount) {
        maxCount = it->second;
        mostCommonBondLengthInt = it->first;
      }
    }
    if (mostCommonBondLengthInt > 0) {
      double mostCommonBondLength =
          static_cast<double>(mostCommonBondLengthInt) * 0.1;
      scaleFactor = RDKIT_BOND_LEN / mostCommonBondLength;
    }
  }
  std::unique_ptr<RDGeom::Transform3D> canonTrans;
  if (canonicalize) {
    auto ctd = MolTransforms::computeCentroid(conf);
    canonTrans.reset(MolTransforms::computeCanonicalTransform(conf, &ctd));
    if (canonicalize < 0) {
      RDGeom::Transform3D rotate90;
      rotate90.SetRotation(0., 1., RDGeom::Point3D(0., 0., 1.));
      *canonTrans *= rotate90;
    }
  }
  bool isScaleFactorSane = (scaleFactor > SCALE_FACTOR_THRESHOLD);
  if (isScaleFactorSane && fabs(scaleFactor - 1.0) > SCALE_FACTOR_THRESHOLD) {
    RDGeom::Transform3D trans;
    trans.setVal(0, 0, scaleFactor);
    trans.setVal(1, 1, scaleFactor);
    if (canonTrans) {
      trans *= *canonTrans;
    }
    MolTransforms::transformConformer(conf, trans);
  } else if (canonTrans) {
    MolTransforms::transformConformer(conf, *canonTrans);
  }
  if (!isScaleFactorSane) {
    scaleFactor = -1.;
  }
  return scaleFactor;
}
}  // namespace RDDepict
