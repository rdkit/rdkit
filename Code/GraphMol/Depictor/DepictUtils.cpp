//
//  Copyright (C) 2003-2022 greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/types.h>
#include <cmath>
#include <Geometry/point.h>
#include "DepictUtils.h"
#include <iostream>
#include <RDGeneral/Invariant.h>
#include <GraphMol/Chirality.h>
#include <algorithm>

namespace {
static const char *FORMER_NBR_INDICES = "__formerNbrIndices";
static const char *FORMER_IDX = "__formerIdx";
}  // end anonymous namespace

namespace RDDepict {
double BOND_LEN = 1.5;
double COLLISION_THRES = 0.70;
double BOND_THRES = 0.50;
double ANGLE_OPEN = 0.1222;  // that is about 7 deg
unsigned int MAX_COLL_ITERS = 15;
double HETEROATOM_COLL_SCALE = 1.3;
unsigned int NUM_BONDS_FLIPS = 3;

RDGeom::INT_POINT2D_MAP embedRing(const RDKit::INT_VECT &ring) {
  // The process here is very straight forward
  // we take the center of the ring to lies at the origin put the first
  // point at the origin and then sweep
  // anticlockwise so by an angle A = 360/n for the next point
  // the length of the arm (l) we want to sweep is easy to compute given the
  // bond length (b) we want to use for each bond in the ring (for now
  // we will assume that this bond length is the same for all bonds in the ring
  //  l = b/sqrt(2*(1 - cos(A))
  // the above formula derives from the triangle formula, where side 'c' is
  // given
  // in terms of sides 'a' and 'b' as
  // c = a^2 + b^2 - 2.a.b.cos(A)
  // where A is the angle between a and b

  // compute the sweep angle
  unsigned int na = ring.size();
  double ang = 2 * M_PI / na;

  // compute the arm length
  double al = BOND_LEN / (sqrt(2 * (1 - cos(ang))));

  RDGeom::INT_POINT2D_MAP res;

  for (unsigned int i = 0; i < na; ++i) {
    auto x = al * cos(i * ang);
    auto y = al * sin(i * ang);
    RDGeom::Point2D loc(x, y);
    res[ring[i]] = loc;
  }

  return res;
}

void transformPoints(RDGeom::INT_POINT2D_MAP &nringCor,
                     const RDGeom::Transform2D &trans) {
  std::for_each(nringCor.begin(), nringCor.end(),
                [&trans](auto &elem) { trans.TransformPoint(elem.second); });
}

RDGeom::Point2D computeBisectPoint(const RDGeom::Point2D &rcr, double ang,
                                   const RDGeom::Point2D &nb1,
                                   const RDGeom::Point2D &nb2) {
  RDGeom::Point2D cloc = nb1;
  cloc += nb2;
  cloc *= 0.5;
  if (ang > M_PI) {
    // invert the cloc
    cloc -= rcr;
    cloc *= -1.0;
    cloc += rcr;
  }
  return cloc;
}

RDGeom::Point2D reflectPoint(const RDGeom::Point2D &point,
                             const RDGeom::Point2D &loc1,
                             const RDGeom::Point2D &loc2) {
  RDGeom::Point2D org(0.0, 0.0);
  RDGeom::Point2D xaxis(1.0, 0.0);
  RDGeom::Point2D cent = (loc1 + loc2);
  cent *= 0.5;

  RDGeom::Transform2D trans;
  trans.SetTransform(org, xaxis, cent, loc1);

  /// reverse transform
  RDGeom::Transform2D itrans;
  itrans.SetTransform(cent, loc1, org, xaxis);

  RDGeom::INT_POINT2D_MAP_I nci;
  RDGeom::Point2D res;
  res = point;
  trans.TransformPoint(res);
  res.y = -res.y;
  itrans.TransformPoint(res);
  return res;
}

void reflectPoints(RDGeom::INT_POINT2D_MAP &coordMap,
                   const RDGeom::Point2D &loc1, const RDGeom::Point2D &loc2) {
  std::for_each(coordMap.begin(), coordMap.end(), [&loc1, &loc2](auto &elem) {
    reflectPoint(elem.second, loc1, loc2);
  });
}

RDKit::INT_VECT setNbrOrder(unsigned int aid, const RDKit::INT_VECT &nbrs,
                            const RDKit::ROMol &mol) {
  PRECONDITION(aid < mol.getNumAtoms(), "");
  PR_QUEUE subsAid;
  int ref = -1;
  // find the neighbor of aid that is not in nbrs i.e. atom A from the comments
  // in the header file and the store the pair <degree, aid> in the order of
  // increasing degree
  for (auto anbr : mol.atomNeighbors(mol.getAtomWithIdx(aid))) {
    // We used to use degree here instead we will start using the CIP rank here
    if (std::find(nbrs.begin(), nbrs.end(), static_cast<int>(anbr->getIdx())) ==
        nbrs.end()) {
      ref = anbr->getIdx();
    }
  }

  RDKit::INT_VECT thold = nbrs;
  if (ref >= 0) {
    thold.push_back(ref);
  }
  // we should be here unless we have more than 3 atoms to worry about
  CHECK_INVARIANT(thold.size() > 3, "");
  thold = rankAtomsByRank(mol, thold);

  // swap the position of the 3rd to last and second to last items in sorted
  // list
  unsigned int ln = thold.size();
  int tint = thold[ln - 3];
  thold[ln - 3] = thold[ln - 2];
  thold[ln - 2] = tint;

  // go clock wise along the list from this position for the arranged neighbor
  // list
  RDKit::INT_VECT res;
  res.reserve(thold.size());
  auto pos = std::find(thold.begin(), thold.end(), ref);
  if (pos != thold.end()) {
    res.insert(res.end(), pos + 1, thold.end());
  }
  if (pos != thold.begin()) {
    res.insert(res.end(), thold.begin(), pos);
  }

  POSTCONDITION(res.size() == nbrs.size(), "");
  return res;
}

int pickFirstRingToEmbed(const RDKit::ROMol &mol,
                         const RDKit::VECT_INT_VECT &fusedRings) {
  // ok this is what we will do here
  // we will pick the ring with the smallest number of substituents
  int res = -1;
  unsigned int maxSize = 0;
  int subs, minsubs = static_cast<int>(1e8);
  int cnt = 0;
  for (const auto &fusedRing : fusedRings) {
    subs = 0;
    for (auto rii : fusedRing) {
      if (mol.getAtomWithIdx(rii)->getDegree() > 2) {
        ++subs;
      }
    }
    if (subs < minsubs) {
      res = cnt;
      minsubs = subs;
      maxSize = fusedRing.size();
    } else if (subs == minsubs) {
      if (fusedRing.size() > maxSize) {
        res = cnt;
        maxSize = fusedRing.size();
      }
    }
    cnt++;
  }
  return res;
}

RDKit::INT_VECT findNextRingToEmbed(const RDKit::INT_VECT &doneRings,
                                    const RDKit::VECT_INT_VECT &fusedRings,
                                    int &nextId) {
  // REVIEW: We are changing this after Issue166
  // Originally the ring that have maximum number of atoms in common with the
  // atoms
  // that have already been embedded will be the ring that will get embedded.
  // But
  // if we can find a ring with two atoms in common with the embedded atoms, we
  // will
  // choose that first before systems with more than 2 atoms in common. Cases
  // with two atoms
  // in common are in general flat systems to start with and can be embedded
  // cleanly.
  // when there are more than 2 atoms in common, these are most likely bridged
  // systems, which are
  // screwed up anyway, might as well screw them up later
  // if we do not have a system with two rings in common then we will return the
  // ring with max,
  // common atoms
  PRECONDITION(doneRings.size() > 0, "");
  PRECONDITION(fusedRings.size() > 1, "");

  RDKit::INT_VECT commonAtoms, res, doneAtoms, notDone;
  for (int i = 0; i < rdcast<int>(fusedRings.size()); i++) {
    if (std::find(doneRings.begin(), doneRings.end(), i) == doneRings.end()) {
      notDone.push_back(i);
    }
  }

  RDKit::Union(fusedRings, doneAtoms, &notDone);

  int maxCommonAtoms = 0;

  int currRingId = 0;
  for (const auto &fusedRing : fusedRings) {
    if (std::find(doneRings.begin(), doneRings.end(), currRingId) !=
        doneRings.end()) {
      currRingId++;
      continue;
    }
    commonAtoms.clear();
    int numCommonAtoms = 0;
    for (auto rii : fusedRing) {
      if (std::find(doneAtoms.begin(), doneAtoms.end(), (rii)) !=
          doneAtoms.end()) {
        commonAtoms.push_back(rii);
        numCommonAtoms++;
      }
    }
    if (numCommonAtoms == 2) {
      // if we found a ring with two atoms in common get out
      nextId = currRingId;
      return commonAtoms;  // FIX: this causes the rendering to be non-canonical
    }
    if (numCommonAtoms > maxCommonAtoms) {
      maxCommonAtoms = numCommonAtoms;
      nextId = currRingId;
      res = commonAtoms;
    }
    ++currRingId;
  }
  // here is an additional constrain we will put on the common atoms it is quite
  // likely that the common atoms form a chain (it is possible we can construct
  // some weird cases where this does not hold true - but for now we will assume
  // this is true. However the IDs in the res may not be in the order of going
  // from one end of the chain to the other -
  //  here is an example C1CCC(CC12)CCC2
  // - two rings here with three atoms in common
  // let ring1:(0,1,2,3,4,5) be a ring that is already embedded, then let
  // ring2:(4,3,6,7,8,5) be the ring that we found to be the next ring we should
  // embed.
  // The commonAtoms are (4,3,5) - note that they will be in this order since
  // the rings are always traversed in order. Now we would like these common
  // atoms to be returned in the order (5,4,3) - then we have a continuous
  // chain, we can do this by simply looking at the original ring order
  // (4,3,6,7,8,5) and observing that 5 need to come to the front

  // find out how many atoms from the end we need to move to the front
  unsigned int cmnLst = 0;
  unsigned int nCmn = res.size();
  for (unsigned int i = 0; i < nCmn; i++) {
    if (res[i] == fusedRings[nextId][i]) {
      cmnLst++;
    } else {
      break;
    }
  }
  // now do the moving if we have to
  if ((cmnLst > 0) && (cmnLst < res.size())) {
    RDKit::INT_VECT tempV = res;

    for (unsigned int i = cmnLst; i < nCmn; i++) {
      res[i - cmnLst] = tempV[i];
    }
    unsigned int nMov = nCmn - cmnLst;
    for (unsigned int i = 0; i < cmnLst; i++) {
      res[nMov + i] = tempV[i];
    }
  }

  POSTCONDITION(res.size() > 0, "");
  return res;
}

RDKit::INT_VECT getAllRotatableBonds(const RDKit::ROMol &mol) {
  RDKit::INT_VECT res;
  for (const auto bond : mol.bonds()) {
    int bid = bond->getIdx();
    if ((bond->getStereo() <= RDKit::Bond::STEREOANY) &&
        (!(mol.getRingInfo()->numBondRings(bid)))) {
      res.push_back(bid);
    }
  }
  return res;
}

RDKit::INT_VECT getRotatableBonds(const RDKit::ROMol &mol, unsigned int aid1,
                                  unsigned int aid2) {
  PRECONDITION(aid1 < mol.getNumAtoms(), "");
  PRECONDITION(aid2 < mol.getNumAtoms(), "");

  RDKit::INT_LIST path = RDKit::MolOps::getShortestPath(mol, aid1, aid2);
  RDKit::INT_VECT res;
  if (path.size() >= 4) {
    // remove the first atom (aid1) and last atom (aid2)
    CHECK_INVARIANT(static_cast<unsigned int>(path.front()) == aid1,
                    "bad first element");
    path.pop_front();
    CHECK_INVARIANT(static_cast<unsigned int>(path.back()) == aid2,
                    "bad last element");
    path.pop_back();

    auto pid = path.front();
    for (auto aid : path) {
      if (aid == pid) {
        continue;
      }
      const RDKit::Bond *bond = mol.getBondBetweenAtoms(pid, aid);
      int bid = bond->getIdx();
      if ((bond->getStereo() <= RDKit::Bond::STEREOANY) &&
          (!(mol.getRingInfo()->numBondRings(bid)))) {
        res.push_back(bid);
      }
      pid = aid;
    }
  }
  return res;
}

void getNbrAtomAndBondIds(unsigned int aid, const RDKit::ROMol *mol,
                          RDKit::INT_VECT &aids, RDKit::INT_VECT &bids) {
  CHECK_INVARIANT(mol, "");
  unsigned int na = mol->getNumAtoms();
  URANGE_CHECK(aid, na);

  for (auto nbr : mol->atomNeighbors(mol->getAtomWithIdx(aid))) {
    auto bi = mol->getBondBetweenAtoms(aid, nbr->getIdx())->getIdx();
    aids.push_back(nbr->getIdx());
    bids.push_back(bi);
  }
}

// find pairs of bonds that can be permuted at a non-ring degree 4
// node. This function will return only those pairs that cannot be
// permuted by flipping a rotatable bond
//
//       D
//       |
//       b3
//       |
//  A-b1-B-b2-C
//       |
//       b4
//       |
//       E
// For example in the above situation on the pairs (b1, b3) and (b1, b4) will be
// returned
// All other permutations can be achieved via a rotatable bond flip.
INT_PAIR_VECT findBondsPairsToPermuteDeg4(const RDGeom::Point2D &center,
                                          const RDKit::INT_VECT &nbrBids,
                                          const VECT_C_POINT &nbrLocs) {
  INT_PAIR_VECT res;

  // make sure there are four of them
  CHECK_INVARIANT(nbrBids.size() == 4, "");
  CHECK_INVARIANT(nbrLocs.size() == 4, "");

  std::vector<RDGeom::Point2D> nbrPts;
  nbrPts.reserve(nbrLocs.size());
  for (const auto &nloc : nbrLocs) {
    RDGeom::Point2D v = (*nloc) - center;
    nbrPts.push_back(v);
  }

  // now find the lay out of the bonds and return the bonds that are 90deg to
  // the
  // the bond to the first neighbor; i.e. we want to find b3 and b4 in the above
  // picture
  double dp1 = nbrPts[0].dotProduct(nbrPts[1]);
  if (fabs(dp1) < 1.e-3) {
    // the first two vectors are perpendicular to each other. We now have b1 and
    // b3 we need to
    // find b4
    INT_PAIR p1(nbrBids[0], nbrBids[1]);
    res.push_back(p1);

    double dp2 = nbrPts[0].dotProduct(nbrPts[2]);
    if (fabs(dp2) < 1.e-3) {
      // now we found b4 as well return the results
      INT_PAIR p2(nbrBids[0], nbrBids[2]);
      res.push_back(p2);
    } else {
      // bids[0] and bids[2] are opposite to each other and we know bids[1] is
      // perpendicular to bids[0]. So bids[3] is also perpendicular to bids[0]
      INT_PAIR p2(nbrBids[0], nbrBids[3]);
      res.push_back(p2);
    }
    return res;
  } else {
    // bids[0] and bids[1] are opposite to each other, so bids[2] and bids[3]
    // must
    // be perpendicular to bids[0]
    INT_PAIR p1(nbrBids[0], nbrBids[2]);
    res.push_back(p1);
    INT_PAIR p2(nbrBids[0], nbrBids[3]);
    res.push_back(p2);
    return res;
  }
}

template <class T>
T rankAtomsByRank(const RDKit::ROMol &mol, const T &commAtms, bool ascending) {
  size_t natms = commAtms.size();
  INT_PAIR_VECT rankAid;
  rankAid.reserve(natms);
  typename T::const_iterator ci;
  for (ci = commAtms.begin(); ci != commAtms.end(); ci++) {
    unsigned int rank;
    const RDKit::Atom *at = mol.getAtomWithIdx(*ci);
    if (at->hasProp(RDKit::common_properties::_CIPRank)) {
      at->getProp(RDKit::common_properties::_CIPRank, rank);
    } else {
      rank = mol.getNumAtoms() * getAtomDepictRank(at) + (*ci);
    }
    rankAid.push_back(std::make_pair(rank, (*ci)));
  }
  if (ascending) {
    std::stable_sort(rankAid.begin(), rankAid.end(),
                     [](const auto &e1, const auto &e2) { return e1 < e2; });
  } else {
    std::stable_sort(rankAid.begin(), rankAid.end(),
                     [](const auto &e1, const auto &e2) { return e1 > e2; });
  }
  T res;
  std::for_each(rankAid.begin(), rankAid.end(),
                [&res](const auto &elem) { res.push_back(elem.second); });
  return res;
}

template RDKit::INT_VECT rankAtomsByRank(const RDKit::ROMol &mol,
                                         const RDKit::INT_VECT &commAtms,
                                         bool ascending);
template RDKit::INT_DEQUE rankAtomsByRank(const RDKit::ROMol &mol,
                                          const RDKit::INT_DEQUE &commAtms,
                                          bool ascending);
template RDKit::INT_LIST rankAtomsByRank(const RDKit::ROMol &mol,
                                         const RDKit::INT_LIST &commAtms,
                                         bool ascending);

bool hasTerminalRGroupOrQueryHydrogen(const RDKit::ROMol &mol) {
  // we do not need the allowRGroups logic if there are no
  // terminal dummy atoms
  auto atoms = mol.atoms();
  return std::any_of(atoms.begin(), atoms.end(),
                     RDKit::isAtomTerminalRGroupOrQueryHydrogen);
}

std::unique_ptr<RDKit::RWMol> prepareTemplateForRGroups(
    RDKit::RWMol &templateMol) {
  auto queryParams = RDKit::MolOps::AdjustQueryParameters::noAdjustments();
  queryParams.adjustSingleBondsToDegreeOneNeighbors = true;
  queryParams.adjustSingleBondsBetweenAromaticAtoms = true;
  RDKit::MolOps::adjustQueryProperties(templateMol, &queryParams);
  std::map<unsigned int, unsigned int> removedIdxToNbrIdx;
  std::unique_ptr<RDKit::RWMol> reducedTemplateMol;
  for (const auto &bond : templateMol.bonds()) {
    int atomIdxToRemove = -1;
    int nbrIdx = -1;
    auto beginAtom = bond->getBeginAtom();
    auto endAtom = bond->getEndAtom();
    if (RDKit::isAtomTerminalRGroupOrQueryHydrogen(beginAtom) &&
        endAtom->hasQuery()) {
      atomIdxToRemove = beginAtom->getIdx();
      nbrIdx = endAtom->getIdx();
    } else if (RDKit::isAtomTerminalRGroupOrQueryHydrogen(endAtom) &&
               beginAtom->hasQuery()) {
      atomIdxToRemove = endAtom->getIdx();
      nbrIdx = beginAtom->getIdx();
    }
    if (atomIdxToRemove != -1) {
      removedIdxToNbrIdx[atomIdxToRemove] = nbrIdx;
    }
  }
  if (!removedIdxToNbrIdx.empty()) {
    reducedTemplateMol.reset(new RDKit::RWMol(templateMol));
    for (auto reducedTemplateAtom : reducedTemplateMol->atoms()) {
      auto formerIdx = reducedTemplateAtom->getIdx();
      reducedTemplateAtom->setProp(FORMER_IDX, formerIdx);
      auto it = removedIdxToNbrIdx.find(formerIdx);
      if (it != removedIdxToNbrIdx.end()) {
        auto otherAtom = reducedTemplateMol->getAtomWithIdx(it->second);
        std::vector<unsigned int> formerNbrIndices;
        otherAtom->getPropIfPresent(FORMER_NBR_INDICES, formerNbrIndices);
        formerNbrIndices.push_back(formerIdx);
        otherAtom->setProp(FORMER_NBR_INDICES, formerNbrIndices);
      }
    }
    reducedTemplateMol->beginBatchEdit();
    for (const auto &pair : removedIdxToNbrIdx) {
      reducedTemplateMol->removeAtom(
          reducedTemplateMol->getAtomWithIdx(pair.first));
    }
    reducedTemplateMol->commitBatchEdit();
  }
  return reducedTemplateMol;
}

void reducedToFullMatches(const RDKit::RWMol &reducedQuery,
                          const RDKit::RWMol &molHs,
                          std::vector<RDKit::MatchVectType> &matches) {
  boost::dynamic_bitset<> molHsMatches(molHs.getNumAtoms());
  for (auto &match : matches) {
    molHsMatches.reset();
    for (const auto &pair : match) {
      molHsMatches.set(pair.second);
    }
    RDKit::MatchVectType newMatch;
    for (auto pairIt = match.begin(); pairIt != match.end(); ++pairIt) {
      const auto reducedQueryAtom = reducedQuery.getAtomWithIdx(pairIt->first);
      const auto molAtom = molHs.getAtomWithIdx(pairIt->second);
      unsigned int formerIdx;
      reducedQueryAtom->getProp(FORMER_IDX, formerIdx);
      pairIt->first = formerIdx;
      std::vector<unsigned int> formerNbrIndices;
      reducedQueryAtom->getPropIfPresent(FORMER_NBR_INDICES, formerNbrIndices);
      for (const auto &molNbr : molHs.atomNeighbors(molAtom)) {
        if (formerNbrIndices.empty()) {
          break;
        }
        auto molNbrIdx = molNbr->getIdx();
        if (!molHsMatches.test(molNbrIdx)) {
          auto formerNbrIdx = formerNbrIndices.back();
          formerNbrIndices.pop_back();
          newMatch.emplace_back(formerNbrIdx, molNbrIdx);
        }
      }
    }
    auto matchSize = match.size();
    match.resize(matchSize + newMatch.size());
    std::move(newMatch.begin(), newMatch.end(), match.begin() + matchSize);
  }
}

bool invertWedgingIfMolHasFlipped(RDKit::ROMol &mol,
                                  const RDGeom::Transform3D &trans) {
  constexpr double FLIP_THRESHOLD = -0.99;
  auto zRot = trans.getVal(2, 2);
  bool shouldFlip = zRot < FLIP_THRESHOLD;
  if (shouldFlip) {
    RDKit::Chirality::invertMolBlockWedgingInfo(mol);
  }
  return shouldFlip;
}
}  // namespace RDDepict
