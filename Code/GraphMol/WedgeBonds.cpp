//
//  Copyright (C) 2023 Greg Landrum and other RDKit contributors
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
#include <GraphMol/QueryOps.h>
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

#include <cstdlib>

namespace RDKit {

namespace Chirality {

namespace {
std::tuple<unsigned int, unsigned int, unsigned int> getDoubleBondPresence(
    const ROMol &mol, const Atom &atom) {
  unsigned int hasDouble = 0;
  unsigned int hasKnownDouble = 0;
  unsigned int hasAnyDouble = 0;
  for (const auto &nbri : boost::make_iterator_range(mol.getAtomBonds(&atom))) {
    const auto bond = mol[nbri];
    if (bond->getBondType() == Bond::BondType::DOUBLE) {
      ++hasDouble;
      if (bond->getStereo() == Bond::BondStereo::STEREOANY) {
        ++hasAnyDouble;
      } else if (bond->getStereo() > Bond::BondStereo::STEREOANY) {
        ++hasKnownDouble;
      }
    }
  }
  return std::make_tuple(hasDouble, hasKnownDouble, hasAnyDouble);
}
}  // namespace

namespace detail {

std::pair<bool, INT_VECT> countChiralNbrs(const ROMol &mol, int noNbrs) {
  // we need ring information; make sure findSSSR has been called before
  // if not call now
  if (!mol.getRingInfo()->isInitialized()) {
    MolOps::findSSSR(mol);
  }

  INT_VECT nChiralNbrs(mol.getNumAtoms(), noNbrs);

  // start by looking for bonds that are already wedged
  for (const auto bond : mol.bonds()) {
    if (bond->getBondDir() == Bond::BEGINWEDGE ||
        bond->getBondDir() == Bond::BEGINDASH ||
        bond->getBondDir() == Bond::UNKNOWN) {
      if (bond->getBeginAtom()->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
          bond->getBeginAtom()->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
        nChiralNbrs[bond->getBeginAtomIdx()] = noNbrs + 1;
      } else if (bond->getEndAtom()->getChiralTag() ==
                     Atom::CHI_TETRAHEDRAL_CW ||
                 bond->getEndAtom()->getChiralTag() ==
                     Atom::CHI_TETRAHEDRAL_CCW) {
        nChiralNbrs[bond->getEndAtomIdx()] = noNbrs + 1;
      }
    }
  }

  // now rank atoms by the number of chiral neighbors or Hs they have:
  bool chiNbrs = false;
  for (const auto at : mol.atoms()) {
    if (nChiralNbrs[at->getIdx()] > noNbrs) {
      // std::cerr << " SKIPPING1: " << at->getIdx() << std::endl;
      continue;
    }
    Atom::ChiralType type = at->getChiralTag();
    if (type != Atom::CHI_TETRAHEDRAL_CW && type != Atom::CHI_TETRAHEDRAL_CCW) {
      continue;
    }
    nChiralNbrs[at->getIdx()] = 0;
    chiNbrs = true;
    for (const auto nat : mol.atomNeighbors(at)) {
      if (nat->getAtomicNum() == 1) {
        // special case: it's an H... we weight these especially high:
        nChiralNbrs[at->getIdx()] -= 10;
        continue;
      }
      type = nat->getChiralTag();
      if (type != Atom::CHI_TETRAHEDRAL_CW &&
          type != Atom::CHI_TETRAHEDRAL_CCW) {
        continue;
      }
      nChiralNbrs[at->getIdx()] -= 1;
    }
  }
  return std::pair<bool, INT_VECT>(chiNbrs, nChiralNbrs);
}

//
// Determine bond wedge state
///
Bond::BondDir determineBondWedgeState(const Bond *bond,
                                      unsigned int fromAtomIdx,
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

  Atom *atom, *bondAtom;  // = bond->getBeginAtom();
  if (bond->getBeginAtom()->getIdx() == fromAtomIdx) {
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
  std::list<int> neighborBondIndices;
  std::list<double> neighborBondAngles;
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
    const Bond *nbrBond = (*mol)[*beg];
    Atom *otherAtom = nbrBond->getOtherAtom(atom);
    if (nbrBond != bond) {
      tmpPt = conf->getAtomPos(otherAtom->getIdx());
      tmpPt.z = 0.0;
      RDGeom::Point3D tmpVect = centerLoc.directionVector(tmpPt);
      double angle = refVect.signedAngleTo(tmpVect);
      if (angle < 0.0) {
        angle += 2. * M_PI;
      }
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
  // here the hydrogen will be between 1 and 2 and we need to add an
  // additional swap
  if (neighborBondAngles.size() == 3) {
    // three coordinated
    auto angleIt = neighborBondAngles.begin();
    ++angleIt;  // the first is the 0 (or reference bond - we will ignoire
                // that
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
Bond::BondDir determineBondWedgeState(const Bond *bond,
                                      const INT_MAP_INT &wedgeBonds,
                                      const Conformer *conf) {
  PRECONDITION(bond, "no bond");
  int bid = bond->getIdx();
  auto wbi = wedgeBonds.find(bid);
  if (wbi == wedgeBonds.end()) {
    return bond->getBondDir();
  }

  unsigned int waid = wbi->second;
  return determineBondWedgeState(bond, waid, conf);
}

// picks a bond for atom that we will wedge when we write the mol file
// returns idx of that bond.
int pickBondToWedge(const Atom *atom, const ROMol &mol,
                    const INT_VECT &nChiralNbrs, const INT_MAP_INT &resSoFar,
                    int noNbrs) {
  // here is what we are going to do
  // - at each chiral center look for a bond that is begins at the atom and
  //   is not yet picked to be wedged for a different chiral center, preferring
  //   bonds to Hs
  // - if we do not find a bond that begins at the chiral center - we will take
  //   the first bond that is not yet picked by any other chiral centers
  // we use the orders calculated above to determine which order to do the
  // wedging

  // If we call wedgeMolBonds() on a fragment, it can happen that we end up with
  // atoms that don't have enough neighbors. Those are going to cause problems,
  // so just bail here.
  // if (atom->getDegree() < 3) {
  //   return -1;
  // }
  std::vector<std::pair<int, int>> nbrScores;
  for (const auto bond : mol.atomBonds(atom)) {
    // can only wedge single bonds:
    if (bond->getBondType() != Bond::SINGLE) {
      continue;
    }

    int bid = bond->getIdx();
    if (resSoFar.find(bid) == resSoFar.end()) {
      // very strong preference for Hs:
      if (bond->getOtherAtom(atom)->getAtomicNum() == 1) {
        nbrScores.emplace_back(-1000000,
                               bid);  // lower than anything else can be
        continue;
      }
      // prefer lower atomic numbers with lower degrees and no specified
      // chirality:
      const Atom *oatom = bond->getOtherAtom(atom);
      int nbrScore = oatom->getAtomicNum() + 100 * oatom->getDegree() +
                     1000 * ((oatom->getChiralTag() != Atom::CHI_UNSPECIFIED));
      // prefer neighbors that are nonchiral or have as few chiral neighbors
      // as possible:
      int oIdx = oatom->getIdx();
      if (nChiralNbrs[oIdx] < noNbrs) {
        // the counts are negative, so we have to subtract them off
        nbrScore -= 100000 * nChiralNbrs[oIdx];
      }
      // prefer bonds to non-ring atoms:
      nbrScore += 10000 * mol.getRingInfo()->numAtomRings(oIdx);
      // prefer non-ring bonds;
      nbrScore += 20000 * mol.getRingInfo()->numBondRings(bid);
      // prefer bonds to atoms which don't have a double bond from them
      unsigned int hasDoubleBond;       // is a double bond there?
      unsigned int hasKnownDoubleBond;  // is specified stereo there?
      unsigned int hasAnyDoubleBond;    // is STEREOANY there?
      std::tie(hasDoubleBond, hasKnownDoubleBond, hasAnyDoubleBond) =
          getDoubleBondPresence(mol, *oatom);
      nbrScore += 11000 * hasDoubleBond;
      nbrScore += 12000 * hasKnownDoubleBond;
      nbrScore += 23000 * hasAnyDoubleBond;

      // std::cerr << "    nrbScore: " << idx << " - " << oIdx << " : "
      //           << nbrScore << " nChiralNbrs: " << nChiralNbrs[oIdx]
      //           << std::endl;
      nbrScores.emplace_back(nbrScore, bid);
    }
  }
  // There's still one situation where this whole thing can fail: an unlucky
  // situation where all neighbors of all neighbors of an atom are chiral and
  // that atom ends up being the last one picked for stereochem assignment. This
  // also happens in cases where the chiral atom doesn't have all of its
  // neighbors (like when working with partially sanitized fragments)
  //
  // We'll bail here by returning -1
  if (nbrScores.empty()) {
    return -1;
  }
  std::sort(nbrScores.begin(), nbrScores.end(), Rankers::pairLess);
  return nbrScores[0].second;
}

}  // namespace detail

INT_MAP_INT pickBondsToWedge(const ROMol &mol, BondWedgingParameters *params) {
  // returns map of bondIdx -> bond begin atom for those bonds that
  // need wedging.
  std::vector<unsigned int> indices(mol.getNumAtoms());
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    indices[i] = i;
  }
  static int noNbrs = 100;
  std::pair<bool, INT_VECT> retVal = detail::countChiralNbrs(mol, noNbrs);
  bool chiNbrs = retVal.first;
  INT_VECT nChiralNbrs = retVal.second;
  if (chiNbrs) {
    std::sort(indices.begin(), indices.end(), [&](auto i1, auto i2) {
      return nChiralNbrs[i1] < nChiralNbrs[i2];
    });
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
  INT_MAP_INT res;
  for (auto idx : indices) {
    if (nChiralNbrs[idx] > noNbrs) {
      // std::cerr << " SKIPPING2: " << idx << std::endl;
      continue;  // already have a wedged bond here
    }
    const Atom *atom = mol.getAtomWithIdx(idx);
    Atom::ChiralType type = atom->getChiralTag();
    // the indices are ordered such that all chiral atoms come first. If
    // this has no chiral flag, we can stop the whole loop:
    if (type != Atom::CHI_TETRAHEDRAL_CW && type != Atom::CHI_TETRAHEDRAL_CCW) {
      break;
    }
    int bnd = detail::pickBondToWedge(atom, mol, nChiralNbrs, res, noNbrs);
    if (bnd >= 0) {
      res[bnd] = idx;
    }
  }
  return res;
}

void wedgeMolBonds(ROMol &mol, const Conformer *conf,
                   const BondWedgingParameters *params) {
  PRECONDITION(conf || mol.getNumConformers(), "no conformer available");
  if (!conf) {
    conf = &mol.getConformer();
  }
  auto wedgeBonds = pickBondsToWedge(mol);

  static BondWedgingParameters defaultWedgingParams;
  if (!params) {
    params = &defaultWedgingParams;
  }
  for (auto bond : mol.bonds()) {
    if (bond->getBondType() == Bond::SINGLE) {
      Bond::BondDir dir =
          detail::determineBondWedgeState(bond, wedgeBonds, conf);
      if (dir == Bond::BEGINWEDGE || dir == Bond::BEGINDASH) {
        bond->setBondDir(dir);

        // it is possible that this
        // wedging was determined by a chiral atom at the end of the
        // bond (instead of at the beginning). In this case we need to
        // reverse the begin and end atoms for the bond
        auto wbi = wedgeBonds.find(bond->getIdx());
        if (wbi != wedgeBonds.end() &&
            static_cast<unsigned int>(wbi->second) != bond->getBeginAtomIdx()) {
          auto tmp = bond->getBeginAtomIdx();
          bond->setBeginAtomIdx(bond->getEndAtomIdx());
          bond->setEndAtomIdx(tmp);
        }
      }
    }
  }
}

void wedgeBond(Bond *bond, unsigned int fromAtomIdx, const Conformer *conf) {
  PRECONDITION(bond, "no bond");
  PRECONDITION(conf, "no conformer");
  PRECONDITION(&conf->getOwningMol() == &bond->getOwningMol(),
               "bond and conformer do not belong to same molecule");
  if (bond->getBondType() != Bond::SINGLE) {
    return;
  }
  Bond::BondDir dir = detail::determineBondWedgeState(bond, fromAtomIdx, conf);
  if (dir == Bond::BEGINWEDGE || dir == Bond::BEGINDASH) {
    bond->setBondDir(dir);
  }
}

}  // namespace Chirality
}  // namespace RDKit
