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

void WedgeBond(Bond *bond, unsigned int fromAtomIdx, const Conformer *conf) {
  PRECONDITION(bond, "no bond");
  PRECONDITION(conf, "no conformer");
  PRECONDITION(&conf->getOwningMol() == &bond->getOwningMol(),
               "bond and conformer do not belong to same molecule");
  if (bond->getBondType() != Bond::SINGLE) {
    return;
  }
  Bond::BondDir dir = DetermineBondWedgeState(bond, fromAtomIdx, conf);
  if (dir == Bond::BEGINWEDGE || dir == Bond::BEGINDASH) {
    bond->setBondDir(dir);
  }
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
  for (ROMol::ConstAtomIterator cai = mol.beginAtoms(); cai != mol.endAtoms();
       ++cai) {
    const Atom *at = *cai;
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
    for(auto *nat : at->nbrs()) {
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
  std::vector<unsigned int> indices(mol.getNumAtoms());
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    indices[i] = i;
  }
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
    if (type != Atom::CHI_TETRAHEDRAL_CW && type != Atom::CHI_TETRAHEDRAL_CCW) {
      break;
    }
    std::vector<std::pair<int, int>> nbrScores;
    for(auto *bond : atom->bonds()) {
      // can only wedge single bonds:
      if (bond->getBondType() != Bond::SINGLE) {
        continue;
      }

      int bid = bond->getIdx();
      if (res.find(bid) == res.end()) {
        // very strong preference for Hs:
        if (bond->getOtherAtom(atom)->getAtomicNum() == 1) {
          nbrScores.emplace_back(
              -1000000, bid);  // lower than anything else can be
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
        nbrScores.emplace_back(nbrScore, bid);
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

  for(auto *nbrBond : atom->bonds()) {
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
Bond::BondDir DetermineBondWedgeState(const Bond *bond,
                                      const INT_MAP_INT &wedgeBonds,
                                      const Conformer *conf) {
  PRECONDITION(bond, "no bond");
  int bid = bond->getIdx();
  auto wbi = wedgeBonds.find(bid);
  if (wbi == wedgeBonds.end()) {
    return bond->getBondDir();
  }

  unsigned int waid = wbi->second;
  return DetermineBondWedgeState(bond, waid, conf);
}

// handles stereochem markers set by the Mol file parser and
// converts them to the RD standard:
void DetectAtomStereoChemistry(RWMol &mol, const Conformer *conf) {
  PRECONDITION(conf, "no conformer");
  PRECONDITION(&(conf->getOwningMol()) == &mol,
               "conformer does not belong to molecule");
  MolOps::assignChiralTypesFromBondDirs(mol, conf->getId(), true);
}

void ClearSingleBondDirFlags(ROMol &mol) {
  for (RWMol::BondIterator bondIt = mol.beginBonds(); bondIt != mol.endBonds();
       ++bondIt) {
    if ((*bondIt)->getBondType() == Bond::SINGLE) {
      if ((*bondIt)->getBondDir() == Bond::UNKNOWN) {
        (*bondIt)->setProp(common_properties::_UnknownStereo, 1);
      }
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
}  // namespace RDKit
