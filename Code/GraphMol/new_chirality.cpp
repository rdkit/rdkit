//
//  Copyright (C) 2020 Greg Landrum and Rational Discovery LLC
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
#include "Chirality.h"

namespace RDKit {
namespace Chirality {

namespace detail {

StereoInfo getStereoInfo(const Bond *bond) {
  PRECONDITION(bond, "bond is null");
  StereoInfo sinfo;
  const auto beginAtom = bond->getBeginAtom();
  const auto endAtom = bond->getEndAtom();
  if (bond->getBondType() == Bond::BondType::DOUBLE) {
    if (beginAtom->getDegree() < 2 || endAtom->getDegree() < 2) {
      throw ValueErrorException("atom degree too low in getStereoInfo(bond)");
    }

    sinfo.type = StereoType::Bond;
    sinfo.centeredOn = bond->getIdx();
    sinfo.controllingAtoms.resize(2);

    const auto mol = bond->getOwningMol();
    for (const auto &nbri :
         boost::make_iterator_range(mol.getAtomBonds(beginAtom))) {
      const auto &nbr = mol[nbri];
      if (nbr->getIdx() != bond->getIdx()) {
        sinfo.controllingAtoms[0] = nbr->getOtherAtomIdx(beginAtom->getIdx());
        break;
      }
    }
    for (const auto &nbri :
         boost::make_iterator_range(mol.getAtomBonds(endAtom))) {
      const auto &nbr = mol[nbri];
      if (nbr->getIdx() != bond->getIdx()) {
        sinfo.controllingAtoms[1] = nbr->getOtherAtomIdx(endAtom->getIdx());
        break;
      }
    }
  } else {
    UNDER_CONSTRUCTION("unsupported bond type in getStereoInfo()");
  }

  return sinfo;
}

StereoInfo getStereoInfo(const Atom *atom) {
  PRECONDITION(atom, "atom is null");
  StereoInfo sinfo;

  sinfo.type = StereoType::Atom;
  sinfo.centeredOn = atom->getIdx();
  sinfo.controllingAtoms.reserve(atom->getDegree());

  const auto mol = atom->getOwningMol();
  for (const auto &nbri :
       boost::make_iterator_range(mol.getAtomNeighbors(atom))) {
    const auto &nbr = mol[nbri];
    sinfo.controllingAtoms.push_back(nbr->getIdx());
  }
  std::sort(sinfo.controllingAtoms.begin(), sinfo.controllingAtoms.end());
  return sinfo;
}

bool isBondPotentialStereoBond(const Bond *bond) {
  PRECONDITION(bond, "bond is null");
  if (bond->getBondType() != Bond::BondType::DOUBLE) {
    return false;
  }

  const auto beginAtom = bond->getBeginAtom();
  const auto endAtom = bond->getEndAtom();
  if ((beginAtom->getTotalDegree() - beginAtom->getTotalNumHs(true)) > 1 &&
      (endAtom->getTotalDegree() - endAtom->getTotalNumHs(true)) > 1) {
    return true;
  } else {
    return false;
  }
}

bool isAtomPotentialTetrahedralCenter(const Atom *atom) {
  PRECONDITION(atom, "atom is null");

  if (atom->getTotalDegree() > 4) {
    return false;
  } else {
    auto mol = atom->getOwningMol();
    auto degree = mol.getAtomDegree(atom);
    if (degree == 4) {
      // chirality is always possible with 4 nbrs
      return true;
    } else if (degree == 1) {
      // chirality is never possible with 1 nbr
      return false;
    } else if (degree < 3 &&
               (atom->getAtomicNum() != 15 && atom->getAtomicNum() != 33)) {
      // less than three neighbors is never stereogenic
      // unless it is a phosphine/arsine with implicit H
      return false;
    } else if (atom->getAtomicNum() == 15 || atom->getAtomicNum() == 33) {
      // from logical flow: degree is 2 or 3 (implicit H)
      // Since InChI Software v. 1.02-standard (2009), phosphines and arsines
      // are always treated as stereogenic even with H atom neighbors.
      // Accept automatically.
      return true;
    } else if (degree == 3) {
      // three-coordinate with a single H we'll accept automatically:
      if (atom->getTotalNumHs() == 1) {
        return true;
      } else {
        // otherwise we default to not being a legal center
        bool legalCenter = false;
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
        return legalCenter;
      }
    } else {
      return false;
    }
  }
}
bool isAtomPotentialStereoAtom(const Atom *atom) {
  return isAtomPotentialTetrahedralCenter(atom);
}
}  // namespace detail

std::vector<StereoInfo> findPotentialStereo(const ROMol &mol) {
  std::vector<unsigned int> atomRanks;
  Canon::chiralRankMolAtoms(mol, atomRanks);

  std::list<StereoInfo> active;
  unsigned int nUndefined = 0;
  for (const auto atom : mol.atoms()) {
    if (detail::isAtomPotentialStereoAtom(atom)) {
      active.push_back(detail::getStereoInfo(atom));
      if (atom->getChiralTag() == Atom::ChiralType::CHI_UNSPECIFIED) {
        ++nUndefined;
      }
    }
  }

  std::vector<StereoInfo> res;
  return res;
}

}  // namespace Chirality
}  // namespace RDKit
