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
const unsigned StereoInfo::NOATOM = std::numeric_limits<unsigned>::max();

namespace detail {

StereoInfo getStereoInfo(const Bond *bond) {
  PRECONDITION(bond, "bond is null");
  StereoInfo sinfo;
  const auto beginAtom = bond->getBeginAtom();
  const auto endAtom = bond->getEndAtom();
  if (bond->getBondType() == Bond::BondType::DOUBLE) {
    if (beginAtom->getDegree() < 2 || endAtom->getDegree() < 2 || beginAtom->getDegree() > 3 || endAtom->getDegree() > 3) {
        throw ValueErrorException("invalid atom degree in getStereoInfo(bond)");
    }

    sinfo.type = StereoType::Bond_Double;
    sinfo.centeredOn = bond->getIdx();
    sinfo.controllingAtoms.reserve(4);

    const auto mol = bond->getOwningMol();
    for (const auto &nbri :
         boost::make_iterator_range(mol.getAtomBonds(beginAtom))) {
      const auto &nbr = mol[nbri];
      if (nbr->getIdx() != bond->getIdx()) {
        sinfo.controllingAtoms.push_back(
            nbr->getOtherAtomIdx(beginAtom->getIdx()));
      }
    }
    if (beginAtom->getDegree() == 2) {
      sinfo.controllingAtoms.push_back(StereoInfo::NOATOM);
    }
    for (const auto &nbri :
         boost::make_iterator_range(mol.getAtomBonds(endAtom))) {
      const auto &nbr = mol[nbri];
      if (nbr->getIdx() != bond->getIdx()) {
        sinfo.controllingAtoms.push_back(
            nbr->getOtherAtomIdx(endAtom->getIdx()));
      }
    }
    if (endAtom->getDegree() == 2) {
      sinfo.controllingAtoms.push_back(StereoInfo::NOATOM);
    }
    Bond::BondStereo stereo = bond->getStereo();
    if(stereo == Bond::BondStereo::STEREONONE){
      // don't need to do anything
    } else if (stereo == Bond::BondStereo::STEREOANY){
      sinfo.specified = Chirality::StereoSpecified::Unknown;
    } else {
      if(stereo == Bond::BondStereo::STEREOE || stereo == Bond::BondStereo::STEREOZ){
        stereo = Chirality::translateEZLabelToCisTrans(stereo);
      }
      sinfo.specified = Chirality::StereoSpecified::Specified;
      const auto satoms = bond->getStereoAtoms();
      if(satoms.size()!=2){
        throw ValueErrorException("only can support 2 stereo neighbors");
      }
      bool firstAtBegin;
      if(satoms[0] == sinfo.controllingAtoms[0]){
        firstAtBegin = true;
      } else if (satoms[0] == sinfo.controllingAtoms[1]) {
        firstAtBegin = false;
      } else {
        throw ValueErrorException("controlling atom mismatch at begin");
      }
      bool firstAtEnd;
      if(satoms[1] == sinfo.controllingAtoms[2]){
        firstAtEnd = true;
      } else if (satoms[1] == sinfo.controllingAtoms[3]) {
        firstAtEnd = false;
      } else {
        throw ValueErrorException("controlling atom mismatch at end");
      }
      auto mismatch = firstAtBegin ^ firstAtEnd;
      if(mismatch) {
        stereo = (stereo == Bond::BondStereo::STEREOCIS ? Bond::BondStereo::STEREOTRANS : Bond::BondStereo::STEREOCIS);
      }
      switch(stereo){
      case Bond::BondStereo::STEREOCIS:
        sinfo.descriptor = Chirality::StereoDescriptor::Bond_Cis;
      break;
        case Bond::BondStereo::STEREOTRANS:
      sinfo.descriptor = Chirality::StereoDescriptor::Bond_Trans;
        break;
      default:
        UNDER_CONSTRUCTION("unrecognized bond stereo type");      
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

  sinfo.type = StereoType::Atom_Tetrahedral;
  sinfo.centeredOn = atom->getIdx();
  sinfo.controllingAtoms.reserve(atom->getDegree());

  const auto mol = atom->getOwningMol();
  for (const auto &nbri :
        boost::make_iterator_range(mol.getAtomBonds(atom))) {
    const auto &bnd = mol[nbri];
    sinfo.controllingAtoms.push_back(bnd->getOtherAtomIdx(atom->getIdx()));
  }
  std::vector<unsigned> origNbrOrder = sinfo.controllingAtoms;
  std::sort(sinfo.controllingAtoms.begin(), sinfo.controllingAtoms.end());

  Atom::ChiralType stereo=atom->getChiralTag();
  if(stereo == Atom::ChiralType::CHI_TETRAHEDRAL_CCW || stereo == Atom::ChiralType::CHI_TETRAHEDRAL_CW){
    sinfo.specified = StereoSpecified::Specified;
    unsigned nSwaps = countSwapsToInterconvert(origNbrOrder,sinfo.controllingAtoms);
    if(nSwaps %2){
      stereo = (stereo==Atom::ChiralType::CHI_TETRAHEDRAL_CCW?Atom::ChiralType::CHI_TETRAHEDRAL_CW:Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    }
    switch(stereo){
    case Atom::ChiralType::CHI_TETRAHEDRAL_CCW:
      sinfo.descriptor = StereoDescriptor::Tet_CCW;
      break;
    case Atom::ChiralType::CHI_TETRAHEDRAL_CW:
      sinfo.descriptor = StereoDescriptor::Tet_CW;
      break;
    default:
      UNDER_CONSTRUCTION("unrecognized chiral flag");      
    }
  }

  return sinfo;
}

bool isBondPotentialStereoBond(const Bond *bond) {
  PRECONDITION(bond, "bond is null");
  if (bond->getBondType() != Bond::BondType::DOUBLE) {
    return false;
  }

  // at the moment the condition for being a potential stereo bond is that
  // each of the beginning and end neighbors must have at least 2 heavy atom neighbors
  // i.e. C/C=N/[H] is not a possible stereo bond
  // but no more than 3 total neighbors.
  const auto beginAtom = bond->getBeginAtom();
  auto begHeavyDegree = beginAtom->getTotalDegree() - beginAtom->getTotalNumHs(true);
  const auto endAtom = bond->getEndAtom();
  auto endHeavyDegree = endAtom->getTotalDegree() - endAtom->getTotalNumHs(true);
  if (begHeavyDegree > 1 && beginAtom->getDegree() < 4 &&
      endHeavyDegree > 1 && endAtom->getDegree() < 4) {
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
  // FIX: need to remove the R/S dependency from this
  Canon::chiralRankMolAtoms(mol, atomRanks);

  std::list<StereoInfo> active;
  unsigned int nUndefined = 0;
  for (const auto atom : mol.atoms()) {
    if (detail::isAtomPotentialStereoAtom(atom)) {
      active.push_back(detail::getStereoInfo(atom));
      if(active.back().specified == Chirality::StereoSpecified::Unknown)
        ++nUndefined;
      } else if(active.back().specified == Chirality::StereoSpecified::Specified){

      }
    }
  }
  for (const auto bond : mol.bonds()) {
    if (detail::isBondPotentialStereoBond(bond)) {
      active.push_back(detail::getStereoInfo(bond));
      if(active.back().specified == Chirality::StereoSpecified::Unknown)
        ++nUndefined;
      }
    }
  }

  std::vector<StereoInfo> res;
  return res;
}

}  // namespace Chirality
}  // namespace RDKit
