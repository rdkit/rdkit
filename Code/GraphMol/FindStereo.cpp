//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
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
#include <algorithm>
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>

#include <boost/dynamic_bitset.hpp>
#include <boost/format.hpp>
#include "Chirality.h"

namespace RDKit {
namespace Chirality {
#ifndef _MSC_VER
const unsigned StereoInfo::NOATOM = std::numeric_limits<unsigned>::max();
#endif
namespace detail {

StereoInfo getStereoInfo(const Bond *bond) {
  PRECONDITION(bond, "bond is null");
  StereoInfo sinfo;
  const auto beginAtom = bond->getBeginAtom();
  const auto endAtom = bond->getEndAtom();
  if (bond->getBondType() == Bond::BondType::DOUBLE) {
    if (beginAtom->getDegree() < 2 || endAtom->getDegree() < 2 ||
        beginAtom->getDegree() > 3 || endAtom->getDegree() > 3) {
      throw ValueErrorException("invalid atom degree in getStereoInfo(bond)");
    }

    sinfo.type = StereoType::Bond_Double;
    sinfo.centeredOn = bond->getIdx();
    sinfo.controllingAtoms.reserve(4);

    bool seenSquiggleBond = false;
    const auto &mol = bond->getOwningMol();
    for (const auto &nbri :
         boost::make_iterator_range(mol.getAtomBonds(beginAtom))) {
      const auto &nbr = mol[nbri];
      if (nbr->getIdx() != bond->getIdx()) {
        if (nbr->getBondDir() == Bond::BondDir::UNKNOWN) {
          seenSquiggleBond = true;
        }
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
        if (nbr->getBondDir() == Bond::BondDir::UNKNOWN) {
          seenSquiggleBond = true;
        }
        sinfo.controllingAtoms.push_back(
            nbr->getOtherAtomIdx(endAtom->getIdx()));
      }
    }
    if (endAtom->getDegree() == 2) {
      sinfo.controllingAtoms.push_back(StereoInfo::NOATOM);
    }
    Bond::BondStereo stereo = bond->getStereo();
    if (stereo == Bond::BondStereo::STEREOANY ||
        bond->getBondDir() == Bond::BondDir::EITHERDOUBLE || seenSquiggleBond) {
      sinfo.specified = Chirality::StereoSpecified::Unknown;
    } else if (stereo != Bond::BondStereo::STEREONONE) {
      if (stereo == Bond::BondStereo::STEREOE ||
          stereo == Bond::BondStereo::STEREOZ) {
        stereo = Chirality::translateEZLabelToCisTrans(stereo);
      }
      sinfo.specified = Chirality::StereoSpecified::Specified;
      const auto satoms = bond->getStereoAtoms();
      if (satoms.size() != 2) {
        throw ValueErrorException("only can support 2 stereo neighbors");
      }
      bool firstAtBegin;
      if (satoms[0] == sinfo.controllingAtoms[0]) {
        firstAtBegin = true;
      } else if (satoms[0] == sinfo.controllingAtoms[1]) {
        firstAtBegin = false;
      } else {
        throw ValueErrorException("controlling atom mismatch at begin");
      }
      bool firstAtEnd;
      if (satoms[1] == sinfo.controllingAtoms[2]) {
        firstAtEnd = true;
      } else if (satoms[1] == sinfo.controllingAtoms[3]) {
        firstAtEnd = false;
      } else {
        throw ValueErrorException("controlling atom mismatch at end");
      }
      auto mismatch = firstAtBegin ^ firstAtEnd;
      if (mismatch) {
        stereo = (stereo == Bond::BondStereo::STEREOCIS
                      ? Bond::BondStereo::STEREOTRANS
                      : Bond::BondStereo::STEREOCIS);
      }
      switch (stereo) {
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

  const auto &mol = atom->getOwningMol();
  int explicitUnknownStereo = 0;
  for (const auto &nbri : boost::make_iterator_range(mol.getAtomBonds(atom))) {
    const auto &bnd = mol[nbri];
    if (bnd->getBondDir() == Bond::UNKNOWN) {
      explicitUnknownStereo = 1;
    } else if (!explicitUnknownStereo) {
      bnd->getPropIfPresent<int>(common_properties::_UnknownStereo,
                                 explicitUnknownStereo);
    }
    sinfo.controllingAtoms.push_back(bnd->getOtherAtomIdx(atom->getIdx()));
  }
  std::vector<unsigned> origNbrOrder = sinfo.controllingAtoms;
  std::sort(sinfo.controllingAtoms.begin(), sinfo.controllingAtoms.end());

  if (explicitUnknownStereo) {
    sinfo.specified = StereoSpecified::Unknown;
  } else {
    Atom::ChiralType stereo = atom->getChiralTag();
    if (stereo == Atom::ChiralType::CHI_TETRAHEDRAL_CCW ||
        stereo == Atom::ChiralType::CHI_TETRAHEDRAL_CW) {
      sinfo.specified = StereoSpecified::Specified;
      unsigned nSwaps =
          countSwapsToInterconvert(origNbrOrder, sinfo.controllingAtoms);
      if (nSwaps % 2) {
        stereo = (stereo == Atom::ChiralType::CHI_TETRAHEDRAL_CCW
                      ? Atom::ChiralType::CHI_TETRAHEDRAL_CW
                      : Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
      }
      switch (stereo) {
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
  }

  return sinfo;
}

bool isBondPotentialStereoBond(const Bond *bond) {
  PRECONDITION(bond, "bond is null");
  if (bond->getBondType() != Bond::BondType::DOUBLE) {
    return false;
  }

  // at the moment the condition for being a potential stereo bond is that
  // each of the beginning and end neighbors must have at least 2 heavy atom
  // neighbors i.e. C/C=N/[H] is not a possible stereo bond but no more than 3
  // total neighbors.
  const auto beginAtom = bond->getBeginAtom();
  auto begHeavyDegree =
      beginAtom->getTotalDegree() - beginAtom->getTotalNumHs(true);
  const auto endAtom = bond->getEndAtom();
  auto endHeavyDegree =
      endAtom->getTotalDegree() - endAtom->getTotalNumHs(true);
  if (begHeavyDegree > 1 && beginAtom->getDegree() < 4 && endHeavyDegree > 1 &&
      endAtom->getDegree() < 4) {
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
    const auto &mol = atom->getOwningMol();
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

std::string getBondSymbol(const Bond *bond) {
  // FIX: this is not complete
  PRECONDITION(bond, "bad bond");
  std::string res;
  if (bond->getIsAromatic()) {
    res = ":";
  } else {
    switch (bond->getBondType()) {
      case Bond::BondType::SINGLE:
        res = "-";
        break;
      case Bond::BondType::DOUBLE:
        res = "=";
        break;
      case Bond::BondType::TRIPLE:
        res = "#";
        break;
      case Bond::BondType::AROMATIC:
        res = ":";
        break;
      default:
        res = "?";
        break;
    }
  }
  return res;
}

std::vector<StereoInfo> findPotentialStereo(ROMol &mol, bool cleanIt,
                                            bool flagPossible) {
  std::map<int, Atom::ChiralType> ochiralTypes;

  // FIX: this never removes stereo

  boost::dynamic_bitset<> knownAtoms(mol.getNumAtoms());
  boost::dynamic_bitset<> possibleAtoms(mol.getNumAtoms());
  std::vector<std::string> atomSymbols(mol.getNumAtoms());
  for (const auto atom : mol.atoms()) {
    auto aidx = atom->getIdx();
    if (detail::isAtomPotentialStereoAtom(atom)) {
      auto sinfo = detail::getStereoInfo(atom);
      switch (sinfo.specified) {
        case Chirality::StereoSpecified::Unknown:
        case Chirality::StereoSpecified::Specified:
          knownAtoms.set(aidx);
          break;
        case Chirality::StereoSpecified::Unspecified:
          break;
        default:
          throw ValueErrorException("bad StereoInfo.specified type");
      }
      if (flagPossible ||
          sinfo.specified != Chirality::StereoSpecified::Unspecified) {
        possibleAtoms.set(aidx);
        // set "fake stereo"
        ochiralTypes[atom->getIdx()] = atom->getChiralTag();
        atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
        atomSymbols[aidx] =
            (boost::format("%s-%d") % atom->getSymbol() % atom->getIdx()).str();
      } else {
        atomSymbols[aidx] = atom->getSymbol();
      }
    } else {
      atomSymbols[aidx] = atom->getSymbol();
    }
  }

  std::vector<std::string> bondSymbols(mol.getNumBonds());
  boost::dynamic_bitset<> knownBonds(mol.getNumBonds());
  boost::dynamic_bitset<> possibleBonds(mol.getNumBonds());
  for (const auto bond : mol.bonds()) {
    auto bidx = bond->getIdx();
    if (detail::isBondPotentialStereoBond(bond)) {
      auto sinfo = detail::getStereoInfo(bond);
      switch (sinfo.specified) {
        case Chirality::StereoSpecified::Unknown:
        case Chirality::StereoSpecified::Specified:
          knownBonds.set(bidx);
          break;
        case Chirality::StereoSpecified::Unspecified:
          break;
        default:
          throw ValueErrorException("bad StereoInfo.specified type");
      }
      if (flagPossible ||
          sinfo.specified != Chirality::StereoSpecified::Unspecified) {
        possibleBonds.set(bidx);
        bondSymbols[bidx] =
            (boost::format("%s-%d") % getBondSymbol(bond) % bond->getIdx())
                .str();
      } else {
        bondSymbols[bidx] = getBondSymbol(bond);
      }
    } else {
      bondSymbols[bidx] = getBondSymbol(bond);
    }
  }

  std::vector<StereoInfo> res;
  while (possibleAtoms.count() || possibleBonds.count()) {
    res.clear();
    bool removedStereo = false;

    // we will use the canonicalization code, pretending that each potential
    // stereo atom and bond is specified and different from all others. After
    // we've done that we can re-examine the potential stereo atoms and bonds
    // and remove any where two controlling atoms have the same rank
    boost::dynamic_bitset<> atomsInPlay(mol.getNumAtoms());
    atomsInPlay.set();
    boost::dynamic_bitset<> bondsInPlay(mol.getNumBonds());
    bondsInPlay.set();
    std::vector<unsigned int> aranks;
    const bool breakTies = false;
    const bool includeChirality = true;
    const bool includeIsotopes = false;
    Canon::rankFragmentAtoms(mol, aranks, atomsInPlay, bondsInPlay,
                             &atomSymbols, &bondSymbols, breakTies,
                             includeChirality, includeIsotopes);
    for (const auto atom : mol.atoms()) {
      auto aidx = atom->getIdx();
      if (ochiralTypes.find(aidx) != ochiralTypes.end()) {
        atom->setChiralTag(ochiralTypes[aidx]);
      }
      if (possibleAtoms[aidx]) {
        auto sinfo = detail::getStereoInfo(atom);
        std::vector<unsigned int> nbrs;
        nbrs.reserve(sinfo.controllingAtoms.size());
        bool haveADupe = false;
        for (auto nbrIdx : sinfo.controllingAtoms) {
          auto rnk = aranks[nbrIdx];
          if (std::find(nbrs.begin(), nbrs.end(), rnk) != nbrs.end()) {
            haveADupe = true;
            break;
          } else {
            nbrs.push_back(rnk);
          }
        }
        if (!haveADupe) {
          res.push_back(std::move(sinfo));
        } else {
          removedStereo = true;
          atomSymbols[aidx] = atom->getSymbol();
          possibleAtoms[aidx] = 0;
          if (cleanIt &&
              sinfo.specified != Chirality::StereoSpecified::Unspecified) {
            atom->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
          }
        }
      }
    }

    for (const auto bond : mol.bonds()) {
      auto bidx = bond->getIdx();
      if (possibleBonds[bidx]) {
        auto sinfo = detail::getStereoInfo(bond);
        ASSERT_INVARIANT(sinfo.controllingAtoms.size() == 4,
                         "bad controlling atoms size");
        bool haveADupe = false;
        if (sinfo.controllingAtoms[0] != Chirality::StereoInfo::NOATOM &&
            sinfo.controllingAtoms[1] != Chirality::StereoInfo::NOATOM &&
            aranks[sinfo.controllingAtoms[0]] ==
                aranks[sinfo.controllingAtoms[1]]) {
          haveADupe = true;
        }
        if (sinfo.controllingAtoms[2] != Chirality::StereoInfo::NOATOM &&
            sinfo.controllingAtoms[3] != Chirality::StereoInfo::NOATOM &&
            aranks[sinfo.controllingAtoms[2]] ==
                aranks[sinfo.controllingAtoms[3]]) {
          haveADupe = true;
        }
        if (!haveADupe) {
          res.push_back(std::move(sinfo));
        } else {
          removedStereo = true;
          bondSymbols[bidx] = getBondSymbol(bond);
          possibleBonds[bidx] = 0;
          if (cleanIt &&
              sinfo.specified != Chirality::StereoSpecified::Unspecified) {
            bond->setStereo(Bond::BondStereo::STEREONONE);
          }
        }
      }
    }
    if (!removedStereo) {
      break;
    }
  }
  return res;
}

// const_casts are always ugly, but we know that findPotentialStereo() doesn't
// modify the molecule if cleanIt is false:
std::vector<StereoInfo> findPotentialStereo(const ROMol &mol) {
  bool cleanIt = false;
  return findPotentialStereo(const_cast<ROMol &>(mol), cleanIt);
}

}  // namespace Chirality
}  // namespace RDKit
