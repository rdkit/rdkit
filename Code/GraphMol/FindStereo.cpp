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
#include <GraphMol/new_canon.h>
#include <RDGeneral/types.h>
#include <algorithm>
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>

#include <boost/dynamic_bitset.hpp>
#include "Chirality.h"
#include <GraphMol/QueryOps.h>

namespace RDKit {
namespace Chirality {
#ifndef _MSC_VER
const unsigned StereoInfo::NOATOM = std::numeric_limits<unsigned>::max();
#endif
namespace detail {

bool isAtomPotentialNontetrahedralCenter(const Atom *atom) {
  PRECONDITION(atom, "atom is null");
  auto nzdegree = Chirality::detail::getAtomNonzeroDegree(atom);
  auto impHDegree = atom->getTotalNumHs();
  auto tnzdegree = nzdegree + impHDegree;
  auto anum = atom->getAtomicNum();
  if (tnzdegree > 6 || tnzdegree < 2 || (anum < 12 && anum != 4)) {
    return false;
  }
  auto chiralType = atom->getChiralTag();
  if (chiralType >= Atom::ChiralType::CHI_SQUAREPLANAR &&
      chiralType <= Atom::ChiralType::CHI_OCTAHEDRAL) {
    return true;
  }

  // with at least four neighbors but nothing specified we can start to imagine
  // that it might be enhanced stereo
  if (chiralType == Atom::ChiralType::CHI_UNSPECIFIED && tnzdegree >= 4) {
    return true;
  }

  return false;
}
bool isAtomPotentialTetrahedralCenter(const Atom *atom) {
  PRECONDITION(atom, "atom is null");
  auto nzDegree = getAtomNonzeroDegree(atom);
  auto tnzDegree = nzDegree + atom->getTotalNumHs();
  if (tnzDegree > 4) {
    return false;
  } else {
    const auto &mol = atom->getOwningMol();
    if (nzDegree == 4) {
      // chirality is always possible with 4 nbrs
      return true;
    } else if (nzDegree <= 1) {
      // chirality is never possible with 0 or 1 nbr
      return false;
    } else if (nzDegree < 3 &&
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
    } else if (nzDegree == 3) {
      // three-coordinate with a single H we'll accept automatically:
      if (atom->getTotalNumHs() == 1) {
        if (detail::has_protium_neighbor(mol, atom)) {
          // more than one H is never stereogenic
          return false;
        }
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
        } else if (atom->getAtomicNum() == 7) {
          // three-coordinate N additional requirements:
          //   in a ring of size 3  (from InChI)
          // OR
          /// is a bridgehead atom (RDKit extension)
          if (mol.getRingInfo()->isAtomInRingOfSize(atom->getIdx(), 3) ||
              queryIsAtomBridgehead(atom)) {
            legalCenter = true;
          }
        }
        return legalCenter;
      }
    } else {
      return false;
    }
  }
}

bool isAtomPotentialStereoAtom(const Atom *atom,
                               bool allowNontetrahehdralStereo) {
  return isAtomPotentialTetrahedralCenter(atom) ||
         (allowNontetrahehdralStereo &&
          isAtomPotentialNontetrahedralCenter(atom));
}

bool isAtomPotentialStereoAtom(const Atom *atom) {
  return isAtomPotentialStereoAtom(atom, getAllowNontetrahedralChirality());
}

StereoInfo getStereoInfo(const Bond *bond) {
  PRECONDITION(bond, "bond is null");
  StereoInfo sinfo;
  const auto beginAtom = bond->getBeginAtom();
  const auto endAtom = bond->getEndAtom();
  if (bond->getBondType() == Bond::BondType::DOUBLE) {
    if (beginAtom->getDegree() < 1 || endAtom->getDegree() < 1 ||
        beginAtom->getDegree() > 3 || endAtom->getDegree() > 3) {
      throw ValueErrorException("invalid atom degree in getStereoInfo(bond)");
    }

    sinfo.type = StereoType::Bond_Double;
    sinfo.centeredOn = bond->getIdx();
    sinfo.controllingAtoms.reserve(4);

    bool seenSquiggleBond = false;
    const auto &mol = bond->getOwningMol();

    auto explore_bond_end = [&mol, &bond, &sinfo,
                             &seenSquiggleBond](const Atom *atom) {
      for (const auto nbr : mol.atomBonds(atom)) {
        if (nbr->getIdx() != bond->getIdx()) {
          if (nbr->getBondDir() == Bond::BondDir::UNKNOWN) {
            seenSquiggleBond = true;
          }
          sinfo.controllingAtoms.push_back(
              nbr->getOtherAtomIdx(atom->getIdx()));
        }
      }

      for (unsigned i = atom->getDegree(); i < 3; ++i) {
        sinfo.controllingAtoms.push_back(StereoInfo::NOATOM);
      }
    };

    explore_bond_end(beginAtom);
    explore_bond_end(endAtom);

    if (!seenSquiggleBond) {
      // check to see if either the begin or end atoms has the _UnknownStereo
      // property set. This happens if there was a squiggle bond to an H
      int explicitUnknownStereo = 0;
      if ((bond->getBeginAtom()->getPropIfPresent<int>(
               common_properties::_UnknownStereo, explicitUnknownStereo) &&
           explicitUnknownStereo) ||
          (bond->getEndAtom()->getPropIfPresent<int>(
               common_properties::_UnknownStereo, explicitUnknownStereo) &&
           explicitUnknownStereo)) {
        seenSquiggleBond = true;
      }
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
      if (satoms[0] == static_cast<int>(sinfo.controllingAtoms[0])) {
        firstAtBegin = true;
      } else if (satoms[0] == static_cast<int>(sinfo.controllingAtoms[1])) {
        firstAtBegin = false;
      } else {
        throw ValueErrorException("controlling atom mismatch at begin");
      }
      bool firstAtEnd;
      if (satoms[1] == static_cast<int>(sinfo.controllingAtoms[2])) {
        firstAtEnd = true;
      } else if (satoms[1] == static_cast<int>(sinfo.controllingAtoms[3])) {
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
    } else {
      sinfo.specified = Chirality::StereoSpecified::Unspecified;
    }
  } else if (bond->getBondType() == Bond::BondType::SINGLE &&
             (bond->getStereo() == Bond::BondStereo::STEREOATROPCCW ||
              bond->getStereo() == Bond::BondStereo::STEREOATROPCW)) {
    if (beginAtom->getDegree() < 2 || endAtom->getDegree() < 2 ||
        beginAtom->getDegree() > 3 || endAtom->getDegree() > 3) {
      throw ValueErrorException("invalid atom degree in getStereoInfo(bond)");
    }

    sinfo.type = StereoType::Bond_Atropisomer;
    sinfo.centeredOn = bond->getIdx();
    sinfo.controllingAtoms.reserve(4);

    const auto &mol = bond->getOwningMol();
    for (const auto nbr : mol.atomBonds(beginAtom)) {
      if (nbr->getIdx() != bond->getIdx()) {
        sinfo.controllingAtoms.push_back(
            nbr->getOtherAtomIdx(beginAtom->getIdx()));
      }
    }
    if (beginAtom->getDegree() == 2) {
      sinfo.controllingAtoms.push_back(StereoInfo::NOATOM);
    }
    for (const auto nbr : mol.atomBonds(endAtom)) {
      if (nbr->getIdx() != bond->getIdx()) {
        sinfo.controllingAtoms.push_back(
            nbr->getOtherAtomIdx(endAtom->getIdx()));
      }
    }
    if (endAtom->getDegree() == 2) {
      sinfo.controllingAtoms.push_back(StereoInfo::NOATOM);
    }

    Bond::BondStereo stereo = bond->getStereo();
    sinfo.specified = Chirality::StereoSpecified::Specified;
    switch (stereo) {
      case Bond::BondStereo::STEREOATROPCW:
        sinfo.descriptor = Chirality::StereoDescriptor::Bond_AtropCW;
        break;
      case Bond::BondStereo::STEREOATROPCCW:
        sinfo.descriptor = Chirality::StereoDescriptor::Bond_AtropCCW;
        break;
      default:
        UNDER_CONSTRUCTION("unrecognized bond stereo type");
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
    } else if (getAllowNontetrahedralChirality() &&
               isAtomPotentialNontetrahedralCenter(atom)) {
      if (stereo == Atom::CHI_UNSPECIFIED) {
        switch (atom->getTotalDegree()) {
          case 4:
            // don't assume non-tetrahedral chirality
            stereo = Atom::ChiralType::CHI_TETRAHEDRAL;
            break;
          case 5:
            stereo = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
            break;
          case 6:
            stereo = Atom::ChiralType::CHI_OCTAHEDRAL;
            break;
          default:
            break;
        }
      }
      sinfo.descriptor = StereoDescriptor::None;
      switch (stereo) {
        case Atom::ChiralType::CHI_TETRAHEDRAL:
          sinfo.type = StereoType::Atom_Tetrahedral;
          break;
        case Atom::ChiralType::CHI_SQUAREPLANAR:
          sinfo.type = StereoType::Atom_SquarePlanar;
          break;
        case Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL:
          sinfo.type = StereoType::Atom_TrigonalBipyramidal;
          break;
        case Atom::ChiralType::CHI_OCTAHEDRAL:
          sinfo.type = StereoType::Atom_Octahedral;
          break;
        default:
          break;
      }
      unsigned int permutation;
      if (atom->getPropIfPresent(common_properties::_chiralPermutation,
                                 permutation)) {
        sinfo.permutation = permutation;
        if (!permutation) {
          // a permutation of zero is an explicit statement that the chirality
          // is unknown
          sinfo.specified = Chirality::StereoSpecified::Unknown;
        } else {
          sinfo.specified = Chirality::StereoSpecified::Specified;
        }
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
  // each of the beginning and end neighbors must have at least 2 explicit
  // neighbors but no more than 3 total neighbors.
  // if it's a ring bond, the smallest ring it's in must have at least 8
  // members
  //  (this is common with InChI)
  const auto beginAtom = bond->getBeginAtom();
  auto begDegree = beginAtom->getTotalDegree();
  const auto endAtom = bond->getEndAtom();
  auto endDegree = endAtom->getTotalDegree();
  if (begDegree > 1 && begDegree < 4 && endDegree > 1 && endDegree < 4 &&
      beginAtom->getTotalNumHs(true) < 2 && endAtom->getTotalNumHs(true) < 2) {
    // check rings
    const auto ri = bond->getOwningMol().getRingInfo();
    for (const auto &bring : ri->bondRings()) {
      if (bring.size() < minRingSizeForDoubleBondStereo &&
          std::find(bring.begin(), bring.end(), bond->getIdx()) !=
              bring.end()) {
        return false;
      }
    }
    return true;
  } else {
    return false;
  }
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

namespace {
inline std::string getAtomCompareSymbol(const Atom &atom) {
  // we originally tried this with boost::format, but it was WAY slower
  return std::to_string(atom.getIsotope()) + atom.getSymbol() +
         std::to_string(atom.getFormalCharge());
}

bool areStereobondControllingAtomsDupes(
    const ROMol &mol, const Bond &bond, unsigned controllingAtom1,
    unsigned controllingAtom2, const std::vector<unsigned> &atomRanks,
    const boost::dynamic_bitset<> &possibleAtoms,
    const boost::dynamic_bitset<> &knownAtoms) {
  PRECONDITION(controllingAtom1 != Chirality::StereoInfo::NOATOM &&
                   controllingAtom2 != Chirality::StereoInfo::NOATOM,
               "Missing a controlling atom");

  if (atomRanks[controllingAtom1] != atomRanks[controllingAtom2]) {
    return false;
  }

  // Now that we know we have 2 neighbors with the same rank, check whether
  // there is a common even sized ring between the controlling atoms in which
  // the atom opposite of the bond may be chiral, which would break the tie
  auto ringInfo = mol.getRingInfo();
  auto atom1Members = ringInfo->atomMembers(controllingAtom1);
  auto atom2Members = ringInfo->atomMembers(controllingAtom2);

  auto it1 = atom1Members.begin();
  auto it2 = atom2Members.begin();
  while (it1 != atom1Members.end() && it2 != atom2Members.end()) {
    if (*it1 < *it2) {
      ++it1;
      continue;
    } else if (*it1 > *it2) {
      ++it2;
      continue;
    }

    auto ring = ringInfo->atomRings().at(*it1);
    ++it1;
    ++it2;

    if (ring.size() % 2) {
      // The common ring is odd-sized, so we can't have a tie-breaking atom
      // directly across the ring, so skip this ring.
      continue;
    }

    for (auto bondEnd : {bond.getBeginAtomIdx(), bond.getEndAtomIdx()}) {
      auto bondEndPosItr = std::find(ring.begin(), ring.end(), bondEnd);
      if (bondEndPosItr != ring.end()) {
        auto bondEndPos = bondEndPosItr - ring.begin();
        auto oppositePos = (bondEndPos + ring.size() / 2) % ring.size();

        auto oppositeIdx = ring[oppositePos];
        if (possibleAtoms[oppositeIdx] || knownAtoms[oppositeIdx]) {
          return false;
        }
      }
    }
  }

  return true;
}

}  // namespace

namespace {
void initAtomInfo(ROMol &mol, bool flagPossible, bool cleanIt,
                  boost::dynamic_bitset<> &knownAtoms,
                  std::vector<std::string> &atomSymbols,
                  boost::dynamic_bitset<> &possibleAtoms) {
  bool allowNontetrahedralStereo = getAllowNontetrahedralChirality();
  for (const auto atom : mol.atoms()) {
    if (atom->needsUpdatePropertyCache()) {
      atom->updatePropertyCache(false);
    }

    auto aidx = atom->getIdx();
    atomSymbols[aidx] = getAtomCompareSymbol(*atom);
    if (detail::isAtomPotentialStereoAtom(atom, allowNontetrahedralStereo)) {
      auto sinfo = detail::getStereoInfo(atom);
      switch (sinfo.specified) {
        case Chirality::StereoSpecified::Unknown:
          knownAtoms.set(aidx);
          atomSymbols[aidx] += std::to_string(aidx);
          break;
        case Chirality::StereoSpecified::Specified:
          knownAtoms.set(aidx);
          if (sinfo.descriptor == StereoDescriptor::Tet_CCW) {
            atomSymbols[aidx] += "_CCW";
          } else if (sinfo.descriptor == StereoDescriptor::Tet_CW) {
            atomSymbols[aidx] += "_CW";
          } else {
            atomSymbols[aidx] += "_STEREO";
          }
          break;
        case Chirality::StereoSpecified::Unspecified:
          if (flagPossible) {
            possibleAtoms.set(aidx);
            if (!cleanIt) {
              atomSymbols[aidx] += "_" + std::to_string(aidx);
            }
          }
          break;
        default:
          throw ValueErrorException("bad StereoInfo.specified type");
      }
    } else if (cleanIt) {
      atom->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
    }
  }
}

void initBondInfo(ROMol &mol, bool flagPossible, bool cleanIt,
                  boost::dynamic_bitset<> &knownBonds,
                  std::vector<std::string> &bondSymbols,
                  boost::dynamic_bitset<> &possibleBonds) {
  for (const auto bond : mol.bonds()) {
    auto bidx = bond->getIdx();
    bondSymbols[bidx] = getBondSymbol(bond);
    if (detail::isBondPotentialStereoBond(bond)) {
      auto sinfo = detail::getStereoInfo(bond);
      switch (sinfo.specified) {
        case Chirality::StereoSpecified::Unknown:
          knownBonds.set(bidx);
          bondSymbols[bidx] += "_" + std::to_string(bidx);
          break;
        case Chirality::StereoSpecified::Specified:
          knownBonds.set(bidx);
          if (sinfo.descriptor == StereoDescriptor::Bond_Cis) {
            bondSymbols[bidx] += "_cis";
          } else if (sinfo.descriptor == StereoDescriptor::Bond_Trans) {
            bondSymbols[bidx] += "_trans";
          } else {
            bondSymbols[bidx] += "_STEREO";
          }
          break;
        case Chirality::StereoSpecified::Unspecified:
          if (flagPossible) {
            possibleBonds.set(bidx);
            if (!cleanIt) {
              bondSymbols[bidx] += "_" + std::to_string(bidx);
            }
          }
          break;
        default:
          throw ValueErrorException("bad StereoInfo.specified type");
      }
    } else {
      auto currentStereo = bond->getStereo();
      if (currentStereo != Bond::BondStereo::STEREOATROPCW &&
          currentStereo != Bond::BondStereo::STEREOATROPCCW) {
        if (cleanIt) {
          bond->setStereo(Bond::BondStereo::STEREONONE);
        }
      } else {
        knownBonds.set(bidx);
        if (currentStereo == Bond::BondStereo::STEREOATROPCW) {
          bondSymbols[bidx] += "_atropcw";
        } else if (currentStereo == Bond::BondStereo::STEREOATROPCCW) {
          bondSymbols[bidx] += "_atropccw";
        }
      }
    }
  }
}
void flagRingStereo(ROMol &mol,
                    std::vector<unsigned int> &possibleRingStereoAtoms,
                    std::vector<unsigned int> &possibleRingStereoBonds,
                    const boost::dynamic_bitset<> &knownAtoms,
                    const boost::dynamic_bitset<> *possibleAtoms,
                    const boost::dynamic_bitset<> &knownBonds,
                    const boost::dynamic_bitset<> *possibleBonds) {
  // flag possible ring stereo cases. The relevant cases here are:
  //    1) even-sized rings with possible (or specified) atoms opposite each
  //       other, like CC1CC(C)C1 or CC1CCC(C)CC1
  //    2) atoms sharing a bond which fuses two or more rings, like the
  //    central
  //       bond in C1CCC2CCCCC2C1

  auto ringInfo = mol.getRingInfo();
  boost::dynamic_bitset<> possibleAtomsInRing(mol.getNumAtoms());
  for (unsigned int ridx = 0; ridx < ringInfo->atomRings().size(); ++ridx) {
    const auto &aring = ringInfo->atomRings()[ridx];
    const auto &bring = ringInfo->bondRings()[ridx];
    unsigned int nHere = 0;
    auto sz = aring.size();
    bool ringIsOddSized = sz % 2;
    auto halfSize = sz / 2 + ringIsOddSized;

    possibleAtomsInRing.reset();
    for (unsigned int ai = 0; ai < sz; ++ai) {
      auto aidx = aring[ai];
      if (!knownAtoms[aidx] && (!possibleAtoms || !possibleAtoms->test(aidx))) {
        continue;
      }
      if (!ringIsOddSized) {
        // find the index of the atom on the opposite side of the even-sized
        // ring
        auto oppositeIdx = aring[(ai + halfSize) % sz];
        bool toAtomOppositePossible = false;
        auto oppositeAtom = mol.getAtomWithIdx(oppositeIdx);
        for (auto bond : mol.atomBonds(oppositeAtom)) {
          auto bidx = bond->getIdx();
          if ((knownBonds[bidx] ||
               (possibleBonds && possibleBonds->test(bidx))) &&
              std::find(bring.begin(), bring.end(), bidx) == bring.end()) {
            toAtomOppositePossible = true;
            break;
          }
        }

        if (knownAtoms[oppositeIdx] ||
            (possibleAtoms && possibleAtoms->test(oppositeIdx)) ||
            toAtomOppositePossible) {
          nHere += 1 + toAtomOppositePossible;
          possibleAtomsInRing.set(aidx);
          possibleAtomsInRing.set(oppositeIdx);
          continue;
        }
      }

      // if the atom is in more than one ring, explore the common edge to see
      // if we can find another potentially chiral atom
      if (ringInfo->numAtomRings(aidx) > 1) {
        auto previousOtherIdx = aidx;
        for (size_t step = 1; step <= halfSize; ++step) {
          auto otherIdx = aring[(ai + step) % sz];
          auto bnd = mol.getBondBetweenAtoms(previousOtherIdx, otherIdx);
          if (ringInfo->numBondRings(bnd->getIdx()) < 2) {
            // We reached the end of the common edge.
            break;
          }
          if (knownAtoms[otherIdx] ||
              (possibleAtoms && possibleAtoms->test(otherIdx))) {
            // We found another chiral atom, no need to keep
            // searching.
            nHere += 2;
            possibleAtomsInRing.set(aidx);
            possibleAtomsInRing.set(otherIdx);
            break;
          }
          previousOtherIdx = otherIdx;
        }
      }
    }
    // if the ring contains at least two atoms with possible stereo,
    // then each of those possibleAtoms should be included for ring stereo
    if (nHere > 1) {
      for (auto aidx : aring) {
        if (possibleAtomsInRing[aidx]) {
          ++possibleRingStereoAtoms[aidx];
        }
      }
      for (auto bidx : bring) {
        ++possibleRingStereoBonds[bidx];
      }
    }
  }
}

bool updateAtoms(ROMol &mol, const std::vector<unsigned int> &aranks,
                 std::vector<std::string> &atomSymbols,
                 boost::dynamic_bitset<> &possibleAtoms,
                 boost::dynamic_bitset<> &knownAtoms,
                 boost::dynamic_bitset<> &fixedAtoms,
                 std::vector<unsigned int> &possibleRingStereoAtoms,
                 std::vector<unsigned int> &possibleRingStereoBonds,
                 std::vector<StereoInfo> &sinfos) {
  bool needAnotherRound = false;
  for (const auto atom : mol.atoms()) {
    auto aidx = atom->getIdx();
    if (knownAtoms[aidx] || possibleAtoms[aidx]) {
      auto sinfo = detail::getStereoInfo(atom);
      if (fixedAtoms[aidx]) {
        sinfos.push_back(std::move(sinfo));
      } else {
        std::vector<unsigned int> nbrs;
        nbrs.reserve(sinfo.controllingAtoms.size());
        bool haveADupe = false;
        if (sinfo.type == StereoType::Atom_Tetrahedral) {
          for (auto nbrIdx : sinfo.controllingAtoms) {
            auto rnk = aranks[nbrIdx];
            if (std::find(nbrs.begin(), nbrs.end(), rnk) != nbrs.end()) {
              // ok, we just hit a duplicate rank. If the atom we're concerned
              // about is a candidate for ring stereo and the bond to the atom
              // with the duplicate rank is a ring bond
              if (possibleRingStereoAtoms[aidx]) {
                auto bnd = mol.getBondBetweenAtoms(aidx, nbrIdx);
                if (!bnd || !possibleRingStereoBonds[bnd->getIdx()]) {
                  haveADupe = true;
                  break;
                }
              } else {
                haveADupe = true;
                break;
              }
            } else {
              nbrs.push_back(rnk);
            }
          }
        }
        if (!haveADupe) {
          // std::cerr << "NBRS from " << aidx << ": ";
          // std::copy(sinfo.controllingAtoms.begin(),
          //           sinfo.controllingAtoms.end(),
          //           std::ostream_iterator<int>(std::cerr, " "));
          // std::cerr << std::endl;
          // std::copy(nbrs.begin(), nbrs.end(),
          //           std::ostream_iterator<int>(std::cerr, " "));
          // std::cerr << std::endl;

          auto acs = atomSymbols[aidx];
          if (!possibleAtoms[aidx]) {
            auto sortednbrs = nbrs;
            std::sort(sortednbrs.begin(), sortednbrs.end());
            // FIX: only works for tetrahedral at the moment
            if (sinfo.type == Chirality::StereoType::Atom_Tetrahedral) {
              auto nSwaps = countSwapsToInterconvert(nbrs, sortednbrs);
              if (nSwaps % 2 &&
                  (sinfo.descriptor == Chirality::StereoDescriptor::Tet_CCW ||
                   sinfo.descriptor == Chirality::StereoDescriptor::Tet_CW)) {
                sinfo.descriptor =
                    sinfo.descriptor == Chirality::StereoDescriptor::Tet_CCW
                        ? Chirality::StereoDescriptor::Tet_CW
                        : Chirality::StereoDescriptor::Tet_CCW;
              }
              if (sinfo.descriptor == Chirality::StereoDescriptor::Tet_CW) {
                acs = getAtomCompareSymbol(*atom) + "_CW";
              } else if (sinfo.descriptor ==
                         Chirality::StereoDescriptor::Tet_CCW) {
                acs = getAtomCompareSymbol(*atom) + "_CCW";
              }
            }
            fixedAtoms.set(aidx);
          }
          if (atomSymbols[aidx] != acs) {
            atomSymbols[aidx] = acs;
            needAnotherRound = true;
          }
          sinfos.push_back(std::move(sinfo));
        } else if (possibleAtoms[aidx]) {
          possibleAtoms[aidx] = 0;
          atomSymbols[aidx] = getAtomCompareSymbol(*atom);
          needAnotherRound = true;

          // if this was creating possible ring stereo, update that info now
          if (possibleRingStereoAtoms[aidx]) {
            --possibleRingStereoAtoms[aidx];
            if (!possibleRingStereoAtoms[aidx]) {
              // we're no longer in any ring with possible ring stereo. Go
              // update all the other atoms/bonds in rings that we're in:
              for (unsigned int ridx = 0;
                   ridx < mol.getRingInfo()->atomRings().size(); ++ridx) {
                const auto &aring = mol.getRingInfo()->atomRings()[ridx];
                unsigned int nHere = 0;
                for (auto raidx : aring) {
                  if (possibleRingStereoAtoms[raidx]) {
                    --possibleRingStereoAtoms[raidx];
                    if (possibleRingStereoAtoms[raidx]) {
                      ++nHere;
                    }
                  }
                }
                if (nHere <= 1) {
                  // update the bondstereo counts too
                  for (auto rbidx : mol.getRingInfo()->bondRings()[ridx]) {
                    if (possibleRingStereoBonds[rbidx]) {
                      --possibleRingStereoBonds[rbidx];
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return needAnotherRound;
}

bool updateBonds(ROMol &mol, const std::vector<unsigned int> &aranks,
                 std::vector<std::string> &bondSymbols,
                 const boost::dynamic_bitset<> &possibleAtoms,
                 boost::dynamic_bitset<> &possibleBonds,
                 const boost::dynamic_bitset<> &knownAtoms,
                 const boost::dynamic_bitset<> &knownBonds,
                 boost::dynamic_bitset<> &fixedBonds,
                 std::vector<StereoInfo> &sinfos) {
  bool needAnotherRound = false;
  for (const auto bond : mol.bonds()) {
    auto bidx = bond->getIdx();
    if (knownBonds[bidx] || possibleBonds[bidx]) {
      auto sinfo = detail::getStereoInfo(bond);
      ASSERT_INVARIANT(sinfo.controllingAtoms.size() == 4,
                       "bad controlling atoms size");

      if ((sinfo.controllingAtoms[0] == Chirality::StereoInfo::NOATOM &&
           sinfo.controllingAtoms[1] == Chirality::StereoInfo::NOATOM) ||
          (sinfo.controllingAtoms[2] == Chirality::StereoInfo::NOATOM &&
           sinfo.controllingAtoms[3] == Chirality::StereoInfo::NOATOM)) {
        // we have a bond with no neighbors on one side, which means it must
        // have a single implicit H on that side. Since the H is implicit,
        // there is no way to know whether it is cis or trans.
        ASSERT_INVARIANT(
            sinfo.specified != StereoSpecified::Specified,
            "stereo bond without neighbors can only be unspecified");
        fixedBonds.set(bidx);
      }

      if (fixedBonds[bidx]) {
        sinfos.push_back(std::move(sinfo));
      } else {
        bool haveADupe = false;
        bool needsSwap = false;
        if (sinfo.controllingAtoms[0] != Chirality::StereoInfo::NOATOM &&
            sinfo.controllingAtoms[1] != Chirality::StereoInfo::NOATOM) {
          if (areStereobondControllingAtomsDupes(
                  mol, *bond, sinfo.controllingAtoms[0],
                  sinfo.controllingAtoms[1], aranks, possibleAtoms,
                  knownAtoms)) {
            haveADupe = true;
          } else if (aranks[sinfo.controllingAtoms[0]] <
                     aranks[sinfo.controllingAtoms[1]]) {
            std::swap(sinfo.controllingAtoms[0], sinfo.controllingAtoms[1]);
            needsSwap = !needsSwap;
          }
        }
        if (sinfo.controllingAtoms[2] != Chirality::StereoInfo::NOATOM &&
            sinfo.controllingAtoms[3] != Chirality::StereoInfo::NOATOM) {
          if (areStereobondControllingAtomsDupes(
                  mol, *bond, sinfo.controllingAtoms[2],
                  sinfo.controllingAtoms[3], aranks, possibleAtoms,
                  knownAtoms)) {
            haveADupe = true;
          } else if (aranks[sinfo.controllingAtoms[2]] <
                     aranks[sinfo.controllingAtoms[3]]) {
            std::swap(sinfo.controllingAtoms[2], sinfo.controllingAtoms[3]);
            needsSwap = !needsSwap;
          }
        }
        if (!haveADupe) {
          if (needsSwap && (sinfo.descriptor == StereoDescriptor::Bond_Cis ||
                            sinfo.descriptor == StereoDescriptor::Bond_Trans)) {
            sinfo.descriptor = sinfo.descriptor == StereoDescriptor::Bond_Cis
                                   ? StereoDescriptor::Bond_Trans
                                   : StereoDescriptor::Bond_Cis;
          }
          auto gbs = bondSymbols[bidx];
          if (sinfo.specified == StereoSpecified::Specified) {
            switch (sinfo.descriptor) {
              case StereoDescriptor::Bond_Cis:
                gbs += "_cis";
                break;
              case StereoDescriptor::Bond_Trans:
                gbs += "_trans";
                break;
              default:
                break;
            }
          } else if (sinfo.specified == StereoSpecified::Unknown) {
            gbs += "_unk";
          }
          if (bondSymbols[bidx] != gbs) {
            bondSymbols[bidx] = gbs;
            needAnotherRound = true;
          }
          if (!possibleBonds[bidx]) {
            fixedBonds.set(bidx);
          }
          sinfos.push_back(std::move(sinfo));
        } else if (possibleBonds[bidx]) {
          possibleBonds[bidx] = 0;
          bondSymbols[bidx] = getBondSymbol(bond);
          needAnotherRound = true;
        }
      }
    }
  }
  return needAnotherRound;
}

void cleanMolStereo(ROMol &mol, const boost::dynamic_bitset<> &fixedAtoms,
                    const boost::dynamic_bitset<> &knownAtoms,
                    const boost::dynamic_bitset<> &fixedBonds,
                    const boost::dynamic_bitset<> &knownBonds) {
  for (auto atom : mol.atoms()) {
    const auto i = atom->getIdx();
    if (!fixedAtoms[i] && knownAtoms[i]) {
      switch (atom->getChiralTag()) {
        case Atom::ChiralType::CHI_TETRAHEDRAL_CCW:
        case Atom::ChiralType::CHI_TETRAHEDRAL_CW:
          atom->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
          for (auto nbrBond : mol.atomBonds(atom)) {
            auto bondDir = nbrBond->getBondDir();
            if (bondDir == Bond::BondDir::BEGINDASH ||
                bondDir == Bond::BondDir::BEGINWEDGE) {
              nbrBond->setBondDir(Bond::BondDir::NONE);
            }
          }
          break;
        case Atom::ChiralType::CHI_TETRAHEDRAL:
        case Atom::ChiralType::CHI_SQUAREPLANAR:
        case Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL:
        case Atom::ChiralType::CHI_OCTAHEDRAL:
          atom->setProp(common_properties::_chiralPermutation, 0);
          break;
        default:
          break;
      }
    }
  }

  bool removedStereo = false;
  for (auto bond : mol.bonds()) {
    const auto i = bond->getIdx();
    if (!fixedBonds[i] && knownBonds[i]) {
      bond->setStereo(Bond::BondStereo::STEREONONE);
      bond->setBondDir(Bond::BondDir::NONE);
      bond->getStereoAtoms().clear();
      removedStereo = true;
    }
  }

  // remove any slash bond dirs that do not have a stereo neighbor

  if (removedStereo) {
    for (auto bond : mol.bonds()) {
      auto bondDir = bond->getBondDir();
      if (bondDir == Bond::BondDir::ENDDOWNRIGHT ||
          bondDir == Bond::BondDir::ENDUPRIGHT) {
        bool dirOk = false;
        for (auto bondEnd : {bond->getBeginAtom(), bond->getEndAtom()}) {
          for (auto nbrBond : mol.atomBonds(bondEnd)) {
            if (nbrBond != bond &&
                nbrBond->getStereo() != Bond::BondStereo::STEREONONE) {
              dirOk = true;
              break;
            }
          }
          if (!dirOk) {
            bond->setBondDir(Bond::BondDir::NONE);
          }
        }
      }
    }
  }
}
}  // namespace

std::vector<StereoInfo> runCleanup(ROMol &mol, bool flagPossible,
                                   bool cleanIt) {
  // This potentially does two passes of "canonicalization" to identify stereo
  // atoms/bonds:
  //   - if cleanIt is true we start with a pass which ignores possible stereo
  //   atoms/bonds and which removes the stereo spec from any atom/bond which
  //   doesn't have unique neighbors
  //   - if flagPossible is true we do a pass where each unspecified possible
  //   stereocenter is treated as if it were different from all others. This
  //   allows us to identify every possible stereo atom/bond

  boost::dynamic_bitset<> knownAtoms(mol.getNumAtoms());
  std::vector<std::string> atomSymbols(mol.getNumAtoms());
  boost::dynamic_bitset<> possibleAtoms(mol.getNumAtoms());
  initAtomInfo(mol, flagPossible, cleanIt, knownAtoms, atomSymbols,
               possibleAtoms);

  std::vector<std::string> bondSymbols(mol.getNumBonds());
  boost::dynamic_bitset<> knownBonds(mol.getNumBonds());
  boost::dynamic_bitset<> possibleBonds(mol.getNumBonds());
  initBondInfo(mol, flagPossible, cleanIt, knownBonds, bondSymbols,
               possibleBonds);

  // copy the original sets of possible atoms/bonds. We need them in the second
  // pass
  auto origPossibleAtoms = possibleAtoms;
  auto origPossibleBonds = possibleBonds;

  // tracks the number of rings with possible ring stereo that the atom is in
  //  (only set for potential stereoatoms)
  std::vector<unsigned int> possibleRingStereoAtoms(mol.getNumAtoms());
  // tracks the number of rings with possible ring stereo that the bond is in
  //  (set for all bonds)
  std::vector<unsigned int> possibleRingStereoBonds(mol.getNumBonds());

  // identify atoms which can be involved in ring stereo
  flagRingStereo(mol, possibleRingStereoAtoms, possibleRingStereoBonds,
                 knownAtoms, cleanIt ? nullptr : &possibleAtoms, knownBonds,
                 cleanIt ? nullptr : &possibleBonds);

  // our return value
  std::vector<StereoInfo> res;

  // these are used to track which atoms/bonds have been altered
  boost::dynamic_bitset<> fixedAtoms(mol.getNumAtoms());
  boost::dynamic_bitset<> fixedBonds(mol.getNumBonds());

  // used to tell rankFragmentAtoms to use all atoms:
  boost::dynamic_bitset<> atomsInPlay(mol.getNumAtoms());
  atomsInPlay.set();
  boost::dynamic_bitset<> bondsInPlay(mol.getNumBonds());
  bondsInPlay.set();

#define LOCAL_CANON 0
#if LOCAL_CANON
  std::vector<Canon::canon_atom> canonAtoms(mol.getNumAtoms());
  Canon::AtomCompareFunctor ftor(&canonAtoms.front(), mol, &atomsInPlay,
                                 &bondsInPlay);
  ftor.df_useIsotopes = false;
  ftor.df_useChirality = false;
  auto atomOrder = new int[mol.getNumAtoms()];
#endif
  std::vector<unsigned int> aranks(mol.getNumAtoms());
  bool needAnotherRound = true;
  while (needAnotherRound) {
    res.clear();
#if LOCAL_CANON
    // find symmetry classes with the canonicalization code
    Canon::detail::initFragmentCanonAtoms(mol, canonAtoms, false, &atomSymbols,
                                          &bondSymbols, atomsInPlay,
                                          bondsInPlay, needsInit);
    needsInit = false;

    const bool includeChirality = false;
    const bool includeIsotopes = false;
    const bool breakTies = false;
    memset(atomOrder, 0, mol.getNumAtoms() * sizeof(int));
    Canon::detail::rankWithFunctor(ftor, breakTies, atomOrder, true,
                                   includeChirality, &atomsInPlay,
                                   &bondsInPlay);
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      aranks[atomOrder[i]] = canonAtoms[atomOrder[i]].index;
    }
#else
    const bool includeChirality = false;
    const bool includeIsotopes = false;
    const bool breakTies = false;
    const bool includeAtomMaps = false;
    // Now apply the canonical atom ranking code with basic connectivity
    // invariants The necessary condition for chirality is that an atom's
    // neighbors must have unique ranks
    Canon::rankFragmentAtoms(
        mol, aranks, atomsInPlay, bondsInPlay, &atomSymbols, &bondSymbols,
        breakTies, includeChirality, includeIsotopes, includeAtomMaps);
#endif
    // check if any new atoms definitely now have stereo; do another loop if so
    needAnotherRound = updateAtoms(
        mol, aranks, atomSymbols, possibleAtoms, knownAtoms, fixedAtoms,
        possibleRingStereoAtoms, possibleRingStereoBonds, res);
    // check if any new bonds definitely now have stereo; do another loop if so
    needAnotherRound |=
        updateBonds(mol, aranks, bondSymbols, possibleAtoms, possibleBonds,
                    knownAtoms, knownBonds, fixedBonds, res);
  }

  if (cleanIt) {
    // remove stereo specs from atoms/bonds which should not have them
    cleanMolStereo(mol, fixedAtoms, knownAtoms, fixedBonds, knownBonds);

    if (flagPossible && (possibleAtoms != origPossibleAtoms ||
                         possibleBonds != origPossibleBonds)) {
      // if we're doing "flagPossible" mode and have done some cleanup, then we
      // need to do another iteration

      possibleAtoms = origPossibleAtoms;
      // flag every center/bond where we removed stereo as possible:
      for (auto i = 0u; i < mol.getNumAtoms(); ++i) {
        if (!fixedAtoms[i] && knownAtoms[i]) {
          possibleAtoms[i] = 1;
          knownAtoms[i] = 0;
        }
        if (possibleAtoms[i]) {
          atomSymbols[i] += "_" + std::to_string(i);
        }
      }
      possibleBonds = origPossibleBonds;
      for (auto i = 0u; i < mol.getNumBonds(); ++i) {
        if (!fixedBonds[i] && knownBonds[i]) {
          possibleBonds[i] = 1;
          knownBonds[i] = 0;
        }
        if (possibleBonds[i]) {
          bondSymbols[i] += "_" + std::to_string(i);
        }
      }

      flagRingStereo(mol, possibleRingStereoAtoms, possibleRingStereoBonds,
                     knownAtoms, &possibleAtoms, knownBonds, &possibleBonds);

      needAnotherRound = true;
      while (needAnotherRound) {
        res.clear();

        // std::copy(atomSymbols.begin(), atomSymbols.end(),
        //           std::ostream_iterator<std::string>(std::cerr, " "));
        // std::cerr << std::endl;
        // std::copy(bondSymbols.begin(), bondSymbols.end(),
        //           std::ostream_iterator<std::string>(std::cerr, " "));
        // std::cerr << std::endl;

#if LOCAL_CANON
        Canon::detail::initFragmentCanonAtoms(mol, canonAtoms, false,
                                              &atomSymbols, &bondSymbols,
                                              atomsInPlay, bondsInPlay, true);
        needsInit = false;

        const bool includeChirality = false;
        const bool breakTies = false;
        memset(atomOrder, 0, mol.getNumAtoms() * sizeof(int));
        Canon::detail::rankWithFunctor(ftor, breakTies, atomOrder, true,
                                       includeChirality, &atomsInPlay,
                                       &bondsInPlay);
        for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
          aranks[atomOrder[i]] = canonAtoms[atomOrder[i]].index;
        }

#else
        // we will use the canonicalization code
        const bool breakTies = false;
        const bool includeChirality = false;
        const bool includeIsotopes = false;
        const bool includeAtomMaps = false;
        Canon::rankFragmentAtoms(
            mol, aranks, atomsInPlay, bondsInPlay, &atomSymbols, &bondSymbols,
            breakTies, includeChirality, includeIsotopes, includeAtomMaps);
#endif
        // fixedAtoms.reset();
        // fixedBonds.reset();
        needAnotherRound = updateAtoms(
            mol, aranks, atomSymbols, possibleAtoms, knownAtoms, fixedAtoms,
            possibleRingStereoAtoms, possibleRingStereoBonds, res);
        needAnotherRound |=
            updateBonds(mol, aranks, bondSymbols, possibleAtoms, possibleBonds,
                        knownAtoms, knownBonds, fixedBonds, res);
      }
    }

    for (const auto atom : mol.atoms()) {
      atom->setProp<unsigned int>(common_properties::_ChiralAtomRank,
                                  aranks[atom->getIdx()]);
    }
  }

#if LOCAL_CANON
  delete[] atomOrder;
  Canon::detail::freeCanonAtoms(canonAtoms);
#endif
  return res;
}

//}  // namespace

std::vector<StereoInfo> findPotentialStereo(ROMol &mol, bool cleanIt,
                                            bool findPossible) {
  if (!mol.getRingInfo()->isSymmSssr()) {
    MolOps::symmetrizeSSSR(mol);
  }

  if (mol.needsUpdatePropertyCache()) {
    mol.updatePropertyCache(false);
  }
  std::vector<StereoInfo> res = runCleanup(mol, findPossible, cleanIt);
  mol.setProp("_potentialStereo", res, true);
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
