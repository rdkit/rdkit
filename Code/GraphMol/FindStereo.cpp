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
#include <boost/format.hpp>
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
  auto tnzdegree =
      Chirality::detail::getAtomNonzeroDegree(atom) + atom->getTotalNumHs();
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
    } else if (nzDegree == 1) {
      // chirality is never possible with 1 nbr
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

bool isAtomPotentialStereoAtom(const Atom *atom) {
  return isAtomPotentialTetrahedralCenter(atom) ||
         isAtomPotentialNontetrahedralCenter(atom);
}

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
    } else if (isAtomPotentialNontetrahedralCenter(atom)) {
      if (stereo == Atom::CHI_UNSPECIFIED) {
        switch (atom->getTotalDegree()) {
          case 4:
            stereo = Atom::ChiralType::CHI_SQUAREPLANAR;
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
  // each of the beginning and end neighbors must have at least 2 heavy atom
  // neighbors i.e. C/C=N/[H] is not a possible stereo bond but no more than 3
  // total neighbors.
  // if it's a ring bond, the smallest ring it's in must have at least 8
  // members
  //  (this is common with InChI)
  const auto beginAtom = bond->getBeginAtom();
  auto begHeavyDegree =
      beginAtom->getTotalDegree() - beginAtom->getTotalNumHs(true);
  const auto endAtom = bond->getEndAtom();
  auto endHeavyDegree =
      endAtom->getTotalDegree() - endAtom->getTotalNumHs(true);
  if (begHeavyDegree > 1 && beginAtom->getDegree() < 4 && endHeavyDegree > 1 &&
      endAtom->getDegree() < 4) {
    // check rings
    const auto ri = bond->getOwningMol().getRingInfo();
    for (const auto &bring : ri->bondRings()) {
      if (bring.size() < 8 && std::find(bring.begin(), bring.end(),
                                        bond->getIdx()) != bring.end()) {
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
std::string getAtomCompareSymbol(const Atom &atom) {
  auto fmt = boost::format("%d%s_%d") % atom.getIsotope() % atom.getSymbol() %
             atom.getFormalCharge();
  return fmt.str();
}
}  // namespace

#if 0
std::vector<StereoInfo> findPotentialStereo(ROMol &mol, bool cleanIt,
                                            bool flagPossible) {
  std::map<int, Atom::ChiralType> ochiralTypes;

  if (!mol.getRingInfo()->isInitialized()) {
    MolOps::symmetrizeSSSR(mol);
  }
  if (mol.needsUpdatePropertyCache()) {
    mol.updatePropertyCache(false);
  }

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
        ochiralTypes[aidx] = atom->getChiralTag();
        atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
        atomSymbols[aidx] = (boost::format("%d%s_%d") % atom->getIsotope() %
                             atom->getSymbol() % aidx)
                                .str();
      } else {
        atomSymbols[aidx] = getAtomCompareSymbol(*atom);
      }
    } else {
      atomSymbols[aidx] = getAtomCompareSymbol(*atom);
      if (cleanIt) {
        atom->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
      }
    }
  }

  // flag possible ring stereo cases. The relevant cases here are:
  //    1) even-sized rings with possible (or specified) atoms opposite each
  //       other, like CC1CC(C)C1 or CC1CCC(C)CC1
  //    2) atoms sharing a bond which fuses two or more rings, like the central
  //       bond in C1CCC2CCCCC2C1

  // tracks the number of rings with possible ring stereo that the atom is in
  //  (only set for potential stereoatoms)
  std::vector<unsigned int> possibleRingStereoAtoms(mol.getNumAtoms());
  // tracks the number of rings with possible ring stereo that the bond is in
  //  (set for all bonds)
  std::vector<unsigned int> possibleRingStereoBonds(mol.getNumBonds());
  boost::dynamic_bitset<> possibleAtomsInRing(mol.getNumAtoms());
  for (unsigned int ridx = 0; ridx < mol.getRingInfo()->atomRings().size();
       ++ridx) {
    const auto &aring = mol.getRingInfo()->atomRings()[ridx];
    unsigned int nHere = 0;
    auto sz = aring.size();
    possibleAtomsInRing.reset();
    for (unsigned int ai = 0; ai < aring.size(); ++ai) {
      auto aidx = aring[ai];
      if (!(aring.size() % 2)) {
        // find the index of the atom on the opposite side of the even-sized
        // ring
        auto oppositeidx = aring[(ai + sz / 2) % sz];
        if ((possibleAtoms[aidx] || knownAtoms[aidx]) &&
            (possibleAtoms[oppositeidx] || knownAtoms[oppositeidx])) {
          ++nHere;
          possibleAtomsInRing.set(aidx);
          continue;
        }
      }
      // if the atom is in more than one bond, see if there's
      // a possible neighbor on a fusion bond
      if (mol.getRingInfo()->numAtomRings(aidx) > 1) {
        auto otheridx = aring[(ai + 1) % aring.size()];
        if (possibleAtoms[otheridx] || knownAtoms[otheridx]) {
          auto bnd = mol.getBondBetweenAtoms(aidx, otheridx);
          CHECK_INVARIANT(bnd, "expected ring bond not found");
          if (mol.getRingInfo()->numBondRings(bnd->getIdx()) > 1) {
            nHere += 2;
            possibleAtomsInRing.set(aidx);
            possibleAtomsInRing.set(otheridx);
          }
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
      for (auto bidx : mol.getRingInfo()->bondRings()[ridx]) {
        ++possibleRingStereoBonds[bidx];
      }
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
  boost::dynamic_bitset<> fixedAtoms(mol.getNumAtoms());
  std::vector<StereoInfo> res;
  while (possibleAtoms.count() || possibleBonds.count()) {
    res.clear();
    bool needAnotherRound = false;

    // std::copy(atomSymbols.begin(), atomSymbols.end(),
    //           std::ostream_iterator<std::string>(std::cerr, " "));
    // std::cerr << std::endl;
    // std::copy(bondSymbols.begin(), bondSymbols.end(),
    //           std::ostream_iterator<std::string>(std::cerr, " "));
    // std::cerr << std::endl;

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
    const bool includeChirality = false;
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
        if (fixedAtoms[aidx]) {
          res.push_back(std::move(sinfo));
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
            if (knownAtoms[aidx]) {
              // std::cerr << "NBRS from " << aidx << ": ";
              // std::copy(sinfo.controllingAtoms.begin(),
              //           sinfo.controllingAtoms.end(),
              //           std::ostream_iterator<int>(std::cerr, " "));
              // std::cerr << std::endl;
              // std::copy(nbrs.begin(), nbrs.end(),
              //           std::ostream_iterator<int>(std::cerr, " "));
              // std::cerr << std::endl;

              auto acs = getAtomCompareSymbol(*atom);
              auto sortednbrs = nbrs;
              std::sort(sortednbrs.begin(), sortednbrs.end());
              // FIX: only works for tetrahedral at the moment
              if (sinfo.type == Chirality::StereoType::Atom_Tetrahedral) {
                auto nSwaps = countSwapsToInterconvert(nbrs, sortednbrs);
                if (nSwaps % 2) {
                  sinfo.descriptor =
                      sinfo.descriptor == Chirality::StereoDescriptor::Tet_CCW
                          ? Chirality::StereoDescriptor::Tet_CW
                          : Chirality::StereoDescriptor::Tet_CCW;
                }
                if (sinfo.descriptor == Chirality::StereoDescriptor::Tet_CW) {
                  acs += "_CW";
                } else if (sinfo.descriptor ==
                           Chirality::StereoDescriptor::Tet_CCW) {
                  acs += "_CCW";
                }
              }
              if (atomSymbols[aidx] != acs) {
                atomSymbols[aidx] = acs;
                needAnotherRound = true;
                fixedAtoms.set(aidx);
              }
            }
            res.push_back(std::move(sinfo));
          } else {
            needAnotherRound = true;
            atomSymbols[aidx] = getAtomCompareSymbol(*atom);
            possibleAtoms[aidx] = 0;
            if (cleanIt &&
                sinfo.specified != Chirality::StereoSpecified::Unspecified) {
              atom->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
            }
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
      } else {
        atom->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
      }
    }
    boost::dynamic_bitset<> fixedBonds(mol.getNumBonds());

    for (const auto bond : mol.bonds()) {
      auto bidx = bond->getIdx();
      if (possibleBonds[bidx]) {
        bond->clearProp(Chirality::_stereoNotPossible);

        auto sinfo = detail::getStereoInfo(bond);
        ASSERT_INVARIANT(sinfo.controllingAtoms.size() == 4,
                         "bad controlling atoms size");
        if (fixedBonds[bidx]) {
          res.push_back(std::move(sinfo));
        } else {
          bool haveADupe = false;
          bool needsSwap = false;
          if (sinfo.controllingAtoms[0] != Chirality::StereoInfo::NOATOM &&
              sinfo.controllingAtoms[1] != Chirality::StereoInfo::NOATOM) {
            if (aranks[sinfo.controllingAtoms[0]] ==
                aranks[sinfo.controllingAtoms[1]]) {
              haveADupe = true;
            } else if (aranks[sinfo.controllingAtoms[0]] <
                       aranks[sinfo.controllingAtoms[1]]) {
              std::swap(sinfo.controllingAtoms[0], sinfo.controllingAtoms[1]);
              needsSwap = !needsSwap;
            }
          }
          if (sinfo.controllingAtoms[2] != Chirality::StereoInfo::NOATOM &&
              sinfo.controllingAtoms[3] != Chirality::StereoInfo::NOATOM) {
            if (aranks[sinfo.controllingAtoms[2]] ==
                aranks[sinfo.controllingAtoms[3]]) {
              haveADupe = true;
            } else if (aranks[sinfo.controllingAtoms[2]] <
                       aranks[sinfo.controllingAtoms[3]]) {
              std::swap(sinfo.controllingAtoms[2], sinfo.controllingAtoms[3]);
              needsSwap = !needsSwap;
            }
          }
          if (!haveADupe) {
            if (knownBonds[bidx]) {
              if (needsSwap) {
                sinfo.descriptor =
                    sinfo.descriptor == StereoDescriptor::Bond_Cis
                        ? StereoDescriptor::Bond_Trans
                        : StereoDescriptor::Bond_Cis;
              }
              auto gbs = getBondSymbol(bond);
              if (bondSymbols[bidx] != gbs) {
                bondSymbols[bidx] = gbs;
                needAnotherRound = true;
              }
            }
            res.push_back(std::move(sinfo));
          } else {
            needAnotherRound = true;
            bondSymbols[bidx] = getBondSymbol(bond);
            possibleBonds[bidx] = 0;
            bond->setProp(Chirality::_stereoNotPossible, 1, true);
            if (cleanIt &&
                sinfo.specified != Chirality::StereoSpecified::Unspecified) {
              bond->setStereo(Bond::BondStereo::STEREONONE);
            }
          }
        }
      } else {
        bond->setProp(Chirality::_stereoNotPossible, 1, true);
      }
    }
    if (!needAnotherRound) {
      break;
    }
  }

  // {
  //   boost::dynamic_bitset<> atomsInPlay(mol.getNumAtoms());
  //   atomsInPlay.set();
  //   boost::dynamic_bitset<> bondsInPlay(mol.getNumBonds());
  //   bondsInPlay.set();
  //   std::vector<unsigned int> aranks;
  //   const bool breakTies = false;
  //   const bool includeChirality = false;
  //   const bool includeIsotopes = false;
  //   Canon::rankFragmentAtoms(mol, aranks, atomsInPlay, bondsInPlay,
  //   nullptr,
  //                            nullptr, breakTies, includeChirality,
  //                            includeIsotopes);
  //   std::cerr << "===============" << std::endl;
  //   std::copy(aranks.begin(), aranks.end(),
  //             std::ostream_iterator<int>(std::cerr, " "));
  //   std::cerr << std::endl;
  // }

  mol.setProp("_potentialStereo", res, true);
  return res;
}

std::vector<StereoInfo> findPotentialStereo(const ROMol &mol) {
  ROMol cp(mol);
  bool cleanIt = false;
  return findPotentialStereo(cp, cleanIt);
}
#else
std::vector<StereoInfo> findPotentialStereo(ROMol &mol, bool cleanIt,
                                            bool flagPossible) {
  std::map<int, Atom::ChiralType> ochiralTypes;

  if (!mol.getRingInfo()->isInitialized()) {
    MolOps::symmetrizeSSSR(mol);
  }

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
        ochiralTypes[aidx] = atom->getChiralTag();
        atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
        atomSymbols[aidx] = (boost::format("%d%s_%d") % atom->getIsotope() %
                             atom->getSymbol() % aidx)
                                .str();
      } else {
        atomSymbols[aidx] = getAtomCompareSymbol(*atom);
      }
    } else {
      atomSymbols[aidx] = getAtomCompareSymbol(*atom);
      if (cleanIt) {
        atom->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
      }
    }
  }

  // flag possible ring stereo cases. The relevant cases here are:
  //    1) even-sized rings with possible (or specified) atoms opposite each
  //       other, like CC1CC(C)C1 or CC1CCC(C)CC1
  //    2) atoms sharing a bond which fuses two or more rings, like the central
  //       bond in C1CCC2CCCCC2C1

  // tracks the number of rings with possible ring stereo that the atom is in
  //  (only set for potential stereoatoms)
  std::vector<unsigned int> possibleRingStereoAtoms(mol.getNumAtoms());
  // tracks the number of rings with possible ring stereo that the bond is in
  //  (set for all bonds)
  std::vector<unsigned int> possibleRingStereoBonds(mol.getNumBonds());
  boost::dynamic_bitset<> possibleAtomsInRing(mol.getNumAtoms());
  for (unsigned int ridx = 0; ridx < mol.getRingInfo()->atomRings().size();
       ++ridx) {
    const auto &aring = mol.getRingInfo()->atomRings()[ridx];
    unsigned int nHere = 0;
    auto sz = aring.size();
    bool ring_is_odd_sized = sz % 2;
    auto half_sz = sz / 2 + ring_is_odd_sized;

    possibleAtomsInRing.reset();
    for (unsigned int ai = 0; ai < sz; ++ai) {
      auto aidx = aring[ai];
      if (!possibleAtoms[aidx] && !knownAtoms[aidx]) {
        continue;
      }

      if (!ring_is_odd_sized) {
        // find the index of the atom on the opposite side of the even-sized
        // ring
        auto oppositeidx = aring[(ai + half_sz) % sz];
        if (possibleAtoms[oppositeidx] || knownAtoms[oppositeidx]) {
          ++nHere;
          possibleAtomsInRing.set(aidx);
          continue;
        }
      }
      // if the atom is in more than one ring, explore the common edge to see if
      // we can find another potentially chiral atom
      if (mol.getRingInfo()->numAtomRings(aidx) > 1) {
        auto previous_otheridx = aidx;
        for (size_t step = 1; step <= half_sz; ++step) {
          auto otheridx = aring[(ai + step) % sz];
          auto bnd = mol.getBondBetweenAtoms(previous_otheridx, otheridx);
          if (mol.getRingInfo()->numBondRings(bnd->getIdx()) < 2) {
            // We reached the end of the common edge.
            break;
          }
          if (possibleAtoms[otheridx] || knownAtoms[otheridx]) {
            // We found another potentially chiral atom, no need to keep
            // searching.
            nHere += 2;
            possibleAtomsInRing.set(aidx);
            possibleAtomsInRing.set(otheridx);
            break;
          }
          previous_otheridx = otheridx;
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
      for (auto bidx : mol.getRingInfo()->bondRings()[ridx]) {
        ++possibleRingStereoBonds[bidx];
      }
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
    const bool includeChirality = false;
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
        if (!haveADupe) {
          res.push_back(std::move(sinfo));
        } else {
          removedStereo = true;
          atomSymbols[aidx] = getAtomCompareSymbol(*atom);
          possibleAtoms[aidx] = 0;
          if (cleanIt &&
              sinfo.specified != Chirality::StereoSpecified::Unspecified) {
            atom->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
          }
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
#endif
std::vector<StereoInfo> cleanExistingStereo(ROMol &mol, bool cleanIt) {
  std::map<int, Atom::ChiralType> ochiralTypes;

  if (!mol.getRingInfo()->isInitialized()) {
    MolOps::symmetrizeSSSR(mol);
  }
  if (mol.needsUpdatePropertyCache()) {
    mol.updatePropertyCache(false);
  }

  boost::dynamic_bitset<> knownAtoms(mol.getNumAtoms());
  std::vector<std::string> atomSymbols(mol.getNumAtoms());
  for (const auto atom : mol.atoms()) {
    auto aidx = atom->getIdx();
    atomSymbols[aidx] = getAtomCompareSymbol(*atom);
    if (detail::isAtomPotentialStereoAtom(atom)) {
      auto sinfo = detail::getStereoInfo(atom);
      switch (sinfo.specified) {
        case Chirality::StereoSpecified::Unknown:
        case Chirality::StereoSpecified::Specified:
          knownAtoms.set(aidx);
          atomSymbols[aidx] += "_@";
          ochiralTypes[aidx] = atom->getChiralTag();
          break;
        case Chirality::StereoSpecified::Unspecified:
          break;
        default:
          throw ValueErrorException("bad StereoInfo.specified type");
      }
    } else if (cleanIt) {
      atom->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
    }
  }

  // flag possible ring stereo cases. The relevant cases here are:
  //    1) even-sized rings with possible (or specified) atoms opposite each
  //       other, like CC1CC(C)C1 or CC1CCC(C)CC1
  //    2) atoms sharing a bond which fuses two or more rings, like the central
  //       bond in C1CCC2CCCCC2C1

  // tracks the number of rings with possible ring stereo that the atom is in
  //  (only set for potential stereoatoms)
  std::vector<unsigned int> possibleRingStereoAtoms(mol.getNumAtoms());
  // tracks the number of rings with possible ring stereo that the bond is in
  //  (set for all bonds)
  std::vector<unsigned int> possibleRingStereoBonds(mol.getNumBonds());
  boost::dynamic_bitset<> possibleAtomsInRing(mol.getNumAtoms());
  for (unsigned int ridx = 0; ridx < mol.getRingInfo()->atomRings().size();
       ++ridx) {
    const auto &aring = mol.getRingInfo()->atomRings()[ridx];
    unsigned int nHere = 0;
    auto sz = aring.size();
    possibleAtomsInRing.reset();
    for (unsigned int ai = 0; ai < aring.size(); ++ai) {
      auto aidx = aring[ai];
      if (!knownAtoms[aidx]) {
        continue;
      }
      if (!(aring.size() % 2) || mol.getRingInfo()->numAtomRings(aidx) > 1) {
        // find the index of the atom on the opposite side of the even-sized
        // ring or the equivalent position for an odd-sized ring if the atom is
        // in more than one ring
        auto oppositeidx = aring[(ai + sz / 2) % sz];
        if (knownAtoms[oppositeidx]) {
          nHere += 2;
          possibleAtomsInRing.set(aidx);
          possibleAtomsInRing.set(oppositeidx);
          continue;
        }
      }
      // if the atom is in more than one bond, see if there's
      // a possible neighbor on a fusion bond
      if (mol.getRingInfo()->numAtomRings(aidx) > 1) {
        auto otheridx = aring[(ai + 1) % aring.size()];
        if (knownAtoms[otheridx]) {
          auto bnd = mol.getBondBetweenAtoms(aidx, otheridx);
          CHECK_INVARIANT(bnd, "expected ring bond not found");
          if (mol.getRingInfo()->numBondRings(bnd->getIdx()) > 1) {
            nHere += 2;
            possibleAtomsInRing.set(aidx);
            possibleAtomsInRing.set(otheridx);
          }
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
      for (auto bidx : mol.getRingInfo()->bondRings()[ridx]) {
        ++possibleRingStereoBonds[bidx];
      }
    }
  }

  std::vector<std::string> bondSymbols(mol.getNumBonds());
  boost::dynamic_bitset<> knownBonds(mol.getNumBonds());
  for (const auto bond : mol.bonds()) {
    auto bidx = bond->getIdx();
    bondSymbols[bidx] = getBondSymbol(bond);
    if (detail::isBondPotentialStereoBond(bond)) {
      auto sinfo = detail::getStereoInfo(bond);
      switch (sinfo.specified) {
        case Chirality::StereoSpecified::Unknown:
        case Chirality::StereoSpecified::Specified:
          knownBonds.set(bidx);
          bondSymbols[bidx] += "_chi";
          break;
        case Chirality::StereoSpecified::Unspecified:
          break;
        default:
          throw ValueErrorException("bad StereoInfo.specified type");
      }
    } else if (cleanIt) {
      bond->setStereo(Bond::BondStereo::STEREONONE);
    }
  }
  boost::dynamic_bitset<> fixedAtoms(mol.getNumAtoms());
  boost::dynamic_bitset<> fixedBonds(mol.getNumBonds());
  std::vector<StereoInfo> res;
  bool needAnotherRound = true;
  while (needAnotherRound) {
    res.clear();
    needAnotherRound = false;

    // std::copy(atomSymbols.begin(), atomSymbols.end(),
    //           std::ostream_iterator<std::string>(std::cerr, " "));
    // std::cerr << std::endl;
    // std::copy(bondSymbols.begin(), bondSymbols.end(),
    //           std::ostream_iterator<std::string>(std::cerr, " "));
    // std::cerr << std::endl;
    // we will use the canonicalization code
    boost::dynamic_bitset<> atomsInPlay(mol.getNumAtoms());
    atomsInPlay.set();
    boost::dynamic_bitset<> bondsInPlay(mol.getNumBonds());
    bondsInPlay.set();
    std::vector<unsigned int> aranks;
    const bool breakTies = false;
    const bool includeChirality = false;
    const bool includeIsotopes = false;
    Canon::rankFragmentAtoms(mol, aranks, atomsInPlay, bondsInPlay,
                             &atomSymbols, &bondSymbols, breakTies,
                             includeChirality, includeIsotopes);

    for (const auto atom : mol.atoms()) {
      auto aidx = atom->getIdx();
      if (knownAtoms[aidx]) {
        auto sinfo = detail::getStereoInfo(atom);
        if (fixedAtoms[aidx]) {
          // FIX: should we be using fixedAtoms here?
          res.push_back(std::move(sinfo));
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

            auto acs = getAtomCompareSymbol(*atom);
            auto sortednbrs = nbrs;
            std::sort(sortednbrs.begin(), sortednbrs.end());
            // FIX: only works for tetrahedral at the moment
            if (sinfo.type == Chirality::StereoType::Atom_Tetrahedral) {
              auto nSwaps = countSwapsToInterconvert(nbrs, sortednbrs);
              if (nSwaps % 2) {
                sinfo.descriptor =
                    sinfo.descriptor == Chirality::StereoDescriptor::Tet_CCW
                        ? Chirality::StereoDescriptor::Tet_CW
                        : Chirality::StereoDescriptor::Tet_CCW;
              }
              if (sinfo.descriptor == Chirality::StereoDescriptor::Tet_CW) {
                acs += "_CW";
              } else if (sinfo.descriptor ==
                         Chirality::StereoDescriptor::Tet_CCW) {
                acs += "_CCW";
              }
            }
            if (atomSymbols[aidx] != acs) {
              atomSymbols[aidx] = acs;
              needAnotherRound = true;
              fixedAtoms.set(aidx);
            }
            res.push_back(std::move(sinfo));
          }
        }
      }
    }

    for (const auto bond : mol.bonds()) {
      auto bidx = bond->getIdx();
      if (knownBonds[bidx]) {
        auto sinfo = detail::getStereoInfo(bond);
        ASSERT_INVARIANT(sinfo.controllingAtoms.size() == 4,
                         "bad controlling atoms size");
        if (fixedBonds[bidx]) {
          res.push_back(std::move(sinfo));
        } else {
          bool haveADupe = false;
          bool needsSwap = false;
          if (sinfo.controllingAtoms[0] != Chirality::StereoInfo::NOATOM &&
              sinfo.controllingAtoms[1] != Chirality::StereoInfo::NOATOM) {
            if (aranks[sinfo.controllingAtoms[0]] ==
                aranks[sinfo.controllingAtoms[1]]) {
              haveADupe = true;
            } else if (aranks[sinfo.controllingAtoms[0]] <
                       aranks[sinfo.controllingAtoms[1]]) {
              std::swap(sinfo.controllingAtoms[0], sinfo.controllingAtoms[1]);
              needsSwap = !needsSwap;
            }
          }
          if (sinfo.controllingAtoms[2] != Chirality::StereoInfo::NOATOM &&
              sinfo.controllingAtoms[3] != Chirality::StereoInfo::NOATOM) {
            if (aranks[sinfo.controllingAtoms[2]] ==
                aranks[sinfo.controllingAtoms[3]]) {
              haveADupe = true;
            } else if (aranks[sinfo.controllingAtoms[2]] <
                       aranks[sinfo.controllingAtoms[3]]) {
              std::swap(sinfo.controllingAtoms[2], sinfo.controllingAtoms[3]);
              needsSwap = !needsSwap;
            }
          }
          if (!haveADupe) {
            if (needsSwap) {
              sinfo.descriptor = sinfo.descriptor == StereoDescriptor::Bond_Cis
                                     ? StereoDescriptor::Bond_Trans
                                     : StereoDescriptor::Bond_Cis;
            }
            auto gbs = getBondSymbol(bond);
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
            if (bondSymbols[bidx] != gbs) {
              bondSymbols[bidx] = gbs;
              needAnotherRound = true;
            }
            res.push_back(std::move(sinfo));
          }
        }
      }
    }
  }

  if (cleanIt) {
    for (auto i = 0u; i < mol.getNumAtoms(); ++i) {
      if (!fixedAtoms[i] && knownAtoms[i]) {
        switch (mol.getAtomWithIdx(i)->getChiralTag()) {
          case Atom::ChiralType::CHI_TETRAHEDRAL_CCW:
          case Atom::ChiralType::CHI_TETRAHEDRAL_CW:
            mol.getAtomWithIdx(i)->setChiralTag(
                Atom::ChiralType::CHI_UNSPECIFIED);
            break;
          case Atom::ChiralType::CHI_TETRAHEDRAL:
          case Atom::ChiralType::CHI_SQUAREPLANAR:
          case Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL:
          case Atom::ChiralType::CHI_OCTAHEDRAL:
            mol.getAtomWithIdx(i)->setProp(
                common_properties::_chiralPermutation, 0);
            break;
          default:
            break;
        }
      }
    }
    for (auto i = 0u; i < mol.getNumBonds(); ++i) {
      if (!fixedBonds[i] && knownBonds[i]) {
        // FIX only does known double bonds
        mol.getBondWithIdx(i)->setStereo(Bond::BondStereo::STEREONONE);
      }
    }
  }

  mol.setProp("_potentialStereo", res, true);
  return res;
}

}  // namespace Chirality
}  // namespace RDKit
