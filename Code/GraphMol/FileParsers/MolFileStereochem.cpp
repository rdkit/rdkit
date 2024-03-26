//
//  Copyright (C) 2004-2021 Greg Landrum and other RDKit contributors
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
#include <Geometry/point.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <algorithm>
#include <RDGeneral/Ranking.h>
#include <RDGeneral/FileParseException.h>

namespace RDKit {

void WedgeBond(Bond *bond, unsigned int fromAtomIdx, const Conformer *conf) {
  Chirality::wedgeBond(bond, fromAtomIdx, conf);
}

void WedgeMolBonds(ROMol &mol, const Conformer *conf) {
  return Chirality::wedgeMolBonds(mol, conf);
}

std::vector<Bond *> getBondNeighbors(ROMol &mol, const Bond &bond) {
  std::vector<Bond *> res;
  for (auto nbri :
       boost::make_iterator_range(mol.getAtomBonds(bond.getBeginAtom()))) {
    auto nbrBond = mol[nbri];
    if (nbrBond == &bond) {
      continue;
    }
    res.push_back(nbrBond);
  }
  for (auto nbri :
       boost::make_iterator_range(mol.getAtomBonds(bond.getEndAtom()))) {
    auto nbrBond = mol[nbri];
    if (nbrBond == &bond) {
      continue;
    }
    res.push_back(nbrBond);
  }
  return res;
}

const Atom *getNonsharedAtom(const Bond &bond1, const Bond &bond2) {
  if (bond1.getBeginAtomIdx() == bond2.getBeginAtomIdx() ||
      bond1.getBeginAtomIdx() == bond2.getEndAtomIdx()) {
    return bond1.getEndAtom();
  } else if (bond1.getEndAtomIdx() == bond2.getBeginAtomIdx() ||
             bond1.getEndAtomIdx() == bond2.getEndAtomIdx()) {
    return bond1.getBeginAtom();
  }
  POSTCONDITION(0, "bonds don't share an atom");
}

const unsigned StereoBondThresholds::DBL_BOND_NO_STEREO;
const unsigned StereoBondThresholds::DBL_BOND_SPECIFIED_STEREO;
const unsigned StereoBondThresholds::CHIRAL_ATOM;
const unsigned StereoBondThresholds::DIRECTION_SET;

// a note on the way the StereoBondThresholds are used:
//  the penalties are all 1/10th of the corresponding threshold, so
//  the penalty for being connected to a chiral atom is
//  StereoBondThresholds::CHIRAL_ATOM/10
//  This allows us to just add up the penalties for a particular
//  single bond and still use one set of thresholds - an individual
//  single bond will never have any particular penalty term applied
//  more than a couple of times
//
void addWavyBondsForStereoAny(ROMol &mol, bool clearDoubleBondFlags,
                              unsigned addWhenImpossible) {
  std::vector<int> singleBondScores(mol.getNumBonds(), 0);
  // used to store the double bond neighbors, if any, of each single bond
  std::map<unsigned, std::vector<unsigned>> singleBondNeighbors;
  boost::dynamic_bitset<> doubleBondsToSet(mol.getNumBonds());
  // mark single bonds adjacent to double bonds
  for (const auto dblBond : mol.bonds()) {
    if (dblBond->getBondType() != Bond::BondType::DOUBLE) {
      continue;
    }
    if (dblBond->getStereo() == Bond::BondStereo::STEREOANY) {
      doubleBondsToSet.set(dblBond->getIdx());
    }
    for (auto singleBond : getBondNeighbors(mol, *dblBond)) {
      if (singleBond->getBondType() != Bond::BondType::SINGLE) {
        continue;
      }
      // NOTE: we could make this canonical by initializing scores to the
      // canonical atom ranks
      int score = singleBondScores[singleBond->getIdx()];
      ++score;

      // penalty for having a direction already set
      if (singleBond->getBondDir() != Bond::BondDir::NONE) {
        score += StereoBondThresholds::DIRECTION_SET / 10;
      }

      // penalties from the double bond itself:

      // penalize being adjacent to a double bond with empty stereo:
      if (dblBond->getStereo() == Bond::BondStereo::STEREONONE) {
        score += StereoBondThresholds::DBL_BOND_NO_STEREO / 10;
      } else if (dblBond->getStereo() > Bond::BondStereo::STEREOANY) {
        // penalize being adjacent to a double bond with specified stereo:
        score += StereoBondThresholds::DBL_BOND_SPECIFIED_STEREO / 10;
      }

      // atom-related penalties
      auto otherAtom = getNonsharedAtom(*singleBond, *dblBond);
      // favor atoms with smaller numbers of neighbors:
      score += 10 * otherAtom->getDegree();
      // penalty for being adjacent to an atom with specified stereo
      if (otherAtom->getChiralTag() != Atom::ChiralType::CHI_UNSPECIFIED &&
          otherAtom->getChiralTag() != Atom::ChiralType::CHI_OTHER) {
        score += StereoBondThresholds::CHIRAL_ATOM / 10;
      }
      singleBondScores[singleBond->getIdx()] = score;
      if (dblBond->getStereo() == Bond::BondStereo::STEREOANY) {
        singleBondNeighbors[singleBond->getIdx()].push_back(dblBond->getIdx());
      }
    }
  }
  std::vector<std::tuple<int, unsigned int, size_t>> sortedScores;
  for (size_t i = 0; i < mol.getNumBonds(); ++i) {
    auto score = singleBondScores[i];
    if (!score) {
      continue;
    }
    sortedScores.push_back(
        std::make_tuple(-1 * singleBondNeighbors[i].size(), score, i));
  }
  std::sort(sortedScores.begin(), sortedScores.end());
  for (const auto &tpl : sortedScores) {
    // FIX: check if dir is already set
    for (auto dblBondIdx : singleBondNeighbors[std::get<2>(tpl)]) {
      if (doubleBondsToSet[dblBondIdx]) {
        if (addWhenImpossible) {
          if (std::get<1>(tpl) > addWhenImpossible) {
            continue;
          }
        } else if (std::get<1>(tpl) >
                   StereoBondThresholds::DBL_BOND_NO_STEREO) {
          BOOST_LOG(rdWarningLog)
              << "Setting wavy bond flag on bond " << std::get<2>(tpl)
              << " which may make other stereo info ambiguous" << std::endl;
        }
        mol.getBondWithIdx(std::get<2>(tpl))
            ->setBondDir(Bond::BondDir::UNKNOWN);
        if (clearDoubleBondFlags) {
          auto dblBond = mol.getBondWithIdx(dblBondIdx);
          if (dblBond->getBondDir() == Bond::BondDir::EITHERDOUBLE) {
            dblBond->setBondDir(Bond::BondDir::NONE);
          }
          dblBond->setStereo(Bond::BondStereo::STEREONONE);
        }
        doubleBondsToSet.reset(dblBondIdx);
      }
    }
  }
  if (addWhenImpossible) {
    if (doubleBondsToSet.count()) {
      std::stringstream sstr;
      sstr << " unable to set wavy bonds for double bonds:";
      for (size_t i = 0; i < mol.getNumBonds(); ++i) {
        if (doubleBondsToSet[i]) {
          sstr << " " << i;
        }
      }
      BOOST_LOG(rdWarningLog) << sstr.str() << std::endl;
    }
  }
}
//
// Determine bond wedge state
///
Bond::BondDir DetermineBondWedgeState(const Bond *bond,
                                      unsigned int fromAtomIdx,
                                      const Conformer *conf) {
  return Chirality::detail::determineBondWedgeState(bond, fromAtomIdx, conf);
}
Bond::BondDir DetermineBondWedgeState(
    const Bond *bond,
    const std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
        &wedgeBonds,
    const Conformer *conf) {
  return Chirality::detail::determineBondWedgeState(bond, wedgeBonds, conf);
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
  MolOps::clearSingleBondDirFlags(mol);
}

void DetectBondStereoChemistry(ROMol &mol, const Conformer *conf) {
  PRECONDITION(conf, "no conformer");
  PRECONDITION(&(conf->getOwningMol()) == &mol,
               "conformer does not belong to molecule");
  MolOps::detectBondStereochemistry(mol, conf->getId());
}
void reapplyMolBlockWedging(RWMol &mol) {
  RDKit::Chirality::reapplyMolBlockWedging(mol);
}

void clearMolBlockWedgingInfo(RWMol &mol) {
  Chirality::clearMolBlockWedgingInfo(mol);
}

void invertMolBlockWedgingInfo(ROMol &mol) {
  Chirality::invertMolBlockWedgingInfo(mol);
}

void markUnspecifiedStereoAsUnknown(ROMol &mol, int confId) {
  const auto conf = mol.getConformer(confId);
  auto wedgeBonds = RDKit::Chirality::pickBondsToWedge(mol, nullptr, &conf);
  for (auto b : mol.bonds()) {
    if (b->getBondType() == Bond::DOUBLE) {
      int dirCode;
      bool reverse;
      RDKit::Chirality::GetMolFileBondStereoInfo(b, wedgeBonds, &conf, dirCode,
                                                 reverse);
      if (dirCode == 3) {
        b->setStereo(Bond::STEREOANY);
      }
    }
  }
  static int noNbrs = 100;
  auto si = Chirality::findPotentialStereo(mol);
  if (si.size()) {
    std::pair<bool, INT_VECT> retVal =
        Chirality::detail::countChiralNbrs(mol, noNbrs);
    INT_VECT nChiralNbrs = retVal.second;
    for (auto i : si) {
      if (i.type == Chirality::StereoType::Atom_Tetrahedral &&
          i.specified == Chirality::StereoSpecified::Unspecified) {
        i.specified = Chirality::StereoSpecified::Unknown;
        auto atom = mol.getAtomWithIdx(i.centeredOn);
        std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
            resSoFar;
        int bndIdx = Chirality::detail::pickBondToWedge(atom, mol, nChiralNbrs,
                                                        resSoFar, noNbrs);
        auto bond = mol.getBondWithIdx(bndIdx);
        bond->setBondDir(Bond::UNKNOWN);
      }
    }
  }
}

void translateChiralFlagToStereoGroups(ROMol &mol,
                                       StereoGroupType zeroFlagGroupType) {
  if (!mol.hasProp(common_properties::_MolFileChiralFlag)) {
    return;
  }
  int flagVal = 0;
  mol.getProp(common_properties::_MolFileChiralFlag, flagVal);
  mol.clearProp(common_properties::_MolFileChiralFlag);

  StereoGroupType sgType =
      flagVal ? StereoGroupType::STEREO_ABSOLUTE : zeroFlagGroupType;

  auto sgs = mol.getStereoGroups();

  boost::dynamic_bitset<> sgAtoms(mol.getNumAtoms());
  boost::dynamic_bitset<> sgBonds(mol.getNumBonds());
  const StereoGroup *absGroup = nullptr;
  for (const auto &sg : sgs) {
    for (const auto aptr : sg.getAtoms()) {
      sgAtoms.set(aptr->getIdx());
    }
    for (const auto bptr : sg.getBonds()) {
      sgBonds.set(bptr->getIdx());
    }
    // if we already have an ABS group, we'll add to it
    if (sgType == StereoGroupType::STEREO_ABSOLUTE && !absGroup &&
        sg.getGroupType() == StereoGroupType::STEREO_ABSOLUTE) {
      absGroup = &sg;
    }
  }
  ROMol::ATOM_PTR_VECT stereoAts;
  ROMol::BOND_PTR_VECT stereoBds;
  for (const auto atom : mol.atoms()) {
    if (!sgAtoms[atom->getIdx()] &&
        (atom->getChiralTag() == Atom::ChiralType::CHI_TETRAHEDRAL_CCW ||
         atom->getChiralTag() == Atom::ChiralType::CHI_TETRAHEDRAL_CW)) {
      stereoAts.push_back(atom);
    }
  }
  for (const auto bond : mol.bonds()) {
    if (!sgBonds[bond->getIdx()] &&
        (bond->getStereo() == Bond::BondStereo::STEREOATROPCCW ||
         bond->getStereo() == Bond::BondStereo::STEREOATROPCW)) {
      stereoBds.push_back(bond);
    }
  }
  if (!stereoAts.empty() || !stereoBds.empty()) {
    if (!absGroup) {
      sgs.emplace_back(sgType, stereoAts, stereoBds, 0);
      mol.setStereoGroups(sgs);
    } else {
      std::vector<StereoGroup> newSgs;
      for (const auto &sg : sgs) {
        if (&sg != absGroup) {
          newSgs.push_back(sg);
        } else {
          ROMol::ATOM_PTR_VECT newStereoAtoms = sg.getAtoms();
          newStereoAtoms.insert(newStereoAtoms.end(), stereoAts.begin(),
                                stereoAts.end());
          ROMol::BOND_PTR_VECT newStereoBonds = sg.getBonds();
          newStereoBonds.insert(newStereoBonds.end(), stereoBds.begin(),
                                stereoBds.end());

          newSgs.emplace_back(StereoGroupType::STEREO_ABSOLUTE, newStereoAtoms,
                              newStereoBonds, 0);
        }
      }
      mol.setStereoGroups(newSgs);
    }
  }
}
}  // namespace RDKit
