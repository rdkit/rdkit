//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Validate.h"
#include "Fragment.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/MolStandardize/FragmentCatalog/FragmentCatalogParams.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/PeriodicTable.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace RDKit {
class RWMol;
class ROMol;

namespace MolStandardize {

std::vector<ValidationErrorInfo> CompositeValidation::validate(
    const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;
  for (const auto &method : validations) {
    auto partial = method->validate(mol, reportAllFailures);
    if (!partial.empty()) {
      std::copy(partial.begin(), partial.end(), std::back_inserter(errors));
      if (!reportAllFailures) {
        break;
      }
    }
  }
  return errors;
}

std::vector<ValidationErrorInfo> RDKitValidation::validate(
    const ROMol &mol, bool reportAllFailures) const {
  ROMol molCopy = mol;
  std::vector<ValidationErrorInfo> errors;

  unsigned int na = mol.getNumAtoms();

  if (!na) {
    errors.emplace_back("ERROR: [NoAtomValidation] Molecule has no atoms");
  }

  // loop over atoms
  for (size_t i = 0; i < na; ++i) {
    if (!reportAllFailures) {
      if (errors.size() >= 1) {
        break;
      }
    }
    Atom *atom = molCopy.getAtomWithIdx(i);
    try {
      atom->calcExplicitValence();
    } catch (const MolSanitizeException &e) {
      errors.push_back("INFO: [ValenceValidation] " + std::string(e.what()));
    }
  }
  return errors;
}

std::vector<ValidationErrorInfo> NoAtomValidation::validate(
    const ROMol &mol, bool /*reportAllFailures*/) const {
  std::vector<ValidationErrorInfo> errors;
  unsigned int na = mol.getNumAtoms();
  if (!na) {
    errors.emplace_back("ERROR: [NoAtomValidation] Molecule has no atoms");
  }
  return errors;
}

std::vector<ValidationErrorInfo> FragmentValidation::validate(
    const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;
  // REVIEW: reportAllFailures is not being used here. is that correct?
  RDUNUSED_PARAM(reportAllFailures);
  std::shared_ptr<FragmentCatalogParams> fparams(new FragmentCatalogParams(""));
  FragmentCatalog fcat(fparams.get());

  const std::vector<std::shared_ptr<ROMol>> &fgrps = fparams->getFuncGroups();
  INT_VECT mapping;
  VECT_INT_VECT atom_mapping;
  std::vector<ROMOL_SPTR> frags =
      MolOps::getMolFrags(mol, true, &mapping, &atom_mapping);

  for (auto &fgrp : fgrps) {
    std::string fname;
    fgrp->getProp(common_properties::_Name, fname);
    std::vector<RDKit::MatchVectType> res;
    unsigned int matches = SubstructMatch(mol, *fgrp, res);
    //		std::cout << fname << " matches " << matches << std::endl;
    if (matches != 0 && frags.size() != 0) {
      VECT_INT_VECT substructmap;  // store idxs of frag from substructmatch
      for (const auto &match : res) {
        std::vector<int> vec;
        for (const auto &pair : match) {
          vec.push_back(pair.second);
        }
        substructmap.push_back(vec);
      }

      // to stop the same fragment being reported many times if present
      // multiple times in molecule
      bool fpresent = false;

      for (auto &molfragidx : atom_mapping) {
        std::sort(molfragidx.begin(), molfragidx.end());
        for (auto &substructidx : substructmap) {
          std::sort(substructidx.begin(), substructidx.end());
          //					// help to debug...
          //					std::cout << "molfragidx: "  <<
          // std::endl; 					for (const auto
          // &i : molfragidx)
          // {
          // std::cout << i; }; std::cout
          // << std::endl; std::cout <<
          //"substructidx: "  << std::endl;
          // for (const auto &i : substructidx) { std::cout << i; };
          // std::cout <<
          // std::endl;
          //					//
          if ((molfragidx == substructidx) && !fpresent) {
            std::string msg = fname + " is present";
            errors.push_back("INFO: [FragmentValidation] " + msg);
            fpresent = true;
          }
        }
      }
    }
  }
  return errors;
}

std::vector<ValidationErrorInfo> NeutralValidation::validate(
    const ROMol &mol, bool /*reportAllFailures*/) const {
  std::vector<ValidationErrorInfo> errors;
  int charge = RDKit::MolOps::getFormalCharge(mol);
  if (charge != 0) {
    std::string charge_str;
    if (charge > 0) {
      charge_str = "+" + std::to_string(charge);
    } else {
      charge_str = std::to_string(charge);
    }
    std::string msg = "Not an overall neutral system (" + charge_str + ')';
    errors.push_back("INFO: [NeutralValidation] " + msg);
  }
  return errors;
}

std::vector<ValidationErrorInfo> IsotopeValidation::validate(
    const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;
  for (auto atom : mol.atoms()) {
    unsigned int isotope = atom->getIsotope();
    if (isotope == 0) {
      continue;
    }

    std::string symbol = atom->getSymbol();
    unsigned int atomicNum = atom->getAtomicNum();
    if (atomicNum && strict) {
      PeriodicTable *periodicTable = PeriodicTable::getTable();
      double mass = periodicTable->getMassForIsotope(atomicNum, isotope);
      if (mass == 0.0) {
        errors.push_back(
            "ERROR: [IsotopeValidation] The molecule contains an unknown isotope: " +
            std::to_string(isotope) + symbol);
      }
    } else {
      errors.push_back("INFO: [IsotopeValidation] Molecule contains isotope " +
                       std::to_string(isotope) + symbol);
    }

    if (!errors.empty() && !reportAllFailures) {
      break;
    }
  }
  return errors;
}

// constructor
MolVSValidation::MolVSValidation()
    : CompositeValidation({std::make_shared<NoAtomValidation>(),
                           std::make_shared<FragmentValidation>(),
                           std::make_shared<NeutralValidation>(),
                           std::make_shared<IsotopeValidation>()}) {}

// overloaded constructor
MolVSValidation::MolVSValidation(
    const std::vector<std::shared_ptr<ValidationMethod>> &validations)
    : CompositeValidation(validations) {}

std::vector<ValidationErrorInfo> AllowedAtomsValidation::validate(
    const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;
  unsigned int na = mol.getNumAtoms();

  for (size_t i = 0; i < na; ++i) {
    if (!reportAllFailures) {
      if (errors.size() >= 1) {
        break;
      }
    }
    const Atom *qatom = mol.getAtomWithIdx(i);
    bool match = false;
    // checks to see qatom matches one of list of allowedAtoms
    for (const auto &allowedAtom : this->d_allowedList) {
      if (allowedAtom->Match(qatom)) {
        match = true;
        break;
      }
    }
    // if no match, append to list of errors.
    if (!match) {
      std::string symbol = qatom->getSymbol();
      errors.push_back("INFO: [AllowedAtomsValidation] Atom " + symbol +
                       " is not in allowedAtoms list");
    }
  }
  return errors;
}

std::vector<ValidationErrorInfo> DisallowedAtomsValidation::validate(
    const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;
  unsigned int na = mol.getNumAtoms();

  for (size_t i = 0; i < na; ++i) {
    if (!reportAllFailures) {
      if (errors.size() >= 1) {
        break;
      }
    }
    const Atom *qatom = mol.getAtomWithIdx(i);
    bool match = false;
    // checks to see qatom matches one of list of allowedAtoms
    for (const auto &disallowedAtom : this->d_disallowedList) {
      if (disallowedAtom->Match(qatom)) {
        match = true;
      }
    }
    // if no match, append to list of errors.
    if (match) {
      std::string symbol = qatom->getSymbol();
      errors.push_back("INFO: [DisallowedAtomsValidation] Atom " + symbol +
                       " is in disallowedAtoms list");
    }
  }
  return errors;
}

std::vector<ValidationErrorInfo> DisallowedRadicalValidation::validate(
    const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;

  for (auto atom : mol.atoms()) {
    unsigned int numRadicalElectrons = atom->getNumRadicalElectrons();
    if (numRadicalElectrons == 0) {
      continue;
    }
    unsigned int atomicNum = atom->getAtomicNum();
    unsigned int degree = atom->getDegree();
    if ((atomicNum == 7 || atomicNum == 8) && numRadicalElectrons == 1 &&
        degree == 1) {
      unsigned int neighborAtomicNum = 0;
      Bond::BondType bondType = Bond::BondType::UNSPECIFIED;
      for (auto neighbor : mol.atomNeighbors(atom)) {
        // only one iteration is performed, because degree == 1
        neighborAtomicNum = neighbor->getAtomicNum();
        bondType = mol.getBondBetweenAtoms(atom->getIdx(), neighbor->getIdx())
                       ->getBondType();
      }
      if (atomicNum == 7 && neighborAtomicNum == 8 &&
          bondType == Bond::BondType::DOUBLE) {
        // nitric oxide
        continue;
      }
      if (atomicNum == 8 && neighborAtomicNum == 7 &&
          bondType == Bond::BondType::SINGLE) {
        // aminoxyl
        continue;
      }
    }
    errors.push_back(
        "ERROR: [DisallowedRadicalValidation] The radical at atom " +
        std::to_string(atom->getIdx()) + " is not allowed");
    if (!reportAllFailures) {
      break;
    }
  }
  return errors;
}

std::vector<ValidationErrorInfo> FeaturesValidation::validate(
    const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;

  // Optionally disallow query and dummy atoms, and aliases
  for (auto atom : mol.atoms()) {
    if (!allowQueries && atom->hasQuery()) {
      errors.push_back("ERROR: [FeaturesValidation] Query atom " +
                       std::to_string(atom->getIdx()) + " is not allowed");
      if (!reportAllFailures) {
        return errors;
      }
    } else if (!allowDummies && isAtomDummy(atom)) {
      errors.push_back("ERROR: [FeaturesValidation] Dummy atom " +
                       std::to_string(atom->getIdx()) + " is not allowed");
      if (!reportAllFailures) {
        return errors;
      }
    }

    if (!allowAtomAliases && atom->hasProp(common_properties::molFileAlias)) {
      errors.push_back(
          "ERROR: [FeaturesValidation] Atom " + std::to_string(atom->getIdx()) +
          " with alias '" +
          atom->getProp<std::string>(common_properties::molFileAlias) +
          "' is not allowed");
      if (!reportAllFailures) {
        return errors;
      }
    }
  }

  // Optionally disallow query, aromatic or dative bonds
  for (auto bond : mol.bonds()) {
    if (!allowQueries && bond->hasQuery()) {
      errors.push_back("ERROR: [FeaturesValidation] Query bond " +
                       std::to_string(bond->getIdx()) + " is not allowed");
      if (!reportAllFailures) {
        return errors;
      }
    }
    if (!allowAromaticBondType &&
        bond->getBondType() == Bond::BondType::AROMATIC) {
      errors.push_back("ERROR: [FeaturesValidation] Bond " +
                       std::to_string(bond->getIdx()) +
                       " of aromatic type is not allowed");
      if (!reportAllFailures) {
        return errors;
      }
    }
    if (!allowDativeBondType && bond->getBondType() == Bond::BondType::DATIVE) {
      errors.push_back("ERROR: [FeaturesValidation] Bond " +
                       std::to_string(bond->getIdx()) +
                       " of dative type is not allowed");
      if (!reportAllFailures) {
        return errors;
      }
    }
  }

  // Optionally disallow using the enahanced stereochemistry
  if (!allowEnhancedStereo && mol.getStereoGroups().size()) {
    errors.emplace_back(
        "ERROR: [FeaturesValidation] Enhanced stereochemistry features are not allowed");
  }

  return errors;
}

std::vector<ValidationErrorInfo> Is2DValidation::validate(
    const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;

  if (!mol.getNumConformers()) {
    errors.emplace_back(
        "ERROR: [Is2DValidation] The molecule has no coordinates");
    return errors;
  }

  const auto &conf = mol.getConformer();

  if (conf.is3D()) {
    errors.emplace_back(
        "ERROR: [Is2DValidation] The molecule includes non-null Z coordinates");
    return errors;
  }

  // conf.is3D() is assigned by the mol format parser based on the input
  // mol block designation, but also taking into account the presence of
  // non-null Z coordinates or stereobonds.
  //
  // the following test is in this sense probably redundant, but it's still
  // implemented in case molecules are built by other means.

  double max_absz{};
  for (const auto &p : conf.getPositions()) {
    max_absz = std::max(std::abs(p.z), max_absz);
  }

  if (max_absz > threshold) {
    errors.emplace_back(
        "ERROR: [Is2DValidation] The molecule includes non-null Z coordinates");
    if (!reportAllFailures) {
      return errors;
    }
  }

  if (conf.getNumAtoms() < 2) {
    // there is nothing else to check here, if there is at most one atom.
    return errors;
  }

  // verify that the atoms are not all in the same position (this often happens
  // because no coordinates were assigned and all atoms appear to be placed in
  // the origin)

  double min_x = std::numeric_limits<double>::max();
  double max_x = std::numeric_limits<double>::min();
  double min_y = std::numeric_limits<double>::max();
  double max_y = std::numeric_limits<double>::min();
  for (const auto &p : conf.getPositions()) {
    min_x = std::min(p.x, min_x);
    max_x = std::max(p.x, max_x);
    min_y = std::min(p.y, min_y);
    max_y = std::max(p.y, max_y);
  }
  auto delta_x = max_x - min_x;
  auto delta_y = max_y - min_y;
  auto max_delta = std::max(delta_x, delta_y);

  if (max_delta < threshold) {
    errors.emplace_back(
        "ERROR: [Is2DValidation] All atoms have the same (x,y) coordinates");
    if (!reportAllFailures) {
      return errors;
    }
  }

  return errors;
}

double Layout2DValidation::squaredMedianBondLength(const ROMol &mol,
                                                   const Conformer &conf) {
  // Compute the squared value of the median bond length, but exclude the bonds
  // of null length.
  double median = 0.0;
  unsigned int numBonds = mol.getNumBonds();
  if (numBonds) {
    std::vector<double> values;
    values.reserve(numBonds);
    for (const auto &bond : mol.bonds()) {
      const auto &p1 = conf.getAtomPos(bond->getBeginAtomIdx());
      const auto &p2 = conf.getAtomPos(bond->getEndAtomIdx());
      auto value = (p1 - p2).lengthSq();
      if (value > 0.) {
        values.push_back(value);
      }
    }
    if (!values.empty()) {
      std::sort(values.begin(), values.end());
      numBonds = values.size();
      if (numBonds % 2) {
        median = values[numBonds / 2];
      } else {
        median = 0.5 * (values[numBonds / 2 - 1] + values[numBonds / 2]);
      }
    }
  }
  return median;
}

std::vector<ValidationErrorInfo> Layout2DValidation::validate(
    const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;

  if (!mol.getNumConformers()) {
    errors.emplace_back(
        "ERROR: [Layout2DValidation] The molecule has no coordinates");
    return errors;
  }

  const auto &conf = mol.getConformer();
  unsigned int natoms = conf.getNumAtoms();

  if (natoms < 2) {
    // there is nothing to check here, if there is only one atom.
    return errors;
  }

  // compute threshold values for the squared atom-atom or atom-bond
  // distance and for the maximum bond length using the median squared
  // bond length as reference.
  auto reference = squaredMedianBondLength(mol, conf);
  if (reference < minMedianBondLength * minMedianBondLength) {
    errors.emplace_back(
        "ERROR: [Layout2DValidation] The median bond length is smaller than the configured limit");
    if (!reportAllFailures) {
      return errors;
    }
  }

  // check for atoms clashing w/ other atoms
  auto atomClashThreshold = clashLimit * clashLimit * reference;
  for (unsigned int i = 0; i < natoms - 1; ++i) {
    const auto &pi = conf.getAtomPos(i);
    for (unsigned int j = i + 1; j < natoms; ++j) {
      const auto &pj = conf.getAtomPos(j);
      auto d2 = (pi - pj).lengthSq();
      if (d2 < atomClashThreshold) {
        errors.push_back("ERROR: [Layout2DValidation] Atom " +
                         std::to_string(i) + " is too close to atom " +
                         std::to_string(j));
        if (!reportAllFailures) {
          return errors;
        }
      }
    }
  }

  // make sure we have the required rings info available
  if (allowLongBondsInRings || allowAtomBondClashExemption) {
    if (!mol.getRingInfo()->isInitialized()) {
      RDKit::MolOps::fastFindRings(mol);
    }
  }

  for (auto bond : mol.bonds()) {
    unsigned int i = bond->getBeginAtomIdx();
    const auto &pi = conf.getAtomPos(i);
    unsigned int j = bond->getEndAtomIdx();
    const auto &pj = conf.getAtomPos(j);

    auto ll = (pi - pj).lengthSq();

    // check for exceedingly long bonds
    auto bondLengthThreshold = bondLengthLimit * bondLengthLimit * reference;
    if (!allowLongBondsInRings ||
        mol.getRingInfo()->numBondRings(bond->getIdx()) == 0) {
      if (ll > bondLengthThreshold) {
        errors.push_back("ERROR: [Layout2DValidation] The length of bond " +
                         std::to_string(bond->getIdx()) + " between atoms " +
                         std::to_string(i) + " and " + std::to_string(j) +
                         " exceeds a configured limit");
        if (!reportAllFailures) {
          return errors;
        }
      }
    }

    if (allowAtomBondClashExemption) {
      // is this bond exempted from atom-bond collision detection?
      if ((ll > 5. * 5. * reference) &&
          mol.getRingInfo()->numBondRings(bond->getIdx()) != 0) {
        continue;
      }
    }

    // check for atoms clashing with this bond
    for (unsigned int k = 0; k < natoms; ++k) {
      if (k == i || k == j) {
        continue;
      }
      const auto &pk = conf.getAtomPos(k);
      /*
               k
              /
            r/
            /
           /
          i---------------j
                  b
       */
      auto vik = pk - pi;
      auto vij = pj - pi;
      auto rr = vik.lengthSq();
      auto bb = vij.lengthSq();
      auto rb = vik.dotProduct(vij);
      static constexpr double EPS{
          1.e-7};  // prevent dividing by zero in extreme cases
      auto kb = (rr * bb - rb * rb) / (bb + EPS);
      if (rb >= 0. &&             /* cos alpha > 0 */
          rb <= bb &&             /* projection of r onto b does not exceed b */
          kb < atomClashThreshold /* distance from bond < limit */
      ) {
        errors.push_back("ERROR: [Layout2DValidation] Atom " +
                         std::to_string(k) + " too close to bond " +
                         std::to_string(bond->getIdx()));
        if (!reportAllFailures) {
          return errors;
        }
      }
    }
  }

  return errors;
}

namespace {
bool hasStereoBond(const ROMol &mol, const Atom *atom) {
  for (auto bond : mol.atomBonds(atom)) {
    if (atom != bond->getBeginAtom()) {
      continue;
    }
    auto bondDir = bond->getBondDir();
    if (bondDir == Bond::BondDir::BEGINDASH ||
        bondDir == Bond::BondDir::BEGINWEDGE ||
        bondDir == Bond::BondDir::UNKNOWN) {
      return true;
    }
  }
  return false;
}

struct BondInfo {
  const Bond *bond = nullptr;
  Bond::BondDir bondDir = Bond::BondDir::NONE;
  double angle = 0.;
};

struct BondDirCount {
  unsigned int wedge = 0;
  unsigned int dash = 0;
  unsigned int unknown = 0;
  unsigned int other = 0;
};

struct NeighborsInfo {
  NeighborsInfo(const ROMol &mol, const Atom *atom);
  std::vector<BondInfo> bonds;
  BondDirCount dirCount;
};

NeighborsInfo::NeighborsInfo(const ROMol &mol, const Atom *atom) {
  for (auto bond : mol.atomBonds(atom)) {
    BondInfo info;
    info.bond = bond;
    if (bond->getBeginAtom() == atom) {
      // do not consider the bond direction
      // settings of bonds that begin from
      // neighboring atoms
      info.bondDir = bond->getBondDir();
    }
    bonds.push_back(info);
  }

  for (const auto &info : bonds) {
    Bond::BondDir dir = info.bondDir;
    switch (dir) {
      case Bond::BondDir::BEGINDASH:
        ++dirCount.dash;
        break;
      case Bond::BondDir::BEGINWEDGE:
        ++dirCount.wedge;
        break;
      case Bond::BondDir::UNKNOWN:
        ++dirCount.unknown;
        break;
      case Bond::BondDir::NONE:
        // ok, bonds with unspecified direction
        // are fine to ignore
      case Bond::ENDUPRIGHT:
      case Bond::ENDDOWNRIGHT:
        // also ignore direction settings that
        // may describe the configuration of an
        // adjacent double bond
        break;
      default:
        ++dirCount.other;
    }
  }

  const auto &conf = mol.getConformer();
  const auto &p = conf.getAtomPos(atom->getIdx());
  const auto bond0 = bonds[0].bond;
  const auto atom0 = bond0->getOtherAtom(atom);
  const auto v0 = conf.getAtomPos(atom0->getIdx()) - p;

  // sort the neighbors based on the angle they form
  // with the first one
  auto degree = bonds.size();
  for (unsigned int n = 1; n < degree; ++n) {
    const auto bondn = bonds[n].bond;
    const auto atomn = bondn->getOtherAtom(atom);
    const auto vn = conf.getAtomPos(atomn->getIdx()) - p;
    bonds[n].angle = v0.signedAngleTo(vn);
  }

  std::sort(
      bonds.begin() + 1, bonds.end(),
      [](const BondInfo &a, const BondInfo &b) { return a.angle < b.angle; });
}

void check3CoordinatedStereo(const ROMol &mol, const Atom *atom,
                             const NeighborsInfo &neighborsInfo,
                             bool /*reportAllFailures*/,
                             std::vector<ValidationErrorInfo> &errors) {
  auto numStereoBonds =
      neighborsInfo.dirCount.dash + neighborsInfo.dirCount.wedge;

  if (numStereoBonds == 1) {
    // identify the stereo bond
    unsigned int i;
    for (i = 0; i < 3; ++i) {
      Bond::BondDir bondDir = neighborsInfo.bonds[i].bondDir;
      if (bondDir == Bond::BondDir::BEGINDASH ||
          bondDir == Bond::BondDir::BEGINWEDGE) {
        break;
      }
    }
    // check for the colinearity of the stereocenter and the other two ligands.
    const auto &conf = mol.getConformer();
    const auto &p = conf.getAtomPos(atom->getIdx());
    const auto atoma =
        neighborsInfo.bonds[(i + 1) % 3].bond->getOtherAtom(atom);
    const auto va = conf.getAtomPos(atoma->getIdx()) - p;
    const auto atomb =
        neighborsInfo.bonds[(i + 2) % 3].bond->getOtherAtom(atom);
    const auto vb = conf.getAtomPos(atomb->getIdx()) - p;

    auto angle = va.angleTo(vb);

    static constexpr auto ANGLE_EPSILON = (M_PI * 5. / 180.);  // 5 degrees
    if (angle < ANGLE_EPSILON || (M_PI - angle) < ANGLE_EPSILON) {
      errors.push_back(
          "ERROR: [StereoValidation] Colinearity of non-stereo bonds at atom " +
          std::to_string(atom->getIdx()));
    }
  } else {
    // configurations with multiple stereo bonds may be formally ambiguous or
    // unambiguos depending on their wedged/dashed direction and relative
    // orientation on the plane. those cases that are formally unambiguous are
    // still most often discouraged or also classified as not acceptable by
    // IUPAC guidelines due to lack of clarity.

    // The AvalonTools' struchk implementation simply doesn't allow multiple
    // stereo bonds on stereo centers with 3 explicit ligands. The validations
    // criteria for this sub-case could be in principle refined, but for now the
    // same policy is implemented.
    errors.push_back("ERROR: [StereoValidation] Atom " +
                     std::to_string(atom->getIdx()) +
                     " has 3 explicit substituents and multiple stereo bonds");
  }
}

void check4CoordinatedStereo(const ROMol &mol, const Atom *atom,
                             const NeighborsInfo &neighborsInfo,
                             bool reportAllFailures,
                             std::vector<ValidationErrorInfo> &errors) {
  if (neighborsInfo.dirCount.dash > 2 || neighborsInfo.dirCount.wedge > 2) {
    // this condition would anyway trigger an "adjacent bonds with like
    // orientation" alert, but this test could be clearer / more explicit.
    errors.push_back("ERROR: [StereoValidation] Atom " +
                     std::to_string(atom->getIdx()) +
                     " has too many stereo bonds with like orientation");
    if (!reportAllFailures) {
      return;
    }
  }

  for (unsigned int i = 0; i < 2; ++i) {
    if ((neighborsInfo.bonds[i].bondDir == Bond::BondDir::BEGINDASH &&
         neighborsInfo.bonds[i + 2].bondDir == Bond::BondDir::BEGINWEDGE) ||
        (neighborsInfo.bonds[i].bondDir == Bond::BondDir::BEGINWEDGE &&
         neighborsInfo.bonds[i + 2].bondDir == Bond::BondDir::BEGINDASH)) {
      errors.push_back(
          "ERROR: [StereoValidation] Atom " + std::to_string(atom->getIdx()) +
          " has opposing stereo bonds with different up/down orientation");
      if (!reportAllFailures) {
        return;
      }
    }
  }

  for (unsigned int i = 0; i < 4; ++i) {
    if ((neighborsInfo.bonds[i].bondDir == Bond::BondDir::BEGINDASH &&
         neighborsInfo.bonds[(i + 1) % 4].bondDir ==
             Bond::BondDir::BEGINDASH) ||
        (neighborsInfo.bonds[i].bondDir == Bond::BondDir::BEGINWEDGE &&
         neighborsInfo.bonds[(i + 1) % 4].bondDir ==
             Bond::BondDir::BEGINWEDGE)) {
      errors.push_back("ERROR: [StereoValidation] Atom " +
                       std::to_string(atom->getIdx()) +
                       " has adjacent stereo bonds with like orientation");
      if (!reportAllFailures) {
        return;
      }
      // it doesn't make sense to output this alert multiple times for the same
      // atom we therefore exit the loop also when reportAllFailures is not set.
      break;
    }
  }

  if (neighborsInfo.dirCount.dash + neighborsInfo.dirCount.wedge == 1) {
    // there is only one wedged/dashed bond. check for 'umbrellas' and
    // other geometric violations. we need the conformation here.
    const auto &conf = mol.getConformer();

    // identify the bond index for the stereo bond with specified direction.
    for (unsigned int i = 0; i < 4; ++i) {
      Bond::BondDir bondDir = neighborsInfo.bonds[i].bondDir;
      if (bondDir == Bond::BondDir::BEGINDASH ||
          bondDir == Bond::BondDir::BEGINWEDGE) {
        // count how many of the other bonds lie on the opposite half-plane,
        // i.e. form an angle > pi/4 with the stereo bond.
        unsigned int opposed = 0;
        const auto &p = conf.getAtomPos(atom->getIdx());
        const auto bondi = neighborsInfo.bonds[i].bond;
        const auto atomi = bondi->getOtherAtom(atom);
        const auto vi = conf.getAtomPos(atomi->getIdx()) - p;
        for (unsigned int j = 0; j < 4; ++j) {
          if (j == i) {
            continue;
          }
          const auto bondj = neighborsInfo.bonds[j].bond;
          const auto atomj = bondj->getOtherAtom(atom);
          const auto vj = conf.getAtomPos(atomj->getIdx()) - p;
          if (vi.angleTo(vj) > 95. * M_PI / 180.) {
            ++opposed;
          }
        }
        if (opposed == 3) {
          errors.push_back(
              "ERROR: [StereoValidation] Atom " +
              std::to_string(atom->getIdx()) +
              " has a potentially ambiguous representation: all non-stereo bonds" +
              " opposite to the only stereo bond");
        }
        if (!reportAllFailures) {
          return;
        }
        // there is only one stereo bond, which means we can exit the
        // outer loop on the first execution of this block.
        break;
      }
    }

    // check for collinearity violations and/or cases where the
    // the middle non-stereo bond is badly positioned (i.e., too short
    // compared to the other two on its sides).
    for (unsigned int i = 0; i < 4; i++) {
      Bond::BondDir bondDir = neighborsInfo.bonds[i].bondDir;
      if (bondDir == Bond::BondDir::BEGINDASH ||
          bondDir == Bond::BondDir::BEGINWEDGE) {
        auto j = (i + 1) % 4;
        auto k = (i + 2) % 4;
        auto l = (i + 3) % 4;
        const auto atomj = neighborsInfo.bonds[j].bond->getOtherAtom(atom);
        const auto atomk = neighborsInfo.bonds[k].bond->getOtherAtom(atom);
        const auto atoml = neighborsInfo.bonds[l].bond->getOtherAtom(atom);
        const auto &pj = conf.getAtomPos(atomj->getIdx());
        const auto &pk = conf.getAtomPos(atomk->getIdx());
        const auto &pl = conf.getAtomPos(atoml->getIdx());
        const auto v1 = pj - pk;
        const auto v2 = pl - pk;
        auto angle = v1.signedAngleTo(v2);
        if (angle < 185. * M_PI / 180.) {
          errors.push_back(
              "ERROR: [StereoValidation] Colinearity or triangle rule violation of "
              "non-stereo bonds at atom " +
              std::to_string(atom->getIdx()) /* +
" due to angle formed by ("  +
std::to_string(atomj->getIdx()+1) + "," +
std::to_string(atomk->getIdx()+1) + "," +
std::to_string(atoml->getIdx()+1) + ")" */
          );
          if (!reportAllFailures) {
            return;
          }
        }
        // there is only one stereo bond, which means we can exit the
        // outer loop on the first execution of this block.
        break;
      }
    }
  }
}

void checkStereo(const ROMol &mol, const Atom *atom, bool reportAllFailures,
                 std::vector<ValidationErrorInfo> &errors) {
  NeighborsInfo neighborsInfo(mol, atom);

  if (neighborsInfo.dirCount.other) {
    errors.push_back(
        "ERROR: [StereoValidation] one or more bonds incident to atom " +
        std::to_string(atom->getIdx()) + " have unexpected direction settings");
    // this is an unlikely condition and it would make little sense to
    // continue the analysis also when reportAllFailures were set.
    return;
  }

  if (neighborsInfo.dirCount.unknown) {
    if (neighborsInfo.dirCount.dash || neighborsInfo.dirCount.wedge) {
      errors.push_back("ERROR: [StereoValidation] Atom " +
                       std::to_string(atom->getIdx()) +
                       " has both unknown and wedged/dashed stereo bonds.");
    }
    // else: if the only stereo bonds have either/unknown direction,
    // we can return here.
    return;
  }

  for (const auto &bondInfo : neighborsInfo.bonds) {
    bool isStereo = bondInfo.bondDir == Bond::BondDir::BEGINDASH ||
                    bondInfo.bondDir == Bond::BondDir::BEGINWEDGE ||
                    bondInfo.bondDir == Bond::BondDir::UNKNOWN;
    if (isStereo && !canHaveDirection(*bondInfo.bond)) {
      errors.push_back("ERROR: [StereoValidation] Bond " +
                       std::to_string(bondInfo.bond->getIdx()) +
                       " has assigned stereo type, but unexpected bond order.");
      if (!reportAllFailures) {
        return;
      }
    }
  }

  // The validation is currently limited to some specific categories of
  // stereocenters
  bool multipleBondFound{}, possibleAllene{};
  for (auto bond : mol.atomBonds(atom)) {
    auto bondType = bond->getBondType();
    if (bondType != Bond::BondType::SINGLE) {
      multipleBondFound = true;
      const Atom *otherAtom = bond->getOtherAtom(atom);
      if (otherAtom->getDegree() == 2) {
        int doubleBondCount{};
        for (auto otherBond : mol.atomBonds(otherAtom)) {
          if (otherBond->getBondType() == Bond::BondType::DOUBLE) {
            ++doubleBondCount;
          }
        }
        if (doubleBondCount == 2) {
          possibleAllene = true;
        }
      }
    }
  }
  auto atomicNum = atom->getAtomicNum();
  if (possibleAllene || (multipleBondFound && atomicNum == 15)) {
    // Allenes and P compounds are not validated at this time.
    return;
  }
  if (multipleBondFound && atomicNum != 16) {
    // A stereo bond was found at an unsaturated atom. This condition used to
    // trigger as error in STRUCHK, but there are valid use cases for it (e.g.,
    // wavy bonds incident to double bonds of undefined/unknown configuration,
    // and atropisomers).
    //
    // Validation of these use cases is not currently implemented.
    return;
  }

  switch (atom->getDegree()) {
    case 1:
    case 2:
      errors.push_back(
          "ERROR: [StereoValidation] Atom " + std::to_string(atom->getIdx()) +
          " has stereo bonds, but less than 3 explicit substituents.");
      break;
    case 3:
      check3CoordinatedStereo(mol, atom, neighborsInfo, reportAllFailures,
                              errors);
      break;
    case 4:
      check4CoordinatedStereo(mol, atom, neighborsInfo, reportAllFailures,
                              errors);
      break;
    default:;
  }
}
}  // namespace

std::vector<ValidationErrorInfo> StereoValidation::validate(
    const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;

  for (auto atom : mol.atoms()) {
    if (hasStereoBond(mol, atom)) {
      checkStereo(mol, atom, reportAllFailures, errors);
    }
    if (!errors.empty() && !reportAllFailures) {
      break;
    }
  }

  return errors;
}

std::vector<ValidationErrorInfo> validateSmiles(const std::string &smiles) {
  RWMOL_SPTR mol(SmilesToMol(smiles));
  if (!mol) {
    std::string message =
        "SMILES Parse Error: syntax error for input: " + smiles;
    throw ValueErrorException(message);
  }

  MolVSValidation vm;
  std::vector<ValidationErrorInfo> errors = vm.validate(*mol, true);

  return errors;
}

}  // namespace MolStandardize
}  // namespace RDKit
