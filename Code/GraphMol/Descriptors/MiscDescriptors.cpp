//
//  Copyright (C) 2026 Glysade LLC and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/GraphMol.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include "MolDescriptors.h"
#include "MiscDescriptors.h"

#include <cmath>
#include <memory>
#include <utility>
#include <vector>

namespace RDKit {
namespace Descriptors {

double calcHeavyAtomMolWt(const ROMol &mol) { return calcAMW(mol, true); }

int calcNumValenceElectrons(const ROMol &mol) {
  const auto *tbl = PeriodicTable::getTable();
  int res = 0;
  for (const auto atom : mol.atoms()) {
    res += tbl->getNouterElecs(atom->getAtomicNum()) - atom->getFormalCharge() +
           static_cast<int>(atom->getTotalNumHs());
  }
  return res;
}

unsigned int calcNumRadicalElectrons(const ROMol &mol) {
  unsigned int res = 0;
  for (const auto atom : mol.atoms()) {
    res += atom->getNumRadicalElectrons();
  }
  return res;
}

namespace {
// returns (minCharge, maxCharge). The comparisons mirror python's builtin
// min()/max() so that NaN charges are handled the same way.
std::pair<double, double> chargeDescriptors(const ROMol &mol) {
  std::vector<double> chgs(mol.getNumAtoms(), 0.0);
  computeGasteigerCharges(mol, chgs);
  double minChg = 500.;
  double maxChg = -500.;
  for (auto chg : chgs) {
    minChg = (minChg < chg) ? minChg : chg;
    maxChg = (maxChg > chg) ? maxChg : chg;
  }
  return std::make_pair(minChg, maxChg);
}
}  // namespace

double calcMaxPartialCharge(const ROMol &mol) {
  return chargeDescriptors(mol).second;
}

double calcMinPartialCharge(const ROMol &mol) {
  return chargeDescriptors(mol).first;
}

double calcMaxAbsPartialCharge(const ROMol &mol) {
  auto chgs = chargeDescriptors(mol);
  double v1 = std::fabs(chgs.first);
  double v2 = std::fabs(chgs.second);
  return (v2 > v1) ? v2 : v1;
}

double calcMinAbsPartialCharge(const ROMol &mol) {
  auto chgs = chargeDescriptors(mol);
  double v1 = std::fabs(chgs.first);
  double v2 = std::fabs(chgs.second);
  return (v2 < v1) ? v2 : v1;
}

double calcFpDensityMorgan(const ROMol &mol, unsigned int radius) {
  auto numHeavyAtoms = mol.getNumHeavyAtoms();
  if (!numHeavyAtoms) {
    return 0.0;
  }
  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpg(
      MorganFingerprint::getMorganGenerator<std::uint64_t>(radius));
  std::unique_ptr<SparseIntVect<std::uint64_t>> fp(
      fpg->getSparseCountFingerprint(mol));
  return static_cast<double>(fp->getNonzeroElements().size()) / numHeavyAtoms;
}

}  // namespace Descriptors
}  // namespace RDKit
