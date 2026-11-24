//
//  Copyright (C) 2026 Brian Kelley
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "NormalizedDescriptors.h"

#include <GraphMol/ROMol.h>
#include <GraphMol/Descriptors/FragmentDescriptors.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace RDKit {
namespace NormalizedDescriptors {

// defined in NormalizedDescriptors.cpp
std::shared_ptr<const HistogramSet> defaultTablesPtr();

namespace {
// The non-fragment descriptors in rdkit.Chem.Descriptors._descList, in that
// order
const std::vector<std::string> descriptorNames = {
    "MaxAbsEStateIndex",
    "MaxEStateIndex",
    "MinAbsEStateIndex",
    "MinEStateIndex",
    "qed",
    "SPS",
    "MolWt",
    "HeavyAtomMolWt",
    "ExactMolWt",
    "NumValenceElectrons",
    "NumRadicalElectrons",
    "MaxPartialCharge",
    "MinPartialCharge",
    "MaxAbsPartialCharge",
    "MinAbsPartialCharge",
    "FpDensityMorgan1",
    "FpDensityMorgan2",
    "FpDensityMorgan3",
    "BCUT2D_MWHI",
    "BCUT2D_MWLOW",
    "BCUT2D_CHGHI",
    "BCUT2D_CHGLO",
    "BCUT2D_LOGPHI",
    "BCUT2D_LOGPLOW",
    "BCUT2D_MRHI",
    "BCUT2D_MRLOW",
    "AvgIpc",
    "BalabanJ",
    "BertzCT",
    "Chi0",
    "Chi0n",
    "Chi0v",
    "Chi1",
    "Chi1n",
    "Chi1v",
    "Chi2n",
    "Chi2v",
    "Chi3n",
    "Chi3v",
    "Chi4n",
    "Chi4v",
    "HallKierAlpha",
    "Ipc",
    "Kappa1",
    "Kappa2",
    "Kappa3",
    "LabuteASA",
    "PEOE_VSA1",
    "PEOE_VSA10",
    "PEOE_VSA11",
    "PEOE_VSA12",
    "PEOE_VSA13",
    "PEOE_VSA14",
    "PEOE_VSA2",
    "PEOE_VSA3",
    "PEOE_VSA4",
    "PEOE_VSA5",
    "PEOE_VSA6",
    "PEOE_VSA7",
    "PEOE_VSA8",
    "PEOE_VSA9",
    "SMR_VSA1",
    "SMR_VSA10",
    "SMR_VSA2",
    "SMR_VSA3",
    "SMR_VSA4",
    "SMR_VSA5",
    "SMR_VSA6",
    "SMR_VSA7",
    "SMR_VSA8",
    "SMR_VSA9",
    "SlogP_VSA1",
    "SlogP_VSA10",
    "SlogP_VSA11",
    "SlogP_VSA12",
    "SlogP_VSA2",
    "SlogP_VSA3",
    "SlogP_VSA4",
    "SlogP_VSA5",
    "SlogP_VSA6",
    "SlogP_VSA7",
    "SlogP_VSA8",
    "SlogP_VSA9",
    "TPSA",
    "EState_VSA1",
    "EState_VSA10",
    "EState_VSA11",
    "EState_VSA2",
    "EState_VSA3",
    "EState_VSA4",
    "EState_VSA5",
    "EState_VSA6",
    "EState_VSA7",
    "EState_VSA8",
    "EState_VSA9",
    "VSA_EState1",
    "VSA_EState10",
    "VSA_EState2",
    "VSA_EState3",
    "VSA_EState4",
    "VSA_EState5",
    "VSA_EState6",
    "VSA_EState7",
    "VSA_EState8",
    "VSA_EState9",
    "FractionCSP3",
    "HeavyAtomCount",
    "NHOHCount",
    "NOCount",
    "NumAliphaticCarbocycles",
    "NumAliphaticHeterocycles",
    "NumAliphaticRings",
    "NumAmideBonds",
    "NumAromaticCarbocycles",
    "NumAromaticHeterocycles",
    "NumAromaticRings",
    "NumAtomStereoCenters",
    "NumBridgeheadAtoms",
    "NumHAcceptors",
    "NumHDonors",
    "NumHeteroatoms",
    "NumHeterocycles",
    "NumRotatableBonds",
    "NumSaturatedCarbocycles",
    "NumSaturatedHeterocycles",
    "NumSaturatedRings",
    "NumSpiroAtoms",
    "NumUnspecifiedAtomStereoCenters",
    "Phi",
    "RingCount",
    "MolLogP",
    "MolMR",
};

const std::vector<std::string> &allNames() {
  static const std::vector<std::string> names = [] {
    auto res = descriptorNames;
    // _descList sorts the fragment descriptors by name
    auto frags = Descriptors::getFragmentDescriptorNames();
    std::sort(frags.begin(), frags.end());
    res.insert(res.end(), frags.begin(), frags.end());
    return res;
  }();
  return names;
}

// _descList names of descriptors the property registry has under other names
const std::unordered_map<std::string, std::string> propertyAliases = {
    {"MolWt", "amw"},
    {"ExactMolWt", "exactmw"},
    {"Chi0n", "chi0n"},
    {"Chi0v", "chi0v"},
    {"Chi1n", "chi1n"},
    {"Chi1v", "chi1v"},
    {"Chi2n", "chi2n"},
    {"Chi2v", "chi2v"},
    {"Chi3n", "chi3n"},
    {"Chi3v", "chi3v"},
    {"Chi4n", "chi4n"},
    {"Chi4v", "chi4v"},
    {"HallKierAlpha", "hallKierAlpha"},
    {"Kappa1", "kappa1"},
    {"Kappa2", "kappa2"},
    {"Kappa3", "kappa3"},
    {"LabuteASA", "labuteASA"},
    {"TPSA", "tpsa"},
    {"HeavyAtomCount", "NumHeavyAtoms"},
    {"NHOHCount", "lipinskiHBD"},
    {"NOCount", "lipinskiHBA"},
    {"NumHAcceptors", "NumHBA"},
    {"NumHDonors", "NumHBD"},
    {"RingCount", "NumRings"},
    {"MolLogP", "CrippenClogP"},
    {"MolMR", "CrippenMR"},
};

std::vector<std::string> toPropertyNames(
    const std::vector<std::string> &names) {
  std::vector<std::string> res;
  res.reserve(names.size());
  for (const auto &name : names) {
    res.push_back(getPropertyName(name));
  }
  return res;
}

const NormalizedProperties &defaultProperties() {
  static const NormalizedProperties props;
  return props;
}
}  // namespace

const std::vector<std::string> &getNormalizedDescriptorNames() {
  static const std::vector<std::string> names = [] {
    auto res = descriptorNames;
    // _descList sorts the fragment descriptors by name
    auto frags = Descriptors::getFragmentDescriptorNames();
    std::sort(frags.begin(), frags.end());
    res.insert(res.end(), frags.begin(), frags.end());
    return res;
  }();
  return names;
}

std::string getPropertyName(const std::string &name) {
  auto it = propertyAliases.find(name);
  return it == propertyAliases.end() ? name : it->second;
}

NormalizedProperties::NormalizedProperties()
    : NormalizedProperties(getNormalizedDescriptorNames(), defaultTablesPtr()) {
}

NormalizedProperties::NormalizedProperties(const HistogramSet &tables)
    : NormalizedProperties(getNormalizedDescriptorNames(), tables) {}

NormalizedProperties::NormalizedProperties(
    const std::vector<std::string> &names, const HistogramSet &tables)
    : NormalizedProperties(names,
                           std::make_shared<const HistogramSet>(tables)) {}

NormalizedProperties::NormalizedProperties(
    const std::vector<std::string> &names,
    std::shared_ptr<const HistogramSet> tables)
    : Properties(toPropertyNames(names)),
      d_descriptorNames(names),
      d_tables(std::move(tables)) {
  // as in descriptastorus, failures normalize to 0.0
  d_tableIndex.reserve(d_descriptorNames.size());
  for (const auto &name : d_descriptorNames) {
    d_tableIndex.push_back(d_tables->getTableIndex(name));
  }
}

std::vector<double> NormalizedProperties::computeRawProperties(
    const ROMol &mol) const {
  return computeValues(mol);
}

std::vector<double> NormalizedProperties::computeProperties(
    const ROMol &mol, bool annotate) const {
  auto res = computeRawProperties(mol);
  for (size_t i = 0; i < res.size(); ++i) {
    auto idx = d_tableIndex[i];
    if (idx < 0 || std::isnan(res[i])) {
      res[i] = 0.0;
    } else {
      res[i] = d_tables->getTable(static_cast<size_t>(idx)).normalize(res[i]);
    }
    if (annotate) {
      mol.setProp<double>(d_descriptorNames[i], res[i]);
    }
  }
  return res;
}

std::vector<double> calcDescriptorValues(const ROMol &mol) {
  return defaultProperties().computeRawProperties(mol);
}

std::vector<double> calcNormalizedDescriptors(const ROMol &mol,
                                              const HistogramSet &tables) {
  // the tables outlive this call, so share them rather than copying
  std::shared_ptr<const HistogramSet> shared(&tables,
                                             [](const HistogramSet *) {});
  return NormalizedProperties(getNormalizedDescriptorNames(), shared)
      .computeProperties(mol);
}

std::vector<double> calcNormalizedDescriptors(const ROMol &mol) {
  return defaultProperties().computeProperties(mol);
}

}  // namespace NormalizedDescriptors
}  // namespace RDKit
