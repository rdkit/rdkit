//
//  Copyright (C) 2004-2026 Rational Discovery LLC and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>

#include <algorithm>

#include <GraphMol/RDKitBase.h>
#include <RDGeneral/Exceptions.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeatureFactory.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeature.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;

namespace {
int getNumMolFeatures(const MolChemicalFeatureFactory &factory,
                      const ROMol &mol,
                      const std::string &includeOnly = "") {
  FeatSPtrList feats = factory.getFeaturesForMol(mol, includeOnly.c_str());
  return feats.size();
}

std::shared_ptr<MolChemicalFeature> getMolFeature(
    const MolChemicalFeatureFactory &factory, const ROMol &mol, int idx,
    const std::string &includeOnly, bool recompute, int confId) {
  static FeatSPtrList feats;
  if (recompute) {
    feats = factory.getFeaturesForMol(mol, includeOnly.c_str(), confId);
  }
  if (idx < 0 || idx >= static_cast<int>(feats.size())) {
    throw IndexErrorException(idx);
  }
  auto fi = feats.begin();
  for (int i = 0; i < idx; ++i) {
    ++fi;
  }
  // Aliasing constructor: wraps the boost::shared_ptr in a std::shared_ptr
  // so nanobind can manage the feature's lifetime correctly.
  auto holder = std::make_shared<FeatSPtr>(*fi);
  return std::shared_ptr<MolChemicalFeature>(holder, holder->get());
}

nb::tuple getFeatureFamilies(const MolChemicalFeatureFactory &factory) {
  std::vector<std::string> fams;
  for (auto iter = factory.beginFeatureDefs();
       iter != factory.endFeatureDefs(); ++iter) {
    std::string fam = (*iter)->getFamily();
    if (std::find(fams.begin(), fams.end(), fam) == fams.end()) {
      fams.push_back(fam);
    }
  }
  nb::list res;
  for (const auto &f : fams) {
    res.append(f);
  }
  return nb::tuple(res);
}

nb::dict getFeatureDefs(const MolChemicalFeatureFactory &factory) {
  nb::dict res;
  for (auto iter = factory.beginFeatureDefs();
       iter != factory.endFeatureDefs(); ++iter) {
    std::string key = (*iter)->getFamily() + "." + (*iter)->getType();
    res[key.c_str()] = (*iter)->getSmarts();
  }
  return res;
}
}  // namespace

void wrap_factory(nb::module_ &m) {
  nb::class_<MolChemicalFeatureFactory>(m, "MolChemicalFeatureFactory",
                                        "Class to featurize a molecule")
      .def("GetNumFeatureDefs",
           &MolChemicalFeatureFactory::getNumFeatureDefs,
           "Get the number of feature definitions")
      .def("GetFeatureFamilies", getFeatureFamilies,
           "Get a tuple of feature types")
      .def("GetFeatureDefs", getFeatureDefs,
           "Get a dictionary with SMARTS definitions for each feature type")
      .def("GetNumMolFeatures", getNumMolFeatures, "mol"_a,
           "includeOnly"_a = std::string(""),
           "Get the number of features the molecule has")
      .def("GetMolFeature", getMolFeature, "mol"_a, "idx"_a,
           "includeOnly"_a = std::string(""), "recompute"_a = true,
           "confId"_a = -1, nb::keep_alive<0, 1>(), nb::keep_alive<0, 2>(),
           "returns a particular feature (by index)");
}
