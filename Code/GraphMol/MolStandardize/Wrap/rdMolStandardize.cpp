//
//  Copyright (C) 2018-2023 Susan H. Leung and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/Wrap.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolStandardize/MolStandardize.h>

#include <vector>

namespace python = boost::python;

namespace {
template <typename FUNCTYPE>
RDKit::ROMol *msHelper(const RDKit::ROMol *mol, python::object params,
                       FUNCTYPE func) {
  if (!mol) {
    throw_value_error("Molecule is None");
  }
  const RDKit::MolStandardize::CleanupParameters *ps =
      &RDKit::MolStandardize::defaultCleanupParameters;
  if (params) {
    ps = python::extract<RDKit::MolStandardize::CleanupParameters *>(params);
  }
  return static_cast<RDKit::ROMol *>(
      func(static_cast<const RDKit::RWMol *>(mol), *ps));
}

RDKit::ROMol *cleanupHelper(const RDKit::ROMol *mol, python::object params) {
  return msHelper(
      mol, params,
      static_cast<
          RDKit::RWMol *(*)(const RDKit::RWMol *,
                            const RDKit::MolStandardize::CleanupParameters &)>(
          RDKit::MolStandardize::cleanup));
}

RDKit::ROMol *normalizeHelper(const RDKit::ROMol *mol, python::object params) {
  return msHelper(mol, params, RDKit::MolStandardize::normalize);
}

RDKit::ROMol *reionizeHelper(const RDKit::ROMol *mol, python::object params) {
  return msHelper(mol, params, RDKit::MolStandardize::reionize);
}

RDKit::ROMol *removeFragsHelper(const RDKit::ROMol *mol,
                                python::object params) {
  return msHelper(mol, params, RDKit::MolStandardize::removeFragments);
}

RDKit::ROMol *canonicalTautomerHelper(const RDKit::ROMol *mol,
                                      python::object params) {
  return msHelper(mol, params, RDKit::MolStandardize::canonicalTautomer);
}

template <typename FUNCTYPE>
void inPlaceHelper(RDKit::ROMol *mol, python::object params, FUNCTYPE func) {
  if (!mol) {
    throw_value_error("Molecule is None");
  }
  const RDKit::MolStandardize::CleanupParameters *ps =
      &RDKit::MolStandardize::defaultCleanupParameters;
  if (params) {
    ps = python::extract<RDKit::MolStandardize::CleanupParameters *>(params);
  }
  func(*static_cast<RDKit::RWMol *>(mol), *ps);
}

template <typename FUNCTYPE>
void inPlaceHelper2(RDKit::ROMol *mol, python::object params,
                    bool skip_standardize, FUNCTYPE func) {
  if (!mol) {
    throw_value_error("Molecule is None");
  }
  const RDKit::MolStandardize::CleanupParameters *ps =
      &RDKit::MolStandardize::defaultCleanupParameters;
  if (params) {
    ps = python::extract<RDKit::MolStandardize::CleanupParameters *>(params);
  }
  func(*static_cast<RDKit::RWMol *>(mol), *ps, skip_standardize);
}
void cleanupInPlaceHelper(RDKit::ROMol *mol, python::object params) {
  inPlaceHelper(
      mol, params,
      static_cast<void (*)(RDKit::RWMol &,
                           const RDKit::MolStandardize::CleanupParameters &)>(
          RDKit::MolStandardize::cleanupInPlace));
}

void normalizeInPlaceHelper(RDKit::ROMol *mol, python::object params) {
  inPlaceHelper(
      mol, params,
      static_cast<void (*)(RDKit::RWMol &,
                           const RDKit::MolStandardize::CleanupParameters &)>(
          RDKit::MolStandardize::normalizeInPlace));
}

void reionizeInPlaceHelper(RDKit::ROMol *mol, python::object params) {
  inPlaceHelper(
      mol, params,
      static_cast<void (*)(RDKit::RWMol &,
                           const RDKit::MolStandardize::CleanupParameters &)>(
          RDKit::MolStandardize::reionizeInPlace));
}

void removeFragmentsInPlaceHelper(RDKit::ROMol *mol, python::object params) {
  inPlaceHelper(
      mol, params,
      static_cast<void (*)(RDKit::RWMol &,
                           const RDKit::MolStandardize::CleanupParameters &)>(
          RDKit::MolStandardize::removeFragmentsInPlace));
}

void fragmentParentInPlaceHelper(RDKit::ROMol *mol, python::object params,
                                 bool skip_standardize) {
  inPlaceHelper2(
      mol, params, skip_standardize,
      static_cast<void (*)(
          RDKit::RWMol &, const RDKit::MolStandardize::CleanupParameters &,
          bool)>(RDKit::MolStandardize::fragmentParentInPlace));
}

void stereoParentInPlaceHelper(RDKit::ROMol *mol, python::object params,
                               bool skip_standardize) {
  inPlaceHelper2(
      mol, params, skip_standardize,
      static_cast<void (*)(RDKit::RWMol &,
                           const RDKit::MolStandardize::CleanupParameters &,
                           bool)>(RDKit::MolStandardize::stereoParentInPlace));
}

void isotopeParentInPlaceHelper(RDKit::ROMol *mol, python::object params,
                                bool skip_standardize) {
  inPlaceHelper2(
      mol, params, skip_standardize,
      static_cast<void (*)(RDKit::RWMol &,
                           const RDKit::MolStandardize::CleanupParameters &,
                           bool)>(RDKit::MolStandardize::isotopeParentInPlace));
}

void chargeParentInPlaceHelper(RDKit::ROMol *mol, python::object params,
                               bool skip_standardize) {
  inPlaceHelper2(
      mol, params, skip_standardize,
      static_cast<void (*)(RDKit::RWMol &,
                           const RDKit::MolStandardize::CleanupParameters &,
                           bool)>(RDKit::MolStandardize::chargeParentInPlace));
}

void superParentInPlaceHelper(RDKit::ROMol *mol, python::object params,
                              bool skip_standardize) {
  inPlaceHelper2(
      mol, params, skip_standardize,
      static_cast<void (*)(RDKit::RWMol &,
                           const RDKit::MolStandardize::CleanupParameters &,
                           bool)>(RDKit::MolStandardize::superParentInPlace));
}

void tautomerParentInPlaceHelper(RDKit::ROMol *mol, python::object params,
                                 bool skip_standardize) {
  inPlaceHelper2(
      mol, params, skip_standardize,
      static_cast<void (*)(
          RDKit::RWMol &, const RDKit::MolStandardize::CleanupParameters &,
          bool)>(RDKit::MolStandardize::tautomerParentInPlace));
}

template <typename FUNCTYPE>
void mtinPlaceHelper(python::object pymols, int numThreads,
                     python::object params, FUNCTYPE func) {
  const RDKit::MolStandardize::CleanupParameters *ps =
      &RDKit::MolStandardize::defaultCleanupParameters;
  if (params) {
    ps = python::extract<RDKit::MolStandardize::CleanupParameters *>(params);
  }
  unsigned int nmols = python::extract<unsigned int>(pymols.attr("__len__")());
  std::vector<RDKit::RWMol *> mols(nmols);
  for (auto i = 0u; i < nmols; ++i) {
    auto mol = static_cast<RDKit::RWMol *>(
        python::extract<RDKit::ROMol *>(pymols[i])());
    mols[i] = mol;
  }
  {
    NOGIL gil;
    func(mols, numThreads, *ps);
  }
}
template <typename FUNCTYPE>
void mtinPlaceHelper2(python::object pymols, int numThreads,
                      python::object params, bool skip_standardize,
                      FUNCTYPE func) {
  const auto *ps = &RDKit::MolStandardize::defaultCleanupParameters;
  if (params) {
    ps = python::extract<RDKit::MolStandardize::CleanupParameters *>(params);
  }
  unsigned int nmols = python::extract<unsigned int>(pymols.attr("__len__")());
  std::vector<RDKit::RWMol *> mols(nmols);
  for (auto i = 0u; i < nmols; ++i) {
    auto mol = static_cast<RDKit::RWMol *>(
        python::extract<RDKit::ROMol *>(pymols[i])());
    mols[i] = mol;
  }
  {
    NOGIL gil;
    func(mols, numThreads, *ps, skip_standardize);
  }
}

void mtcleanupInPlaceHelper(python::object mols, int numThreads,
                            python::object params) {
  mtinPlaceHelper(
      mols, numThreads, params,
      static_cast<void (*)(std::vector<RDKit::RWMol *> &, int,
                           const RDKit::MolStandardize::CleanupParameters &)>(
          RDKit::MolStandardize::cleanupInPlace));
}

void mtnormalizeInPlaceHelper(python::object mols, int numThreads,
                              python::object params) {
  mtinPlaceHelper(
      mols, numThreads, params,
      static_cast<void (*)(std::vector<RDKit::RWMol *> &, int,
                           const RDKit::MolStandardize::CleanupParameters &)>(
          RDKit::MolStandardize::normalizeInPlace));
}

void mtreionizeInPlaceHelper(python::object mols, int numThreads,
                             python::object params) {
  mtinPlaceHelper(
      mols, numThreads, params,
      static_cast<void (*)(std::vector<RDKit::RWMol *> &, int,
                           const RDKit::MolStandardize::CleanupParameters &)>(
          RDKit::MolStandardize::reionizeInPlace));
}

void mtremoveFragmentsInPlaceHelper(python::object mols, int numThreads,
                                    python::object params) {
  mtinPlaceHelper(
      mols, numThreads, params,
      static_cast<void (*)(std::vector<RDKit::RWMol *> &, int,
                           const RDKit::MolStandardize::CleanupParameters &)>(
          RDKit::MolStandardize::removeFragmentsInPlace));
}

void mtfragmentParentInPlaceHelper(python::object mols, int numThreads,
                                   python::object params,
                                   bool skip_standardize) {
  mtinPlaceHelper2(mols, numThreads, params, skip_standardize,
                   static_cast<void (*)(
                       std::vector<RDKit::RWMol *> &, int,
                       const RDKit::MolStandardize::CleanupParameters &, bool)>(
                       RDKit::MolStandardize::fragmentParentInPlace));
}

void mtstereoParentInPlaceHelper(python::object mols, int numThreads,
                                 python::object params, bool skip_standardize) {
  mtinPlaceHelper2(
      mols, numThreads, params, skip_standardize,
      static_cast<void (*)(std::vector<RDKit::RWMol *> &, int,
                           const RDKit::MolStandardize::CleanupParameters &,
                           bool)>(RDKit::MolStandardize::stereoParentInPlace));
}

void mtisotopeParentInPlaceHelper(python::object mols, int numThreads,
                                  python::object params,
                                  bool skip_standardize) {
  mtinPlaceHelper2(
      mols, numThreads, params, skip_standardize,
      static_cast<void (*)(std::vector<RDKit::RWMol *> &, int,
                           const RDKit::MolStandardize::CleanupParameters &,
                           bool)>(RDKit::MolStandardize::isotopeParentInPlace));
}

void mtchargeParentInPlaceHelper(python::object mols, int numThreads,
                                 python::object params, bool skip_standardize) {
  mtinPlaceHelper2(
      mols, numThreads, params, skip_standardize,
      static_cast<void (*)(std::vector<RDKit::RWMol *> &, int,
                           const RDKit::MolStandardize::CleanupParameters &,
                           bool)>(RDKit::MolStandardize::chargeParentInPlace));
}

void mtsuperParentInPlaceHelper(python::object mols, int numThreads,
                                python::object params, bool skip_standardize) {
  mtinPlaceHelper2(
      mols, numThreads, params, skip_standardize,
      static_cast<void (*)(std::vector<RDKit::RWMol *> &, int,
                           const RDKit::MolStandardize::CleanupParameters &,
                           bool)>(RDKit::MolStandardize::superParentInPlace));
}

void mttautomerParentInPlaceHelper(python::object mols, int numThreads,
                                   python::object params,
                                   bool skip_standardize) {
  mtinPlaceHelper2(mols, numThreads, params, skip_standardize,
                   static_cast<void (*)(
                       std::vector<RDKit::RWMol *> &, int,
                       const RDKit::MolStandardize::CleanupParameters &, bool)>(
                       RDKit::MolStandardize::tautomerParentInPlace));
}

template <typename FUNCTYPE>
RDKit::ROMol *parentHelper(const RDKit::ROMol *mol, python::object params,
                           bool skip_standardize, FUNCTYPE func) {
  if (!mol) {
    throw_value_error("Molecule is None");
  }
  const RDKit::MolStandardize::CleanupParameters *ps =
      &RDKit::MolStandardize::defaultCleanupParameters;
  if (params) {
    ps = python::extract<RDKit::MolStandardize::CleanupParameters *>(params);
  }
  return static_cast<RDKit::ROMol *>(
      func(static_cast<const RDKit::RWMol &>(*mol), *ps, skip_standardize));
}

RDKit::ROMol *tautomerParentHelper(const RDKit::ROMol *mol,
                                   python::object params,
                                   bool skip_standardize) {
  return parentHelper(mol, params, skip_standardize,
                      RDKit::MolStandardize::tautomerParent);
}
RDKit::ROMol *fragmentParentHelper(const RDKit::ROMol *mol,
                                   python::object params,
                                   bool skip_standardize) {
  return parentHelper(mol, params, skip_standardize,
                      RDKit::MolStandardize::fragmentParent);
}
RDKit::ROMol *stereoParentHelper(const RDKit::ROMol *mol, python::object params,
                                 bool skip_standardize) {
  return parentHelper(mol, params, skip_standardize,
                      RDKit::MolStandardize::stereoParent);
}
RDKit::ROMol *isotopeParentHelper(const RDKit::ROMol *mol,
                                  python::object params,
                                  bool skip_standardize) {
  return parentHelper(mol, params, skip_standardize,
                      RDKit::MolStandardize::isotopeParent);
}
RDKit::ROMol *chargeParentHelper(const RDKit::ROMol *mol, python::object params,
                                 bool skip_standardize) {
  return parentHelper(mol, params, skip_standardize,
                      RDKit::MolStandardize::chargeParent);
}
RDKit::ROMol *superParentHelper(const RDKit::ROMol *mol, python::object params,
                                bool skip_standardize) {
  return parentHelper(mol, params, skip_standardize,
                      RDKit::MolStandardize::superParent);
}
RDKit::ROMol *disconnectOrganometallicsHelper(RDKit::ROMol &mol,
                                              python::object params) {
  if (params) {
    RDKit::MolStandardize::MetalDisconnectorOptions *mdo =
        python::extract<RDKit::MolStandardize::MetalDisconnectorOptions *>(
            params);
    return RDKit::MolStandardize::disconnectOrganometallics(mol, *mdo);
  } else {
    return RDKit::MolStandardize::disconnectOrganometallics(mol);
  }
}
void disconnectOrganometallicsInPlaceHelper(RDKit::ROMol *mol,
                                            python::object params) {
  if (params) {
    RDKit::MolStandardize::MetalDisconnectorOptions *mdo =
        python::extract<RDKit::MolStandardize::MetalDisconnectorOptions *>(
            params);
    return RDKit::MolStandardize::disconnectOrganometallicsInPlace(
        *static_cast<RDKit::RWMol *>(mol), *mdo);
  } else {
    return RDKit::MolStandardize::disconnectOrganometallicsInPlace(
        *static_cast<RDKit::RWMol *>(mol));
  }
}

}  // namespace

void wrap_validate();
void wrap_charge();
void wrap_metal();
void wrap_fragment();
void wrap_normalize();
void wrap_tautomer();
void wrap_pipeline();

BOOST_PYTHON_MODULE(rdMolStandardize) {
  python::scope().attr("__doc__") =
      "Module containing functions for molecular standardization";

  bool noproxy = true;
  RegisterVectorConverter<RDKit::ROMOL_SPTR>("MOL_SPTR_VECT", noproxy);

  std::string docString = "";

  python::class_<RDKit::MolStandardize::CleanupParameters, boost::noncopyable>(
      "CleanupParameters", "Parameters controlling molecular standardization")
      .def_readwrite("normalizationsFile",
                     &RDKit::MolStandardize::CleanupParameters::normalizations,
                     "file containing the normalization transformations")
      .def_readwrite("acidbaseFile",
                     &RDKit::MolStandardize::CleanupParameters::acidbaseFile,
                     "file containing the acid and base definitions")
      .def_readwrite("fragmentFile",
                     &RDKit::MolStandardize::CleanupParameters::fragmentFile,
                     "file containing the acid and base definitions")
      .def_readwrite(
          "tautomerTransformsFile",
          &RDKit::MolStandardize::CleanupParameters::tautomerTransforms,
          "file containing the tautomer transformations")
      .def_readwrite("maxRestarts",
                     &RDKit::MolStandardize::CleanupParameters::maxRestarts,
                     "maximum number of restarts")
      .def_readwrite("preferOrganic",
                     &RDKit::MolStandardize::CleanupParameters::preferOrganic,
                     "prefer organic fragments to inorganic ones when deciding "
                     "what to keep")
      .def_readwrite("doCanonical",
                     &RDKit::MolStandardize::CleanupParameters::doCanonical,
                     "apply atom-order dependent normalizations (like "
                     "uncharging) in a canonical order")
      .def_readwrite("maxTautomers",
                     &RDKit::MolStandardize::CleanupParameters::maxTautomers,
                     "maximum number of tautomers to generate (defaults to "
                     "1000)")
      .def_readwrite("maxTransforms",
                     &RDKit::MolStandardize::CleanupParameters::maxTransforms,
                     "maximum number of transforms to apply during tautomer "
                     "enumeration (defaults to 1000)")
      .def_readwrite(
          "tautomerRemoveSp3Stereo",
          &RDKit::MolStandardize::CleanupParameters::tautomerRemoveSp3Stereo,
          "remove stereochemistry from sp3 centers involved in "
          "tautomerism (defaults to True)")
      .def_readwrite(
          "tautomerRemoveBondStereo",
          &RDKit::MolStandardize::CleanupParameters::tautomerRemoveBondStereo,
          "remove stereochemistry from double bonds involved in "
          "tautomerism (defaults to True)")
      .def_readwrite(
          "tautomerRemoveIsotopicHs",
          &RDKit::MolStandardize::CleanupParameters::tautomerRemoveIsotopicHs,
          "remove isotopic Hs from centers involved in "
          "tautomerism (defaults to True)")
      .def_readwrite(
          "tautomerReassignStereo",
          &RDKit::MolStandardize::CleanupParameters::tautomerReassignStereo,
          "call AssignStereochemistry on all generated tautomers "
          "(defaults to True)")
      .def_readwrite("largestFragmentChooserUseAtomCount",
                     &RDKit::MolStandardize::CleanupParameters::
                         largestFragmentChooserUseAtomCount,
                     "Whether LargestFragmentChooser should use atom "
                     "count as main criterion before MW (defaults to True)")
      .def_readwrite("largestFragmentChooserCountHeavyAtomsOnly",
                     &RDKit::MolStandardize::CleanupParameters::
                         largestFragmentChooserCountHeavyAtomsOnly,
                     "whether LargestFragmentChooser should only count "
                     "heavy atoms (defaults to False)");

  python::def("UpdateParamsFromJSON",
              &RDKit::MolStandardize::updateCleanupParamsFromJSON,
              python::args("params", "json"),
              "updates the cleanup parameters from the provided JSON string");

  docString = "Standardizes a molecule";
  python::def("Cleanup", cleanupHelper,
              (python::arg("mol"), python::arg("params") = python::object()),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  docString = "Standardizes a molecule in place";
  python::def("CleanupInPlace", cleanupInPlaceHelper,
              (python::arg("mol"), python::arg("params") = python::object()),
              docString.c_str());
  docString = "Standardizes multiple molecules in place";
  python::def("CleanupInPlace", mtcleanupInPlaceHelper,
              (python::arg("mols"), python::arg("numThreads"),
               python::arg("params") = python::object()),
              docString.c_str());
  docString = "Convenience function for standardizing a SMILES";
  python::def("StandardizeSmiles", RDKit::MolStandardize::standardizeSmiles,
              (python::arg("smiles")), docString.c_str());
  docString =
      "Returns the tautomer parent of a given molecule. The fragment parent is "
      "the standardized canonical tautomer of the molecule";
  python::def("TautomerParent", tautomerParentHelper,
              (python::arg("mol"), python::arg("params") = python::object(),
               python::arg("skipStandardize") = false),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  docString = "Generates the tautomer parent in place";
  python::def("TautomerParentInPlace", tautomerParentInPlaceHelper,
              (python::arg("mol"), python::arg("params") = python::object(),
               python::arg("skipStandardize") = false),
              docString.c_str());
  docString = "Generates the tautomer parent in place for multiple molecules";
  python::def("TautomerParentInPlace", mttautomerParentInPlaceHelper,
              (python::arg("mols"), python::arg("numThreads"),
               python::arg("params") = python::object(),
               python::arg("skipStandardize") = false),
              docString.c_str());

  docString = "Returns the largest fragment after doing a cleanup";
  python::def("FragmentParent", fragmentParentHelper,
              (python::arg("mol"), python::arg("params") = python::object(),
               python::arg("skipStandardize") = false),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  docString = "Generates the largest fragment in place";
  python::def("FragmentParentInPlace", fragmentParentInPlaceHelper,
              (python::arg("mol"), python::arg("params") = python::object(),
               python::arg("skipStandardize") = false),
              docString.c_str());
  docString = "Generates the largest fragment in place for multiple molecules";
  python::def("FragmentParentInPlace", mtfragmentParentInPlaceHelper,
              (python::arg("mols"), python::arg("numThreads"),
               python::arg("params") = python::object(),
               python::arg("skipStandardize") = false),
              docString.c_str());

  python::def("StereoParent", stereoParentHelper,
              (python::arg("mol"), python::arg("params") = python::object(),
               python::arg("skipStandardize") = false),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  docString = "Generates the stereo parent in place";
  python::def("StereoParentInPlace", stereoParentInPlaceHelper,
              (python::arg("mol"), python::arg("params") = python::object(),
               python::arg("skipStandardize") = false),
              docString.c_str());
  docString = "Generates the stereo parent in place for multiple molecules";
  python::def("StereoParentInPlace", mtstereoParentInPlaceHelper,
              (python::arg("mols"), python::arg("numThreads"),
               python::arg("params") = python::object(),
               python::arg("skipStandardize") = false),
              docString.c_str());
  docString = "removes all isotopes specifications from the given molecule";
  python::def("IsotopeParent", isotopeParentHelper,
              (python::arg("mol"), python::arg("params") = python::object(),
               python::arg("skipStandardize") = false),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  docString = "Generates the isotope parent in place";
  python::def("IsotopeParentInPlace", isotopeParentInPlaceHelper,
              (python::arg("mol"), python::arg("params") = python::object(),
               python::arg("skipStandardize") = false),
              docString.c_str());
  docString = "Generates the isotope parent in place for multiple molecules";
  python::def("IsotopeParentInPlace", mtisotopeParentInPlaceHelper,
              (python::arg("mols"), python::arg("numThreads"),
               python::arg("params") = python::object(),
               python::arg("skipStandardize") = false),
              docString.c_str());
  docString = "Returns the uncharged version of the largest fragment";
  python::def("ChargeParent", chargeParentHelper,
              (python::arg("mol"), python::arg("params") = python::object(),
               python::arg("skipStandardize") = false),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  docString = "Generates the charge parent in place";
  python::def("ChargeParentInPlace", chargeParentInPlaceHelper,
              (python::arg("mol"), python::arg("params") = python::object(),
               python::arg("skipStandardize") = false),
              docString.c_str());
  docString = "Generates the chargeparent in place for multiple molecules";
  python::def("ChargeParentInPlace", mtchargeParentInPlaceHelper,
              (python::arg("mols"), python::arg("numThreads"),
               python::arg("params") = python::object(),
               python::arg("skipStandardize") = false),
              docString.c_str());

  docString =
      "Returns the super parent. The super parent is the fragment, charge, "
      "isotope, stereo, and tautomer parent of the molecule.";
  python::def("SuperParent", superParentHelper,
              (python::arg("mol"), python::arg("params") = python::object(),
               python::arg("skipStandardize") = false),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  docString = "Generates the super parent in place";
  python::def("SuperParentInPlace", superParentInPlaceHelper,
              (python::arg("mol"), python::arg("params") = python::object(),
               python::arg("skipStandardize") = false),
              docString.c_str());
  docString = "Generates the super parent in place for multiple molecules";
  python::def("SuperParentInPlace", mtsuperParentInPlaceHelper,
              (python::arg("mols"), python::arg("numThreads"),
               python::arg("params") = python::object(),
               python::arg("skipStandardize") = false),
              docString.c_str());

  docString =
      "Applies a series of standard transformations to correct functional "
      "groups and recombine charges";
  python::def("Normalize", normalizeHelper,
              (python::arg("mol"), python::arg("params") = python::object()),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  docString =
      "Applies a series of standard transformations to correct functional "
      "groups and recombine charges, modifies the input molecule";
  python::def("NormalizeInPlace", normalizeInPlaceHelper,
              (python::arg("mol"), python::arg("params") = python::object()),
              docString.c_str());
  docString = "Normalizes multiple molecules in place";
  python::def("NormalizeInPlace", mtnormalizeInPlaceHelper,
              (python::arg("mols"), python::arg("numThreads"),
               python::arg("params") = python::object()),
              docString.c_str());
  docString = "Ensures the strongest acid groups are charged first";
  python::def("Reionize", reionizeHelper,
              (python::arg("mol"), python::arg("params") = python::object()),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  docString =
      "Ensures the strongest acid groups are charged first, modifies the input molecule";
  python::def("ReionizeInPlace", reionizeInPlaceHelper,
              (python::arg("mol"), python::arg("params") = python::object()),
              docString.c_str());
  docString = "Reionizes multiple molecules in place";
  python::def("ReionizeInPlace", mtreionizeInPlaceHelper,
              (python::arg("mols"), python::arg("numThreads"),
               python::arg("params") = python::object()),
              docString.c_str());
  docString = "Removes fragments from the molecule";
  python::def("RemoveFragments", removeFragsHelper,
              (python::arg("mol"), python::arg("params") = python::object()),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  docString =
      "Removes fragments from the molecule, modifies the input molecule";
  python::def("RemoveFragmentsInPlace", removeFragmentsInPlaceHelper,
              (python::arg("mol"), python::arg("params") = python::object()),
              docString.c_str());
  docString = "Removes fragments from multiple molecules in place";
  python::def("RemoveFragmentsInPlace", mtremoveFragmentsInPlaceHelper,
              (python::arg("mols"), python::arg("numThreads"),
               python::arg("params") = python::object()),
              docString.c_str());
  docString = "Returns the canonical tautomer for the molecule";
  python::def("CanonicalTautomer", canonicalTautomerHelper,
              (python::arg("mol"), python::arg("params") = python::object()),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Returns the molecule disconnected using the organometallics"
      " rules.";
  python::def("DisconnectOrganometallics", disconnectOrganometallicsHelper,
              (python::arg("mol"), python::arg("params") = python::object()),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  docString =
      "Disconnects the molecule using the organometallics"
      " rules, modifies the input molecule";
  python::def("DisconnectOrganometallicsInPlace",
              disconnectOrganometallicsInPlaceHelper,
              (python::arg("mol"), python::arg("params") = python::object()),
              docString.c_str());

  wrap_validate();
  wrap_charge();
  wrap_metal();
  wrap_fragment();
  wrap_normalize();
  wrap_tautomer();
  wrap_pipeline();
}
