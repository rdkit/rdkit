//
//  Copyright (C) 2018-2026 Susan H. Leung and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolStandardize/MolStandardize.h>
#include <RDBoost/Wrap_nb.h>

#include <vector>

namespace nb = nanobind;
using namespace nb::literals;

namespace {
template <typename FUNCTYPE>
RDKit::ROMol *msHelper(nb::object pymol, nb::object params, FUNCTYPE func) {
  if (pymol.is_none()) {
    throw nb::value_error("Molecule is None");
  }
  const RDKit::ROMol *mol = nb::cast<RDKit::ROMol *>(pymol);
  const RDKit::MolStandardize::CleanupParameters *ps =
      &RDKit::MolStandardize::defaultCleanupParameters;
  if (!params.is_none()) {
    ps = nb::cast<RDKit::MolStandardize::CleanupParameters *>(params);
  }
  return static_cast<RDKit::ROMol *>(
      func(static_cast<const RDKit::RWMol *>(mol), *ps));
}

RDKit::ROMol *cleanupHelper(nb::object mol, nb::object params) {
  return msHelper(
      mol, params,
      static_cast<
          RDKit::RWMol *(*)(const RDKit::RWMol *,
                            const RDKit::MolStandardize::CleanupParameters &)>(
          RDKit::MolStandardize::cleanup));
}

RDKit::ROMol *normalizeHelper(nb::object mol, nb::object params) {
  return msHelper(mol, params, RDKit::MolStandardize::normalize);
}

RDKit::ROMol *reionizeHelper(nb::object mol, nb::object params) {
  return msHelper(mol, params, RDKit::MolStandardize::reionize);
}

RDKit::ROMol *removeFragsHelper(nb::object mol, nb::object params) {
  return msHelper(mol, params, RDKit::MolStandardize::removeFragments);
}

RDKit::ROMol *canonicalTautomerHelper(nb::object mol, nb::object params) {
  return msHelper(mol, params, RDKit::MolStandardize::canonicalTautomer);
}

template <typename FUNCTYPE>
void inPlaceHelper(RDKit::ROMol *mol, nb::object params, FUNCTYPE func) {
  if (!mol) {
    throw nb::value_error("Molecule is None");
  }
  const RDKit::MolStandardize::CleanupParameters *ps =
      &RDKit::MolStandardize::defaultCleanupParameters;
  if (!params.is_none()) {
    ps = nb::cast<RDKit::MolStandardize::CleanupParameters *>(params);
  }
  func(*static_cast<RDKit::RWMol *>(mol), *ps);
}

template <typename FUNCTYPE>
void inPlaceHelper2(RDKit::ROMol *mol, nb::object params,
                    bool skip_standardize, FUNCTYPE func) {
  if (!mol) {
    throw nb::value_error("Molecule is None");
  }
  const RDKit::MolStandardize::CleanupParameters *ps =
      &RDKit::MolStandardize::defaultCleanupParameters;
  if (!params.is_none()) {
    ps = nb::cast<RDKit::MolStandardize::CleanupParameters *>(params);
  }
  func(*static_cast<RDKit::RWMol *>(mol), *ps, skip_standardize);
}

void cleanupInPlaceHelper(RDKit::ROMol *mol, nb::object params) {
  inPlaceHelper(
      mol, params,
      static_cast<void (*)(RDKit::RWMol &,
                           const RDKit::MolStandardize::CleanupParameters &)>(
          RDKit::MolStandardize::cleanupInPlace));
}

void normalizeInPlaceHelper(RDKit::ROMol *mol, nb::object params) {
  inPlaceHelper(
      mol, params,
      static_cast<void (*)(RDKit::RWMol &,
                           const RDKit::MolStandardize::CleanupParameters &)>(
          RDKit::MolStandardize::normalizeInPlace));
}

void reionizeInPlaceHelper(RDKit::ROMol *mol, nb::object params) {
  inPlaceHelper(
      mol, params,
      static_cast<void (*)(RDKit::RWMol &,
                           const RDKit::MolStandardize::CleanupParameters &)>(
          RDKit::MolStandardize::reionizeInPlace));
}

void removeFragmentsInPlaceHelper(RDKit::ROMol *mol, nb::object params) {
  inPlaceHelper(
      mol, params,
      static_cast<void (*)(RDKit::RWMol &,
                           const RDKit::MolStandardize::CleanupParameters &)>(
          RDKit::MolStandardize::removeFragmentsInPlace));
}

void fragmentParentInPlaceHelper(RDKit::ROMol *mol, nb::object params,
                                 bool skip_standardize) {
  inPlaceHelper2(
      mol, params, skip_standardize,
      static_cast<void (*)(
          RDKit::RWMol &, const RDKit::MolStandardize::CleanupParameters &,
          bool)>(RDKit::MolStandardize::fragmentParentInPlace));
}

void stereoParentInPlaceHelper(RDKit::ROMol *mol, nb::object params,
                               bool skip_standardize) {
  inPlaceHelper2(
      mol, params, skip_standardize,
      static_cast<void (*)(RDKit::RWMol &,
                           const RDKit::MolStandardize::CleanupParameters &,
                           bool)>(RDKit::MolStandardize::stereoParentInPlace));
}

void isotopeParentInPlaceHelper(RDKit::ROMol *mol, nb::object params,
                                bool skip_standardize) {
  inPlaceHelper2(
      mol, params, skip_standardize,
      static_cast<void (*)(RDKit::RWMol &,
                           const RDKit::MolStandardize::CleanupParameters &,
                           bool)>(RDKit::MolStandardize::isotopeParentInPlace));
}

void chargeParentInPlaceHelper(RDKit::ROMol *mol, nb::object params,
                               bool skip_standardize) {
  inPlaceHelper2(
      mol, params, skip_standardize,
      static_cast<void (*)(RDKit::RWMol &,
                           const RDKit::MolStandardize::CleanupParameters &,
                           bool)>(RDKit::MolStandardize::chargeParentInPlace));
}

void superParentInPlaceHelper(RDKit::ROMol *mol, nb::object params,
                              bool skip_standardize) {
  inPlaceHelper2(
      mol, params, skip_standardize,
      static_cast<void (*)(RDKit::RWMol &,
                           const RDKit::MolStandardize::CleanupParameters &,
                           bool)>(RDKit::MolStandardize::superParentInPlace));
}

void tautomerParentInPlaceHelper(RDKit::ROMol *mol, nb::object params,
                                 bool skip_standardize) {
  inPlaceHelper2(
      mol, params, skip_standardize,
      static_cast<void (*)(
          RDKit::RWMol &, const RDKit::MolStandardize::CleanupParameters &,
          bool)>(RDKit::MolStandardize::tautomerParentInPlace));
}

template <typename FUNCTYPE>
void mtinPlaceHelper(nb::object pymols, int numThreads, nb::object params,
                     FUNCTYPE func) {
  const RDKit::MolStandardize::CleanupParameters *ps =
      &RDKit::MolStandardize::defaultCleanupParameters;
  if (!params.is_none()) {
    ps = nb::cast<RDKit::MolStandardize::CleanupParameters *>(params);
  }
  nb::list molList(pymols);
  unsigned int nmols = (unsigned int)nb::len(molList);
  std::vector<RDKit::RWMol *> mols(nmols);
  for (auto i = 0u; i < nmols; ++i) {
    auto mol =
        static_cast<RDKit::RWMol *>(nb::cast<RDKit::ROMol *>(molList[i]));
    mols[i] = mol;
  }
  {
    NOGIL gil;
    func(mols, numThreads, *ps);
  }
}

template <typename FUNCTYPE>
void mtinPlaceHelper2(nb::object pymols, int numThreads, nb::object params,
                      bool skip_standardize, FUNCTYPE func) {
  const auto *ps = &RDKit::MolStandardize::defaultCleanupParameters;
  if (!params.is_none()) {
    ps = nb::cast<RDKit::MolStandardize::CleanupParameters *>(params);
  }
  nb::list molList(pymols);
  unsigned int nmols = (unsigned int)nb::len(molList);
  std::vector<RDKit::RWMol *> mols(nmols);
  for (auto i = 0u; i < nmols; ++i) {
    auto mol =
        static_cast<RDKit::RWMol *>(nb::cast<RDKit::ROMol *>(molList[i]));
    mols[i] = mol;
  }
  {
    NOGIL gil;
    func(mols, numThreads, *ps, skip_standardize);
  }
}

void mtcleanupInPlaceHelper(nb::object mols, int numThreads,
                            nb::object params) {
  mtinPlaceHelper(
      mols, numThreads, params,
      static_cast<void (*)(std::vector<RDKit::RWMol *> &, int,
                           const RDKit::MolStandardize::CleanupParameters &)>(
          RDKit::MolStandardize::cleanupInPlace));
}

void mtnormalizeInPlaceHelper(nb::object mols, int numThreads,
                              nb::object params) {
  mtinPlaceHelper(
      mols, numThreads, params,
      static_cast<void (*)(std::vector<RDKit::RWMol *> &, int,
                           const RDKit::MolStandardize::CleanupParameters &)>(
          RDKit::MolStandardize::normalizeInPlace));
}

void mtreionizeInPlaceHelper(nb::object mols, int numThreads,
                             nb::object params) {
  mtinPlaceHelper(
      mols, numThreads, params,
      static_cast<void (*)(std::vector<RDKit::RWMol *> &, int,
                           const RDKit::MolStandardize::CleanupParameters &)>(
          RDKit::MolStandardize::reionizeInPlace));
}

void mtremoveFragmentsInPlaceHelper(nb::object mols, int numThreads,
                                    nb::object params) {
  mtinPlaceHelper(
      mols, numThreads, params,
      static_cast<void (*)(std::vector<RDKit::RWMol *> &, int,
                           const RDKit::MolStandardize::CleanupParameters &)>(
          RDKit::MolStandardize::removeFragmentsInPlace));
}

void mtfragmentParentInPlaceHelper(nb::object mols, int numThreads,
                                   nb::object params, bool skip_standardize) {
  mtinPlaceHelper2(mols, numThreads, params, skip_standardize,
                   static_cast<void (*)(
                       std::vector<RDKit::RWMol *> &, int,
                       const RDKit::MolStandardize::CleanupParameters &, bool)>(
                       RDKit::MolStandardize::fragmentParentInPlace));
}

void mtstereoParentInPlaceHelper(nb::object mols, int numThreads,
                                 nb::object params, bool skip_standardize) {
  mtinPlaceHelper2(
      mols, numThreads, params, skip_standardize,
      static_cast<void (*)(std::vector<RDKit::RWMol *> &, int,
                           const RDKit::MolStandardize::CleanupParameters &,
                           bool)>(RDKit::MolStandardize::stereoParentInPlace));
}

void mtisotopeParentInPlaceHelper(nb::object mols, int numThreads,
                                  nb::object params, bool skip_standardize) {
  mtinPlaceHelper2(
      mols, numThreads, params, skip_standardize,
      static_cast<void (*)(std::vector<RDKit::RWMol *> &, int,
                           const RDKit::MolStandardize::CleanupParameters &,
                           bool)>(RDKit::MolStandardize::isotopeParentInPlace));
}

void mtchargeParentInPlaceHelper(nb::object mols, int numThreads,
                                 nb::object params, bool skip_standardize) {
  mtinPlaceHelper2(
      mols, numThreads, params, skip_standardize,
      static_cast<void (*)(std::vector<RDKit::RWMol *> &, int,
                           const RDKit::MolStandardize::CleanupParameters &,
                           bool)>(RDKit::MolStandardize::chargeParentInPlace));
}

void mtsuperParentInPlaceHelper(nb::object mols, int numThreads,
                                nb::object params, bool skip_standardize) {
  mtinPlaceHelper2(
      mols, numThreads, params, skip_standardize,
      static_cast<void (*)(std::vector<RDKit::RWMol *> &, int,
                           const RDKit::MolStandardize::CleanupParameters &,
                           bool)>(RDKit::MolStandardize::superParentInPlace));
}

void mttautomerParentInPlaceHelper(nb::object mols, int numThreads,
                                   nb::object params, bool skip_standardize) {
  mtinPlaceHelper2(mols, numThreads, params, skip_standardize,
                   static_cast<void (*)(
                       std::vector<RDKit::RWMol *> &, int,
                       const RDKit::MolStandardize::CleanupParameters &, bool)>(
                       RDKit::MolStandardize::tautomerParentInPlace));
}

template <typename FUNCTYPE>
RDKit::ROMol *parentHelper(nb::object pymol, nb::object params,
                           bool skip_standardize, FUNCTYPE func) {
  if (pymol.is_none()) {
    throw nb::value_error("Molecule is None");
  }
  const RDKit::ROMol *mol = nb::cast<RDKit::ROMol *>(pymol);
  const RDKit::MolStandardize::CleanupParameters *ps =
      &RDKit::MolStandardize::defaultCleanupParameters;
  if (!params.is_none()) {
    ps = nb::cast<RDKit::MolStandardize::CleanupParameters *>(params);
  }
  return static_cast<RDKit::ROMol *>(
      func(static_cast<const RDKit::RWMol &>(*mol), *ps, skip_standardize));
}

RDKit::ROMol *tautomerParentHelper(nb::object mol, nb::object params,
                                   bool skip_standardize) {
  return parentHelper(mol, params, skip_standardize,
                      RDKit::MolStandardize::tautomerParent);
}

RDKit::ROMol *fragmentParentHelper(nb::object mol, nb::object params,
                                   bool skip_standardize) {
  return parentHelper(mol, params, skip_standardize,
                      RDKit::MolStandardize::fragmentParent);
}

RDKit::ROMol *stereoParentHelper(nb::object mol, nb::object params,
                                 bool skip_standardize) {
  return parentHelper(mol, params, skip_standardize,
                      RDKit::MolStandardize::stereoParent);
}

RDKit::ROMol *isotopeParentHelper(nb::object mol, nb::object params,
                                  bool skip_standardize) {
  return parentHelper(mol, params, skip_standardize,
                      RDKit::MolStandardize::isotopeParent);
}

RDKit::ROMol *chargeParentHelper(nb::object mol, nb::object params,
                                 bool skip_standardize) {
  return parentHelper(mol, params, skip_standardize,
                      RDKit::MolStandardize::chargeParent);
}

RDKit::ROMol *superParentHelper(nb::object mol, nb::object params,
                                bool skip_standardize) {
  return parentHelper(mol, params, skip_standardize,
                      RDKit::MolStandardize::superParent);
}

RDKit::ROMol *disconnectOrganometallicsHelper(RDKit::ROMol &mol,
                                              nb::object params) {
  if (!params.is_none()) {
    RDKit::MolStandardize::MetalDisconnectorOptions *mdo =
        nb::cast<RDKit::MolStandardize::MetalDisconnectorOptions *>(params);
    return RDKit::MolStandardize::disconnectOrganometallics(mol, *mdo);
  } else {
    return RDKit::MolStandardize::disconnectOrganometallics(mol);
  }
}

void disconnectOrganometallicsInPlaceHelper(RDKit::ROMol *mol,
                                            nb::object params) {
  if (!params.is_none()) {
    RDKit::MolStandardize::MetalDisconnectorOptions *mdo =
        nb::cast<RDKit::MolStandardize::MetalDisconnectorOptions *>(params);
    return RDKit::MolStandardize::disconnectOrganometallicsInPlace(
        *static_cast<RDKit::RWMol *>(mol), *mdo);
  } else {
    return RDKit::MolStandardize::disconnectOrganometallicsInPlace(
        *static_cast<RDKit::RWMol *>(mol));
  }
}

}  // namespace

void wrap_validate(nb::module_ &m);
void wrap_charge(nb::module_ &m);
void wrap_metal(nb::module_ &m);
void wrap_fragment(nb::module_ &m);
void wrap_normalize(nb::module_ &m);
void wrap_tautomer(nb::module_ &m);
void wrap_pipeline(nb::module_ &m);

NB_MODULE(rdMolStandardize, m) {
  m.doc() = "Module containing functions for molecular standardization";

  nb::class_<RDKit::MolStandardize::CleanupParameters>(
      m, "CleanupParameters",
      "Parameters controlling molecular standardization")
      .def(nb::init<>())
      .def_rw("normalizationsFile",
               &RDKit::MolStandardize::CleanupParameters::normalizations,
               "file containing the normalization transformations")
      .def_rw("acidbaseFile",
               &RDKit::MolStandardize::CleanupParameters::acidbaseFile,
               "file containing the acid and base definitions")
      .def_rw("fragmentFile",
               &RDKit::MolStandardize::CleanupParameters::fragmentFile,
               "file containing the acid and base definitions")
      .def_rw("tautomerTransformsFile",
               &RDKit::MolStandardize::CleanupParameters::tautomerTransforms,
               "file containing the tautomer transformations")
      .def_rw("maxRestarts",
               &RDKit::MolStandardize::CleanupParameters::maxRestarts,
               "maximum number of restarts")
      .def_rw("preferOrganic",
               &RDKit::MolStandardize::CleanupParameters::preferOrganic,
               "prefer organic fragments to inorganic ones when deciding "
               "what to keep")
      .def_rw("doCanonical",
               &RDKit::MolStandardize::CleanupParameters::doCanonical,
               "apply atom-order dependent normalizations (like "
               "uncharging) in a canonical order")
      .def_rw("maxTautomers",
               &RDKit::MolStandardize::CleanupParameters::maxTautomers,
               "maximum number of tautomers to generate (defaults to 1000)")
      .def_rw("maxTransforms",
               &RDKit::MolStandardize::CleanupParameters::maxTransforms,
               "maximum number of transforms to apply during tautomer "
               "enumeration (defaults to 1000)")
      .def_rw("tautomerRemoveSp3Stereo",
               &RDKit::MolStandardize::CleanupParameters::
                   tautomerRemoveSp3Stereo,
               "remove stereochemistry from sp3 centers involved in "
               "tautomerism (defaults to True)")
      .def_rw("tautomerRemoveBondStereo",
               &RDKit::MolStandardize::CleanupParameters::
                   tautomerRemoveBondStereo,
               "remove stereochemistry from double bonds involved in "
               "tautomerism (defaults to True)")
      .def_rw("tautomerRemoveIsotopicHs",
               &RDKit::MolStandardize::CleanupParameters::
                   tautomerRemoveIsotopicHs,
               "remove isotopic Hs from centers involved in "
               "tautomerism (defaults to True)")
      .def_rw("tautomerReassignStereo",
               &RDKit::MolStandardize::CleanupParameters::
                   tautomerReassignStereo,
               "call AssignStereochemistry on all generated tautomers "
               "(defaults to True)")
      .def_rw("largestFragmentChooserUseAtomCount",
               &RDKit::MolStandardize::CleanupParameters::
                   largestFragmentChooserUseAtomCount,
               "Whether LargestFragmentChooser should use atom "
               "count as main criterion before MW (defaults to True)")
      .def_rw("largestFragmentChooserCountHeavyAtomsOnly",
               &RDKit::MolStandardize::CleanupParameters::
                   largestFragmentChooserCountHeavyAtomsOnly,
               "whether LargestFragmentChooser should only count "
               "heavy atoms (defaults to False)")
      .def("__setattr__", &safeSetattr);

  m.def("UpdateParamsFromJSON",
        &RDKit::MolStandardize::updateCleanupParamsFromJSON, "params"_a,
        "json"_a,
        "updates the cleanup parameters from the provided JSON string");

  m.def("Cleanup", cleanupHelper, nb::arg("mol").none(), "params"_a = nb::none(),
        "Standardizes a molecule", nb::rv_policy::take_ownership);
  m.def("CleanupInPlace", cleanupInPlaceHelper, "mol"_a,
        "params"_a = nb::none(), "Standardizes a molecule in place");
  m.def("CleanupInPlace", mtcleanupInPlaceHelper, "mols"_a, "numThreads"_a,
        "params"_a = nb::none(), "Standardizes multiple molecules in place");
  m.def("StandardizeSmiles", RDKit::MolStandardize::standardizeSmiles,
        "smiles"_a, "Convenience function for standardizing a SMILES");
  m.def("TautomerParent", tautomerParentHelper, nb::arg("mol").none(),
        "params"_a = nb::none(), "skipStandardize"_a = false,
        R"DOC(Returns the tautomer parent of a given molecule. The fragment parent is
the standardized canonical tautomer of the molecule)DOC",
        nb::rv_policy::take_ownership);
  m.def("TautomerParentInPlace", tautomerParentInPlaceHelper, "mol"_a,
        "params"_a = nb::none(), "skipStandardize"_a = false,
        "Generates the tautomer parent in place");
  m.def("TautomerParentInPlace", mttautomerParentInPlaceHelper, "mols"_a,
        "numThreads"_a, "params"_a = nb::none(), "skipStandardize"_a = false,
        "Generates the tautomer parent in place for multiple molecules");

  m.def("FragmentParent", fragmentParentHelper, nb::arg("mol").none(),
        "params"_a = nb::none(), "skipStandardize"_a = false,
        "Returns the largest fragment after doing a cleanup",
        nb::rv_policy::take_ownership);
  m.def("FragmentParentInPlace", fragmentParentInPlaceHelper, "mol"_a,
        "params"_a = nb::none(), "skipStandardize"_a = false,
        "Generates the largest fragment in place");
  m.def("FragmentParentInPlace", mtfragmentParentInPlaceHelper, "mols"_a,
        "numThreads"_a, "params"_a = nb::none(), "skipStandardize"_a = false,
        "Generates the largest fragment in place for multiple molecules");

  m.def("StereoParent", stereoParentHelper, nb::arg("mol").none(), "params"_a = nb::none(),
        "skipStandardize"_a = false,
        "Returns the stereo parent of the molecule",
        nb::rv_policy::take_ownership);
  m.def("StereoParentInPlace", stereoParentInPlaceHelper, "mol"_a,
        "params"_a = nb::none(), "skipStandardize"_a = false,
        "Generates the stereo parent in place");
  m.def("StereoParentInPlace", mtstereoParentInPlaceHelper, "mols"_a,
        "numThreads"_a, "params"_a = nb::none(), "skipStandardize"_a = false,
        "Generates the stereo parent in place for multiple molecules");

  m.def("IsotopeParent", isotopeParentHelper, nb::arg("mol").none(),
        "params"_a = nb::none(), "skipStandardize"_a = false,
        "removes all isotopes specifications from the given molecule",
        nb::rv_policy::take_ownership);
  m.def("IsotopeParentInPlace", isotopeParentInPlaceHelper, "mol"_a,
        "params"_a = nb::none(), "skipStandardize"_a = false,
        "Generates the isotope parent in place");
  m.def("IsotopeParentInPlace", mtisotopeParentInPlaceHelper, "mols"_a,
        "numThreads"_a, "params"_a = nb::none(), "skipStandardize"_a = false,
        "Generates the isotope parent in place for multiple molecules");

  m.def("ChargeParent", chargeParentHelper, nb::arg("mol").none(), "params"_a = nb::none(),
        "skipStandardize"_a = false,
        "Returns the uncharged version of the largest fragment",
        nb::rv_policy::take_ownership);
  m.def("ChargeParentInPlace", chargeParentInPlaceHelper, "mol"_a,
        "params"_a = nb::none(), "skipStandardize"_a = false,
        "Generates the charge parent in place");
  m.def("ChargeParentInPlace", mtchargeParentInPlaceHelper, "mols"_a,
        "numThreads"_a, "params"_a = nb::none(), "skipStandardize"_a = false,
        "Generates the chargeparent in place for multiple molecules");

  m.def("SuperParent", superParentHelper, nb::arg("mol").none(), "params"_a = nb::none(),
        "skipStandardize"_a = false,
        R"DOC(Returns the super parent. The super parent is the fragment, charge,
isotope, stereo, and tautomer parent of the molecule.)DOC",
        nb::rv_policy::take_ownership);
  m.def("SuperParentInPlace", superParentInPlaceHelper, "mol"_a,
        "params"_a = nb::none(), "skipStandardize"_a = false,
        "Generates the super parent in place");
  m.def("SuperParentInPlace", mtsuperParentInPlaceHelper, "mols"_a,
        "numThreads"_a, "params"_a = nb::none(), "skipStandardize"_a = false,
        "Generates the super parent in place for multiple molecules");

  m.def("Normalize", normalizeHelper, nb::arg("mol").none(), "params"_a = nb::none(),
        R"DOC(Applies a series of standard transformations to correct functional
groups and recombine charges)DOC",
        nb::rv_policy::take_ownership);
  m.def("NormalizeInPlace", normalizeInPlaceHelper, "mol"_a,
        "params"_a = nb::none(),
        R"DOC(Applies a series of standard transformations to correct functional
groups and recombine charges, modifies the input molecule)DOC");
  m.def("NormalizeInPlace", mtnormalizeInPlaceHelper, "mols"_a, "numThreads"_a,
        "params"_a = nb::none(), "Normalizes multiple molecules in place");

  m.def("Reionize", reionizeHelper, nb::arg("mol").none(), "params"_a = nb::none(),
        "Ensures the strongest acid groups are charged first",
        nb::rv_policy::take_ownership);
  m.def("ReionizeInPlace", reionizeInPlaceHelper, "mol"_a,
        "params"_a = nb::none(),
        "Ensures the strongest acid groups are charged first, modifies the "
        "input molecule");
  m.def("ReionizeInPlace", mtreionizeInPlaceHelper, "mols"_a, "numThreads"_a,
        "params"_a = nb::none(), "Reionizes multiple molecules in place");

  m.def("RemoveFragments", removeFragsHelper, nb::arg("mol").none(),
        "params"_a = nb::none(), "Removes fragments from the molecule",
        nb::rv_policy::take_ownership);
  m.def("RemoveFragmentsInPlace", removeFragmentsInPlaceHelper, "mol"_a,
        "params"_a = nb::none(),
        "Removes fragments from the molecule, modifies the input molecule");
  m.def("RemoveFragmentsInPlace", mtremoveFragmentsInPlaceHelper, "mols"_a,
        "numThreads"_a, "params"_a = nb::none(),
        "Removes fragments from multiple molecules in place");

  m.def("CanonicalTautomer", canonicalTautomerHelper, nb::arg("mol").none(),
        "params"_a = nb::none(),
        "Returns the canonical tautomer for the molecule",
        nb::rv_policy::take_ownership);

  m.def("DisconnectOrganometallics", disconnectOrganometallicsHelper, "mol"_a,
        "params"_a = nb::none(),
        "Returns the molecule disconnected using the organometallics rules.",
        nb::rv_policy::take_ownership);
  m.def("DisconnectOrganometallicsInPlace",
        disconnectOrganometallicsInPlaceHelper, "mol"_a,
        "params"_a = nb::none(),
        "Disconnects the molecule using the organometallics rules, modifies "
        "the input molecule");

  wrap_validate(m);
  wrap_charge(m);
  wrap_metal(m);
  wrap_fragment(m);
  wrap_normalize(m);
  wrap_tautomer(m);
  wrap_pipeline(m);
}
