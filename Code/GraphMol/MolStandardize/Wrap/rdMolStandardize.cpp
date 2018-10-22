//
//  Copyright (C) 2018 Susan H. Leung
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
RDKit::ROMol *cleanupHelper(const RDKit::ROMol *mol, python::object params) {
  const RDKit::MolStandardize::CleanupParameters *ps =
      &RDKit::MolStandardize::defaultCleanupParameters;
  if (params) {
    ps = python::extract<RDKit::MolStandardize::CleanupParameters *>(params);
  }
  return static_cast<RDKit::ROMol *>(RDKit::MolStandardize::cleanup(
      *static_cast<const RDKit::RWMol *>(mol), *ps));
}

RDKit::ROMol *fragmentParentHelper(const RDKit::ROMol *mol,
                                   python::object params,
                                   bool skip_standardize) {
  const RDKit::MolStandardize::CleanupParameters *ps =
      &RDKit::MolStandardize::defaultCleanupParameters;
  if (params) {
    ps = python::extract<RDKit::MolStandardize::CleanupParameters *>(params);
  }
  return static_cast<RDKit::ROMol *>(RDKit::MolStandardize::fragmentParent(
      *static_cast<const RDKit::RWMol *>(mol), *ps, skip_standardize));
}

RDKit::ROMol *chargeParentHelper(const RDKit::ROMol *mol, python::object params,
                                 bool skip_standardize) {
  const RDKit::MolStandardize::CleanupParameters *ps =
      &RDKit::MolStandardize::defaultCleanupParameters;
  if (params) {
    ps = python::extract<RDKit::MolStandardize::CleanupParameters *>(params);
  }
  return static_cast<RDKit::ROMol *>(RDKit::MolStandardize::chargeParent(
      *static_cast<const RDKit::RWMol *>(mol), *ps, skip_standardize));
}

RDKit::ROMol *normalizeHelper(const RDKit::ROMol *mol, python::object params) {
  const RDKit::MolStandardize::CleanupParameters *ps =
      &RDKit::MolStandardize::defaultCleanupParameters;
  if (params) {
    ps = python::extract<RDKit::MolStandardize::CleanupParameters *>(params);
  }
  return static_cast<RDKit::ROMol *>(RDKit::MolStandardize::normalize(
      static_cast<const RDKit::RWMol *>(mol), *ps));
}

RDKit::ROMol *reionizeHelper(const RDKit::ROMol *mol, python::object params) {
  const RDKit::MolStandardize::CleanupParameters *ps =
      &RDKit::MolStandardize::defaultCleanupParameters;
  if (params) {
    ps = python::extract<RDKit::MolStandardize::CleanupParameters *>(params);
  }
  return static_cast<RDKit::ROMol *>(RDKit::MolStandardize::reionize(
      static_cast<const RDKit::RWMol *>(mol), *ps));
}

}  // namespace

void wrap_validate();
void wrap_charge();
void wrap_metal();
void wrap_fragment();
void wrap_normalize();

BOOST_PYTHON_MODULE(rdMolStandardize) {
  python::scope().attr("__doc__") =
      "Module containing functions for molecular standardization";

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
      .def_readwrite("maxTautomers",
                     &RDKit::MolStandardize::CleanupParameters::maxTautomers,
                     "maximum number of tautomers to generate")
      .def_readwrite("preferOrganic",
                     &RDKit::MolStandardize::CleanupParameters::preferOrganic,
                     "prefer organic fragments to inorganic ones when deciding "
                     "what to keep");

  docString = "Standardizes a molecule";
  python::def("Cleanup", cleanupHelper,
              (python::arg("mol"), python::arg("params") = python::object()),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  docString = "Convenience function fo standardizing a SMILES";
  python::def("StandardizeSmiles", RDKit::MolStandardize::standardizeSmiles,
              (python::arg("smiles")), docString.c_str());
  docString = "Returns the largest fragment after doing a cleanup";
  python::def("FragmentParent", fragmentParentHelper,
              (python::arg("mol"), python::arg("params") = python::object(),
               python::arg("skipStandardize") = false),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  docString = "Returns the uncharged version of the largest fragment";
  python::def("ChargeParent", chargeParentHelper,
              (python::arg("mol"), python::arg("params") = python::object(),
               python::arg("skipStandardize") = false),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  docString =
      "Applies a series of standard transformations to correct functional "
      "groups and recombine charges";
  python::def("Normalize", normalizeHelper,
              (python::arg("mol"), python::arg("params") = python::object()),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  docString = "Ensures the strongest acid groups are charged first";
  python::def("Reionize", reionizeHelper,
              (python::arg("mol"), python::arg("params") = python::object()),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  wrap_validate();
  wrap_charge();
  wrap_metal();
  wrap_fragment();
  wrap_normalize();
}
