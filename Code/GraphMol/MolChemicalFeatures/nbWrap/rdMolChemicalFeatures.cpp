//
//  Copyright (C) 2003-2026 Rational Discovery LLC and other RDKit contributors
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

#include <fstream>
#include <sstream>

#include <GraphMol/MolChemicalFeatures/FeatureParser.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeature.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeatureFactory.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;

void wrap_MolChemicalFeat(nb::module_ &m);
void wrap_factory(nb::module_ &m);
void wrap_ChemicalFeatureUtils(nb::module_ &m);

namespace {
MolChemicalFeatureFactory *buildFeatFactory(const std::string &fileName) {
  std::ifstream inStream(fileName);
  if (!inStream.is_open()) {
    PyErr_SetString(PyExc_IOError,
                    ("File: " + fileName + " could not be opened.").c_str());
    throw nb::python_error();
  }
  return buildFeatureFactory(static_cast<std::istream &>(inStream));
}

MolChemicalFeatureFactory *buildFeatFactoryFromString(
    const std::string &fdefString) {
  std::istringstream inStream(fdefString);
  return buildFeatureFactory(static_cast<std::istream &>(inStream));
}
}  // namespace

NB_MODULE(rdMolChemicalFeatures, m) {
  m.doc() =
      R"DOC(Module containing from chemical feature and functions to generate the)DOC";

  nb::exception<FeatureFileParseException>(m, "FeatureFileParseException",
                                           PyExc_ValueError);

  m.def("BuildFeatureFactory", buildFeatFactory, "fileName"_a,
        "Construct a feature factory given a feature definition in a file",
        nb::rv_policy::take_ownership);
  m.def("BuildFeatureFactoryFromString", buildFeatFactoryFromString,
        "fdefString"_a,
        "Construct a feature factory given a feature definition block",
        nb::rv_policy::take_ownership);

  wrap_MolChemicalFeat(m);
  wrap_factory(m);
  wrap_ChemicalFeatureUtils(m);
}
