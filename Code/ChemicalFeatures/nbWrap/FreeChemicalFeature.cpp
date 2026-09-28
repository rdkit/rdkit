//
//  Copyright (C) 2004-2026 Rational Discovery LLC and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>

#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>
#include <RDGeneral/Invariant.h>
#include <ChemicalFeatures/FreeChemicalFeature.h>

namespace nb = nanobind;
using namespace nb::literals;
using ChemicalFeatures::FreeChemicalFeature;

void wrap_freefeat(nb::module_ &m) {
  nb::class_<FreeChemicalFeature>(
      m, "FreeChemicalFeature", nb::dynamic_attr(),
      R"DOC(Class to represent free chemical features.
These chemical features are not associated with a molecule, though they can be matched
to molecular features
)DOC")
      .def(nb::init<>(), "Default Constructor")
      .def(nb::init<const std::string &>(), "pickle"_a,
           "Constructor from a pickle string")
      .def(
          "__init__",
          [](FreeChemicalFeature *self, nb::bytes pkl) {
            new (self) FreeChemicalFeature(
                std::string(static_cast<const char *>(pkl.data()), pkl.size()));
          },
          "pickle"_a, "Constructor from a pickle bytes")
      .def(nb::init<std::string, std::string, const RDGeom::Point3D &, int>(),
           "family"_a, "type"_a, "loc"_a, "id"_a = -1,
           "Constructor with family, type and location specified")
      .def(nb::init<std::string, const RDGeom::Point3D &>(), "family"_a,
           "loc"_a,
           "constructor with family and location specified, empty type and id")
      .def("SetId", &FreeChemicalFeature::setId, "id"_a,
           "Set the id of the feature")
      .def("SetFamily", &FreeChemicalFeature::setFamily, "family"_a,
           "Set the family of the feature")
      .def("SetType", &FreeChemicalFeature::setType, "type"_a,
           "Set the specific type for the feature")
      .def("GetId", &FreeChemicalFeature::getId, "Get the id of the feature")
      .def("GetFamily", &FreeChemicalFeature::getFamily,
           "Get the family of the feature")
      .def("GetType", &FreeChemicalFeature::getType,
           "Get the specific type for the feature")
      .def("SetPos", &FreeChemicalFeature::setPos, "loc"_a,
           "Set the feature position")
      .def("GetPos", &FreeChemicalFeature::getPos,
           "Get the position of the feature")
      .def("__getstate__",
           [](const FreeChemicalFeature &feat) {
             const auto pkl = feat.toString();
             return std::make_tuple(nb::bytes(pkl.c_str(), pkl.size()));
           })
      .def("__setstate__",
           [](FreeChemicalFeature &feat, const std::tuple<nb::bytes> &state) {
             const auto &pkl = std::get<0>(state);
             new (&feat) FreeChemicalFeature(std::string(
                 static_cast<const char *>(pkl.data()), pkl.size()));
           })
      .def("__setstate__",
           [](FreeChemicalFeature &feat, const std::tuple<std::string> &state) {
             new (&feat) FreeChemicalFeature(std::get<0>(state));
           });
}
