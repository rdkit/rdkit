// $Id$
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define NO_IMPORT_ARRAY
#include <RDBoost/python.h>

#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>
#include <RDGeneral/Invariant.h>
#include <RDBoost/PySequenceHolder.h>
#include <ChemicalFeatures/FreeChemicalFeature.h>

namespace ChemicalFeatures {

// support pickling:
struct chemfeat_pickle_suite : rdkit_pickle_suite {
  static python::tuple getinitargs(const FreeChemicalFeature &self) {
    std::string res = self.toString();
    python::object retval = python::object(
        python::handle<>(PyBytes_FromStringAndSize(res.c_str(), res.length())));
    return python::make_tuple(retval);
  };
};

std::string featClassDoc =
    "Class to represent free chemical features.\n\
    These chemical features are not associated with a molecule, though they can be matched \n\
    to molecular features\n";
struct freefeat_wrapper {
  static void wrap() {
    python::class_<FreeChemicalFeature>(
        "FreeChemicalFeature", featClassDoc.c_str(),
        python::init<const std::string &>(python::args("self", "pickle")))
        .def(python::init<>(python::args("self"), "Default Constructor"))
        .def(python::init<std::string, std::string, const RDGeom::Point3D &,
                          int>(
            (python::arg("self"), python::arg("family"), python::arg("type"),
             python::arg("loc"), python::arg("id") = -1),
            "Constructor with family, type and location specified"))
        .def(python::init<std::string, const RDGeom::Point3D &>(
            python::args("self", "family", "loc"),
            "constructor with family and location specified, empty type and "
            "id"))
        .def("SetId", &FreeChemicalFeature::setId, python::args("self", "id"),
             "Set the id of the feature")
        .def("SetFamily", &FreeChemicalFeature::setFamily,
             python::args("self", "family"), "Set the family of the feature")
        .def("SetType", &FreeChemicalFeature::setType,
             python::args("self", "type"),
             "Set the sepcific type for the feature")
        .def("GetId", &FreeChemicalFeature::getId, python::args("self"),
             "Get the id of the feature")
        .def("GetFamily", &FreeChemicalFeature::getFamily,
             "Get the family of the feature",
             python::return_value_policy<python::copy_const_reference>(),
             python::args("self"))
        .def("GetType", &FreeChemicalFeature::getType,
             "Get the sepcific type for the feature",
             python::return_value_policy<python::copy_const_reference>(),
             python::args("self"))
        .def("SetPos", &FreeChemicalFeature::setPos,
             python::args("self", "loc"), "Set the feature position")
        .def("GetPos", &FreeChemicalFeature::getPos, python::args("self"),
             "Get the position of the feature")
        .def_pickle(chemfeat_pickle_suite());
  };
};
}  // namespace ChemicalFeatures

void wrap_freefeat() { ChemicalFeatures::freefeat_wrapper::wrap(); }
