//
//  Copyright (C) 2018 Dan Nealschneider
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define NO_IMPORT_ARRAY
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <string>

// ours
#include <GraphMol/RDKitBase.h>
#include <GraphMol/StereoGroup.h>

namespace python = boost::python;

namespace RDKit {

std::string stereoGroupClassDoc =
    "A collection of atoms with a defined stereochemical relationship.\n\n"
    "Used to help represent a sample with unknown stereochemistry, or that "
    "is a mix\nof diastereomers.\n";

struct stereogroup_wrap {
  static void wrap() {
    python::enum_<RDKit::StereoGroupType>(
        "StereoGroupType")
        .value("STEREO_ABSOLUTE", RDKit::StereoGroupType::STEREO_ABSOLUTE)
        .value("STEREO_OR", RDKit::StereoGroupType::STEREO_OR)
        .value("STEREO_AND", RDKit::StereoGroupType::STEREO_AND)
        .export_values();

    python::class_<ROMol::ATOM_PTR_VECT>("AtomVector")
        .def(python::vector_indexing_suite<ROMol::ATOM_PTR_VECT>());

    python::class_<StereoGroup, boost::shared_ptr<StereoGroup>> (
        "StereoGroup", stereoGroupClassDoc.c_str(), python::init<>())
        .def("GetGroupType", &StereoGroup::getGroupType,
             "Returns the StereoGroupType.\n")
        .def("GetAtoms", &StereoGroup::getAtoms,
             "Access the atoms in the StereoGroup.\n",
             python::return_internal_reference<
                 1, python::with_custodian_and_ward_postcall<0, 1>>());
  }
};
}

void wrap_stereogroup() { RDKit::stereogroup_wrap::wrap(); }
