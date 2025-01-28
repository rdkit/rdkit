//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/python.h>

#include <string>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Chirality.h>

#include <RDBoost/Wrap.h>

namespace python = boost::python;

namespace RDKit {
struct chirality_wrapper {
  static void wrap() {
    python::enum_<Chirality::StereoType>("StereoType")
        .value("Unspecified", Chirality::StereoType::Unspecified)
        .value("Atom_Tetrahedral", Chirality::StereoType::Atom_Tetrahedral)
        .value("Atom_SquarePlanar", Chirality::StereoType::Atom_SquarePlanar)
        .value("Atom_TrigonalBipyramidal",
               Chirality::StereoType::Atom_TrigonalBipyramidal)
        .value("Atom_Octahedral", Chirality::StereoType::Atom_Octahedral)
        .value("Bond_Double", Chirality::StereoType::Bond_Double)
        .value("Bond_Cumulene_Even", Chirality::StereoType::Bond_Cumulene_Even)
        .value("Bond_Atropisomer", Chirality::StereoType::Bond_Atropisomer);
    python::enum_<Chirality::StereoSpecified>("StereoSpecified")
        .value("Unspecified", Chirality::StereoSpecified::Unspecified)
        .value("Specified", Chirality::StereoSpecified::Specified)
        .value("Unknown", Chirality::StereoSpecified::Unknown);
    python::enum_<Chirality::StereoDescriptor>("StereoDescriptor")
        .value("NoValue",
               Chirality::StereoDescriptor::None)  // can't use "None" in Python
        .value("Tet_CW", Chirality::StereoDescriptor::Tet_CW)
        .value("Tet_CCW", Chirality::StereoDescriptor::Tet_CCW)
        .value("Bond_Cis", Chirality::StereoDescriptor::Bond_Cis)
        .value("Bond_Trans", Chirality::StereoDescriptor::Bond_Trans);
    python::class_<Chirality::StereoInfo>("StereoInfo",
                                          "Class describing stereochemistry")
        .def_readonly("NOATOM", &Chirality::StereoInfo::NOATOM,
                      "marker for unspecified int values")
        .def_readwrite("type", &Chirality::StereoInfo::type,
                       "the type of stereo")
        .def_readwrite("specified", &Chirality::StereoInfo::specified,
                       "whether or not it is specified")
        .def_readwrite("centeredOn", &Chirality::StereoInfo::centeredOn,
                       "index of the item the stereo concerns")
        .def_readwrite("descriptor", &Chirality::StereoInfo::descriptor,
                       "stereo descriptor")
        .def_readwrite("permutation", &Chirality::StereoInfo::permutation,
                       "permutation index (used for non-tetrahedral chirality)")
        .def_readonly("controllingAtoms",
                      &Chirality::StereoInfo::controllingAtoms,
                      "indices of the atoms controlling the stereo");
  };
};
}  // namespace RDKit

void wrap_chirality() { RDKit::chirality_wrapper::wrap(); }
