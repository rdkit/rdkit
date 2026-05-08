//
//  Copyright (C) 2026 and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <string>
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Chirality.h>

namespace nb = nanobind;

namespace RDKit {
struct chirality_wrapper {
  static void wrap(nb::module_ &m) {
    nb::enum_<Chirality::StereoType>(m, "StereoType")
        .value("Unspecified", Chirality::StereoType::Unspecified)
        .value("Atom_Tetrahedral", Chirality::StereoType::Atom_Tetrahedral)
        .value("Atom_SquarePlanar", Chirality::StereoType::Atom_SquarePlanar)
        .value("Atom_TrigonalBipyramidal",
               Chirality::StereoType::Atom_TrigonalBipyramidal)
        .value("Atom_Octahedral", Chirality::StereoType::Atom_Octahedral)
        .value("Bond_Double", Chirality::StereoType::Bond_Double)
        .value("Bond_Cumulene_Even", Chirality::StereoType::Bond_Cumulene_Even)
        .value("Bond_Atropisomer", Chirality::StereoType::Bond_Atropisomer)
        .export_values();
    nb::enum_<Chirality::StereoSpecified>(m, "StereoSpecified")
        .value("Unspecified", Chirality::StereoSpecified::Unspecified)
        .value("Specified", Chirality::StereoSpecified::Specified)
        .value("Unknown", Chirality::StereoSpecified::Unknown)
        .export_values();
    nb::enum_<Chirality::StereoDescriptor>(m, "StereoDescriptor")
        .value("NoValue",
               Chirality::StereoDescriptor::None)  // can't use "None" in Python
        .value("Tet_CW", Chirality::StereoDescriptor::Tet_CW)
        .value("Tet_CCW", Chirality::StereoDescriptor::Tet_CCW)
        .value("Bond_Cis", Chirality::StereoDescriptor::Bond_Cis)
        .value("Bond_Trans", Chirality::StereoDescriptor::Bond_Trans)
        .export_values();
    nb::class_<Chirality::StereoInfo>(
        m, "StereoInfo", R"DOC(Class describing stereochemistry)DOC")
        .def_ro_static("NOATOM", &Atom::NOATOM,
                       R"DOC(marker for unspecified int values)DOC")
        .def_rw("type", &Chirality::StereoInfo::type,
                R"DOC(the type of stereo)DOC")
        .def_rw("specified", &Chirality::StereoInfo::specified,
                R"DOC(whether or not it is specified)DOC")
        .def_rw("centeredOn", &Chirality::StereoInfo::centeredOn,
                R"DOC(index of the item the stereo concerns)DOC")
        .def_rw("descriptor", &Chirality::StereoInfo::descriptor,
                R"DOC(stereo descriptor)DOC")
        .def_rw(
            "permutation", &Chirality::StereoInfo::permutation,
            R"DOC(permutation index (used for non-tetrahedral chirality))DOC")
        .def_ro("controllingAtoms", &Chirality::StereoInfo::controllingAtoms,
                R"DOC(indices of the atoms controlling the stereo)DOC");
  };
};
}  // namespace RDKit

void wrap_chirality(nb::module_ &m) { RDKit::chirality_wrapper::wrap(m); }
