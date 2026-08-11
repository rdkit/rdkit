//
//  Copyright (C) 2016-2026 Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <GraphMol/StructChecker/StructChecker.h>
#include <GraphMol/RDKitBase.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;
using namespace RDKit::StructureCheck;

NB_MODULE(rdStructChecker, m) {
  nb::enum_<StructChecker::StructureFlags>(m, "StructureFlags",
                                            nb::is_arithmetic())
      .value("NO_CHANGE", StructChecker::NO_CHANGE)
      .value("BAD_MOLECULE", StructChecker::BAD_MOLECULE)
      .value("ALIAS_CONVERSION_FAILED", StructChecker::ALIAS_CONVERSION_FAILED)
      .value("STEREO_ERROR", StructChecker::STEREO_ERROR)
      .value("STEREO_FORCED_BAD", StructChecker::STEREO_FORCED_BAD)
      .value("ATOM_CLASH", StructChecker::ATOM_CLASH)
      .value("ATOM_CHECK_FAILED", StructChecker::ATOM_CHECK_FAILED)
      .value("SIZE_CHECK_FAILED", StructChecker::SIZE_CHECK_FAILED)
      .value("TRANSFORMED", StructChecker::TRANSFORMED)
      .value("FRAGMENTS_FOUND", StructChecker::FRAGMENTS_FOUND)
      .value("EITHER_WARNING", StructChecker::EITHER_WARNING)
      .value("DUBIOUS_STEREO_REMOVED", StructChecker::DUBIOUS_STEREO_REMOVED)
      .value("RECHARGED", StructChecker::RECHARGED)
      .value("STEREO_TRANSFORMED", StructChecker::STEREO_TRANSFORMED)
      .value("TEMPLATE_TRANSFORMED", StructChecker::TEMPLATE_TRANSFORMED)
      .value("TAUTOMER_TRANSFORMED", StructChecker::TAUTOMER_TRANSFORMED);

  nb::class_<StructCheckerOptions>(m, "StructCheckerOptions")
      .def(nb::init<>())
      .def_rw("AcidityLimit", &StructCheckerOptions::AcidityLimit)
      .def_rw("RemoveMinorFragments", &StructCheckerOptions::RemoveMinorFragments)
      .def_rw("DesiredCharge", &StructCheckerOptions::DesiredCharge)
      .def_rw("CheckCollisions", &StructCheckerOptions::CheckCollisions)
      .def_rw("CollisionLimitPercent", &StructCheckerOptions::CollisionLimitPercent)
      .def_rw("MaxMolSize", &StructCheckerOptions::MaxMolSize)
      .def_rw("ConvertSText", &StructCheckerOptions::ConvertSText)
      .def_rw("StripZeros", &StructCheckerOptions::StripZeros)
      .def_rw("CheckStereo", &StructCheckerOptions::CheckStereo)
      .def_rw("ConvertAtomTexts", &StructCheckerOptions::ConvertAtomTexts)
      .def_rw("GroupsToSGroups", &StructCheckerOptions::GroupsToSGroups)
      .def_rw("Verbose", &StructCheckerOptions::Verbose)
      .def("LoadGoodAugmentedAtoms",
           &StructCheckerOptions::loadGoodAugmentedAtoms, "path"_a,
           "Load the set of good augmented atoms from the specified file path")
      .def("LoadAcidicAugmentedAtoms",
           &StructCheckerOptions::loadAcidicAugmentedAtoms, "path"_a,
           R"DOC(Load the set of acidic augmented atoms from the specified file
path)DOC")
      .def("LoadAugmentedAtomTranslations",
           &StructCheckerOptions::loadAugmentedAtomTranslations, "path"_a,
           R"DOC(Load the set of acidic augmented atoms from the specified file
path)DOC");

  nb::class_<StructChecker>(m, "StructChecker")
      .def(nb::init<>())
      .def(nb::init<const StructCheckerOptions &>())
      .def(
          "CheckMolStructure",
          [](const StructChecker &checker, ROMol &mol) -> unsigned {
            return checker.checkMolStructure(static_cast<RWMol &>(mol));
          },
          "mol"_a, "Check the structure and return a set of structure flags")
      .def_static("StructureFlagsToString",
                  &StructChecker::StructureFlagsToString, "flags"_a,
                  "Return the structure flags as a human readable string")
      .def_static("StringToStructureFlags",
                  &StructChecker::StringToStructureFlags, "str"_a,
                  R"DOC(Convert a comma separated string to the appropriate structure
flags)DOC");
}
