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

NB_MODULE(rdStructChecker, m) {
  nb::enum_<RDKit::StructureCheck::StructChecker::StructureFlags>(
      m, "StructureFlags")
      .value("NO_CHANGE", RDKit::StructureCheck::StructChecker::NO_CHANGE)
      .value("BAD_MOLECULE",
             RDKit::StructureCheck::StructChecker::BAD_MOLECULE)
      .value("ALIAS_CONVERSION_FAILED",
             RDKit::StructureCheck::StructChecker::ALIAS_CONVERSION_FAILED)
      .value("STEREO_ERROR",
             RDKit::StructureCheck::StructChecker::STEREO_ERROR)
      .value("STEREO_FORCED_BAD",
             RDKit::StructureCheck::StructChecker::STEREO_FORCED_BAD)
      .value("ATOM_CLASH", RDKit::StructureCheck::StructChecker::ATOM_CLASH)
      .value("ATOM_CHECK_FAILED",
             RDKit::StructureCheck::StructChecker::ATOM_CHECK_FAILED)
      .value("SIZE_CHECK_FAILED",
             RDKit::StructureCheck::StructChecker::SIZE_CHECK_FAILED)
      .value("TRANSFORMED", RDKit::StructureCheck::StructChecker::TRANSFORMED)
      .value("FRAGMENTS_FOUND",
             RDKit::StructureCheck::StructChecker::FRAGMENTS_FOUND)
      .value("EITHER_WARNING",
             RDKit::StructureCheck::StructChecker::EITHER_WARNING)
      .value("DUBIOUS_STEREO_REMOVED",
             RDKit::StructureCheck::StructChecker::DUBIOUS_STEREO_REMOVED)
      .value("RECHARGED", RDKit::StructureCheck::StructChecker::RECHARGED)
      .value("STEREO_TRANSFORMED",
             RDKit::StructureCheck::StructChecker::STEREO_TRANSFORMED)
      .value("TEMPLATE_TRANSFORMED",
             RDKit::StructureCheck::StructChecker::TEMPLATE_TRANSFORMED)
      .value("TAUTOMER_TRANSFORMED",
             RDKit::StructureCheck::StructChecker::TAUTOMER_TRANSFORMED);

  nb::class_<RDKit::StructureCheck::StructCheckerOptions>(
      m, "StructCheckerOptions")
      .def(nb::init<>())
      .def_rw("AcidityLimit",
              &RDKit::StructureCheck::StructCheckerOptions::AcidityLimit)
      .def_rw("RemoveMinorFragments",
              &RDKit::StructureCheck::StructCheckerOptions::RemoveMinorFragments)
      .def_rw("DesiredCharge",
              &RDKit::StructureCheck::StructCheckerOptions::DesiredCharge)
      .def_rw("CheckCollisions",
              &RDKit::StructureCheck::StructCheckerOptions::CheckCollisions)
      .def_rw(
          "CollisionLimitPercent",
          &RDKit::StructureCheck::StructCheckerOptions::CollisionLimitPercent)
      .def_rw("MaxMolSize",
              &RDKit::StructureCheck::StructCheckerOptions::MaxMolSize)
      .def_rw("ConvertSText",
              &RDKit::StructureCheck::StructCheckerOptions::ConvertSText)
      .def_rw("StripZeros",
              &RDKit::StructureCheck::StructCheckerOptions::StripZeros)
      .def_rw("CheckStereo",
              &RDKit::StructureCheck::StructCheckerOptions::CheckStereo)
      .def_rw("ConvertAtomTexts",
              &RDKit::StructureCheck::StructCheckerOptions::ConvertAtomTexts)
      .def_rw("GroupsToSGroups",
              &RDKit::StructureCheck::StructCheckerOptions::GroupsToSGroups)
      .def_rw("Verbose",
              &RDKit::StructureCheck::StructCheckerOptions::Verbose)
      .def(
          "LoadGoodAugmentedAtoms",
          &RDKit::StructureCheck::StructCheckerOptions::loadGoodAugmentedAtoms,
          "path"_a,
          "Load the set of good augmented atoms from the specified file path")
      .def(
          "LoadAcidicAugmentedAtoms",
          &RDKit::StructureCheck::StructCheckerOptions::loadAcidicAugmentedAtoms,
          "path"_a,
          R"DOC(Load the set of acidic augmented atoms from the specified file
path)DOC")
      .def(
          "LoadAugmentedAtomTranslations",
          &RDKit::StructureCheck::StructCheckerOptions::
              loadAugmentedAtomTranslations,
          "path"_a,
          R"DOC(Load the set of acidic augmented atoms from the specified file
path)DOC");

  nb::class_<RDKit::StructureCheck::StructChecker>(m, "StructChecker")
      .def(nb::init<>())
      .def(nb::init<const RDKit::StructureCheck::StructCheckerOptions &>())
      .def(
          "CheckMolStructure",
          [](const RDKit::StructureCheck::StructChecker &checker,
             RDKit::ROMol &mol) {
            RDKit::RWMol &fixer = static_cast<RDKit::RWMol &>(mol);
            return static_cast<RDKit::StructureCheck::StructChecker::StructureFlags>(
                checker.checkMolStructure(fixer));
          },
          "mol"_a, "Check the structure and return a set of structure flags")
      .def_static(
          "StructureFlagsToString",
          &RDKit::StructureCheck::StructChecker::StructureFlagsToString,
          "flags"_a,
          "Return the structure flags as a human readable string")
      .def_static(
          "StringToStructureFlags",
          &RDKit::StructureCheck::StructChecker::StringToStructureFlags,
          "str"_a,
          R"DOC(Convert a comma separated string to the appropriate structure
flags)DOC");
}
