//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

#include <string>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Chirality.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {
struct chiralityops_wrapper {
  static void wrap(nb::module_ &m) {
    m.def(
        "FindPotentialStereo",
        (std::vector<Chirality::StereoInfo> (*)(
            ROMol &, bool, bool))&Chirality::findPotentialStereo,
        "mol"_a, "cleanIt"_a = false, "flagPossible"_a = true,
        R"DOC(find potential stereo elements in a molecule and returns them as StereoInfo objects
Note that this function is still somewhat experimental and the API
and results may change in a future release.)DOC",
        nb::keep_alive<0, 1>());
    m.def(
        "CleanupStereoGroups", &Chirality::cleanupStereoGroups, "mol"_a,
        R"DOC(removes atoms without specified chirality from stereo groups)DOC");
  };
};
}  // namespace RDKit

void wrap_chiralityops(nb::module_ &m) { RDKit::chiralityops_wrapper::wrap(m); }
