//
//  Copyright (C) 2003-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/FragCatalog/FragCatalogEntry.h>
#include <GraphMol/FragCatalog/FragCatGenerator.h>
#include <GraphMol/FragCatalog/FragCatParams.h>
#include <Catalogs/Catalog.h>
#include <Catalogs/CatalogParams.h>

namespace nb = nanobind;
using namespace nb::literals;

void wrap_fragparams(nb::module_ &m) {
  // this exposed to be read only
  // i.e. one constructed there can be no changes done to the parameter object
  nb::class_<RDKit::FragCatParams>(m, "FragCatParams")
      .def(nb::init<unsigned int, unsigned int, std::string, double>(),
           "lLen"_a, "uLen"_a,
           "fgroupFilename"_a, "tol"_a = 1e-8)
      .def("GetTypeString", &RDKit::FragCatParams::getTypeStr)
      .def("GetUpperFragLength", &RDKit::FragCatParams::getUpperFragLength)
      .def("GetLowerFragLength", &RDKit::FragCatParams::getLowerFragLength)
      .def("GetTolerance", &RDKit::FragCatParams::getTolerance)
      .def("GetNumFuncGroups", &RDKit::FragCatParams::getNumFuncGroups)
      .def("GetFuncGroup",
           (const RDKit::ROMol *(RDKit::FragCatParams::*)(int) const) &
               RDKit::FragCatParams::getFuncGroup,
           nb::rv_policy::reference_internal, "fid"_a)
      .def("Serialize", [](const RDKit::FragCatParams &self) {
        const auto pkl = self.Serialize();
        return nb::bytes(pkl.c_str(), pkl.size());
      });
}
