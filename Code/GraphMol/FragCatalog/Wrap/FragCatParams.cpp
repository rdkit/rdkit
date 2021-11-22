//
//  Copyright (C) 2003-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/python.h>
#include <string>

#include <DataStructs/BitVects.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FragCatalog/FragCatalogEntry.h>
#include <GraphMol/FragCatalog/FragCatGenerator.h>
#include <GraphMol/FragCatalog/FragCatParams.h>
#include <Catalogs/Catalog.h>
#include <Catalogs/CatalogParams.h>

namespace python = boost::python;
namespace RDKit {
struct fragparams_wrapper {
  static void wrap() {
    // this exposed to be read only
    // i.e. one constructed there can be no changes done to the parameter object
    python::class_<FragCatParams>(
        "FragCatParams",
        python::init<int, int, std::string, double>(
            (python::arg("lLen"), python::arg("uLen"),
             python::arg("fgroupFilename"), python::arg("tol") = 1e-8)))
        .def("GetTypeString", &FragCatParams::getTypeStr)
        .def("GetUpperFragLength", &FragCatParams::getUpperFragLength)
        .def("GetLowerFragLength", &FragCatParams::getLowerFragLength)
        .def("GetTolerance", &FragCatParams::getTolerance)
        .def("GetNumFuncGroups", &FragCatParams::getNumFuncGroups)
        .def("GetFuncGroup",
             (const ROMol* (FragCatParams::*)(int) const) &
                 FragCatParams::getFuncGroup,
             python::return_value_policy<python::reference_existing_object>())
        .def("Serialize", &FragCatParams::Serialize);
  };
};
}  // namespace RDKit

void wrap_fragparams() { RDKit::fragparams_wrapper::wrap(); }
