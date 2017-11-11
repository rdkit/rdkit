//
//  Copyright (C) 2017 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define NO_IMPORT_ARRAY
#include <RDBoost/python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <boost/python/list.hpp>
#include <RDGeneral/Exceptions.h>
#include <GraphMol/RDKitBase.h>
#include <CoordGen.h>

namespace python = boost::python;

namespace RDKit {


struct coordgen_wrapper {
  static void wrap() {
    std::string docString = "";

    docString =
        "Add 2D coordinates.\n"
        "ARGUMENTS:\n"
        "   - mol: molecule to modify\n"
        "\n";

    python::def(
        "AddCoords", (void (*)(ROMol &))CoordGen::addCoords,
        (python::arg("mol")),
        docString.c_str());
}
};

} // end of namespace RDKit
BOOST_PYTHON_MODULE(rdCoordGen) {
  python::scope().attr("__doc__") =
      "Module containing interface to the CoordGen library.";
  RDKit::coordgen_wrapper::wrap();
}
