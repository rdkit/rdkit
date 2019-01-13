//
//  Copyright (C) 2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDBoost/python.h>
#include <boost/python/list.hpp>
#include <RDGeneral/Exceptions.h>
#include <GraphMol/RDKitBase.h>
#include <YAeHMOP/EHTTools.h>

namespace python = boost::python;

namespace RDKit {

namespace {}  // end of anonymous namespace

struct EHT_wrapper {
  static void wrap() {
    std::string docString = "";

    docString =
        "Runs the extended Hueckel calculation.\n"
        "ARGUMENTS:\n"
        "   - mol: molecule to use\n"
        "   - confId: (optional) conformation to use\n"
        "\n";
    python::def("RunMol", RDKit::EHTTools::runMol,
                (python::arg("mol"), python::arg("confId") = -1),
                docString.c_str());
  }
};

}  // end of namespace RDKit
BOOST_PYTHON_MODULE(rdEHTTools) {
  python::scope().attr("__doc__") =
      "Module containing interface to the YAeHMOP extended Hueckel library.";
  RDKit::EHT_wrapper::wrap();
}