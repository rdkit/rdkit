//
//  Copyright (C) 2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDBoost/python.h>
#include <GraphMol/GraphMol.h>
#include <RDBoost/Wrap.h>

#include <GraphMol/MolProcessing/Processing.h>

namespace python = boost::python;
using namespace RDKit;

namespace {
python::tuple getFingerprintsHelper(const std::string &fileName) {
  auto fps = MolProccesing::getFingerprintsForMolsInFile(fileName);

  python::list pyFingerprints;
  for (auto &fp : fps) {
    pyFingerprints.append(fp.release());
  }

  return python::tuple(pyFingerprints);
}
}  // namespace

BOOST_PYTHON_MODULE(rdMolProcessing) {
  python::scope().attr("__doc__") =
      "Module containing functions for working with groups of molecules";

  python::def("GetFingerprintsForMolsInFile", &getFingerprintsHelper,
              "returns the fingerprints for the molecules in a file");
}
