//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDBoost/python.h>
#include <RDBoost/Wrap.h>

#include <GraphMol/ROMol.h>
#include <GraphMol/HyperspaceSearch/Hyperspace.h>

namespace python = boost::python;

namespace RDKit {

python::list substructureSearch_helper(HyperspaceSSSearch::Hyperspace &self,
                                       const ROMol &query) {
  auto results = self.substructureSearch(query, 3, -1);
  python::list pyres;
  for (auto &r : results) {
    pyres.append(boost::shared_ptr<ROMol>(r.release()));
  }
  return pyres;
}

void summariseHelper(HyperspaceSSSearch::Hyperspace &self) {
  self.summarise(std::cout);
}

BOOST_PYTHON_MODULE(rdHyperspaceSearch) {
  python::scope().attr("__doc__") =
      "Module containing implementation of Hyperspace search of"
      " Synthon-based chemical libraries such as Enamine REAL.";

  std::string docString = "HyperspaceSearch object.";
  python::class_<RDKit::HyperspaceSSSearch::Hyperspace, boost::noncopyable>(
      "Hyperspace", docString.c_str(), python::init<>())
      .def("ReadTextFile", &RDKit::HyperspaceSSSearch::Hyperspace::readTextFile,
           (python::arg("self"), python::arg("inFile")),
           "Reads text file of the sort used by ChemSpace/Enamine.")
      .def("ReadDBFile", &RDKit::HyperspaceSSSearch::Hyperspace::readDBFile,
           (python::arg("self"), python::arg("inFile")),
           "Reads binary database file.")
      .def("WriteDBFile", &RDKit::HyperspaceSSSearch::Hyperspace::writeDBFile,
           (python::arg("self"), python::arg("outFile")),
           "Writes binary database file.")
      .def("Summarise", &RDKit::summariseHelper, (python::arg("self")),
           "Writes a summary of the Hyperspace to stdout.")
      .def("SubstructureSearch", &RDKit::substructureSearch_helper,
           (python::arg("self"), python::arg("query")),
           "Does a substructure search in the Hyperspace.");
}

}  // namespace RDKit