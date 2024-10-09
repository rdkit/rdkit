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

python::list substructureSearch_helper(HyperspaceSearch::Hyperspace &self,
                                       const ROMol &query,
                                       const python::object &py_params) {
  HyperspaceSearch::HyperspaceSearchParams params;
  if (!py_params.is_none()) {
    params =
        python::extract<HyperspaceSearch::HyperspaceSearchParams>(py_params);
  }
  auto results = self.substructureSearch(query, params);
  python::list pyres;
  for (auto &r : results) {
    pyres.append(boost::shared_ptr<ROMol>(r.release()));
  }
  return pyres;
}

void summariseHelper(HyperspaceSearch::Hyperspace &self) {
  self.summarise(std::cout);
}

BOOST_PYTHON_MODULE(rdHyperspaceSearch) {
  python::scope().attr("__doc__") =
      "Module containing implementation of Hyperspace search of"
      " Synthon-based chemical libraries such as Enamine REAL.";
  std::string docString = "HyperspaceSearch parameters.";
  python::class_<RDKit::HyperspaceSearch::HyperspaceSearchParams,
                 boost::noncopyable>("HyperspaceSearchParams",
                                     docString.c_str())
      .def_readwrite(
          "maxBondSplits",
          &RDKit::HyperspaceSearch::HyperspaceSearchParams::maxBondSplits,
          "The maximum number of bonds to break in the query."
          "  It should be no more than 1 less than the maximum"
          " number of Synthon sets in Hyperspace.  More than"
          " that doesn't matter, but will slow the search down"
          " to no good effect.  Default=3.")
      .def_readwrite("maxHits",
                     &RDKit::HyperspaceSearch::HyperspaceSearchParams::maxHits,
                     "The maximum number of hits to return.  Default=1000."
                     "Use -1 for no maximum.")
      .def_readwrite(
          "buildHits",
          &RDKit::HyperspaceSearch::HyperspaceSearchParams::buildHits,
          "If false, reports the maximum number of hits that"
          " the search could produce, but doesn't return them.");

  docString = "HyperspaceSearch object.";
  python::class_<RDKit::HyperspaceSearch::Hyperspace, boost::noncopyable>(
      "Hyperspace", docString.c_str(), python::init<>())
      .def("ReadTextFile", &RDKit::HyperspaceSearch::Hyperspace::readTextFile,
           (python::arg("self"), python::arg("inFile")),
           "Reads text file of the sort used by ChemSpace/Enamine.")
      .def("ReadDBFile", &RDKit::HyperspaceSearch::Hyperspace::readDBFile,
           (python::arg("self"), python::arg("inFile")),
           "Reads binary database file.")
      .def("WriteDBFile", &RDKit::HyperspaceSearch::Hyperspace::writeDBFile,
           (python::arg("self"), python::arg("outFile")),
           "Writes binary database file.")
      .def("Summarise", &RDKit::summariseHelper, (python::arg("self")),
           "Writes a summary of the Hyperspace to stdout.")
      .def("SubstructureSearch", &RDKit::substructureSearch_helper,
           (python::arg("self"), python::arg("query"),
            python::arg("params") = python::object()),
           "Does a substructure search in the Hyperspace.");
}

}  // namespace RDKit