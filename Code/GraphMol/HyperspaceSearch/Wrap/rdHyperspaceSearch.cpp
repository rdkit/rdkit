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

python::list hitMolecules_helper(
    const HyperspaceSearch::SubstructureResults &res) {
  python::list pyres;
  for (auto &r : res.hitMolecules()) {
    pyres.append(boost::shared_ptr<ROMol>(new ROMol(*r)));
  }
  return pyres;
}

struct HyperspaceResults_wrapper {
  static void wrap() {
    std::string docString = "Used to return results of Hyperspace searches.";
    python::class_<RDKit::HyperspaceSearch::SubstructureResults>(
        "HyperspaceResult", docString.c_str(), python::no_init)
        .def("hitMolecules", hitMolecules_helper, python::args("self"),
             "A function returning hits from the search")
        .def_readonly("maxNumResults",
                      &HyperspaceSearch::SubstructureResults::maxNumResults,
                      "The upper bound on number of results possible.  There"
                      " may be fewer than this in practice for several reasons"
                      " such as duplicate reagent sets being removed or the"
                      " final product not matching the query even though the"
                      " synthons suggested they would.");
  }
};

HyperspaceSearch::SubstructureResults substructureSearch_helper(
    HyperspaceSearch::Hyperspace &self, const ROMol &query,
    const python::object &py_params) {
  HyperspaceSearch::HyperspaceSearchParams params;
  if (!py_params.is_none()) {
    params =
        python::extract<HyperspaceSearch::HyperspaceSearchParams>(py_params);
  }
  auto results = self.substructureSearch(query, params);
  return results;
}

void summariseHelper(HyperspaceSearch::Hyperspace &self) {
  self.summarise(std::cout);
}

BOOST_PYTHON_MODULE(rdHyperspaceSearch) {
  python::scope().attr("__doc__") =
      "Module containing implementation of Hyperspace search of"
      " Synthon-based chemical libraries such as Enamine REAL.";

  HyperspaceResults_wrapper::wrap();

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