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

#include <GraphMol/ROMol.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>

namespace python = boost::python;

namespace RDKit {

python::list hitMolecules_helper(
    const SynthonSpaceSearch::SubstructureResults &res) {
  python::list pyres;
  for (auto &r : res.getHitMolecules()) {
    pyres.append(boost::shared_ptr<ROMol>(new ROMol(*r)));
  }
  return pyres;
}

struct SubstructureResults_wrapper {
  static void wrap() {
    std::string docString = "Used to return results of SynthonSpace searches.";
    python::class_<RDKit::SynthonSpaceSearch::SubstructureResults>(
        "SubstructureResult", docString.c_str(), python::no_init)
        .def("GetHitMolecules", hitMolecules_helper, python::args("self"),
             "A function returning hits from the search")
        .def_readonly(
            "GetMaxNumResults",
            &SynthonSpaceSearch::SubstructureResults::getMaxNumResults,
            "The upper bound on number of results possible.  There"
            " may be fewer than this in practice for several reasons"
            " such as duplicate reagent sets being removed or the"
            " final product not matching the query even though the"
            " synthons suggested they would.");
  }
};

SynthonSpaceSearch::SubstructureResults substructureSearch_helper(
    SynthonSpaceSearch::SynthonSpace &self, const ROMol &query,
    const python::object &py_params) {
  SynthonSpaceSearch::SynthonSpaceSearchParams params;
  if (!py_params.is_none()) {
    params = python::extract<SynthonSpaceSearch::SynthonSpaceSearchParams>(
        py_params);
  }
  auto results = self.substructureSearch(query, params);
  return results;
}

void summariseHelper(SynthonSpaceSearch::SynthonSpace &self) {
  self.summarise(std::cout);
}

BOOST_PYTHON_MODULE(rdSynthonSpaceSearch) {
  python::scope().attr("__doc__") =
      "Module containing implementation of SynthonSpace search of"
      " Synthon-based chemical libraries such as Enamine REAL.";

  SubstructureResults_wrapper::wrap();

  std::string docString = "SynthonSpaceSearch parameters.";
  python::class_<RDKit::SynthonSpaceSearch::SynthonSpaceSearchParams,
                 boost::noncopyable>("SynthonSpaceSearchParams",
                                     docString.c_str())
      .def_readwrite(
          "maxBondSplits",
          &RDKit::SynthonSpaceSearch::SynthonSpaceSearchParams::maxBondSplits,
          "The maximum number of bonds to break in the query."
          "  It should be between 1 and 4 and will be constrained to be so."
          "  Default=4.")
      .def_readwrite(
          "maxHits",
          &RDKit::SynthonSpaceSearch::SynthonSpaceSearchParams::maxHits,
          "The maximum number of hits to return.  Default=1000."
          "Use -1 for no maximum.")
      .def_readwrite(
          "hitStart",
          &RDKit::SynthonSpaceSearch::SynthonSpaceSearchParams::hitStart,
          "The sequence number of the hit to start from.  So that you"
          " can return the next N hits of a search having already"
          " obtained N-1.  Default=0")
      .def_readwrite(
          "randomSample",
          &RDKit::SynthonSpaceSearch::SynthonSpaceSearchParams::randomSample,
          "If True, returns a random sample of the hits, up to maxHits"
          " in number.  Default=False.")
      .def_readwrite(
          "randomSeed",
          &RDKit::SynthonSpaceSearch::SynthonSpaceSearchParams::randomSeed,
          "If using randomSample, this seeds the random number"
          " generator so as to give reproducible results.  Default=-1"
          " means use a random seed.")
      .def_readwrite(
          "buildHits",
          &RDKit::SynthonSpaceSearch::SynthonSpaceSearchParams::buildHits,
          "If false, reports the maximum number of hits that"
          " the search could produce, but doesn't return them.");

  docString = "SynthonSpaceSearch object.";
  python::class_<RDKit::SynthonSpaceSearch::SynthonSpace, boost::noncopyable>(
      "SynthonSpace", docString.c_str(), python::init<>())
      .def("ReadTextFile",
           &RDKit::SynthonSpaceSearch::SynthonSpace::readTextFile,
           (python::arg("self"), python::arg("inFile")),
           "Reads text file of the sort used by ChemSpace/Enamine.")
      .def("ReadDBFile", &RDKit::SynthonSpaceSearch::SynthonSpace::readDBFile,
           (python::arg("self"), python::arg("inFile")),
           "Reads binary database file.")
      .def("WriteDBFile", &RDKit::SynthonSpaceSearch::SynthonSpace::writeDBFile,
           (python::arg("self"), python::arg("outFile")),
           "Writes binary database file.")
      .def("GetNumReactions",
           &RDKit::SynthonSpaceSearch::SynthonSpace::getNumReactions,
           (python::arg("self")),
           "Returns number of reactions in the SynthonSpace.")
      .def("GetNumProducts",
           &RDKit::SynthonSpaceSearch::SynthonSpace::getNumProducts,
           (python::arg("self")),
           "Returns number of products in the SynthonSpace, with multiple"
           " counting of any duplicates.")
      .def("Summarise", &RDKit::summariseHelper, (python::arg("self")),
           "Writes a summary of the SynthonSpace to stdout.")
      .def("SubstructureSearch", &RDKit::substructureSearch_helper,
           (python::arg("self"), python::arg("query"),
            python::arg("params") = python::object()),
           "Does a substructure search in the SynthonSpace.");
}

}  // namespace RDKit