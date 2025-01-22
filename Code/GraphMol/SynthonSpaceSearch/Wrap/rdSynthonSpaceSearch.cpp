//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <csignal>

#include <RDBoost/python.h>
#include <RDBoost/Wrap.h>

#include <GraphMol/ROMol.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>

namespace python = boost::python;

namespace RDKit {

python::list hitMolecules_helper(const SynthonSpaceSearch::SearchResults &res) {
  python::list pyres;
  for (auto &r : res.getHitMolecules()) {
    pyres.append(boost::make_shared<ROMol>(*r));
  }
  return pyres;
}

struct SearchResults_wrapper {
  static void wrap() {
    const std::string docString =
        "Used to return results of SynthonSpace searches.";
    python::class_<SynthonSpaceSearch::SearchResults>(
        "SubstructureResult", docString.c_str(), python::no_init)
        .def("GetHitMolecules", hitMolecules_helper, python::args("self"),
             "A function returning hits from the search")
        .def("GetMaxNumResults",
             &SynthonSpaceSearch::SearchResults::getMaxNumResults,
             "The upper bound on number of results possible.  There"
             " may be fewer than this in practice for several reasons"
             " such as duplicate reagent sets being removed or the"
             " final product not matching the query even though the"
             " synthons suggested they would.")
        .def("GetTimedOut", &SynthonSpaceSearch::SearchResults::getTimedOut,
             "Returns whether the search timed out or not.")
        .def("GetCancelled", &SynthonSpaceSearch::SearchResults::getCancelled,
             "Returns whether the search was cancelled or not.");
  }
};

SynthonSpaceSearch::SearchResults substructureSearch_helper(
    SynthonSpaceSearch::SynthonSpace &self, const ROMol &query,
    const python::object &py_params) {
  SynthonSpaceSearch::SynthonSpaceSearchParams params;
  if (!py_params.is_none()) {
    params = python::extract<SynthonSpaceSearch::SynthonSpaceSearchParams>(
        py_params);
  }
  SynthonSpaceSearch::SearchResults results;
  {
    NOGIL gil;
    results = self.substructureSearch(query, params);
  }
  if (results.getCancelled()) {
    throw_runtime_error("SubstructureSearch cancelled");
  }
  return results;
}

SynthonSpaceSearch::SearchResults fingerprintSearch_helper(
    SynthonSpaceSearch::SynthonSpace &self, const ROMol &query,
    const python::object &fingerprintGenerator,
    const python::object &py_params) {
  SynthonSpaceSearch::SynthonSpaceSearchParams params;
  if (!py_params.is_none()) {
    params = python::extract<SynthonSpaceSearch::SynthonSpaceSearchParams>(
        py_params);
  }
  SynthonSpaceSearch::SearchResults results;
  {
    NOGIL gil;
    const FingerprintGenerator<std::uint64_t> *fpGen =
        python::extract<FingerprintGenerator<std::uint64_t> *>(
            fingerprintGenerator);
    results = self.fingerprintSearch(query, *fpGen, params);
  }
  if (results.getCancelled()) {
    throw_runtime_error("FingerprintSearch cancelled");
  }
  return results;
}

void summariseHelper(const SynthonSpaceSearch::SynthonSpace &self) {
  self.summarise(std::cout);
}

BOOST_PYTHON_MODULE(rdSynthonSpaceSearch) {
  python::scope().attr("__doc__") =
      "Module containing implementation of SynthonSpace search of"
      " Synthon-based chemical libraries such as Enamine REAL."
      "  NOTE: This functionality is experimental and the API"
      " and/or results may change in future releases.";

  SearchResults_wrapper::wrap();

  std::string docString = "SynthonSpaceSearch parameters.";
  python::class_<SynthonSpaceSearch::SynthonSpaceSearchParams,
                 boost::noncopyable>("SynthonSpaceSearchParams",
                                     docString.c_str())
      .def_readwrite(
          "maxBondSplits",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::maxBondSplits,
          "The maximum number of bonds to break in the query."
          "  It should be between 1 and 4 and will be constrained to be so."
          "  Default=4.")
      .def_readwrite("maxHits",
                     &SynthonSpaceSearch::SynthonSpaceSearchParams::maxHits,
                     "The maximum number of hits to return.  Default=1000."
                     "Use -1 for no maximum.")
      .def_readwrite(
          "maxNumFrags",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::maxNumFrags,
          "The maximum number of fragments the query can be broken into."
          "  Big molecules will create huge numbers of fragments that may cause"
          " excessive memory use.  If the number of fragments hits this number,"
          " fragmentation stops and the search results will likely be incomplete."
          "  Default=100000.")
      .def_readwrite(
          "hitStart", &SynthonSpaceSearch::SynthonSpaceSearchParams::hitStart,
          "The sequence number of the hit to start from.  So that you"
          " can return the next N hits of a search having already"
          " obtained N-1.  Default=0")
      .def_readwrite(
          "randomSample",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::randomSample,
          "If True, returns a random sample of the hits, up to maxHits"
          " in number.  Default=False.")
      .def_readwrite(
          "randomSeed",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::randomSeed,
          "If using randomSample, this seeds the random number"
          " generator so as to give reproducible results.  Default=-1"
          " means use a random seed.")
      .def_readwrite("buildHits",
                     &SynthonSpaceSearch::SynthonSpaceSearchParams::buildHits,
                     "If false, reports the maximum number of hits that"
                     " the search could produce, but doesn't return them.")
      .def_readwrite(
          "numRandomSweeps",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::numRandomSweeps,
          "The random sampling doesn't always produce the"
          " required number of hits in 1 go.  This parameter"
          " controls how many loops it makes to try and get"
          " the hits before giving up.  Default=10.")
      .def_readwrite(
          "similarityCutoff",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::similarityCutoff,
          "Similarity cutoff for returning hits by fingerprint similarity."
          "  At present the fp is hard-coded to be Morgan, bits, radius=2."
          "  Default=0.5.")
      .def_readwrite(
          "fragSimilarityAdjuster",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::fragSimilarityAdjuster,
          "Similarities of fragments are generally low due to low bit"
          " densities.  For the fragment matching, reduce the similarity cutoff"
          " off by this amount.  Default=0.1.")
      .def_readwrite(
          "approxSimilarityAdjuster",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::
              approxSimilarityAdjuster,
          "The fingerprint search uses an approximate similarity method"
          " before building a product and doing a final check.  The"
          " similarityCutoff is reduced by this value for the approximate"
          " check.  A lower value will give faster run times at the"
          " risk of missing some hits.  The value you use should have a"
          " positive correlation with your FOMO.  The default of 0.1 is"
          " appropriate for Morgan fingerprints.  With RDKit fingerprints,"
          " 0.05 is adequate, and higher than that has been seen to"
          " produce long run times.")
      .def_readwrite(
          "timeOut", &SynthonSpaceSearch::SynthonSpaceSearchParams::timeOut,
          "Time limit for search, in seconds.  Default is 600s, 0 means no"
          " timeout.  Requires an integer");

  docString = "SynthonSpaceSearch object.";
  python::class_<SynthonSpaceSearch::SynthonSpace, boost::noncopyable>(
      "SynthonSpace", docString.c_str(), python::init<>())
      .def("ReadTextFile", &SynthonSpaceSearch::SynthonSpace::readTextFile,
           (python::arg("self"), python::arg("inFile")),
           "Reads text file of the sort used by ChemSpace/Enamine.")
      .def("ReadDBFile", &SynthonSpaceSearch::SynthonSpace::readDBFile,
           (python::arg("self"), python::arg("inFile")),
           "Reads binary database file.")
      .def("WriteDBFile", &SynthonSpaceSearch::SynthonSpace::writeDBFile,
           (python::arg("self"), python::arg("outFile")),
           "Writes binary database file.")
      .def("WriteEnumeratedFile",
           &SynthonSpaceSearch::SynthonSpace::writeEnumeratedFile,
           (python::arg("self"), python::arg("outFile")),
           "Writes enumerated library to file.")
      .def("GetNumReactions",
           &SynthonSpaceSearch::SynthonSpace::getNumReactions,
           python::arg("self"),
           "Returns number of reactions in the SynthonSpace.")
      .def("GetNumProducts", &SynthonSpaceSearch::SynthonSpace::getNumProducts,
           python::arg("self"),
           "Returns number of products in the SynthonSpace, with multiple"
           " counting of any duplicates.")
      .def("Summarise", &summariseHelper, python::arg("self"),
           "Writes a summary of the SynthonSpace to stdout.")
      .def("SubstructureSearch", &substructureSearch_helper,
           (python::arg("self"), python::arg("query"),
            python::arg("params") = python::object()),
           "Does a substructure search in the SynthonSpace.")
      .def("FingerprintSearch", &fingerprintSearch_helper,
           (python::arg("self"), python::arg("query"),
            python::arg("fingerprintGenerator"),
            python::arg("params") = python::object()),
           "Does a fingerprint search in the SynthonSpace using the"
           " FingerprintGenerator passed in.")
      .def(
          "BuildSynthonFingerprints",
          &SynthonSpaceSearch::SynthonSpace::buildSynthonFingerprints,
          (python::arg("self"), python::arg("fingerprintGenerator")),
          "Build the synthon fingerprints ready for similarity searching.  This"
          " is done automatically when the first similarity search is done, but if"
          " converting a text file to binary format it might need to be done"
          " explicitly.");
}

}  // namespace RDKit