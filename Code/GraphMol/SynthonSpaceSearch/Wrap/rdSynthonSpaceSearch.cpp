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
#include <GraphMol/RascalMCES/RascalOptions.h>
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

SynthonSpaceSearch::SearchResults substructureSearch_helper1(
    SynthonSpaceSearch::SynthonSpace &self, const ROMol &query,
    const python::object &py_smParams, const python::object &py_params) {
  SynthonSpaceSearch::SynthonSpaceSearchParams params;
  SubstructMatchParameters smParams;
  if (!py_smParams.is_none()) {
    smParams = python::extract<SubstructMatchParameters>(py_smParams);
  }
  if (!py_params.is_none()) {
    params = python::extract<SynthonSpaceSearch::SynthonSpaceSearchParams>(
        py_params);
  }

  SynthonSpaceSearch::SearchResults results;
  {
    NOGIL gil;
    results = self.substructureSearch(query, smParams, params);
  }
  if (results.getCancelled()) {
    throw_runtime_error("SubstructureSearch cancelled");
  }
  return results;
}

SynthonSpaceSearch::SearchResults substructureSearch_helper2(
    SynthonSpaceSearch::SynthonSpace &self,
    const GeneralizedSubstruct::ExtendedQueryMol &query,
    const python::object &py_smParams, const python::object &py_params) {
  SubstructMatchParameters smParams;
  if (!py_smParams.is_none()) {
    smParams = python::extract<SubstructMatchParameters>(py_smParams);
  }
  SynthonSpaceSearch::SynthonSpaceSearchParams params;
  if (!py_params.is_none()) {
    params = python::extract<SynthonSpaceSearch::SynthonSpaceSearchParams>(
        py_params);
  }
  SynthonSpaceSearch::SearchResults results;
  {
    NOGIL gil;
    results = self.substructureSearch(query, smParams, params);
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

SynthonSpaceSearch::SearchResults rascalSearch_helper(
    SynthonSpaceSearch::SynthonSpace &self, const ROMol &query,
    const python::object &py_rascalOptions, const python::object &py_params) {
  RascalMCES::RascalOptions rascalOptions;
  rascalOptions = python::extract<RascalMCES::RascalOptions>(py_rascalOptions);
  SynthonSpaceSearch::SynthonSpaceSearchParams params;
  if (!py_params.is_none()) {
    params = python::extract<SynthonSpaceSearch::SynthonSpaceSearchParams>(
        py_params);
  }
  {
    NOGIL gil;
    return self.rascalSearch(query, rascalOptions, params);
  }
}

void summariseHelper(SynthonSpaceSearch::SynthonSpace &self) {
  self.summarise(std::cout);
}

void convertTextToDBFile_helper(const std::string &inFilename,
                                const std::string &outFilename,
                                python::object fpGen) {
  const FingerprintGenerator<std::uint64_t> *fpGenCpp = nullptr;
  if (fpGen) {
    fpGenCpp = python::extract<FingerprintGenerator<std::uint64_t> *>(fpGen);
  }
  bool cancelled = false;
  SynthonSpaceSearch::convertTextToDBFile(inFilename, outFilename, cancelled,
                                          fpGenCpp);
  if (cancelled) {
    throw_runtime_error("Database conversion cancelled");
  }
}

void readTextFile_helper(SynthonSpaceSearch::SynthonSpace &self,
                         const std::string &inFilename) {
  bool cancelled = false;
  {
    NOGIL gil;
    self.readTextFile(inFilename, cancelled);
  }
  if (cancelled) {
    throw_runtime_error("Database read cancelled.");
  }
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
      .def_readwrite("maxHits",
                     &SynthonSpaceSearch::SynthonSpaceSearchParams::maxHits,
                     "The maximum number of hits to return.  Default=1000."
                     "Use -1 for no maximum.")
      .def_readwrite(
          "maxNumFrags",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::maxNumFragSets,
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
          " timeout.  Requires an integer")
      .def_readwrite(
          "numThreads",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::numThreads,
          "The number of threads to use for search.  If > 0, will use that"
          " number.  If <= 0, will use the number of hardware"
          " threads plus this number.  So if the number of"
          " hardware threads is 8, and numThreads is -1, it will"
          " use 7 threads.  Default=1.");

  docString = "SynthonSpaceSearch object.";
  python::class_<SynthonSpaceSearch::SynthonSpace, boost::noncopyable>(
      "SynthonSpace", docString.c_str(), python::init<>())
      .def("ReadTextFile", &readTextFile_helper,
           (python::arg("self"), python::arg("inFile")),
           "Reads text file of the sort used by ChemSpace/Enamine.")
      .def("ReadDBFile", &SynthonSpaceSearch::SynthonSpace::readDBFile,
           (python::arg("self"), python::arg("inFile"),
            python::arg("numThreads") = 1),
           "Reads binary database file.  Takes optional number of threads,"
           "default=1.")
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
      .def("GetSynthonFingerprintType",
           &SynthonSpaceSearch::SynthonSpace::getSynthonFingerprintType,
           python::arg("self"),
           "Returns the information string for the fingerprint generator"
           " used to create this space.")
      .def("SubstructureSearch", &substructureSearch_helper1,
           (python::arg("self"), python::arg("query"),
            python::arg("substructMatchParams") = python::object(),
            python::arg("params") = python::object()),
           "Does a substructure search in the SynthonSpace.")
      .def("SubstructureSearch", &substructureSearch_helper2,
           (python::arg("self"), python::arg("query"),
            python::arg("substructMatchParams") = python::object(),
            python::arg("params") = python::object()),
           "Does a substructure search in the SynthonSpace using an"
           " extended query.")
      .def("FingerprintSearch", &fingerprintSearch_helper,
           (python::arg("self"), python::arg("query"),
            python::arg("fingerprintGenerator"),
            python::arg("params") = python::object()),
           "Does a fingerprint search in the SynthonSpace using the"
           " FingerprintGenerator passed in.")
      .def("RascalSearch", &rascalSearch_helper,
           (python::arg("self"), python::arg("query"),
            python::arg("rascalOptions"),
            python::arg("params") = python::object()),
           "Does a search using the Rascal similarity score.  The similarity"
           " threshold used is provided by rascalOptions, and the one in"
           " params is ignored.")
      .def(
          "BuildSynthonFingerprints",
          &SynthonSpaceSearch::SynthonSpace::buildSynthonFingerprints,
          (python::arg("self"), python::arg("fingerprintGenerator")),
          "Build the synthon fingerprints ready for similarity searching.  This"
          " is done automatically when the first similarity search is done, but if"
          " converting a text file to binary format it might need to be done"
          " explicitly.");

  docString =
      "Convert the text file into the binary DB file in our format."
      "  Assumes that all synthons from a reaction are contiguous in the input file."
      "  This uses a lot less memory than using ReadTextFile() followed by"
      "  WriteDBFile()."
      "- inFilename the name of the text file"
      "- outFilename the name of the binary file"
      "- optional fingerprint generator";
  python::def("ConvertTextToDBFile", &RDKit::convertTextToDBFile_helper,
              (python::arg("inFilename"), python::arg("outFilename"),
               python::arg("fpGen") = python::object()),
              docString.c_str());

  docString =
      "Format an integer with spaces every 3 digits for ease of reading";
  python::def("FormattedIntegerString",
              &RDKit::SynthonSpaceSearch::formattedIntegerString,
              python::arg("value"), docString.c_str());
}

}  // namespace RDKit