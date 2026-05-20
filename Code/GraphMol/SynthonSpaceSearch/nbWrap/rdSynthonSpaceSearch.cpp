//
// Copyright (C) 2024-2026 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <csignal>
#include <stdexcept>

#include <nanobind/nanobind.h>
#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/unique_ptr.h>

#include <RDBoost/Wrap_nb.h>

#include <GraphMol/ROMol.h>
#include <GraphMol/RascalMCES/RascalOptions.h>
#include <GraphMol/GeneralizedSubstruct/XQMol.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;

namespace {

nb::list hitMolecules_helper(const SynthonSpaceSearch::SearchResults &res) {
  nb::list pyres;
  for (const auto &r : res.getHitMolecules()) {
    pyres.append(nb::cast(new ROMol(*r), nb::rv_policy::take_ownership));
  }
  return pyres;
}

SynthonSpaceSearch::SearchResults substructureSearch_helper1(
    SynthonSpaceSearch::SynthonSpace &self, const ROMol &query,
    const nb::object &py_smParams, const nb::object &py_params) {
  SynthonSpaceSearch::SynthonSpaceSearchParams params;
  SubstructMatchParameters smParams;
  if (!py_smParams.is_none()) {
    smParams = nb::cast<SubstructMatchParameters>(py_smParams);
  }
  if (!py_params.is_none()) {
    params = nb::cast<SynthonSpaceSearch::SynthonSpaceSearchParams>(py_params);
  }

  SynthonSpaceSearch::SearchResults results;
  {
    NOGIL gil;
    results = self.substructureSearch(query, smParams, params);
  }
  if (results.getCancelled()) {
    throw std::runtime_error("SubstructureSearch cancelled");
  }
  return results;
}

SynthonSpaceSearch::SearchResults substructureSearch_helper2(
    SynthonSpaceSearch::SynthonSpace &self,
    const GeneralizedSubstruct::ExtendedQueryMol &query,
    const nb::object &py_smParams, const nb::object &py_params) {
  SubstructMatchParameters smParams;
  if (!py_smParams.is_none()) {
    smParams = nb::cast<SubstructMatchParameters>(py_smParams);
  }
  SynthonSpaceSearch::SynthonSpaceSearchParams params;
  if (!py_params.is_none()) {
    params = nb::cast<SynthonSpaceSearch::SynthonSpaceSearchParams>(py_params);
  }
  SynthonSpaceSearch::SearchResults results;
  {
    NOGIL gil;
    results = self.substructureSearch(query, smParams, params);
  }
  if (results.getCancelled()) {
    throw std::runtime_error("SubstructureSearch cancelled");
  }
  return results;
}

struct CallbackAdapter {
  nb::object py_callable;

  bool operator()(std::vector<std::unique_ptr<ROMol>> &results) const {
    nb::list pyres;
    for (auto &mol : results) {
      pyres.append(nb::cast(mol.release(), nb::rv_policy::take_ownership));
    }
    nb::object ret = py_callable(pyres);
    if (ret.is_none()) {
      return false;
    }
    return nb::cast<bool>(ret);
  }
};

void substructureSearch_helper3(SynthonSpaceSearch::SynthonSpace &self,
                                const ROMol &query, nb::object py_callable,
                                const nb::object &py_smParams,
                                const nb::object &py_params) {
  SynthonSpaceSearch::SynthonSpaceSearchParams params;
  SubstructMatchParameters smParams;
  if (!py_smParams.is_none()) {
    smParams = nb::cast<SubstructMatchParameters>(py_smParams);
  }
  if (!py_params.is_none()) {
    params = nb::cast<SynthonSpaceSearch::SynthonSpaceSearchParams>(py_params);
  }
  CallbackAdapter callback{py_callable};
  self.substructureSearch(query, callback, smParams, params);
}

SynthonSpaceSearch::SearchResults fingerprintSearch_helper(
    SynthonSpaceSearch::SynthonSpace &self, const ROMol &query,
    const nb::object &fingerprintGenerator, const nb::object &py_params) {
  SynthonSpaceSearch::SynthonSpaceSearchParams params;
  if (!py_params.is_none()) {
    params = nb::cast<SynthonSpaceSearch::SynthonSpaceSearchParams>(py_params);
  }
  SynthonSpaceSearch::SearchResults results;
  {
    NOGIL gil;
    const FingerprintGenerator<std::uint64_t> *fpGen =
        nb::cast<FingerprintGenerator<std::uint64_t> *>(fingerprintGenerator);
    results = self.fingerprintSearch(query, *fpGen, params);
  }
  if (results.getCancelled()) {
    throw std::runtime_error("FingerprintSearch cancelled");
  }
  return results;
}

void fingerprintSearch_helper2(SynthonSpaceSearch::SynthonSpace &self,
                               const ROMol &query,
                               const nb::object &fingerprintGenerator,
                               nb::object py_callable,
                               const nb::object &py_params) {
  SynthonSpaceSearch::SynthonSpaceSearchParams params;
  if (!py_params.is_none()) {
    params = nb::cast<SynthonSpaceSearch::SynthonSpaceSearchParams>(py_params);
  }
  const FingerprintGenerator<std::uint64_t> *fpGen =
      nb::cast<FingerprintGenerator<std::uint64_t> *>(fingerprintGenerator);
  CallbackAdapter callback{py_callable};
  self.fingerprintSearch(query, *fpGen, callback, params);
}

SynthonSpaceSearch::SearchResults rascalSearch_helper(
    SynthonSpaceSearch::SynthonSpace &self, const ROMol &query,
    const nb::object &py_rascalOptions, const nb::object &py_params) {
  RascalMCES::RascalOptions rascalOptions =
      nb::cast<RascalMCES::RascalOptions>(py_rascalOptions);
  SynthonSpaceSearch::SynthonSpaceSearchParams params;
  if (!py_params.is_none()) {
    params = nb::cast<SynthonSpaceSearch::SynthonSpaceSearchParams>(py_params);
  }
  {
    NOGIL gil;
    return self.rascalSearch(query, rascalOptions, params);
  }
}

void rascalSearch_helper2(SynthonSpaceSearch::SynthonSpace &self,
                          const ROMol &query,
                          const nb::object &py_rascalOptions,
                          nb::object py_callable,
                          const nb::object &py_params) {
  RascalMCES::RascalOptions rascalOptions =
      nb::cast<RascalMCES::RascalOptions>(py_rascalOptions);
  SynthonSpaceSearch::SynthonSpaceSearchParams params;
  if (!py_params.is_none()) {
    params = nb::cast<SynthonSpaceSearch::SynthonSpaceSearchParams>(py_params);
  }
  CallbackAdapter callback{py_callable};
  self.rascalSearch(query, rascalOptions, callback, params);
}

void summariseHelper(SynthonSpaceSearch::SynthonSpace &self) {
  self.summarise(std::cout);
}

void convertTextToDBFile_helper(const std::filesystem::path &inFilename,
                                const std::filesystem::path &outFilename,
                                nb::object fpGen) {
  const FingerprintGenerator<std::uint64_t> *fpGenCpp = nullptr;
  if (!fpGen.is_none()) {
    fpGenCpp =
        nb::cast<FingerprintGenerator<std::uint64_t> *>(fpGen);
  }
  bool cancelled = false;
  SynthonSpaceSearch::convertTextToDBFile(inFilename.string(),
                                          outFilename.string(), cancelled,
                                          fpGenCpp);
  if (cancelled) {
    throw std::runtime_error("Database conversion cancelled");
  }
}

void readTextFile_helper(SynthonSpaceSearch::SynthonSpace &self,
                         const std::filesystem::path &inFilename) {
  bool cancelled = false;
  {
    NOGIL gil;
    self.readTextFile(inFilename.string(), cancelled);
  }
  if (cancelled) {
    throw std::runtime_error("Database read cancelled.");
  }
}

}  // namespace

NB_MODULE(rdSynthonSpaceSearch, m) {
  m.doc() =
      R"DOC(Module containing implementation of SynthonSpace search of
Synthon-based chemical libraries such as Enamine REAL.
  NOTE: This functionality is experimental and the API
and/or results may change in future releases.)DOC";

  nb::class_<SynthonSpaceSearch::SearchResults>(
      m, "SubstructureResult",
      "Used to return results of SynthonSpace searches.")
      .def("GetHitMolecules", hitMolecules_helper,
           "A function returning hits from the search")
      .def("GetMaxNumResults",
           &SynthonSpaceSearch::SearchResults::getMaxNumResults,
           R"DOC(The upper bound on number of results possible.  There
may be fewer than this in practice for several reasons
such as duplicate reagent sets being removed or the
final product not matching the query even though the
synthons suggested they would.)DOC")
      .def("GetTimedOut", &SynthonSpaceSearch::SearchResults::getTimedOut,
           "Returns whether the search timed out or not.")
      .def("GetCancelled", &SynthonSpaceSearch::SearchResults::getCancelled,
           "Returns whether the search was cancelled or not.");

  nb::class_<SynthonSpaceSearch::SynthonSpaceSearchParams>(
      m, "SynthonSpaceSearchParams", "SynthonSpaceSearch parameters.")
      .def(nb::init<>())
      .def_rw("maxHits",
              &SynthonSpaceSearch::SynthonSpaceSearchParams::maxHits,
              R"DOC(The maximum number of hits to return.  Default=1000.Use -1 for no maximum.)DOC")
      .def_rw(
          "maxNumFrags",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::maxNumFragSets,
          R"DOC(The maximum number of fragments the query can be broken into.
  Big molecules will create huge numbers of fragments that may cause
excessive memory use.  If the number of fragments hits this number,
fragmentation stops and the search results will likely be incomplete.
  Default=100000.)DOC")
      .def_rw(
          "hitStart",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::hitStart,
          R"DOC(The sequence number of the hit to start from.  So that you
can return the next N hits of a search having already
obtained N-1.  Default=0)DOC")
      .def_rw(
          "randomSample",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::randomSample,
          R"DOC(If True, returns a random sample of the hits, up to maxHits
in number.  Default=False.)DOC")
      .def_rw(
          "randomSeed",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::randomSeed,
          R"DOC(If using randomSample, this seeds the random number
generator so as to give reproducible results.  Default=-1
means use a random seed.)DOC")
      .def_rw("buildHits",
              &SynthonSpaceSearch::SynthonSpaceSearchParams::buildHits,
              R"DOC(If false, reports the maximum number of hits that
the search could produce, but doesn't return them.)DOC")
      .def_rw(
          "numRandomSweeps",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::numRandomSweeps,
          R"DOC(The random sampling doesn't always produce the
required number of hits in 1 go.  This parameter
controls how many loops it makes to try and get
the hits before giving up.  Default=10.)DOC")
      .def_rw(
          "similarityCutoff",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::similarityCutoff,
          R"DOC(Similarity cutoff for returning hits by fingerprint similarity.
  At present the fp is hard-coded to be Morgan, bits, radius=2.
  Default=0.5.)DOC")
      .def_rw(
          "fragSimilarityAdjuster",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::fragSimilarityAdjuster,
          R"DOC(Similarities of fragments are generally low due to low bit
densities.  For the fragment matching, reduce the similarity cutoff
off by this amount.  Default=0.1.)DOC")
      .def_rw(
          "approxSimilarityAdjuster",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::
              approxSimilarityAdjuster,
          R"DOC(The fingerprint search uses an approximate similarity method
before building a product and doing a final check.  The
similarityCutoff is reduced by this value for the approximate
check.  A lower value will give faster run times at the
risk of missing some hits.  The value you use should have a
positive correlation with your FOMO.  The default of 0.1 is
appropriate for Morgan fingerprints.  With RDKit fingerprints,
0.05 is adequate, and higher than that has been seen to
produce long run times.)DOC")
      .def_rw(
          "minHitHeavyAtoms",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::minHitHeavyAtoms,
          "Minimum number of heavy atoms in a hit.  Default=0.")
      .def_rw(
          "maxHitHeavyAtoms",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::maxHitHeavyAtoms,
          "Maximum number of heavy atoms in a hit.  Default=-1 means no maximum.")
      .def_rw("minHitMolWt",
              &SynthonSpaceSearch::SynthonSpaceSearchParams::minHitMolWt,
              "Minimum molecular weight for a hit.  Default=0.0.")
      .def_rw(
          "maxHitMolWt",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::maxHitMolWt,
          "Maximum molecular weight for a hit.  Default=0.0 mean no maximum.")
      .def_rw(
          "minHitChiralAtoms",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::minHitChiralAtoms,
          "Minimum number of chiral atoms in a hit.  Default=0.")
      .def_rw(
          "maxHitChiralAtoms",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::maxHitChiralAtoms,
          "Maximum number of chiral atoms in a hit.  Default=-1 means no maximum.")
      .def_rw(
          "timeOut",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::timeOut,
          R"DOC(Time limit for search, in seconds.  Default is 600s, 0 means no
timeout.  Requires an integer)DOC")
      .def_rw(
          "toTryChunkSize",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::toTryChunkSize,
          "Process possible hits using the given chunk size")
      .def_rw(
          "numThreads",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::numThreads,
          R"DOC(The number of threads to use for search.  If > 0, will use that
number.  If <= 0, will use the number of hardware
threads plus this number.  So if the number of
hardware threads is 8, and numThreads is -1, it will
use 7 threads.  Default=1.)DOC")
      .def("__setattr__", &safeSetattr);

  nb::class_<SynthonSpaceSearch::SynthonSpace>(m, "SynthonSpace",
                                               "SynthonSpaceSearch object.")
      .def(nb::init<>())
      .def("ReadTextFile", &readTextFile_helper, "inFile"_a,
           "Reads text file of the sort used by ChemSpace/Enamine.")
      .def(
          "ReadDBFile",
          [](SynthonSpaceSearch::SynthonSpace &self,
             const std::filesystem::path &inFile, int numThreads) {
            self.readDBFile(inFile.string(), numThreads);
          },
          "inFile"_a, "numThreads"_a = 1,
          R"DOC(Reads binary database file.  Takes optional number of threads,default=1.)DOC")
      .def(
          "WriteDBFile",
          [](const SynthonSpaceSearch::SynthonSpace &self,
             const std::filesystem::path &outFile) {
            self.writeDBFile(outFile.string());
          },
          "outFile"_a, "Writes binary database file.")
      .def(
          "WriteEnumeratedFile",
          [](const SynthonSpaceSearch::SynthonSpace &self,
             const std::filesystem::path &outFile) {
            self.writeEnumeratedFile(outFile.string());
          },
          "outFile"_a, "Writes enumerated library to file.")
      .def("GetNumReactions",
           &SynthonSpaceSearch::SynthonSpace::getNumReactions,
           "Returns number of reactions in the SynthonSpace.")
      .def("GetNumProducts", &SynthonSpaceSearch::SynthonSpace::getNumProducts,
           R"DOC(Returns number of products in the SynthonSpace, with multiple
counting of any duplicates.)DOC")
      .def("Summarise", &summariseHelper,
           "Writes a summary of the SynthonSpace to stdout.")
      .def("GetSynthonFingerprintType",
           &SynthonSpaceSearch::SynthonSpace::getSynthonFingerprintType,
           R"DOC(Returns the information string for the fingerprint generator
used to create this space.)DOC")
      .def("SubstructureSearch", &substructureSearch_helper1, "query"_a,
           "substructMatchParams"_a = nb::none(), "params"_a = nb::none(),
           "Does a substructure search in the SynthonSpace.")
      .def("SubstructureSearch", &substructureSearch_helper2, "query"_a,
           "substructMatchParams"_a = nb::none(), "params"_a = nb::none(),
           R"DOC(Does a substructure search in the SynthonSpace using an
extended query.)DOC")
      .def("SubstructureSearchIncremental", &substructureSearch_helper3,
           "query"_a, "callback"_a, "substructMatchParams"_a = nb::none(),
           "params"_a = nb::none(),
           "Does a substructure search in the SynthonSpace returning results in the callback.")
      .def("FingerprintSearch", &fingerprintSearch_helper, "query"_a,
           "fingerprintGenerator"_a, "params"_a = nb::none(),
           R"DOC(Does a fingerprint search in the SynthonSpace using the
FingerprintGenerator passed in.)DOC")
      .def("FingerprintSearchIncremental", &fingerprintSearch_helper2,
           "query"_a, "fingerprintGenerator"_a, "callback"_a,
           "params"_a = nb::none(),
           R"DOC(Does a fingerprint search in the SynthonSpace using the
FingerprintGenerator passed in, returning results the callback.)DOC")
      .def("RascalSearch", &rascalSearch_helper, "query"_a, "rascalOptions"_a,
           "params"_a = nb::none(),
           R"DOC(Does a search using the Rascal similarity score.  The similarity
threshold used is provided by rascalOptions, and the one in
params is ignored.)DOC")
      .def("RascalSearchIncremental", &rascalSearch_helper2, "query"_a,
           "rascalOptions"_a, "callback"_a, "params"_a = nb::none(),
           R"DOC(Does a search using the Rascal similarity score.  The similarity
threshold used is provided by rascalOptions, and the one in
params is ignored.  Returns results iteratively in the callback.)DOC")
      .def(
          "BuildSynthonFingerprints",
          &SynthonSpaceSearch::SynthonSpace::buildSynthonFingerprints,
          "fingerprintGenerator"_a,
          R"DOC(Build the synthon fingerprints ready for similarity searching.  This
is done automatically when the first similarity search is done, but if
converting a text file to binary format it might need to be done
explicitly.)DOC");

  m.def("ConvertTextToDBFile", &convertTextToDBFile_helper, "inFilename"_a,
        "outFilename"_a, "fpGen"_a = nb::none(),
        R"DOC(Convert the text file into the binary DB file in our format.
  Assumes that all synthons from a reaction are contiguous in the input file.
  This uses a lot less memory than using ReadTextFile() followed by
  WriteDBFile().
- inFilename the name of the text file
- outFilename the name of the binary file
- optional fingerprint generator)DOC");

  m.def("FormattedIntegerString",
        &RDKit::SynthonSpaceSearch::formattedIntegerString, "value"_a,
        "Format an integer with spaces every 3 digits for ease of reading");
}
