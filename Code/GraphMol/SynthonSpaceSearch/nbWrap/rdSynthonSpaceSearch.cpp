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
#include <utility>

#include <nanobind/nanobind.h>
#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/unique_ptr.h>
#include <nanobind/stl/vector.h>

#include <RDBoost/Wrap_nb.h>

#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/EnumerateStereoisomers/EnumerateStereoisomers.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/GaussianShape/ShapeOverlayOptions.h>
#include <GraphMol/GeneralizedSubstruct/XQMol.h>
#include <GraphMol/RascalMCES/RascalOptions.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;

namespace {

void throwCancelled(const char *msg) {
  PyErr_SetString(PyExc_KeyboardInterrupt, msg);
  throw nb::python_error();
}

nb::list hitMolecules_helper(const SynthonSpaceSearch::SearchResults &res) {
  nb::list pyres;
  for (const auto &r : res.getHitMolecules()) {
    pyres.append(nb::cast(new ROMol(*r), nb::rv_policy::take_ownership));
  }
  return pyres;
}

nb::object bestHitMolecule_helper(
    const SynthonSpaceSearch::SearchResults &res) {
  const auto &bh = res.getBestHit();
  if (bh) {
    return nb::cast(new ROMol(*bh), nb::rv_policy::take_ownership);
  }
  return nb::none();
}

nb::object get_excludedVolume(
    const SynthonSpaceSearch::SynthonSpaceSearchParams &params) {
  if (!params.excludedVolume) {
    return nb::none();
  }
  return nb::cast(params.excludedVolume, nb::rv_policy::reference);
}

void set_excludedVolume(SynthonSpaceSearch::SynthonSpaceSearchParams &params,
                        const nb::object &pyExcVol) {
  if (pyExcVol.is_none()) {
    params.excludedVolume = nullptr;
    return;
  }
  params.excludedVolume = new GaussianShape::ShapeInput(
      nb::cast<GaussianShape::ShapeInput>(pyExcVol));
}

SynthonSpaceSearch::SearchResults substructureSearch_helper1(
    SynthonSpaceSearch::SynthonSpace &self, const ROMol &query,
    const std::optional<SubstructMatchParameters> &py_smParams,
    const std::optional<SynthonSpaceSearch::SynthonSpaceSearchParams>
        &py_params) {
  SynthonSpaceSearch::SearchResults results;
  {
    NOGIL gil;
    results = self.substructureSearch(
        query, py_smParams.value_or(SubstructMatchParameters()),
        py_params.value_or(SynthonSpaceSearch::SynthonSpaceSearchParams()));
  }
  if (results.getCancelled()) {
    throwCancelled("SubstructureSearch cancelled");
  }
  return results;
}

SynthonSpaceSearch::SearchResults substructureSearch_helper2(
    SynthonSpaceSearch::SynthonSpace &self,
    const GeneralizedSubstruct::ExtendedQueryMol &query,
    const std::optional<SubstructMatchParameters> &py_smParams,
    const std::optional<SynthonSpaceSearch::SynthonSpaceSearchParams>
        &py_params) {
  SynthonSpaceSearch::SearchResults results;
  {
    NOGIL gil;
    results = self.substructureSearch(
        query, py_smParams.value_or(SubstructMatchParameters()),
        py_params.value_or(SynthonSpaceSearch::SynthonSpaceSearchParams()));
  }
  if (results.getCancelled()) {
    throwCancelled("SubstructureSearch cancelled");
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

void substructureSearch_helper3(
    SynthonSpaceSearch::SynthonSpace &self, const ROMol &query,
    nb::object py_callable,
    const std::optional<SubstructMatchParameters> &py_smParams,
    const std::optional<SynthonSpaceSearch::SynthonSpaceSearchParams>
        &py_params) {
  CallbackAdapter callback{py_callable};
  self.substructureSearch(
      query, callback, py_smParams.value_or(SubstructMatchParameters()),
      py_params.value_or(SynthonSpaceSearch::SynthonSpaceSearchParams()));
}

SynthonSpaceSearch::SearchResults substructureSearch_helper4(
    SynthonSpaceSearch::SynthonSpace &self, const ROMol &query,
    const std::optional<SubstructMatchParameters> &py_smParams,
    const std::optional<SynthonSpaceSearch::SynthonSpaceSearchParams>
        &py_params,
    std::uint64_t startLine, std::uint64_t finishLine) {
  SynthonSpaceSearch::SearchResults results;
  {
    NOGIL gil;
    results = self.substructureSearch(
        query, py_smParams.value_or(SubstructMatchParameters()),
        py_params.value_or(SynthonSpaceSearch::SynthonSpaceSearchParams()),
        startLine, finishLine);
  }
  if (results.getCancelled()) {
    throwCancelled("SubstructureSearch cancelled");
  }
  return results;
}

SynthonSpaceSearch::SearchResults substructureSearch_helper5(
    SynthonSpaceSearch::SynthonSpace &self,
    const GeneralizedSubstruct::ExtendedQueryMol &query,
    const std::optional<SubstructMatchParameters> &py_smParams,
    const std::optional<SynthonSpaceSearch::SynthonSpaceSearchParams>
        &py_params,
    std::uint64_t startLine, std::uint64_t finishLine) {
  SynthonSpaceSearch::SearchResults results;
  {
    NOGIL gil;
    results = self.substructureSearch(
        query, py_smParams.value_or(SubstructMatchParameters()),
        py_params.value_or(SynthonSpaceSearch::SynthonSpaceSearchParams()),
        startLine, finishLine);
  }
  if (results.getCancelled()) {
    throwCancelled("SubstructureSearch cancelled");
  }
  return results;
}

SynthonSpaceSearch::SearchResults fingerprintSearch_helper(
    SynthonSpaceSearch::SynthonSpace &self, const ROMol &query,
    const nb::object &fingerprintGenerator,
    const std::optional<SynthonSpaceSearch::SynthonSpaceSearchParams>
        &py_params) {
  SynthonSpaceSearch::SearchResults results;
  {
    NOGIL gil;
    const FingerprintGenerator<std::uint64_t> *fpGen =
        nb::cast<FingerprintGenerator<std::uint64_t> *>(fingerprintGenerator);
    results = self.fingerprintSearch(
        query, *fpGen,
        py_params.value_or(SynthonSpaceSearch::SynthonSpaceSearchParams()));
  }
  if (results.getCancelled()) {
    throwCancelled("FingerprintSearch cancelled");
  }
  return results;
}

void fingerprintSearch_helper2(
    SynthonSpaceSearch::SynthonSpace &self, const ROMol &query,
    const nb::object &fingerprintGenerator, nb::object py_callable,
    const std::optional<SynthonSpaceSearch::SynthonSpaceSearchParams>
        &py_params) {
  const FingerprintGenerator<std::uint64_t> *fpGen =
      nb::cast<FingerprintGenerator<std::uint64_t> *>(fingerprintGenerator);
  CallbackAdapter callback{py_callable};
  self.fingerprintSearch(
      query, *fpGen, callback,
      py_params.value_or(SynthonSpaceSearch::SynthonSpaceSearchParams()));
}

SynthonSpaceSearch::SearchResults fingerprintSearch_helper3(
    SynthonSpaceSearch::SynthonSpace &self, const ROMol &query,
    const nb::object &fingerprintGenerator,
    const std::optional<SynthonSpaceSearch::SynthonSpaceSearchParams>
        &py_params,
    std::uint64_t startLine, std::uint64_t finishLine) {
  const FingerprintGenerator<std::uint64_t> *fpGen =
      nb::cast<FingerprintGenerator<std::uint64_t> *>(fingerprintGenerator);
  SynthonSpaceSearch::SearchResults results;
  {
    NOGIL gil;
    results = self.fingerprintSearch(
        query, *fpGen,
        py_params.value_or(SynthonSpaceSearch::SynthonSpaceSearchParams()),
        startLine, finishLine);
  }
  if (results.getCancelled()) {
    throwCancelled("FingerprintSearch cancelled");
  }
  return results;
}

SynthonSpaceSearch::SearchResults rascalSearch_helper(
    SynthonSpaceSearch::SynthonSpace &self, const ROMol &query,
    const nb::object &py_rascalOptions,
    const std::optional<SynthonSpaceSearch::SynthonSpaceSearchParams>
        &py_params) {
  RascalMCES::RascalOptions rascalOptions =
      nb::cast<RascalMCES::RascalOptions>(py_rascalOptions);
  SynthonSpaceSearch::SearchResults results;
  {
    NOGIL gil;
    results = self.rascalSearch(
        query, rascalOptions,
        py_params.value_or(SynthonSpaceSearch::SynthonSpaceSearchParams()));
  }
  if (results.getCancelled()) {
    throwCancelled("RascalSearch cancelled");
  }
  return results;
}

void rascalSearch_helper2(
    SynthonSpaceSearch::SynthonSpace &self, const ROMol &query,
    const std::optional<RascalMCES::RascalOptions> &py_rascalOptions,
    nb::object py_callable,
    const std::optional<SynthonSpaceSearch::SynthonSpaceSearchParams>
        &py_params) {
  CallbackAdapter callback{py_callable};
  self.rascalSearch(
      query, py_rascalOptions.value_or(RascalMCES::RascalOptions()), callback,
      py_params.value_or(SynthonSpaceSearch::SynthonSpaceSearchParams()));
}

SynthonSpaceSearch::SearchResults rascalSearch_helper3(
    SynthonSpaceSearch::SynthonSpace &self, const ROMol &query,
    const nb::object &py_rascalOptions,
    const std::optional<SynthonSpaceSearch::SynthonSpaceSearchParams>
        &py_params,
    std::uint64_t startLine, std::uint64_t finishLine) {
  RascalMCES::RascalOptions rascalOptions =
      nb::cast<RascalMCES::RascalOptions>(py_rascalOptions);
  SynthonSpaceSearch::SearchResults results;
  {
    NOGIL gil;
    results = self.rascalSearch(
        query, rascalOptions,
        py_params.value_or(SynthonSpaceSearch::SynthonSpaceSearchParams()),
        startLine, finishLine);
  }
  if (results.getCancelled()) {
    throwCancelled("RascalSearch cancelled");
  }
  return results;
}

SynthonSpaceSearch::SearchResults shapeSearch_helper(
    SynthonSpaceSearch::SynthonSpace &self, const ROMol &query,
    const std::optional<SynthonSpaceSearch::SynthonSpaceSearchParams>
        &py_params) {
  SynthonSpaceSearch::SearchResults results;
  {
    NOGIL gil;
    results = self.shapeSearch(
        query,
        py_params.value_or(SynthonSpaceSearch::SynthonSpaceSearchParams()));
  }
  if (results.getCancelled()) {
    throwCancelled("ShapeSearch cancelled");
  }
  return results;
}

SynthonSpaceSearch::SearchResults shapeSearch_helper3(
    SynthonSpaceSearch::SynthonSpace &self, const ROMol &query,
    std::uint64_t startLine, std::uint64_t finishLine,
    const SynthonSpaceSearch::SynthonSpaceSearchParams &params) {
  SynthonSpaceSearch::SearchResults results;
  {
    NOGIL gil;
    results = self.shapeSearch(query, params, startLine, finishLine);
  }
  if (results.getCancelled()) {
    throwCancelled("ShapeSearch cancelled");
  }
  return results;
}

void summariseHelper(SynthonSpaceSearch::SynthonSpace &self) {
  self.summarise(std::cout);
}

void reportSynthonUsage_helper(const SynthonSpaceSearch::SynthonSpace &self) {
  self.reportSynthonUsage(std::cout);
}

void convertTextToDBFile_helper(
    const std::filesystem::path &inFilename,
    const std::filesystem::path &outFilename, nb::object fpGen,
    const std::optional<SynthonSpaceSearch::ShapeBuildParams> &py_shapeParams) {
  const FingerprintGenerator<std::uint64_t> *fpGenCpp = nullptr;
  if (!fpGen.is_none()) {
    fpGenCpp = nb::cast<FingerprintGenerator<std::uint64_t> *>(fpGen);
  }
  std::optional<SynthonSpaceSearch::ShapeBuildParams> shapeParams;
  if (py_shapeParams) {
    shapeParams = *py_shapeParams;
  }

  bool cancelled = false;
  SynthonSpaceSearch::convertTextToDBFile(
      inFilename.string(), outFilename.string(), cancelled, fpGenCpp,
      shapeParams ? &(*shapeParams) : nullptr);
  if (cancelled) {
    throwCancelled("Database conversion cancelled");
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
    throwCancelled("Database read cancelled");
  }
}

void buildShapes_helper(
    SynthonSpaceSearch::SynthonSpace &spc,
    const std::optional<SynthonSpaceSearch::ShapeBuildParams> &py_params) {
  bool cancelled = false;
  SynthonSpaceSearch::ShapeBuildParams params;
  if (py_params) {
    params = *py_params;
  }
  {
    NOGIL gil;
    spc.buildSynthonShapes(cancelled, params);
  }
  if (cancelled) {
    throwCancelled("Shape building cancelled");
  }
}

void setUserConfGen_helper(SynthonSpaceSearch::ShapeBuildParams &ps,
                           nb::object func) {
  ps.userConformerGenerator =
      [func = std::move(func)](
          const std::string &smiles,
          unsigned int numConformers) -> std::unique_ptr<RWMol> {
    nb::gil_scoped_acquire gil;
    auto res = func(smiles, numConformers);
    if (res.is_none()) {
      return nullptr;
    }
    auto *m = nb::cast<ROMol *>(res);
    if (!m) {
      return nullptr;
    }
    return std::make_unique<RWMol>(*m);
  };
}

void setUserConfGen_helper2(SynthonSpaceSearch::SynthonSpaceSearchParams &ps,
                            nb::object func) {
  ps.userConformerGenerator =
      [func = std::move(func)](
          const std::string &smiles,
          unsigned int numConformers) -> std::unique_ptr<RWMol> {
    nb::gil_scoped_acquire gil;
    auto res = func(smiles, numConformers);
    if (res.is_none()) {
      return nullptr;
    }
    auto *m = nb::cast<ROMol *>(res);
    if (!m) {
      return nullptr;
    }
    return std::make_unique<RWMol>(*m);
  };
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
           "Returns whether the search was cancelled or not.")
      .def("GetBestHit", &bestHitMolecule_helper,
           "Returns the best hit found in the similarity search, even"
           " when none were under the search threshold.  May be empty if not"
           " a similarity search or nothing came close to a match.");

  nb::class_<SynthonSpaceSearch::SynthonSpaceSearchParams>(
      m, "SynthonSpaceSearchParams", "SynthonSpaceSearch parameters.")
      .def(nb::init<>())
      .def_rw(
          "maxHits", &SynthonSpaceSearch::SynthonSpaceSearchParams::maxHits,
          R"DOC(The maximum number of hits to return.  Default=1000.Use -1 for no maximum.)DOC")
      .def_rw(
          "maxNumFrags",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::maxNumFragSets,
          R"DOC(The maximum number of fragments the query can be broken into.
  Big molecules will create huge numbers of fragments that may cause
excessive memory use.  If the number of fragments hits this number,
fragmentation stops and the search results will likely be incomplete.
  Default=100000.)DOC")
      .def_rw("hitStart",
              &SynthonSpaceSearch::SynthonSpaceSearchParams::hitStart,
              R"DOC(The sequence number of the hit to start from.  So that you
can return the next N hits of a search having already
obtained N-1.  Default=0)DOC")
      .def_rw("randomSample",
              &SynthonSpaceSearch::SynthonSpaceSearchParams::randomSample,
              R"DOC(If True, returns a random sample of the hits, up to maxHits
in number.  Default=False.)DOC")
      .def_rw("randomSeed",
              &SynthonSpaceSearch::SynthonSpaceSearchParams::randomSeed,
              R"DOC(If using randomSample, this seeds the random number
generator so as to give reproducible results.  Default=-1
means use a random seed.)DOC")
      .def_rw("buildHits",
              &SynthonSpaceSearch::SynthonSpaceSearchParams::buildHits,
              R"DOC(If false, reports the maximum number of hits that
the search could produce, but doesn't return them.)DOC")
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
      .def_rw("minHitHeavyAtoms",
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
      .def_rw("minHitChiralAtoms",
              &SynthonSpaceSearch::SynthonSpaceSearchParams::minHitChiralAtoms,
              "Minimum number of chiral atoms in a hit.  Default=0.")
      .def_rw(
          "maxHitChiralAtoms",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::maxHitChiralAtoms,
          "Maximum number of chiral atoms in a hit.  Default=-1 means no maximum.")
      .def_rw("numConformers",
              &SynthonSpaceSearch::SynthonSpaceSearchParams::numConformers,
              "When doing a shape search, the number of conformers to generate"
              " for molecules.  Default=100.")
      .def_rw("confRMSThreshold",
              &SynthonSpaceSearch::SynthonSpaceSearchParams::confRMSThreshold,
              "When doing a shape search, the RMS threshold to use when pruning"
              " conformers.  Default=1.0.")
      .def_rw(
          "shapeOverlayOptions",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::shapeOverlayOptions,
          "Options for the shape overlays.")
      .def_rw(
          "bestHit", &SynthonSpaceSearch::SynthonSpaceSearchParams::bestHit,
          "If True, when doing a shape search it will return the hit conformer"
          " with the best shape match to the query conformer.  If False, it just"
          " returns the first hit conformer that exceeds the similarity cutoff."
          "  The latter will be faster but the returned hit conformations are likely"
          " to be less relevant.")
      .def_rw("enumerateUnspecifiedStereo",
              &SynthonSpaceSearch::SynthonSpaceSearchParams::
                  enumerateUnspecifiedStereo,
              "When doing a shape search, if there is"
              " unspecified stereochemistry in either"
              " the query or potential hit, enumerate"
              " test all possibilities.  Default=False.")
      .def_rw("stereoEnumOpts",
              &SynthonSpaceSearch::SynthonSpaceSearchParams::stereoEnumOpts,
              "Options for stereoisomer enumeration.")
      .def_rw(
          "timeOut", &SynthonSpaceSearch::SynthonSpaceSearchParams::timeOut,
          R"DOC(Time limit for search, in seconds.  Default is 600s, 0 means no
timeout.  Requires an integer)DOC")
      .def_rw("toTryChunkSize",
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
      .def_rw(
          "useProgressBar",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::useProgressBar,
          "Makes a progress bar of given width.  The number given is the"
          " number of '*' characters in a full bar.  There will"
          " be about another 35 characters or so depending on the size of the"
          " job.  Default=0 means no bar.")
      .def_prop_rw(
          "excludedVolume", &get_excludedVolume, &set_excludedVolume,
          "  Add an excluded volume to use in the shape search.  The volume"
          " overlap and mean overlap over clashing atoms will be reported.")
      .def_rw("maxExcludedVolume",
              &SynthonSpaceSearch::SynthonSpaceSearchParams::maxExcludedVolume,
              "Maximum allowed excluded volume for a hit to be accepted."
              "  Default -1.0 means no maximum.")
      .def_rw(
          "maxMeanExcludedVolume",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::maxMeanExcludedVolume,
          "Maximum mean excluded volume for a hit to be accepted.  The"
          " mean is the total excluded volume divided by the number"
          " of clashing atoms (within 2 CARBON_RAD of an excluded volume"
          " atom).  To try and distinguish between a mild clash"
          " over the whole hit and a few atoms having a really bad clash.")
      .def_rw("possibleHitsFile",
              &SynthonSpaceSearch::SynthonSpaceSearchParams::possibleHitsFile,
              "Name of a file to save the possible hits to."
              " These are the combinations of synthons that"
              " might match the query but need building and"
              " final checking.  Each line has a"
              " space-separated list of the synthons and the"
              " hit's name.  The file will be emptied and"
              " re-filled if it already exists.")
      .def_rw(
          "maxPossibleHitsToWrite",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::maxPossibleHitsToWrite,
          "Maximum number of lines to write to possibleHitsFile.  When dealing"
          " with huge synthon spaces it's very easy to fill a disk.  Default=10M.")
      .def_rw("writePossibleHitsAndStop",
              &SynthonSpaceSearch::SynthonSpaceSearchParams::
                  writePossibleHitsAndStop,
              "If True, creates the possibleHitsFile and stops without doing"
              " the final building and checking.  Default is False.")
      .def(
          "setUserConformerGenerator", &setUserConfGen_helper2,
          nb::keep_alive<1, 2>(), "func"_a,
          R"DOC(Allows you to provide a function that will be called instead of the default
 conformer generator to generate conformers for the synthons.  The function should
 take a SMILES string and the maximum number of conformers to generated and
 return a molecule object.)DOC")
      .def("__setattr__", &safeSetattr);

  nb::class_<SynthonSpaceSearch::ShapeBuildParams>(
      m, "ShapeBuildParams",
      "Parameters for building shape objects for SynthonSpaceSearch.")
      .def(nb::init<>())
      .def_rw("numConfs", &SynthonSpaceSearch::ShapeBuildParams::numConfs,
              "Maximum number of conformers per synthon or query.  Default=10")
      .def_rw("rmsThreshold",
              &SynthonSpaceSearch::ShapeBuildParams::rmsThreshold,
              "RMS threshold to use when pruning conformations.  Default=1.0.")
      .def_rw(
          "shapeSimThreshold",
          &SynthonSpaceSearch::ShapeBuildParams::shapeSimThreshold,
          "When generating shapes, similarity threshold for pruning.  No 2 shapes"
          " for each synthon or query will be more similar than this threshold.  Default=1.9.")
      .def_rw(
          "numThreads", &SynthonSpaceSearch::ShapeBuildParams::numThreads,
          "The number of threads to use for shape building.  If > 0, will use that"
          " number.  If <= 0, will use the number of hardware"
          " threads plus this number.Default=1.")
      .def_rw(
          "randomSeed", &SynthonSpaceSearch::ShapeBuildParams::randomSeed,
          "Seed for random number generator.  Default=-1 means use system random seed.")
      .def_rw("stereoEnumOpts",
              &SynthonSpaceSearch::ShapeBuildParams::stereoEnumOpts,
              "Options for stereoisomer enumeration.")
      .def_rw(
          "useProgressBar",
          &SynthonSpaceSearch::ShapeBuildParams::useProgressBar,
          "Makes a progress bar of given width.  The number given is the"
          " number of '*' characters in a full bar.  There will"
          " be about another 35 characters or so depending on the size of the"
          " job.  Default=0 means no bar.")
      .def_rw("maxSynthonAtoms",
              &SynthonSpaceSearch::ShapeBuildParams::maxSynthonAtoms,
              "If >0, sets a maximum number of heavy atoms, excluding dummies,"
              " for synthon to have a shape made.  Default=0.")
      .def_rw(
          "maxEmbedAttempts",
          &SynthonSpaceSearch::ShapeBuildParams::maxEmbedAttempts,
          "Maximum number of attempts for embedding a single synthon.  Default=10.")
      .def_rw("timeOut", &SynthonSpaceSearch::ShapeBuildParams::timeOut,
              "Maximum time in seconds to spend on each synthon when generating"
              " conformers.  Default=600 means no timeout.")
      .def_rw(
          "interimFile", &SynthonSpaceSearch::ShapeBuildParams::interimFile,
          "Interim file to write the SynthonSpace to during shape generation.  In the"
          " event of a failure, a restart from this file will be possible.")
      .def_rw(
          "interimWrites", &SynthonSpaceSearch::ShapeBuildParams::interimWrites,
          "If an interim file has been given, every this many shapes write a"
          " new version of the file.  Default=1000.")
      .def(
          "setUserConformerGenerator", &setUserConfGen_helper,
          nb::keep_alive<1, 2>(), "func"_a,
          R"DOC(Allows you to provide a function that will be called instead of the default
 conformer generator to generate conformers for the synthons.  The function should
 take a SMILES string and the maximum number of conformers to generated and
 return a molecule object.)DOC")
      .def("__setattr__", &safeSetattr);

  nb::class_<SynthonSpaceSearch::SynthonSpace>(m, "SynthonSpace",
                                               "SynthonSpaceSearch object.")
      .def(nb::init<>())
      .def("ReadTextFile", &readTextFile_helper, "inFile"_a,
           "Reads text file of the sort used by ChemSpace/Enamine.")
      .def(
          "ReadDBFile",
          [](SynthonSpaceSearch::SynthonSpace &self,
             const std::filesystem::path &inFile,
             int numThreads) { self.readDBFile(inFile.string(), numThreads); },
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
      .def("GetNumSynthons", &SynthonSpaceSearch::SynthonSpace::getNumSynthons,
           "Returns number of synthons in the SynthonSpace.")
      .def(
          "GetNumSynthonsWithShapes",
          &SynthonSpaceSearch::SynthonSpace::getNumSynthonsWithShapes,
          "Returns the number of synthons in the SynthonSpace that have a shape.")
      .def("Summarise", &summariseHelper,
           "Writes a summary of the SynthonSpace to stdout.")
      .def(
          "ReportSynthonUsage", &reportSynthonUsage_helper,
          "Writes a summary of the synthon usage in the SynthonSpace to stdout.")
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
      .def(
          "SubstructureSearch", &substructureSearch_helper4, "query"_a,
          "substructMatchParams"_a = nb::none(), "params"_a = nb::none(),
          "startLine"_a, "finishLine"_a,
          R"DOC(Take the contents of params.possibleHitsFile, which is assumed to have
been written by an earlier search, and extract those that are indeed
hits.  It makes sense that params is the same as the one used to
generate the possible hits, but this is not essential.  You could search
at a higher similarity threshold than used to create the possible hits,
for example.)DOC")
      .def(
          "SubstructureSearch", &substructureSearch_helper5, "query"_a,
          "substructMatchParams"_a = nb::none(), "params"_a = nb::none(),
          "startLine"_a, "finishLine"_a,
          R"DOC(Take the contents of params.possibleHitsFile, which is assumed to have
been written by an earlier search, and extract those that are indeed
hits.  It makes sense that params is the same as the one used to
generate the possible hits, but this is not essential.  You could search
at a higher similarity threshold than used to create the possible hits,
for example.)DOC")
      .def(
          "SubstructureSearchIncremental", &substructureSearch_helper3,
          "query"_a, "callback"_a, "substructMatchParams"_a = nb::none(),
          "params"_a = nb::none(),
          "Does a substructure search in the SynthonSpace returning results in the callback.")
      .def("FingerprintSearch", &fingerprintSearch_helper, "query"_a,
           "fingerprintGenerator"_a, "params"_a = nb::none(),
           R"DOC(Does a fingerprint search in the SynthonSpace using the
FingerprintGenerator passed in.)DOC")
      .def(
          "FingerprintSearch", &fingerprintSearch_helper3, "query"_a,
          "fingerprintGenerator"_a, "params"_a, "startLine"_a, "finishLine"_a,
          R"DOC(Take the contents of params.possibleHitsFile, which is assumed to have
been written by an earlier search, and extract those that are indeed
hits.  It makes sense that params is the same as the one used to
generate the possible hits, but this is not essential.  You could search
at a higher similarity threshold than used to create the possible hits,
for example.
Duplicate SMILES strings produced by different reactions will
be returned.)DOC")
      .def("FingerprintSearchIncremental", &fingerprintSearch_helper2,
           "query"_a, "fingerprintGenerator"_a, "callback"_a,
           "params"_a = nb::none(),
           R"DOC(Does a fingerprint search in the SynthonSpace using the
FingerprintGenerator passed in, returning results the callback.)DOC")
      .def(
          "RascalSearch", &rascalSearch_helper, "query"_a, "rascalOptions"_a,
          "params"_a = nb::none(),
          R"DOC(Does a search using the Rascal similarity score.  The similarity
threshold used is provided by rascalOptions, and the one in
params is ignored.)DOC")
      .def(
          "RascalSearch", &rascalSearch_helper3, "query"_a, "rascalOptions"_a,
          "params"_a, "startLine"_a, "finishLine"_a,
          R"DOC(Take the contents of params.possibleHitsFile, which is assumed to have
been written by an earlier search, and extract those that are indeed
hits.  It makes sense that params is the same as the one used to
generate the possible hits, but this is not essential.  You could search
at a higher similarity threshold than used to create the possible hits,
for example.
Duplicate SMILES strings produced by different reactions will
be returned.)DOC")
      .def(
          "RascalSearchIncremental", &rascalSearch_helper2, "query"_a,
          "rascalOptions"_a, "callback"_a, "params"_a = nb::none(),
          R"DOC(Does a search using the Rascal similarity score.  The similarity
threshold used is provided by rascalOptions, and the one in
params is ignored.  Returns results iteratively in the callback.)DOC")
      .def(
          "ShapeSearch", &shapeSearch_helper, "query"_a,
          "params"_a = nb::none(),
          R"DOC(Perform a shape similarity search with the given query molecule
across the synthonspace library.  Duplicate SMILES strings produced by
different reactions will be returned.  Requires a query with at least
1 3D conformer.  Only the first conformer will be used in the search.)DOC")
      .def(
          "ShapeSearch", &shapeSearch_helper3, "query"_a, "startLine"_a,
          "finishLine"_a, "params"_a,
          R"DOC(Take the contents of params.possibleHitsFile, which is assumed to have
been written by an earlier search, and extract those that are indeed
hits.  It makes sense that params is the same as the one used to
generate the possible hits, but this is not essential.  You could search
at a higher similarity threshold than used to create the possible hits,
for example.
Duplicate SMILES strings produced by different reactions will
be returned.  Requires a query with at least 1 3D conformer.  Only
the first conformer will be used in the search.)DOC")
      .def(
          "BuildSynthonFingerprints",
          &SynthonSpaceSearch::SynthonSpace::buildSynthonFingerprints,
          "fingerprintGenerator"_a, "progressBarWidth"_a = 0,
          R"DOC(Build the synthon fingerprints ready for similarity searching.  This
is done automatically when the first similarity search is done, but if
converting a text file to binary format it might need to be done
explicitly.  If progressBarWidth is > 0, a progress bar of that width
plus about 35 characters is displayed.)DOC")
      .def(
          "BuildSynthonShapes", &buildShapes_helper, "py_params"_a = nb::none(),
          "Build shapes for the synthons.  The conformations are generated, pruned"
          " with the given threshold, which is passed directly to EmbedMultipleConfs.");

  m.def("ConvertTextToDBFile", &convertTextToDBFile_helper, "inFilename"_a,
        "outFilename"_a, "fpGen"_a = nb::none(),
        "py_shapeParams"_a = nb::none(),
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
