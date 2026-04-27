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
#include <GraphMol/EnumerateStereoisomers/EnumerateStereoisomers.h>
#include <GraphMol/RascalMCES/RascalOptions.h>
#include <GraphMol/GaussianShape/ShapeOverlayOptions.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>

namespace python = boost::python;

namespace RDKit {
namespace helpers {

python::list hitMolecules_helper(const SynthonSpaceSearch::SearchResults &res) {
  python::list pyres;
  for (auto &r : res.getHitMolecules()) {
    pyres.append(boost::make_shared<ROMol>(*r));
  }
  return pyres;
}

boost::shared_ptr<ROMol> bestHitMolecule_helper(
    const SynthonSpaceSearch::SearchResults &res) {
  const auto &bh = res.getBestHit();
  if (bh) {
    return boost::make_shared<ROMol>(*bh);
  }
  return boost::shared_ptr<ROMol>();
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
             "Returns whether the search was cancelled or not.")
        .def("GetBestHit", &bestHitMolecule_helper,
             "Returns the best hit found in the similarity search, even"
             " when none were under the search threshold.  May be empty if not"
             " a similarity search or nothing came close to a match.");
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
    PyErr_SetString(PyExc_KeyboardInterrupt, "SubstructureSearch cancelled");
    python::throw_error_already_set();
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
    PyErr_SetString(PyExc_KeyboardInterrupt, "SubstructureSearch cancelled");
    python::throw_error_already_set();
  }
  return results;
}

struct CallbackAdapter {
  python::object py_callable;

  bool operator()(std::vector<std::unique_ptr<ROMol>> &results) const {
    python::list pyres;
    for (auto &mol : results) {
      pyres.append(boost::shared_ptr<ROMol>(mol.release()));
    }
    return static_cast<bool>(py_callable(pyres));
  }
};

static void substructureSearch_helper3(SynthonSpaceSearch::SynthonSpace &self,
                                       const ROMol &query,
                                       python::object py_callable,
                                       const python::object &py_smParams,
                                       const python::object &py_params) {
  SynthonSpaceSearch::SynthonSpaceSearchParams params;
  SubstructMatchParameters smParams;
  if (!py_smParams.is_none()) {
    smParams = python::extract<SubstructMatchParameters>(py_smParams);
  }
  if (!py_params.is_none()) {
    params = python::extract<SynthonSpaceSearch::SynthonSpaceSearchParams>(
        py_params);
  }
  CallbackAdapter callback{py_callable};
  self.substructureSearch(query, callback, smParams, params);
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
    PyErr_SetString(PyExc_KeyboardInterrupt, "FingerprintSearch cancelled");
    python::throw_error_already_set();
  }
  return results;
}

static void fingerprintSearch_helper_2(
    SynthonSpaceSearch::SynthonSpace &self, const ROMol &query,
    const python::object &fingerprintGenerator, python::object py_callable,
    const python::object &py_params) {
  SynthonSpaceSearch::SynthonSpaceSearchParams params;
  if (!py_params.is_none()) {
    params = python::extract<SynthonSpaceSearch::SynthonSpaceSearchParams>(
        py_params);
  }
  const FingerprintGenerator<std::uint64_t> *fpGen =
      python::extract<FingerprintGenerator<std::uint64_t> *>(
          fingerprintGenerator);

  CallbackAdapter callback{py_callable};
  self.fingerprintSearch(query, *fpGen, callback, params);
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
  SynthonSpaceSearch::SearchResults results;
  {
    NOGIL gil;
    results = self.rascalSearch(query, rascalOptions, params);
  }
  if (results.getCancelled()) {
    PyErr_SetString(PyExc_KeyboardInterrupt, "RascalSearch cancelled");
    python::throw_error_already_set();
  }
  return results;
}

SynthonSpaceSearch::SearchResults shapeSearch_helper(
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
    results = self.shapeSearch(query, params);
  }
  if (results.getCancelled()) {
    PyErr_SetString(PyExc_KeyboardInterrupt, "ShapeSearch cancelled");
    python::throw_error_already_set();
  }
  return results;
}

static void rascalSearch_helper_2(SynthonSpaceSearch::SynthonSpace &self,
                                  const ROMol &query,
                                  const python::object &py_rascalOptions,
                                  python::object py_callable,
                                  const python::object &py_params) {
  RascalMCES::RascalOptions rascalOptions;
  rascalOptions = python::extract<RascalMCES::RascalOptions>(py_rascalOptions);
  SynthonSpaceSearch::SynthonSpaceSearchParams params;
  if (!py_params.is_none()) {
    params = python::extract<SynthonSpaceSearch::SynthonSpaceSearchParams>(
        py_params);
  }
  CallbackAdapter callback{py_callable};
  self.rascalSearch(query, rascalOptions, callback, params);
}

void summariseHelper(SynthonSpaceSearch::SynthonSpace &self) {
  self.summarise(std::cout);
}

void reportSynthonUsage_helper(const SynthonSpaceSearch::SynthonSpace &self) {
  self.reportSynthonUsage(std::cout);
}

void convertTextToDBFile_helper(const std::string &inFilename,
                                const std::string &outFilename,
                                python::object fpGen,
                                python::object py_params) {
  const FingerprintGenerator<std::uint64_t> *fpGenCpp = nullptr;
  if (fpGen) {
    fpGenCpp = python::extract<FingerprintGenerator<std::uint64_t> *>(fpGen);
  }
  std::unique_ptr<SynthonSpaceSearch::ShapeBuildParams> shapeParams;
  if (!py_params.is_none()) {
    shapeParams.reset(new SynthonSpaceSearch::ShapeBuildParams());
    *shapeParams =
        python::extract<SynthonSpaceSearch::ShapeBuildParams>(py_params);
  }
  bool cancelled = false;
  SynthonSpaceSearch::convertTextToDBFile(inFilename, outFilename, cancelled,
                                          fpGenCpp, shapeParams.get());
  if (cancelled) {
    PyErr_SetString(PyExc_KeyboardInterrupt, "Database conversion cancelled");
    python::throw_error_already_set();
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
    PyErr_SetString(PyExc_KeyboardInterrupt, "Database read cancelled");
    python::throw_error_already_set();
  }
}

void buildShapes_helper(SynthonSpaceSearch::SynthonSpace &spc,
                        const python::object &py_params) {
  SynthonSpaceSearch::ShapeBuildParams params;
  if (!py_params.is_none()) {
    params = python::extract<SynthonSpaceSearch::ShapeBuildParams>(py_params);
  }
  bool cancelled = false;
  {
    NOGIL gil;
    spc.buildSynthonShapes(cancelled, params);
  }
  if (cancelled) {
    PyErr_SetString(PyExc_KeyboardInterrupt, "Shape building cancelled");
    python::throw_error_already_set();
  }
}

class pyUserConfGenFunctor {
 public:
  pyUserConfGenFunctor(python::object obj) : dp_obj(std::move(obj)) {}
  ~pyUserConfGenFunctor() = default;
  std::unique_ptr<RWMol> operator()(const std::string &smiles,
                                    unsigned int numConformers) {
    // grab the GIL
    PyGILStateHolder h;
    auto res = dp_obj(smiles, numConformers);
    ROMol *m = python::extract<ROMol *>(res);
    return std::make_unique<RWMol>(*m);
  }

 private:
  python::object dp_obj;
};

void setUserConfGen_helper(SynthonSpaceSearch::ShapeBuildParams &ps,
                           python::object func) {
  ps.userConformerGenerator = pyUserConfGenFunctor(func);
}

void setUserConfGen_helper2(SynthonSpaceSearch::SynthonSpaceSearchParams &ps,
                            python::object func) {
  ps.userConformerGenerator = pyUserConfGenFunctor(func);
}
}  // namespace helpers

BOOST_PYTHON_MODULE(rdSynthonSpaceSearch) {
  python::scope().attr("__doc__") =
      "Module containing implementation of SynthonSpace search of"
      " Synthon-based chemical libraries such as Enamine REAL."
      "  NOTE: This functionality is experimental and the API"
      " and/or results may change in future releases.";

  helpers::SearchResults_wrapper::wrap();

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
          "minHitHeavyAtoms",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::minHitHeavyAtoms,
          "Minimum number of heavy atoms in a hit.  Default=0.")
      .def_readwrite(
          "maxHitHeavyAtoms",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::maxHitHeavyAtoms,
          "Maximum number of heavy atoms in a hit.  Default=-1 means no maximum.")
      .def_readwrite("minHitMolWt",
                     &SynthonSpaceSearch::SynthonSpaceSearchParams::minHitMolWt,
                     "Minimum molecular weight for a hit.  Default=0.0.")
      .def_readwrite(
          "maxHitMolWt",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::maxHitMolWt,
          "Maximum molecular weight for a hit.  Default=0.0 mean no maximum.")
      .def_readwrite(
          "minHitChiralAtoms",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::minHitChiralAtoms,
          "Minimum number of chiral atoms in a hit.  Default=0.")
      .def_readwrite(
          "maxHitChiralAtoms",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::maxHitChiralAtoms,
          "Maximum number of chiral atoms in a hit.  Default=-1 means no maximum.")
      .def_readwrite(
          "numConformers",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::numConformers,
          "When doing a shape search, the number of conformers to generate"
          " for molecules.  Default=100.")
      .def_readwrite(
          "confRMSThreshold",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::confRMSThreshold,
          "When doing a shape search, the RMS threshold to use when pruning"
          " conformers.  Default=1.0.")
      .def_readwrite(
          "shapeOverlayOptions",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::shapeOverlayOptions,
          "Options for the shape overlays.")
      .def_readwrite(
          "bestHit", &SynthonSpaceSearch::SynthonSpaceSearchParams::bestHit,
          "If True, when doing a shape search it will return the hit conformer"
          " with the best shape match to the query conformer.  If False, it just"
          " returns the first hit conformer that exceeds the similarity cutoff."
          "  The latter will be faster but the returned hit conformations are likely"
          " to be less relevant.")
      .def_readwrite("enumerateUnspecifiedStereo",
                     &SynthonSpaceSearch::SynthonSpaceSearchParams::
                         enumerateUnspecifiedStereo,
                     "When doing a shape search, if there is"
                     " unspecified stereochemistry in either"
                     " the query or potential hit, enumerate"
                     " test all possibilities.  Default=False.")
      .def_readwrite(
          "stereoEnumOpts",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::stereoEnumOpts,
          "Options for stereoisomer enumeration.")
      .def_readwrite(
          "timeOut", &SynthonSpaceSearch::SynthonSpaceSearchParams::timeOut,
          "Time limit for search, in seconds.  Default is 600s, 0 means no"
          " timeout.  Requires an integer")
      .def_readwrite(
          "toTryChunkSize",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::toTryChunkSize,
          "Process possible hits using the given chunk size")
      .def_readwrite(
          "numThreads",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::numThreads,
          "The number of threads to use for search.  If > 0, will use that"
          " number.  If <= 0, will use the number of hardware"
          " threads plus this number.  So if the number of"
          " hardware threads is 8, and numThreads is -1, it will"
          " use 7 threads.  Default=1.")
      .def_readwrite(
          "useProgressBar",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::useProgressBar,
          "Makes a progress bar of given width.  The number given is the"
          " number of '*' characters in a full bar.  There will"
          " be about another 35 characters or so depending on the size of the"
          " job.  Default=0 means no bar.")
      .def_readwrite(
          "excludedVolume",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::excludedVolume,
          "  An excluded volume to use in the shape search.  The volume"
          " overlap and mean overlap over clashing atoms will be reported.")
      .def_readwrite(
          "maxExcludedVolume",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::maxExcludedVolume,
          "Maximum allowed excluded volume for a hit to be accepted."
          "  Default -1.0 means no maximum.")
      .def_readwrite(
          "maxMeanExcludedVolume",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::maxMeanExcludedVolume,
          "Maximum mean excluded volume for a hit to be accepted.  The"
          " mean is the total excluded volume divided by the number"
          " of clashing atoms (within 2 CARBON_RAD of an excluded volume"
          " atom).  To try and distinguish between a mild clash"
          " over the whole hit and a few atoms having a really bad clash.")
      .def_readwrite(
          "possibleHitsFile",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::possibleHitsFile,
          "Name of a file to save the possible hits to."
          " These are the combinations of synthons that"
          " might match the query but need building and"
          " final checking.  Each line has a"
          " space-separated list of the synthons and the"
          " hit's name.  The file will be emptied and"
          " re-filled if it already exists.")
      .def_readwrite(
          "writePossibleHitsAndStop",
          &SynthonSpaceSearch::SynthonSpaceSearchParams::
              writePossibleHitsAndStop,
          "If True, creates the possibleHitsFile and stops without doing"
          " the final building and checking.  Default is False.")
      .def(
          "setUserConformerGenerator", helpers::setUserConfGen_helper2,
          python::with_custodian_and_ward<1, 2>(), python::args("self", "func"),
          R"DOC(Allows you to provide a function that will be called instead of the default
 conformer generator to generate conformers for the synthons.  The function should
 take a SMILES string and the maximum number of conformers to generated and
 return a molecule object.)DOC")
      .def("__setattr__", &safeSetattr);

  docString = "Parameters for building shape objects for SynthonSpaceSearch.";
  python::class_<SynthonSpaceSearch::ShapeBuildParams, boost::noncopyable>(
      "ShapeBuildParams", docString.c_str())
      .def_readwrite(
          "numConfs", &SynthonSpaceSearch::ShapeBuildParams::numConfs,
          "Maximum number of conformers per synthon or query.  Default=10")
      .def_readwrite(
          "rmsThreshold", &SynthonSpaceSearch::ShapeBuildParams::rmsThreshold,
          "RMS threshold to use when pruning conformations.  Default=1.0.")
      .def_readwrite(
          "shapeSimThreshold",
          &SynthonSpaceSearch::ShapeBuildParams::shapeSimThreshold,
          "When generating shapes, similarity threshold for pruning.  No 2 shapes"
          " for each synthon or query will be more similar than this threshold.  Default=1.9.")
      .def_readwrite(
          "numThreads", &SynthonSpaceSearch::ShapeBuildParams::numThreads,
          "The number of threads to use for shape building.  If > 0, will use that"
          " number.  If <= 0, will use the number of hardware"
          " threads plus this number.Default=1.")
      .def_readwrite(
          "randomSeed", &SynthonSpaceSearch::ShapeBuildParams::randomSeed,
          "Seed for random number generator.  Default=-1 means use system random seed.")
      .def_readwrite("stereoEnumOpts",
                     &SynthonSpaceSearch::ShapeBuildParams::stereoEnumOpts,
                     "Options for stereoisomer enumeration.")
      .def_readwrite(
          "useProgressBar",
          &SynthonSpaceSearch::ShapeBuildParams::useProgressBar,
          "Makes a progress bar of given width.  The number given is the"
          " number of '*' characters in a full bar.  There will"
          " be about another 35 characters or so depending on the size of the"
          " job.  Default=0 means no bar.")
      .def_readwrite(
          "maxSynthonAtoms",
          &SynthonSpaceSearch::ShapeBuildParams::maxSynthonAtoms,
          "If >0, sets a maximum number of heavy atoms, excluding dummies,"
          " for synthon to have a shape made.  Default=0.")
      .def_readwrite(
          "maxEmbedAttempts",
          &SynthonSpaceSearch::ShapeBuildParams::maxEmbedAttempts,
          "Maximum number of attempts for embedding a single synthon.  Default=10.")
      .def_readwrite(
          "timeOut", &SynthonSpaceSearch::ShapeBuildParams::timeOut,
          "Maximum time in seconds to spend on each synthon when generating"
          " conformers.  Default=600 means no timeout.")
      .def_readwrite(
          "interimFile", &SynthonSpaceSearch::ShapeBuildParams::interimFile,
          "Interim file to write the SynthonSpace to during shape generation.  In the"
          " event of a failure, a restart from this file will be possible.")
      .def_readwrite(
          "interimWrites", &SynthonSpaceSearch::ShapeBuildParams::interimWrites,
          "If an interim file has been given, every this many shapes write a"
          " new version of the file.  Default=1000.")
      .def(
          "setUserConformerGenerator", helpers::setUserConfGen_helper,
          python::with_custodian_and_ward<1, 2>(), python::args("self", "func"),
          R"DOC(Allows you to provide a function that will be called instead of the default
 conformer generator to generate conformers for the synthons.  The function should
 take a SMILES string and the maximum number of conformers to generated and
 return a molecule object.)DOC")
      .def("__setattr__", &safeSetattr);

  docString = "SynthonSpaceSearch object.";
  python::class_<SynthonSpaceSearch::SynthonSpace, boost::noncopyable>(
      "SynthonSpace", docString.c_str(), python::init<>())
      .def("ReadTextFile", &helpers::readTextFile_helper,
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
      .def("GetNumSynthons", &SynthonSpaceSearch::SynthonSpace::getNumSynthons,
           python::arg("self"),
           "Returns number of synthons in the SynthonSpace.")
      .def("Summarise", &helpers::summariseHelper, python::arg("self"),
           "Writes a summary of the SynthonSpace to stdout.")
      .def(
          "ReportSynthonUsage", &helpers::reportSynthonUsage_helper,
          python::arg("self"),
          "Writes a summary of the synthon usage in the SynthonSpace to stdout.")
      .def("GetSynthonFingerprintType",
           &SynthonSpaceSearch::SynthonSpace::getSynthonFingerprintType,
           python::arg("self"),
           "Returns the information string for the fingerprint generator"
           " used to create this space.")
      .def("SubstructureSearch", &helpers::substructureSearch_helper1,
           (python::arg("self"), python::arg("query"),
            python::arg("substructMatchParams") = python::object(),
            python::arg("params") = python::object()),
           "Does a substructure search in the SynthonSpace.")
      .def("SubstructureSearch", &helpers::substructureSearch_helper2,
           (python::arg("self"), python::arg("query"),
            python::arg("substructMatchParams") = python::object(),
            python::arg("params") = python::object()),
           "Does a substructure search in the SynthonSpace using an"
           " extended query.")
      .def(
          "SubstructureSearchIncremental", &helpers::substructureSearch_helper3,
          (python::arg("self"), python::arg("query"), python::arg("callback"),
           python::arg("substructMatchParams") = python::object(),
           python::arg("params") = python::object()),
          "Does a substructure search in the SynthonSpace returning results in the callback.")
      .def("FingerprintSearch", &helpers::fingerprintSearch_helper,
           (python::arg("self"), python::arg("query"),
            python::arg("fingerprintGenerator"),
            python::arg("params") = python::object()),
           "Does a fingerprint search in the SynthonSpace using the"
           " FingerprintGenerator passed in.")
      .def("FingerprintSearchIncremental", &helpers::fingerprintSearch_helper_2,
           (python::arg("self"), python::arg("query"),
            python::arg("fingerprintGenerator"), python::arg("callback"),
            python::arg("params") = python::object()),
           "Does a fingerprint search in the SynthonSpace using the"
           " FingerprintGenerator passed in, returning results the callback.")
      .def("RascalSearch", &helpers::rascalSearch_helper,
           (python::arg("self"), python::arg("query"),
            python::arg("rascalOptions"),
            python::arg("params") = python::object()),
           "Does a search using the Rascal similarity score.  The similarity"
           " threshold used is provided by rascalOptions, and the one in"
           " params is ignored.")
      .def("RascalSearchIncremental", &helpers::rascalSearch_helper_2,
           (python::arg("self"), python::arg("query"),
            python::arg("rascalOptions"), python::arg("callback"),
            python::arg("params") = python::object()),
           "Does a search using the Rascal similarity score.  The similarity"
           " threshold used is provided by rascalOptions, and the one in"
           " params is ignored.  Returns results iteratively in the callback.")
      .def("ShapeSearch", &helpers::shapeSearch_helper,
           (python::arg("self"), python::arg("query"),
            python::arg("params") = python::object()),
           "Does a search using the pubchem-align3d shape similarity metric.")
      .def(
          "BuildSynthonFingerprints",
          &SynthonSpaceSearch::SynthonSpace::buildSynthonFingerprints,
          (python::arg("self"), python::arg("fingerprintGenerator"),
           python::arg("progressBarWidth") = 0),
          "Build the synthon fingerprints ready for similarity searching.  This"
          " is done automatically when the first similarity search is done, but if"
          " converting a text file to binary format it might need to be done"
          " explicitly.  If progressBarWidth is > 0, a progress bar of that width"
          " plus about 35 characters is displayed.")
      .def(
          "BuildSynthonShapes", &helpers::buildShapes_helper,
          (python::arg("self"), python::arg("py_params") = python::object()),
          "Build shapes for the synthons.  The conformations are generated, pruned"
          " with the given threshold, which is passed directly to EmbedMultipleConfs.");

  docString =
      "Convert the text file into the binary DB file in our format."
      "  Assumes that all synthons from a reaction are contiguous in the input file."
      "  This uses a lot less memory than using ReadTextFile() followed by"
      "  WriteDBFile()."
      "- inFilename the name of the text file"
      "- outFilename the name of the binary file"
      "- optional fingerprint generator";
  python::def("ConvertTextToDBFile", &helpers::convertTextToDBFile_helper,
              (python::arg("inFilename"), python::arg("outFilename"),
               python::arg("fpGen") = python::object(),
               python::arg("py_shapeParams") = python::object()),
              docString.c_str());

  docString =
      "Format an integer with spaces every 3 digits for ease of reading";
  python::def("FormattedIntegerString",
              &SynthonSpaceSearch::formattedIntegerString, python::arg("value"),
              docString.c_str());
}

}  // namespace RDKit