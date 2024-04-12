//
// Copyright (C) David Cosgrove 2023
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
#include <GraphMol/RascalMCES/RascalMCES.h>
#include <GraphMol/RascalMCES/RascalClusterOptions.h>
#include <GraphMol/RascalMCES/RascalOptions.h>
#include <GraphMol/RascalMCES/RascalResult.h>

namespace python = boost::python;

namespace {

python::list convertVecPairInt(const std::vector<std::pair<int, int>> &vec) {
  python::list pyres;
  for (const auto &p : vec) {
    python::tuple tup = python::make_tuple(p.first, p.second);
    pyres.append(tup);
  }
  return pyres;
}

python::list bondMatches(const RDKit::RascalMCES::RascalResult &self) {
  return convertVecPairInt(self.getBondMatches());
}
python::list atomMatches(const RDKit::RascalMCES::RascalResult &self) {
  return convertVecPairInt(self.getAtomMatches());
}

void largestFragmentOnly(RDKit::RascalMCES::RascalResult &self) {
  self.largestFragOnly();
}

struct RascalResult_wrapper {
  static void wrap() {
    std::string docString = "Used to return RASCAL MCES results.";
    python::class_<RDKit::RascalMCES::RascalResult>(
        "RascalResult", docString.c_str(), python::no_init)
        .def_readonly("smartsString",
                      &RDKit::RascalMCES::RascalResult::getSmarts,
                      "SMARTS string defining the MCES.")
        .def("bondMatches", bondMatches, python::args("self"),
             "A function returning a list of list "
             "of tuples, each inner list containing the matching bonds in the "
             "MCES as tuples of bond indices from mol1 and mol2")
        .def("atomMatches", atomMatches, python::args("self"),
             "Likewise for atoms.")
        .def(
            "largestFragmentOnly", largestFragmentOnly, python::args("self"),
            "Function that cuts the MCES down to the single largest frag.  This cannot be undone.")
        .def_readonly("similarity",
                      &RDKit::RascalMCES::RascalResult::getSimilarity,
                      "Johnson similarity between 2 molecules.")
        .def_readonly("numFragments",
                      &RDKit::RascalMCES::RascalResult::getNumFrags,
                      "Number of fragments in MCES.")
        .def_readonly("largestFragmentSize",
                      &RDKit::RascalMCES::RascalResult::getLargestFragSize,
                      "Number of atoms in largest fragment.")
        .def_readonly("tier1Sim", &RDKit::RascalMCES::RascalResult::getTier1Sim,
                      "The tier 1 similarity estimate.")
        .def_readonly("tier2Sim", &RDKit::RascalMCES::RascalResult::getTier2Sim,
                      "The tier 2 similarity estimate.")
        .def_readonly("timedOut", &RDKit::RascalMCES::RascalResult::getTimedOut,
                      "Whether it timed out.");
  }
};
}  // namespace

namespace RDKit {

python::list findMCESWrapper(const ROMol &mol1, const ROMol &mol2,
                             const python::object &py_opts) {
  RascalMCES::RascalOptions opts;
  if (!py_opts.is_none()) {
    opts = python::extract<RascalMCES::RascalOptions>(py_opts);
  }
  std::vector<RDKit::RascalMCES::RascalResult> results;
  {
    NOGIL gil;
    results = RascalMCES::rascalMCES(mol1, mol2, opts);
  }
  python::list pyres;
  for (auto &res : results) {
    pyres.append(res);
  }
  return pyres;
}

std::vector<std::shared_ptr<ROMol>> extractMols(python::object mols) {
  std::vector<std::shared_ptr<ROMol>> cmols;
  unsigned int nElems = python::extract<unsigned int>(mols.attr("__len__")());
  cmols.resize(nElems);
  for (unsigned int i = 0; i < nElems; ++i) {
    if (!mols[i]) {
      throw_value_error("molecule is None");
    }
    cmols[i] = python::extract<std::shared_ptr<ROMol>>(mols[i]);
  }
  return cmols;
}

python::list packOutputMols(
    const std::vector<std::vector<unsigned int>> &clusters) {
  python::list pyres;
  for (auto &clus : clusters) {
    python::list mols;
    for (auto &m : clus) {
      mols.append(m);
    }
    pyres.append(mols);
  }
  return pyres;
}

python::list rascalClusterWrapper(python::object mols,
                                  const python::object &py_opts) {
  RascalMCES::RascalClusterOptions opts;
  if (!py_opts.is_none()) {
    opts = python::extract<RascalMCES::RascalClusterOptions>(py_opts);
  }
  auto cmols = extractMols(mols);
  std::vector<RDKit::UINT_VECT> clusters;
  {
    NOGIL gil;
    clusters = RascalMCES::rascalCluster(cmols, opts);
  }
  return packOutputMols(clusters);
}

python::list rascalButinaClusterWrapper(python::object mols,
                                        const python::object &py_opts) {
  RascalMCES::RascalClusterOptions opts;
  if (!py_opts.is_none()) {
    opts = python::extract<RascalMCES::RascalClusterOptions>(py_opts);
  }
  auto cmols = extractMols(mols);
  std::vector<RDKit::UINT_VECT> clusters;
  {
    NOGIL gil;
    clusters = RascalMCES::rascalButinaCluster(cmols, opts);
  }
  return packOutputMols(clusters);
}

BOOST_PYTHON_MODULE(rdRascalMCES) {
  python::scope().attr("__doc__") =
      "Module containing implementation of RASCAL Maximum Common Edge Substructure algorithm.";
  RascalResult_wrapper::wrap();

  std::string docString = "RASCAL Options";
  python::class_<RDKit::RascalMCES::RascalOptions, boost::noncopyable>(
      "RascalOptions", docString.c_str())
      .def_readwrite(
          "similarityThreshold",
          &RDKit::RascalMCES::RascalOptions::similarityThreshold,
          "Threshold below which MCES won't be run.  Between 0.0 and 1.0, default=0.7.")
      .def_readwrite(
          "singleLargestFrag",
          &RDKit::RascalMCES::RascalOptions::singleLargestFrag,
          "Return the just single largest fragment of the MCES.  This is equivalent to running with allBestMCEs=True, finding the result with the largest largestFragmentSize, and calling its largestFragmentOnly method.")
      .def_readwrite(
          "completeAromaticRings",
          &RDKit::RascalMCES::RascalOptions::completeAromaticRings,
          "If True (default), partial aromatic rings won't be returned.")
      .def_readwrite("ringMatchesRingOnly",
                     &RDKit::RascalMCES::RascalOptions::ringMatchesRingOnly,
                     "If True (default), ring bonds won't match ring bonds.")
      .def_readwrite(
          "exactConnectionsMatch",
          &RDKit::RascalMCES::RascalOptions::exactConnectionsMatch,
          "If True (default is False), atoms will only match atoms if they have the same\n"
          " number of explicit connections.  E.g. the central atom of\n"
          " C(C)(C) won't match either atom in CC")
      .def_readwrite(
          "minFragSize", &RDKit::RascalMCES::RascalOptions::minFragSize,
          "Imposes a minimum on the number of atoms in a fragment that may be part of the MCES.  Default -1 means no minimum.")
      .def_readwrite(
          "maxFragSeparation",
          &RDKit::RascalMCES::RascalOptions::maxFragSeparation,
          "Maximum number of bonds between fragments in the MCES for both to be reported.  Default -1 means no maximum.  If exceeded, the smaller fragment will be removed.")
      .def_readwrite(
          "allBestMCESs", &RDKit::RascalMCES::RascalOptions::allBestMCESs,
          "If True, reports all MCESs found of the same maximum size.  Default False means just report the first found.")
      .def_readwrite(
          "returnEmptyMCES", &RDKit::RascalMCES::RascalOptions::returnEmptyMCES,
          "If the estimated similarity between the 2 molecules doesn't meet the similarityThreshold, no results are returned.  If you want to know what the"
          " estimates were, set this to True, and examine the tier1Sim and tier2Sim properties of the result then returned.")
      .def_readwrite(
          "timeout", &RDKit::RascalMCES::RascalOptions::timeout,
          "Maximum time (in seconds) to spend on an individual MCESs determination.  Default 60, -1 means no limit.")
      .def_readwrite(
          "maxBondMatchPairs",
          &RDKit::RascalMCES::RascalOptions::maxBondMatchPairs,
          "Too many matching bond (vertex) pairs can cause the process to run out of memory."
          "  The default of 1000 is fairly safe.  Increase with caution, as memory use increases"
          " with the square of this number.  ");

  docString =
      "Find one or more MCESs between the 2 molecules given.  Returns a list of "
      "RascalResult objects."
      "- mol1"
      "- mol2 The two molecules for which to find the MCES"
      "- opts Optional RascalOptions object changing the default run mode."
      "";
  python::def("FindMCES", &RDKit::findMCESWrapper,
              (python::arg("mol1"), python::arg("mol2"),
               python::arg("opts") = python::object()),
              docString.c_str());

  docString =
      "RASCAL Cluster Options.  Most of these pertain to RascalCluster calculations.  Only similarityCutoff is used by RascalButinaCluster.";
  python::class_<RDKit::RascalMCES::RascalClusterOptions, boost::noncopyable>(
      "RascalClusterOptions", docString.c_str())
      .def_readwrite(
          "similarityCutoff",
          &RDKit::RascalMCES::RascalClusterOptions::similarityCutoff,
          "Similarity cutoff for molecules to be in the same cluster.  Between 0.0 and 1.0, default=0.7.")
      .def_readwrite(
          "minFragSize", &RDKit::RascalMCES::RascalClusterOptions::minFragSize,
          "The minimum number of atoms in a fragment for it to be included in the MCES.  Default=3.")
      .def_readwrite(
          "maxNumFrags", &RDKit::RascalMCES::RascalClusterOptions::maxNumFrags,
          "The maximum number of fragments allowed in the MCES for each pair of molecules. Default=2.  So that the MCES"
          " isn't a lot of small fragments scattered around the molecules giving an inflated estimate of similarity.")
      .def_readwrite(
          "numThreads", &RDKit::RascalMCES::RascalClusterOptions::numThreads,
          "Number of threads to use during clustering.  Default=-1 means all the hardware threads less one.")
      .def_readwrite(
          "a", &RDKit::RascalMCES::RascalClusterOptions::a,
          "The penalty score for each unconnected component in the MCES. Default=0.05.")
      .def_readwrite(
          "b", &RDKit::RascalMCES::RascalClusterOptions::a,
          "The weight of matched bonds over matched atoms. Default=2.")
      .def_readwrite(
          "minIntraClusterSim",
          &RDKit::RascalMCES::RascalClusterOptions::minIntraClusterSim,
          "Two pairs of molecules are included in the same cluster if the similarity between"
          " their MCESs is greater than this.  Default=0.9.")
      .def_readwrite(
          "clusterMergeSim",
          &RDKit::RascalMCES::RascalClusterOptions::clusterMergeSim,
          "Two clusters are merged if the fraction of molecules they have in common is greater than this.  Default=0.6.");

  docString =
      "Use the RASCAL MCES similarity metric to do fuzzy clustering.  Returns a list of lists "
      "of molecules, each inner list being a cluster.  The last cluster is all the "
      "molecules that didn't fit into another cluster (the singletons)."
      "- mols List of molecules to be clustered"
      "- opts Optional RascalOptions object changing the default run mode."
      "";
  python::def("RascalCluster", &RDKit::rascalClusterWrapper,
              (python::arg("mols"), python::arg("opts") = python::object()),
              docString.c_str());
  docString =
      "Use the RASCAL MCES similarity metric to do Butina clustering"
      " (Butina JCICS 39 747-750 (1999)).  Returns a list of lists of molecules,"
      " each inner list being a cluster.  The last cluster is all the"
      " molecules that didn't fit into another cluster (the singletons)."
      "- mols List of molecules to be clustered"
      "- opts Optional RascalOptions object changing the default run mode."
      "";
  python::def("RascalButinaCluster", &RDKit::rascalButinaClusterWrapper,
              (python::arg("mols"), python::arg("opts") = python::object()),
              docString.c_str());
}

}  // namespace RDKit