//
// Copyright (C) 2023-2026 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/shared_ptr.h>

#include <RDBoost/Wrap_nb.h>

#include <GraphMol/ROMol.h>
#include <GraphMol/RascalMCES/RascalMCES.h>
#include <GraphMol/RascalMCES/RascalClusterOptions.h>
#include <GraphMol/RascalMCES/RascalOptions.h>
#include <GraphMol/RascalMCES/RascalResult.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace {

nb::list convertVecPairInt(const std::vector<std::pair<int, int>> &vec) {
  nb::list pyres;
  for (const auto &p : vec) {
    pyres.append(nb::make_tuple(p.first, p.second));
  }
  return pyres;
}

std::vector<std::shared_ptr<RDKit::ROMol>> extractMols(nb::object mols) {
  std::vector<std::shared_ptr<RDKit::ROMol>> cmols;
  unsigned int nElems = nb::len(mols);
  cmols.resize(nElems);
  for (unsigned int i = 0; i < nElems; ++i) {
    nb::object mol = mols[i];
    if (mol.is_none()) {
      throw nb::value_error("molecule is None");
    }
    cmols[i] = nb::cast<std::shared_ptr<RDKit::ROMol>>(mol);
  }
  return cmols;
}

nb::list packOutputMols(
    const std::vector<std::vector<unsigned int>> &clusters) {
  nb::list pyres;
  for (const auto &clus : clusters) {
    nb::list mols;
    for (auto m : clus) {
      mols.append(m);
    }
    pyres.append(mols);
  }
  return pyres;
}

}  // namespace

NB_MODULE(rdRascalMCES, m) {
  m.doc() = "Module containing implementation of RASCAL Maximum Common Edge Substructure algorithm.";

  nb::class_<RDKit::RascalMCES::RascalResult>(m, "RascalResult",
      "Used to return RASCAL MCES results.")
      .def_prop_ro("smartsString",
                   &RDKit::RascalMCES::RascalResult::getSmarts,
                   "SMARTS string defining the MCES.")
      .def("bondMatches",
           [](const RDKit::RascalMCES::RascalResult &self) {
             return convertVecPairInt(self.getBondMatches());
           },
           R"DOC(A function returning a list of list
of tuples, each inner list containing the matching bonds in the
MCES as tuples of bond indices from mol1 and mol2)DOC")
      .def("atomMatches",
           [](const RDKit::RascalMCES::RascalResult &self) {
             return convertVecPairInt(self.getAtomMatches());
           },
           "Likewise for atoms.")
      .def("largestFragmentOnly",
           [](RDKit::RascalMCES::RascalResult &self) {
             self.largestFragOnly();
           },
           "Function that cuts the MCES down to the single largest frag.  This cannot be undone.")
      .def_prop_ro("similarity",
                   &RDKit::RascalMCES::RascalResult::getSimilarity,
                   "Johnson similarity between 2 molecules.")
      .def_prop_ro("numFragments",
                   &RDKit::RascalMCES::RascalResult::getNumFrags,
                   "Number of fragments in MCES.")
      .def_prop_ro("largestFragmentSize",
                   &RDKit::RascalMCES::RascalResult::getLargestFragSize,
                   "Number of atoms in largest fragment.")
      .def_prop_ro("tier1Sim",
                   &RDKit::RascalMCES::RascalResult::getTier1Sim,
                   "The tier 1 similarity estimate.")
      .def_prop_ro("tier2Sim",
                   &RDKit::RascalMCES::RascalResult::getTier2Sim,
                   "The tier 2 similarity estimate.")
      .def_prop_ro("timedOut",
                   &RDKit::RascalMCES::RascalResult::getTimedOut,
                   "Whether it timed out.");

  nb::class_<RDKit::RascalMCES::RascalOptions>(m, "RascalOptions",
      "RASCAL Options")
      .def(nb::init<>())
      .def_rw("similarityThreshold",
              &RDKit::RascalMCES::RascalOptions::similarityThreshold,
              "Threshold below which MCES won't be run.  Between 0.0 and 1.0, default=0.7.")
      .def_rw("singleLargestFrag",
              &RDKit::RascalMCES::RascalOptions::singleLargestFrag,
              R"DOC(Return the just single largest fragment of the MCES. It is
equivalent to running with allBestMCEs=True, finding the result
with the largest largestFragmentSize, and calling its
largestFragmentOnly method.  This option may not produce the largest
possible single fragment that the molecules have in common. If you
definitely want that you may be better off using rdFMCS.)DOC")
      .def_rw("completeAromaticRings",
              &RDKit::RascalMCES::RascalOptions::completeAromaticRings,
              "If True (default), partial aromatic rings won't be returned.")
      .def_rw("ringMatchesRingOnly",
              &RDKit::RascalMCES::RascalOptions::ringMatchesRingOnly,
              "If True (default is False), ring bonds won't match non-ring bonds.")
      .def_rw("completeSmallestRings",
              &RDKit::RascalMCES::RascalOptions::completeSmallestRings,
              "If True (default is False), only complete rings present in both input molecule's RingInfo will be returned. Implies completeAromaticRings and ringMatchesRingOnly.")
      .def_rw("exactConnectionsMatch",
              &RDKit::RascalMCES::RascalOptions::exactConnectionsMatch,
              R"DOC(If True (default is False), atoms will only match atoms if they have the same
number of explicit connections.  E.g. the central atom of
C(C)(C) won't match either atom in CC)DOC")
      .def_rw("minFragSize",
              &RDKit::RascalMCES::RascalOptions::minFragSize,
              "Imposes a minimum on the number of atoms in a fragment that may be part of the MCES.  Default -1 means no minimum.")
      .def_rw("maxFragSeparation",
              &RDKit::RascalMCES::RascalOptions::maxFragSeparation,
              "Maximum number of bonds between fragments in the MCES for both to be reported.  Default -1 means no maximum.  If exceeded, the smaller fragment will be removed.")
      .def_rw("allBestMCESs",
              &RDKit::RascalMCES::RascalOptions::allBestMCESs,
              "If True, reports all MCESs found of the same maximum size.  Default False means just report the first found.")
      .def_rw("maxBestMCESs",
              &RDKit::RascalMCES::RascalOptions::maxBestMCESs,
              R"DOC(Some pathological cases produce huge numbers of equivalent solutions that can crash
the program due to memory depletion.  This caps the number of such solutions to prevent
this happening.  Default=10000.)DOC")
      .def_rw("returnEmptyMCES",
              &RDKit::RascalMCES::RascalOptions::returnEmptyMCES,
              R"DOC(If the estimated similarity between the 2 molecules doesn't meet the similarityThreshold, no results are returned.  If you want to know what the
estimates were, set this to True, and examine the tier1Sim and tier2Sim properties of the result then returned.)DOC")
      .def_rw("timeout",
              &RDKit::RascalMCES::RascalOptions::timeout,
              "Maximum time (in seconds) to spend on an individual MCESs determination.  Default 60, -1 means no limit.")
      .def_rw("maxBondMatchPairs",
              &RDKit::RascalMCES::RascalOptions::maxBondMatchPairs,
              R"DOC(Too many matching bond (vertex) pairs can cause the process to run out of memory.
The default of 1000 is fairly safe.  Increase with caution, as memory use increases
with the square of this number.)DOC")
      .def_rw("equivalentAtoms",
              &RDKit::RascalMCES::RascalOptions::equivalentAtoms,
              R"DOC(SMARTS strings defining atoms that should
be considered equivalent. e.g.
[F,Cl,Br,I] so all halogens will match each other.
Space-separated list allowing more than 1
class of equivalent atoms.)DOC")
      .def_rw("ignoreBondOrders",
              &RDKit::RascalMCES::RascalOptions::ignoreBondOrders,
              R"DOC(If True, will treat all bonds as the same,
irrespective of order.  Default=False.)DOC")
      .def_rw("ignoreAtomAromaticity",
              &RDKit::RascalMCES::RascalOptions::ignoreAtomAromaticity,
              R"DOC(If True, matches atoms solely on atomic number.
If False, will treat aromatic and aliphatic atoms
as different.  Default=True.)DOC")
      .def_rw("minCliqueSize",
              &RDKit::RascalMCES::RascalOptions::minCliqueSize,
              R"DOC(Normally, the minimum clique size is specified
via the similarityThreshold.  Sometimes it's
more convenient to specify it directly.  If this
is > 0, it will over-ride the
similarityThreshold.
Note that this refers to the
minimum number of BONDS in the MCES. Default=0.)DOC")
      .def("__setattr__", &safeSetattr);

  m.def("FindMCES",
        [](const RDKit::ROMol &mol1, const RDKit::ROMol &mol2,
           nb::object py_opts) {
          RDKit::RascalMCES::RascalOptions opts;
          if (!py_opts.is_none()) {
            opts = nb::cast<RDKit::RascalMCES::RascalOptions>(py_opts);
          }
          std::vector<RDKit::RascalMCES::RascalResult> results;
          {
            NOGIL gil;
            results = RDKit::RascalMCES::rascalMCES(mol1, mol2, opts);
          }
          nb::list pyres;
          for (auto &res : results) {
            pyres.append(res);
          }
          return pyres;
        },
        "mol1"_a, "mol2"_a, "opts"_a = nb::none(),
        R"DOC(Find one or more MCESs between the 2 molecules given.  Returns a list of
RascalResult objects.
- mol1
- mol2 The two molecules for which to find the MCES
- opts Optional RascalOptions object changing the default run mode.)DOC");

  nb::class_<RDKit::RascalMCES::RascalClusterOptions>(m, "RascalClusterOptions",
      "RASCAL Cluster Options.  Most of these pertain to RascalCluster calculations.  Only similarityCutoff is used by RascalButinaCluster.")
      .def(nb::init<>())
      .def_rw("similarityCutoff",
              &RDKit::RascalMCES::RascalClusterOptions::similarityCutoff,
              "Similarity cutoff for molecules to be in the same cluster.  Between 0.0 and 1.0, default=0.7.")
      .def_rw("minFragSize",
              &RDKit::RascalMCES::RascalClusterOptions::minFragSize,
              "The minimum number of atoms in a fragment for it to be included in the MCES.  Default=3.")
      .def_rw("maxNumFrags",
              &RDKit::RascalMCES::RascalClusterOptions::maxNumFrags,
              R"DOC(The maximum number of fragments allowed in the MCES for each pair of molecules. Default=2.  So that the MCES
isn't a lot of small fragments scattered around the molecules giving an inflated estimate of similarity.)DOC")
      .def_rw("numThreads",
              &RDKit::RascalMCES::RascalClusterOptions::numThreads,
              "Number of threads to use during clustering.  Default=-1 means all the hardware threads less one.")
      .def_rw("a",
              &RDKit::RascalMCES::RascalClusterOptions::a,
              "The penalty score for each unconnected component in the MCES. Default=0.05.")
      .def_rw("b",
              &RDKit::RascalMCES::RascalClusterOptions::b,
              "The weight of matched bonds over matched atoms. Default=2.")
      .def_rw("minIntraClusterSim",
              &RDKit::RascalMCES::RascalClusterOptions::minIntraClusterSim,
              R"DOC(Two pairs of molecules are included in the same cluster if the similarity between
their MCESs is greater than this.  Default=0.9.)DOC")
      .def_rw("clusterMergeSim",
              &RDKit::RascalMCES::RascalClusterOptions::clusterMergeSim,
              "Two clusters are merged if the fraction of molecules they have in common is greater than this.  Default=0.6.")
      .def("__setattr__", &safeSetattr);

  m.def("RascalCluster",
        [](nb::object mols, nb::object py_opts) {
          RDKit::RascalMCES::RascalClusterOptions opts;
          if (!py_opts.is_none()) {
            opts = nb::cast<RDKit::RascalMCES::RascalClusterOptions>(py_opts);
          }
          auto cmols = extractMols(mols);
          std::vector<RDKit::UINT_VECT> clusters;
          {
            NOGIL gil;
            clusters = RDKit::RascalMCES::rascalCluster(cmols, opts);
          }
          return packOutputMols(clusters);
        },
        "mols"_a, "opts"_a = nb::none(),
        R"DOC(Use the RASCAL MCES similarity metric to do fuzzy clustering.  Returns a list of lists
of molecules, each inner list being a cluster.  The last cluster is all the
molecules that didn't fit into another cluster (the singletons).
- mols List of molecules to be clustered
- opts Optional RascalOptions object changing the default run mode.)DOC");

  m.def("RascalButinaCluster",
        [](nb::object mols, nb::object py_opts) {
          RDKit::RascalMCES::RascalClusterOptions opts;
          if (!py_opts.is_none()) {
            opts = nb::cast<RDKit::RascalMCES::RascalClusterOptions>(py_opts);
          }
          auto cmols = extractMols(mols);
          std::vector<RDKit::UINT_VECT> clusters;
          {
            NOGIL gil;
            clusters = RDKit::RascalMCES::rascalButinaCluster(cmols, opts);
          }
          return packOutputMols(clusters);
        },
        "mols"_a, "opts"_a = nb::none(),
        R"DOC(Use the RASCAL MCES similarity metric to do Butina clustering
(Butina JCICS 39 747-750 (1999)).  Returns a list of lists of molecules,
each inner list being a cluster.  The last cluster is all the
molecules that didn't fit into another cluster (the singletons).
- mols List of molecules to be clustered
- opts Optional RascalOptions object changing the default run mode.)DOC");
}