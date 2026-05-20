//
//  Copyright (C) 2019-2026 Greg Landrum, T5 Informatics GmbH,
//                           and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <sstream>

#include <GraphMol/GraphMol.h>
#include <RDBoost/Wrap_nb.h>

#include <GraphMol/ScaffoldNetwork/ScaffoldNetwork.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;

namespace {

ScaffoldNetwork::ScaffoldNetwork createNetworkHelper(
    const std::vector<std::shared_ptr<ROMol>> &mols,
    const ScaffoldNetwork::ScaffoldNetworkParams &params) {
  ScaffoldNetwork::ScaffoldNetwork res;
  {
    NOGIL gil;
    ScaffoldNetwork::updateScaffoldNetwork(mols, res, params);
  }
  return res;
}

void updateNetworkHelper(const std::vector<std::shared_ptr<ROMol>> &mols,
                         ScaffoldNetwork::ScaffoldNetwork &net,
                         const ScaffoldNetwork::ScaffoldNetworkParams &params) {
  NOGIL gil;
  ScaffoldNetwork::updateScaffoldNetwork(mols, net, params);
}

}  // namespace

NB_MODULE(rdScaffoldNetwork, m) {
  m.doc() = "Module containing functions for creating a Scaffold Network";

  nb::class_<ScaffoldNetwork::ScaffoldNetworkParams>(m, "ScaffoldNetworkParams")
      .def(nb::init<>(), "Default constructor")
      .def(nb::init<const std::vector<std::string> &>(),
           "bondBreakerSmartsList"_a,
           "Constructor taking a list of Reaction SMARTS for the fragmentation "
           "reactions")
      .def_rw("includeGenericScaffolds",
              &ScaffoldNetwork::ScaffoldNetworkParams::includeGenericScaffolds,
              "include scaffolds with all atoms replaced by dummies")
      .def_rw(
          "includeGenericBondScaffolds",
          &ScaffoldNetwork::ScaffoldNetworkParams::includeGenericBondScaffolds,
          "include scaffolds with all bonds replaced by single bonds")
      .def_rw(
          "includeScaffoldsWithoutAttachments",
          &ScaffoldNetwork::ScaffoldNetworkParams::
              includeScaffoldsWithoutAttachments,
          "remove attachment points from scaffolds and include the result")
      .def_rw(
          "includeScaffoldsWithAttachments",
          &ScaffoldNetwork::ScaffoldNetworkParams::
              includeScaffoldsWithAttachments,
          "Include the version of the scaffold with attachment points")
      .def_rw("includeNames",
              &ScaffoldNetwork::ScaffoldNetworkParams::includeNames,
              "Include molecules names of the input molecules")
      .def_rw(
          "keepOnlyFirstFragment",
          &ScaffoldNetwork::ScaffoldNetworkParams::keepOnlyFirstFragment,
          "keep only the first fragment from the bond breaking rule")
      .def_rw(
          "pruneBeforeFragmenting",
          &ScaffoldNetwork::ScaffoldNetworkParams::pruneBeforeFragmenting,
          "Do a pruning/flattening step before starting fragmenting")
      .def_rw("flattenIsotopes",
              &ScaffoldNetwork::ScaffoldNetworkParams::flattenIsotopes,
              "remove isotopes when flattening")
      .def_rw("flattenChirality",
              &ScaffoldNetwork::ScaffoldNetworkParams::flattenChirality,
              "remove chirality and bond stereo when flattening")
      .def_rw(
          "flattenKeepLargest",
          &ScaffoldNetwork::ScaffoldNetworkParams::flattenKeepLargest,
          "keep only the largest fragment when doing flattening")
      .def_rw("collectMolCounts",
              &ScaffoldNetwork::ScaffoldNetworkParams::collectMolCounts,
              "keep track of the number of molecules each scaffold was "
              "found in");

  nb::enum_<ScaffoldNetwork::EdgeType>(m, "EdgeType")
      .value("Fragment", ScaffoldNetwork::EdgeType::Fragment)
      .value("Generic", ScaffoldNetwork::EdgeType::Generic)
      .value("GenericBond", ScaffoldNetwork::EdgeType::GenericBond)
      .value("RemoveAttachment", ScaffoldNetwork::EdgeType::RemoveAttachment)
      .value("Initialize", ScaffoldNetwork::EdgeType::Initialize);

  nb::class_<ScaffoldNetwork::NetworkEdge>(m, "NetworkEdge")
      .def_ro("beginIdx", &ScaffoldNetwork::NetworkEdge::beginIdx,
              "index of the begin node in node list")
      .def_ro("endIdx", &ScaffoldNetwork::NetworkEdge::endIdx,
              "index of the end node in node list")
      .def_ro("type", &ScaffoldNetwork::NetworkEdge::type, "type of the edge")
      .def("__str__", [](const ScaffoldNetwork::NetworkEdge &self) {
        std::ostringstream oss;
        oss << self;
        return oss.str();
      });

  nb::class_<ScaffoldNetwork::ScaffoldNetwork>(m, "ScaffoldNetwork")
      .def(nb::init<>(), "Default constructor")
#ifdef RDK_USE_BOOST_SERIALIZATION
      .def("__init__",
           [](ScaffoldNetwork::ScaffoldNetwork &self, nb::bytes pkl) {
             std::string s(static_cast<const char *>(pkl.data()), pkl.size());
             new (&self) ScaffoldNetwork::ScaffoldNetwork(s);
           },
           "pkl"_a)
#endif
      .def_ro("nodes", &ScaffoldNetwork::ScaffoldNetwork::nodes,
              "the sequence of SMILES defining the nodes")
      .def_ro("counts", &ScaffoldNetwork::ScaffoldNetwork::counts,
              "the number of times each node was encountered while "
              "building the network.")
      .def_ro("molCounts", &ScaffoldNetwork::ScaffoldNetwork::molCounts,
              "the number of moleclues each node was found in.")
      .def_ro("edges", &ScaffoldNetwork::ScaffoldNetwork::edges,
              "the sequence of network edges")
#ifdef RDK_USE_BOOST_SERIALIZATION
      .def("__getstate__",
           [](const ScaffoldNetwork::ScaffoldNetwork &self) {
             std::stringstream oss;
             boost::archive::text_oarchive oa(oss);
             oa << self;
             const std::string res = oss.str();
             return std::make_tuple(nb::bytes(res.c_str(), res.size()));
           })
      .def("__setstate__",
           [](ScaffoldNetwork::ScaffoldNetwork &self, nb::tuple state) {
             if (nb::len(state) != 1) {
               throw nb::value_error("invalid ScaffoldNetwork pickle state");
             }
             auto b = nb::cast<nb::bytes>(state[0]);
             std::string pkl(static_cast<const char *>(b.data()), b.size());
             new (&self) ScaffoldNetwork::ScaffoldNetwork(pkl);
           });
#else
      .def("__getstate__",
           [](const ScaffoldNetwork::ScaffoldNetwork &self) {
             throw nb::runtime_error(
                 "Pickling of ScaffoldNetwork instances is not "
                 "enabled in this build");
           })
      .def("__setstate__",
           [](ScaffoldNetwork::ScaffoldNetwork &self, nb::tuple state) {
             throw nb::runtime_error(
                 "Pickling of ScaffoldNetwork instances is not "
                 "enabled in this build");
           });
#endif

  m.def("CreateScaffoldNetwork", &createNetworkHelper, "mols"_a, "params"_a,
        "create (and return) a new network from a sequence of molecules");
  m.def("UpdateScaffoldNetwork", &updateNetworkHelper, "mols"_a, "network"_a,
        "params"_a, "update an existing network by adding molecules");
  m.def("BRICSScaffoldParams", &ScaffoldNetwork::getBRICSNetworkParams,
        "Returns parameters for generating scaffolds using BRICS "
        "fragmentation rules");
}
