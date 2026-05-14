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
    nb::object pmols,
    const ScaffoldNetwork::ScaffoldNetworkParams &params) {
  auto mols = pythonObjectToVect<std::shared_ptr<ROMol>>(pmols);
  ScaffoldNetwork::ScaffoldNetwork res;
  if (mols) {
    NOGIL gil;
    updateScaffoldNetwork(*mols, res, params);
  }
  return res;
}

void updateNetworkHelper(nb::object pmols,
                         ScaffoldNetwork::ScaffoldNetwork &net,
                         const ScaffoldNetwork::ScaffoldNetworkParams &params) {
  auto mols = pythonObjectToVect<std::shared_ptr<ROMol>>(pmols);
  if (mols) {
    NOGIL gil;
    updateScaffoldNetwork(*mols, net, params);
  }
}

ScaffoldNetwork::ScaffoldNetworkParams getBRICSParams() {
  return ScaffoldNetwork::getBRICSNetworkParams();
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
      .def_rw(
          "includeGenericScaffolds",
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
              "found in")
      .def("__setattr__", &safeSetattr);

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
      .def_ro("type", &ScaffoldNetwork::NetworkEdge::type,
              "type of the edge")
      .def("__str__", [](const ScaffoldNetwork::NetworkEdge &self) {
        std::ostringstream oss;
        oss << self;
        return oss.str();
      });

  nb::class_<ScaffoldNetwork::ScaffoldNetwork>(m, "ScaffoldNetwork")
      .def(nb::init<>(), "Default constructor")
      .def_ro("nodes", &ScaffoldNetwork::ScaffoldNetwork::nodes,
              "the sequence of SMILES defining the nodes")
      .def_ro("counts", &ScaffoldNetwork::ScaffoldNetwork::counts,
              "the number of times each node was encountered while "
              "building the network.")
      .def_ro("molCounts", &ScaffoldNetwork::ScaffoldNetwork::molCounts,
              "the number of moleclues each node was found in.")
      .def_ro("edges", &ScaffoldNetwork::ScaffoldNetwork::edges,
              "the sequence of network edges")
      .def("__getstate__",
           [](const ScaffoldNetwork::ScaffoldNetwork &self) {
             nb::list edges;
             for (const auto &e : self.edges) {
               edges.append(
                   nb::make_tuple(e.beginIdx, e.endIdx, (int)e.type));
             }
             return nb::make_tuple(self.nodes, self.counts, self.molCounts,
                                   edges);
           })
      .def("__setstate__",
           [](ScaffoldNetwork::ScaffoldNetwork &self, nb::tuple state) {
             if (nb::len(state) != 4) {
               throw std::runtime_error("invalid ScaffoldNetwork pickle state");
             }
             new (&self) ScaffoldNetwork::ScaffoldNetwork();
             self.nodes = nb::cast<std::vector<std::string>>(state[0]);
             self.counts = nb::cast<std::vector<unsigned>>(state[1]);
             self.molCounts = nb::cast<std::vector<unsigned>>(state[2]);
             for (auto e : nb::cast<nb::list>(state[3])) {
               auto et = nb::cast<nb::tuple>(e);
               self.edges.emplace_back(
                   nb::cast<size_t>(et[0]), nb::cast<size_t>(et[1]),
                   static_cast<ScaffoldNetwork::EdgeType>(nb::cast<int>(et[2])));
             }
           });

  m.def("CreateScaffoldNetwork", &createNetworkHelper, "mols"_a, "params"_a,
        "create (and return) a new network from a sequence of molecules");
  m.def("UpdateScaffoldNetwork", &updateNetworkHelper, "mols"_a, "network"_a,
        "params"_a, "update an existing network by adding molecules");
  m.def("BRICSScaffoldParams", &getBRICSParams,
        "Returns parameters for generating scaffolds using BRICS "
        "fragmentation rules");
}
