//
//  Copyright (C) 2019 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDBoost/python.h>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <GraphMol/GraphMol.h>
#include <RDBoost/Wrap.h>

#include <GraphMol/ScaffoldNetwork/ScaffoldNetwork.h>

namespace python = boost::python;
using namespace RDKit;

namespace {
ScaffoldNetwork::ScaffoldNetwork *createNetworkHelper(
    python::object pmols,
    const ScaffoldNetwork::ScaffoldNetworkParams &params) {
  auto mols = pythonObjectToVect<ROMOL_SPTR>(pmols);
  ScaffoldNetwork::ScaffoldNetwork *res = new ScaffoldNetwork::ScaffoldNetwork;
  updateScaffoldNetwork(*mols, *res, params);
  return res;
}
void updateNetworkHelper(python::object pmols,
                         ScaffoldNetwork::ScaffoldNetwork &net,
                         const ScaffoldNetwork::ScaffoldNetworkParams &params) {
  auto mols = pythonObjectToVect<ROMOL_SPTR>(pmols);
  updateScaffoldNetwork(*mols, net, params);
}

ScaffoldNetwork::ScaffoldNetworkParams *getBRICSParams() {
  return new ScaffoldNetwork::ScaffoldNetworkParams(
      ScaffoldNetwork::getBRICSNetworkParams());
}

}  // namespace

#ifdef RDK_USE_BOOST_SERIALIZATION
struct scaffoldnetwork_pickle_suite : python::pickle_suite {
  static python::tuple getinitargs(
      const RDKit::ScaffoldNetwork::ScaffoldNetwork &self) {
    std::stringstream oss;
    boost::archive::text_oarchive oa(oss);
    oa << self;
    std::string res = oss.str();
    return python::make_tuple(python::object(python::handle<>(
        PyBytes_FromStringAndSize(res.c_str(), res.length()))));
  };
};
#else
struct scaffoldnetwork_pickle_suite : python::pickle_suite {
  static python::tuple getinitargs(
      const RDKit::ScaffoldNetwork::ScaffoldNetwork &self) {
    throw_runtime_error("Pickling of ScaffoldNetwork instances is not enabled");
  };
};
#endif

BOOST_PYTHON_MODULE(rdScaffoldNetwork) {
  python::scope().attr("__doc__") =
      "Module containing functions for creating a Scaffold Network";

  // register the vector_indexing_suite for NetworkEdges
  // if it hasn't already been done.
  // logic from https://stackoverflow.com/a/13017303
  boost::python::type_info info =
      boost::python::type_id<std::vector<ScaffoldNetwork::NetworkEdge>>();
  const boost::python::converter::registration *reg =
      boost::python::converter::registry::query(info);
  if (reg == NULL || (*reg).m_to_python == NULL) {
    python::class_<std::vector<ScaffoldNetwork::NetworkEdge>>(
        "NetworkEdge_VECT")
        .def(python::vector_indexing_suite<
             std::vector<ScaffoldNetwork::NetworkEdge>>());
  }

  iterable_converter().from_python<std::vector<std::string>>();

  python::class_<ScaffoldNetwork::ScaffoldNetworkParams>(
      "ScaffoldNetworkParams", "Scaffold network parameters", python::init<>())
      .def(python::init<const std::vector<std::string> &>(
          (python::arg("bondBreakerSmartsList")),
          "Constructor taking a list of Reaction SMARTS for the fragmentation "
          "reactions"))
      .def_readwrite(
          "includeGenericScaffolds",
          &ScaffoldNetwork::ScaffoldNetworkParams::includeGenericScaffolds,
          "include scaffolds with all atoms replaced by dummies")
      .def_readwrite(
          "includeGenericBondScaffolds",
          &ScaffoldNetwork::ScaffoldNetworkParams::includeGenericBondScaffolds,
          "include scaffolds with all bonds replaced by single bonds")
      .def_readwrite(
          "includeScaffoldsWithoutAttachments",
          &ScaffoldNetwork::ScaffoldNetworkParams::
              includeScaffoldsWithoutAttachments,
          "remove attachment points from scaffolds and include the result")
      .def_readwrite(
          "includeScaffoldsWithAttachments",
          &ScaffoldNetwork::ScaffoldNetworkParams::
              includeScaffoldsWithAttachments,
          "Include the version of the scaffold with attachment points")
      .def_readwrite(
          "keepOnlyFirstFragment",
          &ScaffoldNetwork::ScaffoldNetworkParams::keepOnlyFirstFragment,
          "keep only the first fragment from the bond breaking rule")
      .def_readwrite(
          "pruneBeforeFragmenting",
          &ScaffoldNetwork::ScaffoldNetworkParams::pruneBeforeFragmenting,
          "Do a pruning/flattening step before starting fragmenting")
      .def_readwrite("flattenIsotopes",
                     &ScaffoldNetwork::ScaffoldNetworkParams::flattenIsotopes,
                     "remove isotopes when flattening")
      .def_readwrite("flattenChirality",
                     &ScaffoldNetwork::ScaffoldNetworkParams::flattenChirality,
                     "remove chirality and bond stereo when flattening")
      .def_readwrite(
          "flattenKeepLargest",
          &ScaffoldNetwork::ScaffoldNetworkParams::flattenKeepLargest,
          "keep only the largest fragment when doing flattening");

  python::enum_<ScaffoldNetwork::EdgeType>("EdgeType")
      .value("Fragment", ScaffoldNetwork::EdgeType::Fragment)
      .value("Generic", ScaffoldNetwork::EdgeType::Generic)
      .value("GenericBond", ScaffoldNetwork::EdgeType::GenericBond)
      .value("RemoveAttachment", ScaffoldNetwork::EdgeType::RemoveAttachment)
      .value("Initialize", ScaffoldNetwork::EdgeType::Initialize);

  python::class_<ScaffoldNetwork::NetworkEdge>(
      "NetworkEdge", "A scaffold network edge", python::no_init)
      .def_readonly("beginIdx", &ScaffoldNetwork::NetworkEdge::beginIdx,
                    "index of the begin node in node list")
      .def_readonly("endIdx", &ScaffoldNetwork::NetworkEdge::endIdx,
                    "index of the end node in node list")
      .def_readonly("type", &ScaffoldNetwork::NetworkEdge::type,
                    "type of the edge")
      .def(python::self_ns::str(python::self_ns::self));

  python::class_<ScaffoldNetwork::ScaffoldNetwork>(
      "ScaffoldNetwork", "A scaffold network", python::init<>())
#ifdef RDK_USE_BOOST_SERIALIZATION
      .def(python::init<const std::string &>())
      // enable pickle support
      .def_pickle(scaffoldnetwork_pickle_suite())
#endif
      .def_readonly("nodes", &ScaffoldNetwork::ScaffoldNetwork::nodes,
                    "the sequence of SMILES defining the nodes")
      .def_readonly("counts", &ScaffoldNetwork::ScaffoldNetwork::counts,
                    "the number of times each node was encountered while "
                    "building the network.")
      .def_readonly("edges", &ScaffoldNetwork::ScaffoldNetwork::edges,
                    "the sequence of network edges");

  python::def("CreateScaffoldNetwork", &createNetworkHelper,
              (python::arg("mols"), python::arg("params")),
              "create (and return) a new network from a sequence of molecules",
              python::return_value_policy<python::manage_new_object>());
  python::def(
      "UpdateScaffoldNetwork", &updateNetworkHelper,
      (python::arg("mols"), python::arg("network"), python::arg("params")),
      "update an existing network by adding molecules");

  python::def("BRICSScaffoldParams", &getBRICSParams,
              "Returns parameters for generating scaffolds using BRICS "
              "fragmentation rules",
              python::return_value_policy<python::manage_new_object>());
}
