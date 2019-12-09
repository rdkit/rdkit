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
    python::object pmols, const ScaffoldNetwork::ScaffoldNetworkParams &params) {
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

BOOST_PYTHON_MODULE(rdScaffoldNetwork) {
  python::scope().attr("__doc__") =
      "Module containing functions for creating a Scaffold Network";

  // register the vector_indexing_suite for SubstanceGroups
  // if it hasn't already been done.
  // logic from https://stackoverflow.com/a/13017303
  boost::python::type_info info =
      boost::python::type_id<std::vector<ScaffoldNetwork::NetworkEdge>>();
  const boost::python::converter::registration *reg =
      boost::python::converter::registry::query(info);
  if (reg == NULL || (*reg).m_to_python == NULL) {
    python::class_<std::vector<ScaffoldNetwork::NetworkEdge>>("NetworkEdge_VECT")
        .def(python::vector_indexing_suite<std::vector<ScaffoldNetwork::NetworkEdge>>());
  }

  std::string docString = R"DOC()DOC";

  python::class_<ScaffoldNetwork::ScaffoldNetworkParams>(
      "ScaffoldNetworkParams", "Scaffold network parameters", python::init<>())
      .def_readwrite("includeGenericScaffolds",
                     &ScaffoldNetwork::ScaffoldNetworkParams::includeGenericScaffolds)
      .def_readwrite("includeGenericBondScaffolds",
                     &ScaffoldNetwork::ScaffoldNetworkParams::includeGenericScaffolds)
      .def_readwrite("includeScaffoldsWithoutAttachments",
                     &ScaffoldNetwork::ScaffoldNetworkParams::includeScaffoldsWithoutAttachments)
      .def_readwrite("includeScaffoldsWithAttachments",
                     &ScaffoldNetwork::ScaffoldNetworkParams::includeScaffoldsWithAttachments)
      .def_readwrite("keepOnlyFirstFragment",
                     &ScaffoldNetwork::ScaffoldNetworkParams::keepOnlyFirstFragment)
      .def_readwrite("pruneBeforeFragmenting",
                     &ScaffoldNetwork::ScaffoldNetworkParams::pruneBeforeFragmenting)
      .def_readwrite("flattenIsotopes", &ScaffoldNetwork::ScaffoldNetworkParams::flattenIsotopes)
      .def_readwrite("flattenChirality",
                     &ScaffoldNetwork::ScaffoldNetworkParams::flattenChirality)
      .def_readwrite("flattenKeepLargest",
                     &ScaffoldNetwork::ScaffoldNetworkParams::flattenKeepLargest);

  python::enum_<ScaffoldNetwork::EdgeType>("EdgeType")
      .value("Fragment", ScaffoldNetwork::EdgeType::Fragment)
      .value("Generic", ScaffoldNetwork::EdgeType::Generic)
      .value("GenericBond", ScaffoldNetwork::EdgeType::GenericBond)
      .value("RemoveAttachment", ScaffoldNetwork::EdgeType::RemoveAttachment)
      .value("Initialize", ScaffoldNetwork::EdgeType::Initialize);

  python::class_<ScaffoldNetwork::NetworkEdge>("NetworkEdge", "A scaffold network edge",
                              python::no_init)
      .def_readonly("beginIdx", &ScaffoldNetwork::NetworkEdge::beginIdx)
      .def_readonly("endIdx", &ScaffoldNetwork::NetworkEdge::endIdx)
      .def_readonly("type", &ScaffoldNetwork::NetworkEdge::type)
      .def(python::self_ns::str(python::self_ns::self));

  python::class_<ScaffoldNetwork::ScaffoldNetwork>(
      "ScaffoldNetwork", "A scaffold network", python::init<>())
      .def_readonly("nodes", &ScaffoldNetwork::ScaffoldNetwork::nodes)
      .def_readonly("counts", &ScaffoldNetwork::ScaffoldNetwork::counts)
      .def_readonly("edges", &ScaffoldNetwork::ScaffoldNetwork::edges);

  python::def("CreateScaffoldNetwork", &createNetworkHelper,
              (python::arg("mols"), python::arg("params")),
              "create (and return) a new network",
              python::return_value_policy<python::manage_new_object>());
  python::def(
      "UpdateScaffoldNetwork", &updateNetworkHelper,
      (python::arg("mols"), python::arg("network"), python::arg("params")),
      "update an existing network");

  python::def("BRICSScaffoldParams", &getBRICSParams,
              "Returns parameters for generating BRICS scaffolds.",
              python::return_value_policy<python::manage_new_object>());
}
