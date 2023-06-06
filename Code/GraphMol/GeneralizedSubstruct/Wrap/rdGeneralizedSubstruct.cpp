//
//  Copyright (C) 2023 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDBoost/python.h>
#include <GraphMol/GraphMol.h>
#include <RDBoost/Wrap.h>

#include <GraphMol/GeneralizedSubstruct/XQMol.h>
#include <GraphMol/Wrap/substructmethods.h>

namespace python = boost::python;
using namespace RDKit;
using namespace RDKit::GeneralizedSubstruct;

namespace {
bool hasSubstructHelper(const ROMol &mol, const ExtendedQueryMol &query,
                        SubstructMatchParameters *iparams) {
  SubstructMatchParameters params;
  if (iparams) {
    params = *iparams;
  }
  return hasSubstructMatch(mol, query, params);
}

PyObject *getSubstructHelper(const ROMol &mol, const ExtendedQueryMol &query,
                             SubstructMatchParameters *iparams) {
  std::vector<MatchVectType> matches;
  SubstructMatchParameters params;
  if (iparams) {
    params = *iparams;
  }
  {
    NOGIL gil;
    params.maxMatches = 1;
    matches = SubstructMatch(mol, query, params);
  }
  if (matches.empty()) {
    matches.push_back(MatchVectType());
  }
  return convertMatches(matches[0]);
}
PyObject *getSubstructsHelper(const ROMol &mol, const ExtendedQueryMol &query,
                              SubstructMatchParameters *iparams) {
  std::vector<MatchVectType> matches;
  SubstructMatchParameters params;
  if (iparams) {
    params = *iparams;
  }
  {
    NOGIL gil;
    matches = SubstructMatch(mol, query, params);
  }

  PyObject *res = PyTuple_New(matches.size());
  for (int idx = 0; idx < matches.size(); idx++) {
    PyTuple_SetItem(res, idx, convertMatches(matches[idx]));
  }
  return res;
}

ExtendedQueryMol *createExtendedQueryMolHelper(const ROMol &mol) {
  return new ExtendedQueryMol(std::move(createExtendedQueryMol(mol)));
}

}  // namespace

BOOST_PYTHON_MODULE(rdGeneralizedSubstruct) {
  python::scope().attr("__doc__") =
      "Module containing functions for generalized substructure searching";

  python::class_<ExtendedQueryMol, boost::noncopyable>(
      "ExtendedQueryMol", "Extended query molecule",
      python::init<const std::string &, bool>(
          (python::args("text"), python::args("isJSON") = false),
          "constructor from either a binary string (from ToBinary()) or a JSON string."))
      .def("InitFromBinary", &ExtendedQueryMol::initFromBinary)
      .def("InitFromJSON", &ExtendedQueryMol::initFromJSON)
      .def("ToBinary", &ExtendedQueryMol::toBinary)
      .def("ToJSON", &ExtendedQueryMol::toJSON);

  python::def(
      "MolHasXQMSubstructMatch", &hasSubstructHelper,
      (python::arg("mol"), python::arg("query"),
       python::arg("params") = python::object()),
      "determines whether or not a molecule is a match to a generalized substructure query");
  python::def(
      "MolGetXQMSubstructMatch", &getSubstructHelper,
      (python::arg("mol"), python::arg("query"),
       python::arg("params") = python::object()),
      "returns first match (if any) of a molecule to a generalized substructure query");
  python::def(
      "MolGetXQMSubstructMatches", &getSubstructsHelper,
      (python::arg("mol"), python::arg("query"),
       python::arg("params") = python::object()),
      "returns all matches (if any) of a molecule to a generalized substructure query");

  /*
  RDKIT_GENERALIZEDSUBSTRUCT_EXPORT std::vector<MatchVectType>
  SubstructMatch( const ROMol &mol, const ExtendedQueryMol &query, const
  SubstructMatchParameters &params = SubstructMatchParameters());

  inline bool hasSubstructMatch(

  */

  python::def(
      "CreateExtendedQueryMol", createExtendedQueryMolHelper,
      (python::arg("mol")),
      python::return_value_policy<python::manage_new_object>(),
      R"DOC(Converts a molecule into an extended query molecule for generalized substructure searching. 
  This includes enumerating link nodes, varaible attachment points, SRUs, and constructing TautomerQueries)DOC");
}
