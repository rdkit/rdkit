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

python::object XQMolToBinary(const ExtendedQueryMol &self) {
  std::string res;
  {
    NOGIL gil;
    res = self.toBinary();
  }
  python::object retval = python::object(
      python::handle<>(PyBytes_FromStringAndSize(res.c_str(), res.length())));
  return retval;
}

bool hasSubstructHelper(const ROMol &mol, const ExtendedQueryMol &query,
                        SubstructMatchParameters *iparams) {
  SubstructMatchParameters params;
  if (iparams) {
    params = *iparams;
  }
  bool res;
  {
    NOGIL gil;
    res = hasSubstructMatch(mol, query, params);
  }
  return res;
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
  for (auto idx = 0u; idx < matches.size(); idx++) {
    PyTuple_SetItem(res, idx, convertMatches(matches[idx]));
  }
  return res;
}

ExtendedQueryMol *createExtendedQueryMolHelper(
    const ROMol &mol, bool doEnumeration, bool doTautomers,
    bool adjustQueryProperties, MolOps::AdjustQueryParameters *ps) {
  MolOps::AdjustQueryParameters defaults;
  if (!ps) {
    ps = &defaults;
  }
  return new ExtendedQueryMol(createExtendedQueryMol(
      mol, doEnumeration, doTautomers, adjustQueryProperties, *ps));
}

}  // namespace

BOOST_PYTHON_MODULE(rdGeneralizedSubstruct) {
  python::scope().attr("__doc__") =
      "Module containing functions for generalized substructure searching";

  python::class_<ExtendedQueryMol, boost::noncopyable>(
      "ExtendedQueryMol",
      "Extended query molecule for use in generalized substructure searching.",
      python::init<const std::string &, bool>(
          (python::args("self"), python::args("text"),
           python::args("isJSON") = false),
          "constructor from either a binary string (from ToBinary()) or a JSON string."))
      .def("InitFromBinary", &ExtendedQueryMol::initFromBinary,
           python::args("self", "pkl"))
      .def("InitFromJSON", &ExtendedQueryMol::initFromJSON,
           python::args("self", "text"))
      .def("ToBinary", XQMolToBinary, python::args("self"))
      .def("ToJSON", &ExtendedQueryMol::toJSON, python::args("self"))
      .def("PatternFingerprintQuery",
           &ExtendedQueryMol::patternFingerprintQuery,
           (python::arg("self"), python::arg("fingerprintSize") = 2048));

  python::def(
      "MolHasSubstructMatch", &hasSubstructHelper,
      (python::arg("mol"), python::arg("query"),
       python::arg("params") = python::object()),
      "determines whether or not a molecule is a match to a generalized substructure query");
  python::def(
      "MolGetSubstructMatch", &getSubstructHelper,
      (python::arg("mol"), python::arg("query"),
       python::arg("params") = python::object()),
      "returns first match (if any) of a molecule to a generalized substructure query");
  python::def(
      "MolGetSubstructMatches", &getSubstructsHelper,
      (python::arg("mol"), python::arg("query"),
       python::arg("params") = python::object()),
      "returns all matches (if any) of a molecule to a generalized substructure query");

  python::def(
      "PatternFingerprintTarget", &patternFingerprintTargetMol,
      (python::arg("target"), python::arg("fingerprintSize") = 2048),
      "Creates a pattern fingerprint for a target molecule that is compatible with an extended query");

  python::def("CreateExtendedQueryMol", createExtendedQueryMolHelper,
              (python::arg("mol"), python::arg("doEnumeration") = true,
               python::arg("doTautomers") = true,
               python::arg("adjustQueryProperties") = false,
               python::arg("adjustQueryParameters") = python::object()),
              python::return_value_policy<python::manage_new_object>(),
              R"DOC(Creates an ExtendedQueryMol from the input molecule

  This takes a query molecule and, conceptually, performs the following steps to
  produce an ExtendedQueryMol:

    1. Enumerates features like Link Nodes and SRUs
    2. Converts everything into TautomerQueries
    3. Runs adjustQueryProperties()

  Each step is optional
)DOC");
}
