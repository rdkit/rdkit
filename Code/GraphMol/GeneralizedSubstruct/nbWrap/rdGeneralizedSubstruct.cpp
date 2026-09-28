//
//  Copyright (C) 2023-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/unique_ptr.h>

#include <GraphMol/GraphMol.h>
#include <RDBoost/Wrap_nb.h>

#include <GraphMol/GeneralizedSubstruct/XQMol.h>
#include <GraphMol/nbWrap/substructmethods.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;
using namespace RDKit::GeneralizedSubstruct;

namespace {

nb::bytes XQMolToBinary(const ExtendedQueryMol &self) {
  std::string res;
  {
    NOGIL gil;
    res = self.toBinary();
  }
  return nb::bytes(res.c_str(), res.length());
}

bool hasSubstructHelper(const ROMol &mol, const ExtendedQueryMol &query,
                        const SubstructMatchParameters &params) {
  bool res;
  {
    NOGIL gil;
    res = hasSubstructMatch(mol, query, params);
  }
  return res;
}

std::vector<int> getSubstructHelper(const ROMol &mol,
                                    const ExtendedQueryMol &query,
                                    SubstructMatchParameters params) {
  std::vector<MatchVectType> matches;
  {
    NOGIL gil;
    params.maxMatches = 1;
    matches = SubstructMatch(mol, query, params);
  }
  if (matches.empty()) {
    matches.push_back(MatchVectType());
  }
  return convertMatch(matches[0]);
}

std::vector<std::vector<int>> getSubstructsHelper(
    const ROMol &mol, const ExtendedQueryMol &query,
    const SubstructMatchParameters &params) {
  std::vector<MatchVectType> matches;
  {
    NOGIL gil;
    matches = SubstructMatch(mol, query, params);
  }
  std::vector<std::vector<int>> res;
  res.reserve(matches.size());
  for (const auto &match : matches) {
    res.push_back(convertMatch(match));
  }
  return res;
}

}  // namespace

NB_MODULE(rdGeneralizedSubstruct, m) {
  m.doc() =
      R"DOC(Module containing functions for generalized substructure searching)DOC";

  nb::class_<ExtendedQueryMol>(
      m, "ExtendedQueryMol",
      R"DOC(Extended query molecule for use in generalized substructure searching.)DOC")
      .def(nb::init<const std::string &, bool>(), "text"_a, "isJSON"_a = false,
           R"DOC(constructor from either a binary string (from ToBinary()) or a JSON string.)DOC")
      .def(
          "__init__",
          [](ExtendedQueryMol &self, nb::bytes data) {
            std::string s(data.c_str(), data.size());
            new (&self) ExtendedQueryMol(s, false);
          },
          "data"_a,
          R"DOC(constructor from binary data returned by ToBinary().)DOC")
      .def("InitFromBinary", &ExtendedQueryMol::initFromBinary, "pkl"_a)
      .def("InitFromJSON", &ExtendedQueryMol::initFromJSON, "text"_a)
      .def("ToBinary", XQMolToBinary)
      .def("ToJSON", &ExtendedQueryMol::toJSON)
      .def("PatternFingerprintQuery",
           &ExtendedQueryMol::patternFingerprintQuery,
           "fingerprintSize"_a = 2048);

  m.def(
      "MolHasSubstructMatch", &hasSubstructHelper, "mol"_a, "query"_a,
      "params"_a = SubstructMatchParameters(),
      R"DOC(determines whether or not a molecule is a match to a generalized substructure query)DOC");
  m.def(
      "MolGetSubstructMatch", &getSubstructHelper, "mol"_a, "query"_a,
      "params"_a = SubstructMatchParameters(),
      R"DOC(returns first match (if any) of a molecule to a generalized substructure query)DOC");
  m.def(
      "MolGetSubstructMatches", &getSubstructsHelper, "mol"_a, "query"_a,
      "params"_a = SubstructMatchParameters(),
      R"DOC(returns all matches (if any) of a molecule to a generalized substructure query)DOC");

  m.def(
      "PatternFingerprintTarget", &patternFingerprintTargetMol, "target"_a,
      "fingerprintSize"_a = 2048,
      R"DOC(Creates a pattern fingerprint for a target molecule that is compatible with an extended query)DOC");

  m.def(
      "CreateExtendedQueryMol",
      [](const ROMol &mol, bool doEnumeration, bool doTautomers,
         bool adjustQueryProperties,
         const MolOps::AdjustQueryParameters &ps) {
        return createExtendedQueryMol(mol, doEnumeration, doTautomers,
                                      adjustQueryProperties, ps);
      },
      "mol"_a, "doEnumeration"_a = true, "doTautomers"_a = true,
      "adjustQueryProperties"_a = false,
      "adjustQueryParameters"_a = MolOps::AdjustQueryParameters(),
      R"DOC(Creates an ExtendedQueryMol from the input molecule

This takes a query molecule and, conceptually, performs the following steps to
produce an ExtendedQueryMol:

  1. Enumerates features like Link Nodes and SRUs
  2. Converts everything into TautomerQueries
  3. Runs adjustQueryProperties()

Each step is optional
)DOC");
}
