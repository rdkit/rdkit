//
// Created by Gareth Jones on 5/30/2020.
//
// Copyright 2020-2026 Schrodinger, Inc and other RDKit contributors
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/pair.h>

#include <GraphMol/TautomerQuery/TautomerQuery.h>
#include <GraphMol/nbWrap/substructmethods.h>
#include <RDBoost/python_streambuf_nb.h>

namespace nb = nanobind;
using namespace nb::literals;
using boost_adaptbx::python::streambuf;
using namespace RDKit;

namespace {

bool tautomerIsSubstructOf(const TautomerQuery &self, const ROMol &target,
                           bool recursionPossible = true,
                           bool useChirality = false,
                           bool useQueryQueryMatches = false) {
  return HasSubstructMatch(target, self, recursionPossible, useChirality,
                           useQueryQueryMatches);
}

bool tautomerIsSubstructOfWithParams(const TautomerQuery &self,
                                     const ROMol &target,
                                     const SubstructMatchParameters &params) {
  return helpHasSubstructMatch(target, self, params);
}

std::vector<int> tautomerGetSubstructMatch(const TautomerQuery &self,
                                           const ROMol &target,
                                           bool useChirality = false,
                                           bool useQueryQueryMatches = false) {
  return GetSubstructMatch(target, self, useChirality, useQueryQueryMatches);
}

std::vector<int> tautomerGetSubstructMatchWithParams(
    const TautomerQuery &self, const ROMol &target,
    const SubstructMatchParameters &params) {
  return helpGetSubstructMatch(target, self, params);
}

std::vector<std::vector<int>> tautomerGetSubstructMatchesWithParams(
    const TautomerQuery &self, const ROMol &target,
    const SubstructMatchParameters &params) {
  return helpGetSubstructMatches(target, self, params);
}

std::vector<std::vector<int>> tautomerGetSubstructMatches(
    const TautomerQuery &self, const ROMol &target, bool uniquify = true,
    bool useChirality = false, bool useQueryQueryMatches = false,
    unsigned int maxMatches = 1000) {
  return GetSubstructMatches(target, self, uniquify, useChirality,
                             useQueryQueryMatches, maxMatches);
}

std::vector<std::shared_ptr<RDKit::ROMol>> getTautomers(
    const TautomerQuery &self) {
  std::vector<std::shared_ptr<RDKit::ROMol>> out;
  for (const auto &sptr : self.getTautomers()) out.push_back(toStd(sptr));
  return out;
}

// Returns a list of (match, tautomer) pairs
std::vector<std::tuple<std::vector<int>, std::shared_ptr<RDKit::ROMol>>>
matchesWithTautomersToList(std::vector<MatchVectType> &matches,
                           const MOL_SPTR_VECT &matchingTautomers) {
  std::vector<std::tuple<std::vector<int>, std::shared_ptr<RDKit::ROMol>>> res;
  res.reserve(matches.size());
  for (size_t idx = 0; idx < matches.size(); ++idx) {
    res.emplace_back(convertMatch(matches[idx]), toStd(matchingTautomers[idx]));
  }
  return res;
}

std::vector<std::tuple<std::vector<int>, std::shared_ptr<RDKit::ROMol>>>
tautomerGetSubstructMatchesWithTautomersWithParams(
    const TautomerQuery &self, const ROMol &target,
    const SubstructMatchParameters &params) {
  std::vector<MatchVectType> matches;
  std::vector<ROMOL_SPTR> matchingTautomers;

  if (params.extraFinalCheck) {
    // NOTE: Because we are going into/out of python here, we can't
    // run with NOGIL
    matches = self.substructOf(target, params, &matchingTautomers);
  } else {
    NOGIL gil;
    matches = self.substructOf(target, params, &matchingTautomers);
  }

  return matchesWithTautomersToList(matches, matchingTautomers);
}

std::vector<std::tuple<std::vector<int>, std::shared_ptr<RDKit::ROMol>>>
tautomerGetSubstructMatchesWithTautomers(const TautomerQuery &self,
                                         const ROMol &target,
                                         bool uniquify = true,
                                         bool useChirality = false,
                                         bool useQueryQueryMatches = false,
                                         unsigned int maxMatches = 1000) {
  std::vector<MatchVectType> matches;
  std::vector<ROMOL_SPTR> matchingTautomers;
  SubstructMatchParameters params;
  params.uniquify = uniquify;
  params.useChirality = useChirality;
  params.useQueryQueryMatches = useQueryQueryMatches;
  params.maxMatches = maxMatches;

  {
    NOGIL gil;
    matches = self.substructOf(target, params, &matchingTautomers);
  }
  return matchesWithTautomersToList(matches, matchingTautomers);
}

nb::bytes TQToBinary(const TautomerQuery &tq) {
  auto res = tq.serialize();
  return nb::bytes(res.c_str(), res.length());
}
std::string TQToString(const TautomerQuery &tq) { return tq.serialize(); }

void toStream(const TautomerQuery &tq, nb::object fileobj) {
  streambuf ss(fileobj, 't');
  streambuf::ostream ost(ss);
  tq.toStream(ost);
}

void initFromStream(TautomerQuery &tq, nb::object fileobj) {
  streambuf ss(fileobj,
               'b');  // python StringIO can't seek, so need binary data
  streambuf::istream is(ss);
  tq.initFromStream(is);
}

}  // namespace

NB_MODULE(rdTautomerQuery, m) {
  m.doc() = R"DOC(Module for tautomer-aware substructure searching.

Provides the TautomerQuery class which enables substructure searching
that accounts for tautomeric forms of the query molecule.
)DOC";

  nb::class_<TautomerQuery>(m, "TautomerQuery",
                            R"DOC(The Tautomer Query Class.
Creates a query that enables structure search accounting for matching of
Tautomeric forms
)DOC")
      .def(nb::new_([](nb::bytes b) {
             return new TautomerQuery(
                 std::string(static_cast<const char *>(b.data()), b.size()));
           }),
           "pickle"_a, "Construct a TautomerQuery from a pickle bytes.")
      .def(nb::new_([](std::string pkl) { return new TautomerQuery(pkl); }),
           "pickle"_a, "Construct a TautomerQuery from a pickle string.")
      .def(nb::new_(
               [](const ROMol &mol, const std::string &tautomerTransformFile) {
                 return TautomerQuery::fromMol(mol, tautomerTransformFile);
               }),
           "mol"_a, "tautomerTransformFile"_a = std::string(),
           "Construct a TautomerQuery from a molecule.")
      .def(
          "IsSubstructOf", tautomerIsSubstructOf, "target"_a,
          "recursionPossible"_a = true, "useChirality"_a = false,
          "useQueryQueryMatches"_a = false,
          R"DOC(Check if this tautomer query is a substructure of the target molecule.

ARGUMENTS:
 - target: the target molecule
 - recursionPossible: (optional) allow recursive queries (default True)
 - useChirality: (optional) use chirality in matching (default False)
 - useQueryQueryMatches: (optional) use query-query matching logic (default False)

RETURNS: True or False
)DOC")
      .def(
          "IsSubstructOf", tautomerIsSubstructOfWithParams, "target"_a,
          "params"_a,
          R"DOC(Check if this tautomer query is a substructure of the target molecule.

ARGUMENTS:
 - target: the target molecule
 - params: SubstructMatchParameters object

RETURNS: True or False
)DOC")
      .def(
          "GetSubstructMatch", tautomerGetSubstructMatch, "target"_a,
          "useChirality"_a = false, "useQueryQueryMatches"_a = false,
          R"DOC(Return the first substructure match of this tautomer query in the target.

ARGUMENTS:
 - target: the target molecule
 - useChirality: (optional) use chirality in matching (default False)
 - useQueryQueryMatches: (optional) use query-query matching logic (default False)

RETURNS: a tuple of atom indices on match, or empty tuple on no match
)DOC")
      .def(
          "GetSubstructMatch", tautomerGetSubstructMatchWithParams, "target"_a,
          "params"_a,
          R"DOC(Return the first substructure match of this tautomer query in the target.

ARGUMENTS:
 - target: the target molecule
 - params: SubstructMatchParameters object

RETURNS: a tuple of atom indices on match, or empty tuple on no match
)DOC")
      .def(
          "GetSubstructMatches", tautomerGetSubstructMatches, "target"_a,
          "uniquify"_a = true, "useChirality"_a = false,
          "useQueryQueryMatches"_a = false, "maxMatches"_a = 1000,
          R"DOC(Return all substructure matches of this tautomer query in the target.

ARGUMENTS:
 - target: the target molecule
 - uniquify: (optional) only return unique matches (default True)
 - useChirality: (optional) use chirality in matching (default False)
 - useQueryQueryMatches: (optional) use query-query matching logic (default False)
 - maxMatches: (optional) maximum number of matches to return (default 1000)

RETURNS: a tuple of tuples of atom indices
)DOC")
      .def(
          "GetSubstructMatches", tautomerGetSubstructMatchesWithParams,
          "target"_a, "params"_a,
          R"DOC(Return all substructure matches of this tautomer query in the target.

ARGUMENTS:
 - target: the target molecule
 - params: SubstructMatchParameters object

RETURNS: a tuple of tuples of atom indices
)DOC")
      .def("GetSubstructMatchesWithTautomers",
           tautomerGetSubstructMatchesWithTautomers, "target"_a,
           "uniquify"_a = true, "useChirality"_a = false,
           "useQueryQueryMatches"_a = false, "maxMatches"_a = 1000,
           R"DOC(Return all substructure matches with their matching tautomers.

ARGUMENTS:
 - target: the target molecule
 - uniquify: (optional) only return unique matches (default True)
 - useChirality: (optional) use chirality in matching (default False)
 - useQueryQueryMatches: (optional) use query-query matching logic (default False)
 - maxMatches: (optional) maximum number of matches to return (default 1000)

RETURNS: a list of (match, tautomer) pairs
)DOC")
      .def("GetSubstructMatchesWithTautomers",
           tautomerGetSubstructMatchesWithTautomersWithParams, "target"_a,
           "params"_a,
           R"DOC(Return all substructure matches with their matching tautomers.

ARGUMENTS:
 - target: the target molecule
 - params: SubstructMatchParameters object

RETURNS: a list of (match, tautomer) pairs
)DOC")
      .def("PatternFingerprintTemplate",
           &TautomerQuery::patternFingerprintTemplate,
           "fingerprintSize"_a = 2048, nb::rv_policy::take_ownership,
           R"DOC(Return the pattern fingerprint of the template molecule.

ARGUMENTS:
 - fingerprintSize: (optional) size of the fingerprint (default 2048)

RETURNS: an ExplicitBitVect fingerprint
)DOC")
      .def("GetTemplateMolecule", &TautomerQuery::getTemplateMolecule,
           nb::rv_policy::reference_internal,
           R"DOC(Return the template molecule used for substructure searching.
)DOC")
      .def("GetModifiedAtoms", &TautomerQuery::getModifiedAtoms,
           R"DOC(Return the indices of tautomeric atoms.
)DOC")
      .def("GetModifiedBonds", &TautomerQuery::getModifiedBonds,
           R"DOC(Return the indices of tautomeric bonds.
)DOC")
      .def("GetTautomers", getTautomers,
           R"DOC(Return the list of tautomers of the query molecule.
)DOC")
      .def(
          "ToBinary", TQToBinary,
          R"DOC(Return a binary string (pickle) representation of this TautomerQuery.
)DOC")
      .def("ToStream", toStream, "fileobj"_a,
           R"DOC(Serialize this TautomerQuery to a file-like object.
)DOC")
      .def("InitFromStream", initFromStream, "fileobj"_a,
           R"DOC(Initialize this TautomerQuery from a file-like object.
)DOC")
      .def("__getstate__", getObjectState<TautomerQuery, TQToString>)
      .def("__setstate__", setObjectState<TautomerQuery>);
  m.def(
      "PatternFingerprintTautomerTarget",
      &TautomerQuery::patternFingerprintTarget, "target"_a,
      "fingerprintSize"_a = 2048, nb::rv_policy::take_ownership,
      R"DOC(Return the pattern fingerprint of a target molecule for tautomer searching.

ARGUMENTS:
 - target: the target molecule
 - fingerprintSize: (optional) size of the fingerprint (default 2048)

RETURNS: an ExplicitBitVect fingerprint
)DOC");

  m.def("TautomerQueryCanSerialize", TautomerQueryCanSerialize,
        R"DOC(Returns True if the TautomerQuery is serializable
(requires that the RDKit was built with boost::serialization)
)DOC");
}
