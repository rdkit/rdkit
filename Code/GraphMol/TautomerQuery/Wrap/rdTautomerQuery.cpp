//
// Created by Gareth Jones on 5/30/2020.
//
// Copyright 2020-2022 Schrodinger, Inc and other RDKit contributors
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include <RDBoost/Wrap.h>
#include <RDBoost/python.h>
#include <GraphMol/TautomerQuery/TautomerQuery.h>
#include <GraphMol/Wrap/substructmethods.h>
#include <RDBoost/python_streambuf.h>

#define PY_ARRAY_UNIQUE_SYMBOL rdtautomerquery_array_API

namespace python = boost::python;
using boost_adaptbx::python::streambuf;
using namespace RDKit;

namespace {

TautomerQuery *createDefaultTautomerQuery(const ROMol &mol) {
  return TautomerQuery::fromMol(mol);
}

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

PyObject *tautomerGetSubstructMatch(const TautomerQuery &self,
                                    const ROMol &target,
                                    bool useChirality = false,
                                    bool useQueryQueryMatches = false) {
  return GetSubstructMatch(target, self, useChirality, useQueryQueryMatches);
}

PyObject *tautomerGetSubstructMatchWithParams(
    const TautomerQuery &self, const ROMol &target,
    const SubstructMatchParameters &params) {
  return helpGetSubstructMatch(target, self, params);
}

PyObject *tautomerGetSubstructMatchesWithParams(
    const TautomerQuery &self, const ROMol &target,
    const SubstructMatchParameters &params) {
  return helpGetSubstructMatches(target, self, params);
}

PyObject *tautomerGetSubstructMatches(const TautomerQuery &self,
                                      const ROMol &target, bool uniquify = true,
                                      bool useChirality = false,
                                      bool useQueryQueryMatches = false,
                                      unsigned int maxMatches = 1000) {
  return GetSubstructMatches(target, self, uniquify, useChirality,
                             useQueryQueryMatches, maxMatches);
}

PyObject *getTautomers(const TautomerQuery &self) {
  auto tauts = self.getTautomers();
  auto nTauts = tauts.size();
  auto tuple = PyTuple_New(nTauts);
  for (size_t i = 0; i < nTauts; i++) {
    PyTuple_SetItem(tuple, i,
                    python::converter::shared_ptr_to_python(tauts[i]));
  }
  return tuple;
}

PyObject *matchesWithTautomersToTuple(std::vector<MatchVectType> &matches,
                                      const MOL_SPTR_VECT &matchingTautomers) {
  int const numberMatches = matches.size();
  PyObject *res = PyTuple_New(numberMatches);
  for (int idx = 0; idx < numberMatches; idx++) {
    PyObject *pair = PyTuple_New(2);
    PyTuple_SetItem(pair, 0, convertMatches(matches[idx]));
    PyTuple_SetItem(
        pair, 1,
        python::converter::shared_ptr_to_python(matchingTautomers[idx]));
    PyTuple_SetItem(res, idx, pair);
  }
  return res;
}

PyObject *tautomerGetSubstructMatchesWithTautomersWithParams(
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

  return matchesWithTautomersToTuple(matches, matchingTautomers);
}

PyObject *tautomerGetSubstructMatchesWithTautomers(
    const TautomerQuery &self, const ROMol &target, bool uniquify = true,
    bool useChirality = false, bool useQueryQueryMatches = false,
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
  return matchesWithTautomersToTuple(matches, matchingTautomers);
}

}  // namespace

namespace RDKit {
struct tautomerquery_pickle_suite : python::pickle_suite {
  static python::tuple getinitargs(const TautomerQuery &self) {
    if (!TautomerQueryCanSerialize()) {
      throw_runtime_error("Pickling of TautomerQuery instances is not enabled");
    }
    auto res = self.serialize();
    return python::make_tuple(python::object(python::handle<>(
        PyBytes_FromStringAndSize(res.c_str(), res.length()))));
  };
};

void toStream(const TautomerQuery &tq, python::object &fileobj) {
  streambuf ss(fileobj, 't');
  streambuf::ostream ost(ss);
  tq.toStream(ost);
}

void initFromStream(TautomerQuery &tq, python::object &fileobj) {
  streambuf ss(fileobj,
               'b');  // python StringIO can't seek, so need binary data
  streambuf::istream is(ss);
  tq.initFromStream(is);
}

python::object TQToBinary(const TautomerQuery &tq) {
  auto res = tq.serialize();
  return python::object(
      python::handle<>(PyBytes_FromStringAndSize(res.c_str(), res.length())));
}
}  // namespace RDKit

struct TautomerQuery_wrapper {
  static void wrap() {
    RegisterVectorConverter<size_t>("UnsignedLong_Vect");

    auto docString =
        "The Tautomer Query Class.\n\
  Creates a query that enables structure search accounting for matching of\n\
  Tautomeric forms\n";

    python::class_<TautomerQuery, boost::noncopyable>(
        "TautomerQuery", docString, python::no_init)
        .def("__init__", python::make_constructor(createDefaultTautomerQuery))
        .def("__init__", python::make_constructor(&TautomerQuery::fromMol))
        .def(python::init<std::string>())
        .def("IsSubstructOf", tautomerIsSubstructOf,
             (python::arg("self"), python::arg("target"),
              python::arg("recursionPossible") = true,
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false))
        .def(
            "IsSubstructOf", tautomerIsSubstructOfWithParams,
            (python::arg("self"), python::arg("target"), python::arg("params")))
        .def("GetSubstructMatch", tautomerGetSubstructMatch,
             (python::arg("self"), python::arg("target"),
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false))
        .def(
            "GetSubstructMatch", tautomerGetSubstructMatchWithParams,
            (python::arg("self"), python::arg("target"), python::arg("params")))
        .def("GetSubstructMatches", tautomerGetSubstructMatches,
             (python::arg("self"), python::arg("target"),
              python::arg("uniquify") = true,
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false,
              python::arg("maxMatches") = 1000))
        .def(
            "GetSubstructMatches", tautomerGetSubstructMatchesWithParams,
            (python::arg("self"), python::arg("target"), python::arg("params")))
        .def("GetSubstructMatchesWithTautomers",
             tautomerGetSubstructMatchesWithTautomers,
             (python::arg("self"), python::arg("target"),
              python::arg("uniquify") = true,
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false,
              python::arg("maxMatches") = 1000))
        .def(
            "GetSubstructMatchesWithTautomers",
            tautomerGetSubstructMatchesWithTautomersWithParams,
            (python::arg("self"), python::arg("target"), python::arg("params")))
        .def("PatternFingerprintTemplate",
             &TautomerQuery::patternFingerprintTemplate,
             (python::arg("fingerprintSize") = 2048),
             python::return_value_policy<python::manage_new_object>())
        .def("GetTemplateMolecule", &TautomerQuery::getTemplateMolecule,
             python::return_internal_reference<>())
        .def("GetModifiedAtoms", &TautomerQuery::getModifiedAtoms)
        .def("GetModifiedBonds", &TautomerQuery::getModifiedBonds)
        .def("GetTautomers", getTautomers)
        .def("ToBinary", TQToBinary)
        .def_pickle(tautomerquery_pickle_suite());

    python::def("PatternFingerprintTautomerTarget",
                &TautomerQuery::patternFingerprintTarget,
                (python::arg("target"), python::arg("fingerprintSize") = 2048),
                python::return_value_policy<python::manage_new_object>());

    python::def(
        "TautomerQueryCanSerialize", TautomerQueryCanSerialize,
        "Returns True if the TautomerQuery is serializable "
        "(requires that the RDKit was built with boost::serialization)");
  }
};

void wrap_TautomerQuery() { TautomerQuery_wrapper::wrap(); }

BOOST_PYTHON_MODULE(rdTautomerQuery) { wrap_TautomerQuery(); }
