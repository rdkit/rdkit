//
// Created by Gareth Jones on 5/30/2020.
//
// Copyright 2020 Schrodinger, Inc
//

#include <GraphMol/TautomerQuery/TautomerQuery.h>
#include <RDBoost/python.h>
#include <GraphMol/Wrap/substructmethods.h>

namespace python = boost::python;

namespace RDKit {

boost::shared_ptr<TautomerQuery> createDefaultTautomerQuery(const ROMol &mol) {
  return TautomerQuery::fromMol(mol);
}

bool tautomerIsSubstructOf(const TautomerQuery &self, const ROMol &target,
                           bool recursionPossible = true,
                           bool useChirality = false,
                           bool useQueryQueryMatches = false) {
  return HasSubstructMatch(target, self, recursionPossible, useChirality,
                           useQueryQueryMatches);
}

PyObject *tautomerGetSubstructMatch(const TautomerQuery &self,
                                    const ROMol &target,
                                    bool useChirality = false,
                                    bool useQueryQueryMatches = false) {
  return GetSubstructMatch(target, self, useChirality, useQueryQueryMatches);
}

PyObject *tautomerGetSubstructMatches(const TautomerQuery &self,
                                      const ROMol &target, bool uniquify = true,
                                      bool useChirality = false,
                                      bool useQueryQueryMatches = false,
                                      unsigned int maxMatches = 1000) {
  return GetSubstructMatches(target, self, uniquify, useChirality,
                             useQueryQueryMatches, maxMatches);
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

    SubstructMatchParameters matchParamters;
    matches = self.substructOf(target, params, &matchingTautomers);
  }
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

struct TautomerQuery_wrapper {
  static void wrap() {
    auto docString =
        "The Tautomer Query Class.\n\
  Creates a query that enables structure search accounting for matching of\n\
  Tautomeric forms\n";

    python::class_<TautomerQuery, boost::noncopyable>(
        "TautomerQuery", docString, python::no_init)
        .def("__init__", python::make_constructor(createDefaultTautomerQuery))
        .def("IsSubstructOf",
             (bool (*)(const TautomerQuery &self, const ROMol &target, bool,
                       bool, bool))tautomerIsSubstructOf,
             (python::arg("self"), python::arg("target"),
              python::arg("recursionPossible") = true,
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false))
        .def("GetSubstructMatch",
             (PyObject * (*)(const TautomerQuery &self, const ROMol &target,
                             bool, bool)) tautomerGetSubstructMatch,
             (python::arg("self"), python::arg("target"),
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false))
        .def("GetSubstructMatches",
             (PyObject * (*)(const TautomerQuery self, const ROMol &mol, bool,
                             bool, bool, unsigned int))
                 tautomerGetSubstructMatches,
             (python::arg("self"), python::arg("target"),
              python::arg("uniquify") = true,
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false,
              python::arg("maxMatches") = 1000))
        .def("GetSubstructMatchesWithTautomers",
             (PyObject * (*)(const TautomerQuery self, const ROMol &mol, bool,
                             bool, bool, unsigned int))
                 tautomerGetSubstructMatchesWithTautomers,
             (python::arg("self"), python::arg("target"),
              python::arg("uniquify") = true,
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false,
              python::arg("maxMatches") = 1000))
        .def("PatternFingerprintTemplate",
             &TautomerQuery::patternFingerprintTemplate,
             (python::arg("fingerprintSize") = 2048),
             python::return_value_policy<python::manage_new_object>());

    python::def("PatternFingerprintTautomerTarget",
                &TautomerQuery::patternFingerprintTarget,
                (python::arg("target"), python::arg("fingerprintSize") = 2048),
                python::return_value_policy<python::manage_new_object>());
  };
};

}  // namespace RDKit

void wrap_TautomerQuery() { RDKit::TautomerQuery_wrapper::wrap(); }
