//
//  Copyright (C) 2017-2025 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RDKIT_SUBSTRUCT_METHODS_H
#define RDKIT_SUBSTRUCT_METHODS_H
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <RDBoost/Wrap_nb.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <nanobind/nanobind.h>

namespace RDKit {

class pyFunctor {
 public:
  virtual ~pyFunctor() = default;
};

class pyFinalMatchFunctor : public pyFunctor {
 public:
  pyFinalMatchFunctor(nb::object obj) : dp_obj(std::move(obj)) {}
  ~pyFinalMatchFunctor() = default;
  bool operator()(const ROMol &m, std::span<const unsigned int> match) {
    // grab the GIL
    PyGILStateHolder h;
    // boost::python doesn't handle std::span, so we need to convert the span to
    // a vector before calling into python:
    std::vector<unsigned int> matchVec(match.begin(), match.end());
    return nb::cast<bool>(dp_obj(m, matchVec));
  }

 private:
  nb::object dp_obj;
};

//   template <typename T>
// class pyMatchFunctor : public pyFunctor {
//  public:
//   pyMatchFunctor(python::object obj) : dp_obj(std::move(obj)) {}
//   ~pyMatchFunctor() = default;
//   bool operator()(const T &a1, const T &a2) {
//     // grab the GIL
//     PyGILStateHolder h;
//     return python::extract<bool>(dp_obj(boost::ref(a1), boost::ref(a2)));
//   }

//  private:
//   python::object dp_obj;
// };

inline PyObject *convertMatches(const MatchVectType &matches) {
  PyObject *res = PyTuple_New(matches.size());
  std::for_each(matches.begin(), matches.end(), [res](const auto &pair) {
    PyTuple_SetItem(res, pair.first, PyLong_FromLong(pair.second));
  });
  return res;
}

inline PyObject *convertMatchesToTupleOfPairs(const MatchVectType &matches) {
  PyObject *res = PyTuple_New(matches.size());
  std::for_each(matches.begin(), matches.end(),
                [res, &matches](const auto &pair) {
                  PyObject *pyPair = PyTuple_New(2);
                  PyTuple_SetItem(pyPair, 0, PyLong_FromLong(pair.first));
                  PyTuple_SetItem(pyPair, 1, PyLong_FromLong(pair.second));
                  PyTuple_SetItem(res, &pair - &matches.front(), pyPair);
                });
  return res;
}

template <typename T1, typename T2>
bool HasSubstructMatch(T1 &mol, T2 &query, bool recursionPossible = true,
                       bool useChirality = false,
                       bool useQueryQueryMatches = false) {
  NOGIL gil;
  MatchVectType res;
  return SubstructMatch(mol, query, res, recursionPossible, useChirality,
                        useQueryQueryMatches);
}

template <typename T1, typename T2>
nb::object GetSubstructMatch(T1 &mol, T2 &query, bool useChirality = false,
                             bool useQueryQueryMatches = false) {
  MatchVectType matches;
  {
    NOGIL gil;
    SubstructMatch(mol, query, matches, true, useChirality,
                   useQueryQueryMatches);
  }
  return nb::steal<nb::object>(convertMatches(matches));
}

template <typename T1, typename T2>
nb::object GetSubstructMatches(T1 &mol, T2 &query, bool uniquify = true,
                               bool useChirality = false,
                               bool useQueryQueryMatches = false,
                               unsigned int maxMatches = 1000) {
  std::vector<MatchVectType> matches;
  int matched;
  {
    NOGIL gil;
    matched = SubstructMatch(mol, query, matches, uniquify, true, useChirality,
                             useQueryQueryMatches, maxMatches);
  }
  PyObject *res = PyTuple_New(matched);
  for (int idx = 0; idx < matched; idx++) {
    PyTuple_SetItem(res, idx, convertMatches(matches[idx]));
  }
  return nb::steal<nb::object>(res);
}

template <typename T1, typename T2>
void pySubstructHelper(T1 &mol, T2 &query,
                       const SubstructMatchParameters &params,
                       std::vector<MatchVectType> &matches) {
  // it's safe to release the GIL here since the functors wrapping the python
  // callbacks will reacquire it as needed
  NOGIL gil;
  matches = SubstructMatch(mol, query, params);
}
template <typename T1, typename T2>
bool helpHasSubstructMatch(T1 &mol, T2 &query,
                           const SubstructMatchParameters &params) {
  SubstructMatchParameters ps = params;
  ps.maxMatches = 1;
  std::vector<MatchVectType> matches;
  pySubstructHelper(mol, query, params, matches);
  return matches.size() != 0;
}

template <typename T1, typename T2>
nb::object helpGetSubstructMatch(T1 &mol, T2 &query,
                                 const SubstructMatchParameters &params) {
  SubstructMatchParameters ps = params;
  ps.maxMatches = 1;
  std::vector<MatchVectType> matches;
  pySubstructHelper(mol, query, params, matches);
  MatchVectType match;
  if (matches.size()) {
    match = matches[0];
  }
  return nb::steal<nb::object>(convertMatches(match));
}

template <typename T1, typename T2>
nb::object helpGetSubstructMatches(T1 &mol, T2 &query,
                                   const SubstructMatchParameters &params) {
  std::vector<MatchVectType> matches;
  pySubstructHelper(mol, query, params, matches);
  PyObject *res = PyTuple_New(matches.size());
  for (size_t idx = 0; idx < matches.size(); idx++) {
    PyTuple_SetItem(res, idx, convertMatches(matches[idx]));
  }
  return nb::steal<nb::object>(res);
}

}  // namespace RDKit
#endif
