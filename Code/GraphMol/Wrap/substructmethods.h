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
#include <boost/python.hpp>
#include <RDBoost/Wrap.h>
#include <GraphMol/Substruct/SubstructMatch.h>

namespace RDKit {

class pyFunctor {
 public:
  virtual ~pyFunctor() = default;
};

class pyFinalMatchFunctor : public pyFunctor {
 public:
  pyFinalMatchFunctor(python::object obj) : dp_obj(std::move(obj)) {}
  ~pyFinalMatchFunctor() = default;
  bool operator()(const ROMol &m, std::span<const unsigned int> match) {
    // grab the GIL
    PyGILStateHolder h;
    // boost::python doesn't handle std::span, so we need to convert the span to
    // a vector before calling into python:
    std::vector<unsigned int> matchVec(match.begin(), match.end());
    return python::extract<bool>(dp_obj(boost::ref(m), boost::ref(matchVec)));
  }

 private:
  python::object dp_obj;
};
template <typename T>
class pyMatchFunctor : public pyFunctor {
 public:
  pyMatchFunctor(python::object obj) : dp_obj(std::move(obj)) {}
  ~pyMatchFunctor() = default;
  bool operator()(const T &a1, const T &a2) {
    // grab the GIL
    PyGILStateHolder h;
    return python::extract<bool>(dp_obj(boost::ref(a1), boost::ref(a2)));
  }

 private:
  python::object dp_obj;
};

inline PyObject *convertMatches(const MatchVectType &matches) {
  PyObject *res = PyTuple_New(matches.size());
  std::for_each(matches.begin(), matches.end(), [res](const auto &pair) {
    PyTuple_SetItem(res, pair.first, PyInt_FromLong(pair.second));
  });
  return res;
}

inline PyObject *convertMatchesToTupleOfPairs(const MatchVectType &matches) {
  PyObject *res = PyTuple_New(matches.size());
  std::for_each(matches.begin(), matches.end(),
                [res, &matches](const auto &pair) {
                  PyObject *pyPair = PyTuple_New(2);
                  PyTuple_SetItem(pyPair, 0, PyInt_FromLong(pair.first));
                  PyTuple_SetItem(pyPair, 1, PyInt_FromLong(pair.second));
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
PyObject *GetSubstructMatch(T1 &mol, T2 &query, bool useChirality = false,
                            bool useQueryQueryMatches = false) {
  MatchVectType matches;
  {
    NOGIL gil;
    SubstructMatch(mol, query, matches, true, useChirality,
                   useQueryQueryMatches);
  }
  return convertMatches(matches);
}

template <typename T1, typename T2>
PyObject *GetSubstructMatches(T1 &mol, T2 &query, bool uniquify = true,
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
  return res;
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
PyObject *helpGetSubstructMatch(T1 &mol, T2 &query,
                                const SubstructMatchParameters &params) {
  SubstructMatchParameters ps = params;
  ps.maxMatches = 1;
  std::vector<MatchVectType> matches;
  pySubstructHelper(mol, query, params, matches);
  MatchVectType match;
  if (matches.size()) {
    match = matches[0];
  }
  return convertMatches(match);
}

template <typename T1, typename T2>
PyObject *helpGetSubstructMatches(T1 &mol, T2 &query,
                                  const SubstructMatchParameters &params) {
  std::vector<MatchVectType> matches;
  pySubstructHelper(mol, query, params, matches);
  PyObject *res = PyTuple_New(matches.size());
  for (size_t idx = 0; idx < matches.size(); idx++) {
    PyTuple_SetItem(res, idx, convertMatches(matches[idx]));
  }
  return res;
}

}  // namespace RDKit
#endif
