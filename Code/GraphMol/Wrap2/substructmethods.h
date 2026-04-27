//
//  Copyright (C) 2026 Greg Landrum
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
#include <RDBoost/Wrap_nb.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>
#include <GraphMol/Substruct/SubstructMatch.h>

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
    return nb::cast<bool>(dp_obj(boost::ref(m), boost::ref(matchVec)));
  }

 private:
  nb::object dp_obj;
};
template <typename T>
class pyMatchFunctor : public pyFunctor {
 public:
  pyMatchFunctor(nb::object obj) : dp_obj(std::move(obj)) {}
  ~pyMatchFunctor() = default;
  bool operator()(const T &a1, const T &a2) {
    // grab the GIL
    PyGILStateHolder h;
    return nb::cast<bool>(dp_obj(boost::ref(a1), boost::ref(a2)));
  }

 private:
  nb::object dp_obj;
};

inline std::vector<int> convertMatch(const MatchVectType &match) {
  std::vector<int> res(match.size());
  std::for_each(match.begin(), match.end(),
                [&res](const auto &pair) { res[pair.first] = pair.second; });
  return res;
}

#if 0
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
#endif

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
bool HasSubstructMatch(T1 &mol, T2 &query, bool recursionPossible = true,
                       bool useChirality = false,
                       bool useQueryQueryMatches = false) {
  NOGIL gil;
  MatchVectType res;
  return SubstructMatch(mol, query, res, recursionPossible, useChirality,
                        useQueryQueryMatches);
}
template <typename T1, typename T2>
bool helpHasSubstructMatch(
    T1 &mol, T2 &query, const std::optional<SubstructMatchParameters> params) {
  SubstructMatchParameters ps = params.value_or(SubstructMatchParameters());
  ps.maxMatches = 1;
  std::vector<MatchVectType> matches;
  pySubstructHelper(mol, query, ps, matches);
  return matches.size() != 0;
}

template <typename T1, typename T2>
std::vector<int> GetSubstructMatch(T1 &mol, T2 &query,
                                   bool useChirality = false,
                                   bool useQueryQueryMatches = false) {
  MatchVectType matches;
  {
    NOGIL gil;
    SubstructMatch(mol, query, matches, true, useChirality,
                   useQueryQueryMatches);
  }
  std::vector<int> res = convertMatch(matches);
  return res;
}

template <typename T1, typename T2>
std::vector<int> helpGetSubstructMatch(
    T1 &mol, T2 &query, const std::optional<SubstructMatchParameters> params) {
  SubstructMatchParameters ps = params.value_or(SubstructMatchParameters());
  ps.maxMatches = 1;
  std::vector<MatchVectType> matches;
  pySubstructHelper(mol, query, ps, matches);
  MatchVectType match;
  if (matches.size()) {
    match = matches[0];
  }
  return convertMatch(match);
}

template <typename T1, typename T2>
std::vector<std::vector<int>> GetSubstructMatches(
    T1 &mol, T2 &query, bool uniquify = true, bool useChirality = false,
    bool useQueryQueryMatches = false, unsigned int maxMatches = 1000) {
  std::vector<MatchVectType> matches;
  int matched;
  {
    NOGIL gil;
    matched = SubstructMatch(mol, query, matches, uniquify, true, useChirality,
                             useQueryQueryMatches, maxMatches);
  }
  std::vector<std::vector<int>> res;
  res.reserve(matched);
  std::for_each(matches.begin(), matches.end(), [&res](const auto &match) {
    res.push_back(convertMatch(match));
  });
  return res;
}

template <typename T1, typename T2>
std::vector<std::vector<int>> helpGetSubstructMatches(
    T1 &mol, T2 &query, const std::optional<SubstructMatchParameters> params) {
  SubstructMatchParameters ps = params.value_or(SubstructMatchParameters());
  std::vector<MatchVectType> matches;
  pySubstructHelper(mol, query, ps, matches);
  std::vector<std::vector<int>> res;
  res.reserve(matches.size());
  std::for_each(matches.begin(), matches.end(), [&res](const auto &match) {
    res.push_back(convertMatch(match));
  });
  return res;
}

}  // namespace RDKit
#endif
