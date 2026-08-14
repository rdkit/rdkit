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
#include <nanobind/stl/optional.h>
#include <nanobind/stl/tuple.h>
#include <GraphMol/Substruct/SubstructMatch.h>

namespace RDKit {

template <typename T, typename U>
class pyMatchFunctor {
 public:
  pyMatchFunctor(nb::object obj) : dp_callable(nb::borrow<nb::callable>(obj)) {
    if (!dp_callable.is_valid()) {
      throw std::runtime_error("Matcher object is not callable");
    }

    // Special case to create a shortcut for AtomCoordsMatchFunctors:
    if constexpr (std::is_same_v<T, Atom> && std::is_same_v<U, Atom>) {
      nb::try_cast<AtomCoordsMatchFunctor *>(obj, dp_coordsMatchedFunc);
    }
  }

  ~pyMatchFunctor() = default;

  bool operator()(const T &a1, const U &a2) const {
    // grab the GIL
    PyGILStateHolder h;

    if constexpr (std::is_same_v<T, ROMol> &&
                  std::is_same_v<U, std::span<const unsigned int>>) {
      // nanobind doesn't support std::span, so we need to convert the span to
      // a vector before calling into python. This might be dependent
      // on the nanobind version.
      std::vector<unsigned int> matchVec(a2.begin(), a2.end());
      return nb::cast<bool>(dp_callable(&a1, &matchVec));
    }

    if constexpr (std::is_same_v<T, Atom> && std::is_same_v<U, Atom>) {
      // If the callable is a subclass of AtomCoordsMatchFunctor,
      // we can take a shortcut and avoid passing the args through
      // the nanobind wrappers twice.
      if (dp_coordsMatchedFunc) {
        return (*dp_coordsMatchedFunc)(a1, a2);
      }
    }

    // We want to pass pointers to the callable: the args will
    // be passed through the nanobind wrappers (twice?), and
    // it seems there are copies made along the way. If we pass
    // Atom or Bond references, they seem to lose Mol ownership
    // info, resulting in an exception during the fn call.
    return nb::cast<bool>(dp_callable(&a1, &a2));
  }

 private:
  // The callable is borrowed, so ownership of the Python function
  // is shared with instances.
  nb::callable dp_callable;

  // This function is guaranteed to exist for the lifetime of the callable.
  AtomCoordsMatchFunctor *dp_coordsMatchedFunc = nullptr;
};

inline std::vector<int> convertMatch(const MatchVectType &match) {
  std::vector<int> res(match.size());
  std::for_each(match.begin(), match.end(),
                [&res](const auto &pair) { res[pair.first] = pair.second; });
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
