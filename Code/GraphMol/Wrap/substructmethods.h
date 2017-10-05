//
//  Copyright (C) 2017 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RDKIT_SUBSTRUCT_METHODS_H
#define RDKIT_SUBSTRUCT_METHODS_H
#include <boost/python.hpp>
#include <RDBoost/Wrap.h>
#include <GraphMol/Substruct/SubstructMatch.h>

namespace RDKit {

inline PyObject *convertMatches(MatchVectType &matches) {
  PyObject *res = PyTuple_New(matches.size());
  MatchVectType::const_iterator i;
  for (i = matches.begin(); i != matches.end(); i++) {
    PyTuple_SetItem(res, i->first, PyInt_FromLong(i->second));
  }
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
}  // end of RDKit namespace
#endif
