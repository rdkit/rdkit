//
//  Copyright (C) 2019  Greg Landrum
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RDKIT_PICKERHELPERS_H
#define RDKIT_PICKERHELPERS_H

#include <vector>
#include <DataStructs/BitOps.h>

// NOTE: TANIMOTO and DICE provably return the same results for the diversity
// picking this is still here just in case we ever later want to support other
//    methods.
typedef enum { TANIMOTO = 1, DICE } DistanceMethod;

template <typename BV>
class pyBVFunctor {
 public:
  pyBVFunctor(const std::vector<const BV *> &obj, DistanceMethod method)
      : d_obj(obj), d_method(method) {}
  ~pyBVFunctor() = default;
  double operator()(unsigned int i, unsigned int j) {
    double res = 0.0;
    switch (d_method) {
      case TANIMOTO:
        res = 1. - TanimotoSimilarity(*d_obj[i], *d_obj[j]);
        break;
      case DICE:
        res = 1. - DiceSimilarity(*d_obj[i], *d_obj[j]);
        break;
      default:
        throw_value_error("unsupported similarity value");
    }
    return res;
  }

 private:
  const std::vector<const BV *> &d_obj;
  DistanceMethod d_method;
};

class pyobjFunctor {
 public:
  pyobjFunctor(python::object obj) : dp_obj(std::move(obj)) {}
  ~pyobjFunctor() = default;
  double operator()(unsigned int i, unsigned int j) {
    return python::extract<double>(dp_obj(i, j));
  }

 private:
  python::object dp_obj;
};

#endif  // RDKIT_PICKERHELPERS_H
