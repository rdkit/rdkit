// $Id$
//
//                 Copyright 2001-2006
//                   Rational Discovery LLC 
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//


#include "types.h"

namespace RDKit{
//  template <typename T>
//  T larger_of(T arg1,T arg2) { return arg1>arg2 ? arg1 : arg2; };

double round(double num){
  double floorVal = floor(num);
  double ceilVal = ceil(num);
  return num-floorVal > ceilVal-num ? ceilVal : floorVal;
};


  void Union(const INT_VECT &r1, const INT_VECT &r2, INT_VECT &res) {
    res.resize(0);
    res = r1;
    INT_VECT_CI ri;
    for (ri = r2.begin(); ri != r2.end(); ri++) {
      if (std::find(res.begin(), res.end(), (*ri)) == res.end()){
	res.push_back(*ri);
      }
    }
  }

  void Intersect(const INT_VECT &r1, const INT_VECT &r2, INT_VECT &res) {
    res.resize(0);
    INT_VECT_CI ri;
    for (ri = r1.begin(); ri != r1.end(); ri++) {
      if (std::find(r2.begin(), r2.end(), (*ri)) != r2.end()){
	res.push_back(*ri);
      }
    }
  }

  void Union(const VECT_INT_VECT &rings, INT_VECT &res, const INT_VECT *exclude) {
    res.resize(0);
    INT_VECT ring;
    unsigned int id;
    unsigned int nrings = static_cast<unsigned int>(rings.size());
    INT_VECT_CI ri;
    
    for (id = 0; id < nrings; id++) {
      if (exclude) {
	if (std::find(exclude->begin(), exclude->end(), static_cast<int>(id)) != exclude->end()) {
	  continue;
	}
      }
      ring = rings[id];
      for (ri = ring.begin(); ri != ring.end(); ri++) {
	if (std::find(res.begin(), res.end(), (*ri)) == res.end()) {
	  res.push_back(*ri);
	}
      }
    }
  }

  int nextCombination(INT_VECT &comb, int tot) {
    int nelem = static_cast<int>(comb.size());
    int celem = nelem - 1;
    
    while (comb[celem] == (tot - nelem + celem)) {
      celem--;
      if (celem < 0) {
	return -1;
      }
    }

    unsigned int i;
    comb[celem] += 1;
    for (i = celem+1; i < comb.size(); i++) {
      comb[i] = comb[i-1] + 1;
    }
    return celem;
  }

  

}
