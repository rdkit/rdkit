//
//  Copyright (C) 2003-2007 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_METRICFUNCS_H__
#define __RD_METRICFUNCS_H__
#include <cmath>
#include <DataStructs/BitOps.h>

namespace RDDataManip {
  //! return the Euclidean distance between two vectors
  template <typename T1, typename T2>
  double EuclideanDistanceMetric(const T1 &v1, const T2 &v2, unsigned int dim) {
    double dist = 0.0;
    for (unsigned int i = 0; i < dim; i++) {
      double diff = static_cast<double>(v1[i]) - static_cast<double>(v2[i]);
      dist += (diff*diff);
    }
    return sqrt(dist);
  };


  // FIX: there's no reason to have this tied to TanimotoSimilarity... could include
  // a different sim function as a template param
  //! return the Tanimoto distance (1-TanimotoSimilarity) between two bit vectors
  template <typename T1, typename T2>
  double TanimotoDistanceMetric(const T1 &bv1, const T2 &bv2, unsigned int dim) {
    // the dim parameter is actually irrelevant here but we have to include it to deal with 
    // template version of setMetricFunc in MetricMatricCalc
    return (1.0 - SimilarityWrapper(bv1, bv2,(double (*)(const T1&,const T2&))TanimotoSimilarity));
  };

  //! return the Tanimoto similarity between two bit vectors
  template <typename T1, typename T2>
  double TanimotoSimilarityMetric(const T1 &bv1, const T2 &bv2, unsigned int dim) {
    return SimilarityWrapper(bv1,bv2,(double (*)(const T1&,const T2&))TanimotoSimilarity);
  };
}

#endif
    
    
