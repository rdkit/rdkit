//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#ifndef __RD_METRICFUNCS_H__
#define __RD_METRICFUNCS_H__
#include <cmath>
#include <DataStructs/BitOps.h>

namespace RDDataManip {
  
  template <typename T1, typename T2> double EuclideanDistance(const T1 &v1, const T2 &v2, int dim) 
    {
      int i;
      double diff, dist = 0.0;
      for (i = 0; i < dim; i++) {
        diff = (double)v1[i] - (double)v2[i];
        dist += (diff*diff);
      }
      return sqrt(dist);
    };

  template <typename T1, typename T2> double TanimotoDistance(const T1 &bv1, const T2 &bv2, int dim) 
    {
      // the dim parameter is actually irrelevant here but we have to include it to deal with 
      // template version of setMetricFunc in MetricMatricCalc
      return (1.0 - TanimotoSimilarity(bv1, bv2));
    }

  template <typename T1, typename T2> double TanimotoSimilarity(const T1 &bv1, const T2 &bv2, int dim)
    {
      double res;
      if(bv1.GetNumBits()!=bv2.GetNumBits()){
        // fold down automatically:
        unsigned int nb1=bv1.GetNumBits();
        unsigned int nb2=bv2.GetNumBits();
    	
    	if(nb1 < nb2){
    	  T2 *bv3=FoldFingerprint(bv2,nb2/nb1);
    	  res = TanimotoSimilarity(bv1,*bv3);
    	  delete bv3;
    	} else {
    	  T2 *bv3=FoldFingerprint(bv1,nb1/nb2);
    	  res = TanimotoSimilarity(*bv3,bv2);
    	  delete bv3;
    	}
          } else {
    	res = TanimotoSimilarity(bv1, bv2);
      }
      return res;
    }

}

#endif
    
    
