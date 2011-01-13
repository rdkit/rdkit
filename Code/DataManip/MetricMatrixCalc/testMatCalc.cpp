// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MetricFuncs.h"
#include "MetricMatrixCalc.h"

#include <cstdlib>
#include <time.h>

using namespace RDDataManip;
int main() {
  
  int n = 10;
  int m = 3;
  int dlen = n*(n-1)/2;
  int i, j;
  double *desc = new double[n*m];
  double **desc2D = new double*[n];

  for (i = 0; i < n; i++) {
    desc2D[i] = desc;
    desc += m;
  }
  desc = desc2D[0];

  for (i = 0; i < n; i++) {
    for (j = 0; j < m ; j++) {
      desc[i*m + j] = ((double)rand())/10;
    }
  }

  //double x = EuclideanDistanceMetric(desc2D[0], desc2D[1], m);
  double *dmat = new double[dlen];
  MetricMatrixCalc<double**, double*> mmCalc;
  mmCalc.setMetricFunc(&EuclideanDistanceMetric<double *, double *>);
  mmCalc.calcMetricMatrix(desc2D, n, m, dmat);

  for (i = 0; i < dlen; i++) {
    std::cout << dmat[i] << "\n";
  }

  delete [] desc2D;
  delete [] desc;
  delete [] dmat;
  
  exit(0);
}
