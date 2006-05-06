//
//              LAPACK++ 1.1 Linear Algebra Package 1.1
//               University of Tennessee, Knoxvilee, TN.
//            Oak Ridge National Laboratory, Oak Ridge, TN.
//        Authors: J. J. Dongarra, E. Greaser, R. Pozo, D. Walker
//                 (C) 1992-1996 All Rights Reserved
//
//                             NOTICE
//
// Permission to use, copy, modify, and distribute this software and
// its documentation for any purpose and without fee is hereby granted
// provided that the above copyright notice appear in all copies and
// that both the copyright notice and this permission notice appear in
// supporting documentation.
//
// Neither the Institutions (University of Tennessee, and Oak Ridge National
// Laboratory) nor the Authors make any representations about the suitability 
// of this software for any purpose.  This software is provided ``as is'' 
// without express or implied warranty.
//
// LAPACK++ was funded in part by the U.S. Department of Energy, the
// National Science Foundation and the State of Tennessee.


#include <stdlib.h>
#include "lafnames"
//#include LA_SYMM_MAT_DOUBLE_H // changed of VC++
#include "symd.h"
//#include LA_SYMM_FACT_DOUBLE_H // changed of VC++
#include "syfd.h"


int main(int argc, char *argv[])
{

    int m = 5, n = 5;

  LaSymmFactDouble SF;
  LaSymmMatDouble SM(m,n);
  LaVectorInt V(m);

    
  SF.ref(SM);
  cout << "SF.S(): \n";
  cout << SF.S() << endl;

  SF.pivot().ref(V);
  cout << "SF.pivot(): \n";
  cout << SF.pivot() << endl;

  cout << "SF.info(): \n";
  cout << SF.info() << endl;

  exit(0);
}
