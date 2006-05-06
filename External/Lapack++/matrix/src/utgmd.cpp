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

#include "lafnames.h"
//#include LA_UPPER_TRIANG_MAT_DOUBLE_H // changed of VC++
#include "utgmd.h"


double LaUpperTriangMatDouble::outofbounds_ = 0; // initialize outobounds

int LaUpperTriangMatDouble::debug_ = 0;  // turn off global debug 
                                // use A.debug(1) to turn on/off
                                // and A.debug() to check current status.

int* LaUpperTriangMatDouble::info_= new int;  // turn off info print flag.

LaUpperTriangMatDouble& LaUpperTriangMatDouble::copy(LaUpperTriangMatDouble &ob)
{

  if (debug())
  {
	  std::cout << " ob: " << ob.info() << std::endl;
  }


  int M = ob.size(0);
  int N = ob.size(1);
  int i,j;

  resize(ob);

  for (i=0; i<M; i++)
    for (j=0; j<N; j++)
      if (j>=i)
        (*this)(i,j) = ob(i,j);

  if (debug())
  {
	  std::cout << " *this: " << this->info() << std::endl;
  }

  return *this;
}

LaUpperTriangMatDouble& LaUpperTriangMatDouble::operator=(const double &s)
{
    int M = size(0);
    int N = size(1);
    int i, j;

    for (i=0; i<M; i++)
        for (j=i; j<N; j++)
            (*this)(i,j) = s;

    return *this;
}


std::ostream &operator<<(std::ostream &s, const LaUpperTriangMatDouble &ob)
{
  if (*(ob.info_))     // print out only matrix info, not actual values
  {
      *(ob.info_) = 0; // reset the flag
      s << "(" << ob.size(0) << "x" << ob.size(1) << ") " ;
      s << "Indices: " << ob.index(0) << " " << ob.index(1);
      s << " #ref: " << ob.ref_count() ;
      s << " sa:" << ob.shallow();
  }
  else
  {
    int M = ob.size(0);
    int N = ob.size(1);
    int i,j;

    for (i=0; i<M; i++)
    {
     for (j=0;j<N;j++)
     {
       if (j>=i)
       s << ob(i,j) << "  ";
     }
     s << std::endl;
    }
  }
  return s;
}
