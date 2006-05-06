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
//#include LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H // changed of VC++
#include "ultgmd.h"


double LaUnitLowerTriangMatDouble::outofbounds_ = 0; // set outofbounds_. 

int LaUnitLowerTriangMatDouble::debug_ = 0; // set debug to 0 initially.

int* LaUnitLowerTriangMatDouble::info_= new int;  // turn off info print flag.

LaUnitLowerTriangMatDouble& LaUnitLowerTriangMatDouble::copy(LaUnitLowerTriangMatDouble &ob)
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
      if (i>j)
        (*this)(i,j) = ob(i,j);

  if (debug())
  {
	  std::cout << " *this: " << this->info() << std::endl;
  }

  return *this;
}

 LaUnitLowerTriangMatDouble& LaUnitLowerTriangMatDouble::operator=(double s)
{

  int M = (*this).size(0);
  int N = (*this).size(1);
  int i,j;
  
  for (j=0; j<N; j++)
    for (i=0; i<M; i++)
      if (i>j)
        (*this)(i,j) = s;

  return *this;
}

 

std::ostream &operator<<(std::ostream &s, const LaUnitLowerTriangMatDouble &ob)
{
  if (*(ob.info_))     // print out only matrix info, not actual values
  {
      ob.info_ = 0; // reset the flag
      s << "(" << ob.size(0) << "x" << ob.size(1) << ") " ;
      s << "Indices: " << ob.index(0) << " " << ob.index(1);
      s << " #ref: " << ob.ref_count()  ;
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
       if (i>j)
       s << ob(i,j) << "  ";
     }
     s << std::endl;
    }
  }
  return s;
}
