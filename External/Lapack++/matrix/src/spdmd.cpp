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
//#include LA_SPD_MAT_DOUBLE_H
#include "spdmd.h"


int LaSpdMatDouble::debug_ = 0; // set debug to 0 initially.

int* LaSpdMatDouble::info_= new int;  // turn off info print flag.


LaSpdMatDouble& LaSpdMatDouble::copy(const LaSpdMatDouble &ob)
{

  if (debug())
  {
	  std::cout << " ob: " << ob.info() << std::endl;
  }

  int M = ob.size(0);
  int N = ob.size(1);
  int i,j;

  // current scheme in copy() is to detach the left-hand-side
  // from whatever it was pointing to.
  //

  resize(ob);

  for (j=0; j<N; j++)
    for (i=0; i<M; i++)
      if (i>=j)
        (*this)(i,j) = ob(i,j);

  if (debug())
  {
	  std::cout << " *this: " << this->info() << std::endl;
  }

  return *this;
}


LaSpdMatDouble& LaSpdMatDouble::operator=(const double &s)
{
  int M = size(0);
  int N = size(1);
  int i,j;

  for (j=0; j<N; j++)
    for (i=0; i<M; i++)
      if (j<=i)
        (*this)(i,j) = s;

  return *this;
}


std::ostream& operator<<(std::ostream& s, const LaSpdMatDouble& G)
{
    
  if (*(G.info_))     // print out only matrix info, not actual values
  {
      *(G.info_) = 0; // reset the flag
      s << "(" << G.size(0) << "x" << G.size(1) << ") " ;
      s << "Indices: " << G.index(0) << " " << G.index(1);
      s << " #ref: " << G.ref_count() ;
      s << " sa:" << G.shallow();
  }
  else
  {
    int M = G.size(0);
    int N = G.size(1);
    int i,j;

    for (j=0; j<N; j++)
    {
      for (i=0;i<M;i++)
        if (j>i)
          s << G(j,i) << " ";
        else
          s << G(i,j) << " ";
        s << std::endl;
    }
  }
  return s;
}


LaSpdMatDouble::operator LaGenMatDouble()
{
  int M = (*this).size(0);
  int N = (*this).size(1);
  int i,j;

  LaGenMatDouble G(M,N);

  for (j=0; j<N; j++)
    for (i=0; i<M; i++)
        G(i,j) = (*this)(i,j); 

  return G;
}


LaSpdMatDouble::operator LaLowerTriangMatDouble()
{
  int M = (*this).size(0);
  int N = (*this).size(1);

  LaLowerTriangMatDouble Lower(M,N);

  Lower.copy((*this).lower_data_);
  
  return Lower;
} 
