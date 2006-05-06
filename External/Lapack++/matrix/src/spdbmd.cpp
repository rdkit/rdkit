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
//#include LA_SPD_BAND_MAT_DOUBLE_H
#include "spdbmd.h"

double LaSpdBandMatDouble::outofbounds_ = 0; // initialize outofbounds_. 

int LaSpdBandMatDouble::debug_ = 0; // initialize debug to off.

int* LaSpdBandMatDouble::info_ = new int; // turn off info print flag.

LaSpdBandMatDouble& LaSpdBandMatDouble::operator=(double scalar)
{

  int i,j;

  for (i=0; i<N_; i++)
    for (j=0; j<N_; j++)
    {
      if(((i>=j)&&(i-j<=kl_)))
        (*this)(i,j) = scalar;
    }

  return *this;
}

LaSpdBandMatDouble& LaSpdBandMatDouble::copy(const LaSpdBandMatDouble &ob)
{

  int i,j;

  resize(ob);

  for (i=0; i<ob.N_; i++)
    for (j=0; j<ob.N_; j++)
        (*this)(i,j) = ob(i,j);

  return *this;
}

std::ostream& operator<<(std::ostream &s, const LaSpdBandMatDouble &ob)
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
    int i,j;
    int N_ = ob.N_;
    int kl_ = ob.kl_;

    for (i=0; i<N_; i++)
    {
      for (j=0; j<N_; j++)
        {
          if(((i>=j)&&(i-j<=kl_)))
            s << ob(i,j) << ' ';
          else if (((j>=i)&&(j-i<=kl_)))
            s << ob(j,i) << ' ';
        }
      s << "\n";
    }
  }
  return s;
}

double& LaSpdBandMatDouble::operator()(int i, int j)
{

#ifdef LA_BOUNDS_CHK
  if ((i<0||i>=N_)||(j<0||j>=N_))
   {
     cerr << "Index to SPD Banded Matrix out of range!\n";
     exit (1);
   }
#endif

  if (i>=j)
    if (i-j<=kl_)
        return data_(kl_+i-j,j);
    else
        return outofbounds_;

  else // (j>i)
    if (j-i<=kl_)
        return data_(kl_+j-i,i);
    else
        return outofbounds_;

}


double& LaSpdBandMatDouble::operator()(int i, int j) const
{

#ifdef LA_BOUNDS_CHK
  if ((i<0||i>=N_)||(j<0||j>=N_))
   {
     cerr << "Index to SPD Banded Matrix out of range!\n";
     exit (1);
   }
#endif

  if (i>=j)
    if (i-j<=kl_)
        return data_(kl_+i-j,j);
    else
        return outofbounds_;

  else  // (j>i)
    if (j-i<=kl_)
        return data_(kl_+j-i,i);
    else
        return outofbounds_;

}

