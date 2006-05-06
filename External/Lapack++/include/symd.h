//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.


#ifndef _LA_SYMM_MAT_DOUBLE_H_
#define _LA_SYMM_MAT_DOUBLE_H_

//#include LA_LOWER_TRIANG_MAT_DOUBLE_H // changed of VC++
#include "ltgmd.h"

class LaSymmMatDouble
{
  LaLowerTriangMatDouble lower_data_;
  static int debug_;        // print out debug info.
  static int* info_;        // print out information on matrix.

public:

        /*::::::::::::::::::::::::::*/

        /* Constructors/Destructors */

        /*::::::::::::::::::::::::::*/


  inline LaSymmMatDouble();
  inline LaSymmMatDouble(int,int);
  inline LaSymmMatDouble(double*,int,int);
  inline LaSymmMatDouble(const LaSymmMatDouble &);
  inline ~LaSymmMatDouble();

  // operators

  inline LaSymmMatDouble& ref(LaSymmMatDouble &);
  LaSymmMatDouble& copy(const LaSymmMatDouble &);
  LaSymmMatDouble& operator=(const double &);
  inline LaSymmMatDouble& operator=(const LaSymmMatDouble &);
  inline double& operator()(int,int);
  inline double& operator()(int,int) const;

  operator LaGenMatDouble();
  operator LaLowerTriangMatDouble();

  inline int size(int) const ;           // submatrix size
  inline int inc(int d) const ;          // explicit increment
  inline int gdim(int d) const ;         // global dimensions
  inline double* addr() const {        // return address of data.
        return lower_data_.addr();}
  inline int ref_count() const {        // return ref_count of matrix.
        return lower_data_.ref_count();}
  inline LaIndex index(int d) const {     // return indices of matrix.
        return lower_data_.index(d);}
  inline int shallow() const {      // return indices of matrix.
        return lower_data_.shallow();}

  inline int debug() const {    // return debug flag.
        return debug_;}
  inline int debug(int d) { // set debug flag.
        return debug_ = d;}

  inline LaSymmMatDouble& resize(const LaSymmMatDouble&);

  inline const LaSymmMatDouble& info() const {
        int *t = info_;
        *t = 1;
        return *this;}

  //* I/O *//

  friend std::ostream& operator<<(std::ostream&, const LaSymmMatDouble&);

};


  // constructor functions

inline LaSymmMatDouble::LaSymmMatDouble() : lower_data_()
{
}

inline LaSymmMatDouble::LaSymmMatDouble(int i,int j) : lower_data_(i,j)
{
}

inline LaSymmMatDouble::LaSymmMatDouble(double *d,int i,int j):lower_data_(d,i,j)
{
}

inline LaSymmMatDouble::LaSymmMatDouble(const LaSymmMatDouble &S)
{
  lower_data_.copy(S.lower_data_);
}

  // destructor function

inline LaSymmMatDouble::~LaSymmMatDouble()
{
        // automatically calls the destructor for LaLowerTriangMatDouble
}

  // operator functions

inline double& LaSymmMatDouble::operator()(int i, int j)
{

  if (j>i)
     return (lower_data_(j,i));
  else
     return (lower_data_(i,j));

}

inline double& LaSymmMatDouble::operator()(int i, int j) const
{

  if (j>i)
     return (lower_data_(j,i));
  else
     return (lower_data_(i,j));

}

inline LaSymmMatDouble& LaSymmMatDouble::ref(LaSymmMatDouble &S)
{
  lower_data_.ref(S.lower_data_);

  return *this;
}

inline LaSymmMatDouble& LaSymmMatDouble::resize(const LaSymmMatDouble &S)
{
    lower_data_.resize(S.lower_data_);

    return *this;
}

inline LaSymmMatDouble& LaSymmMatDouble::operator=(const LaSymmMatDouble &S)
{
        return copy(S);
}

inline int LaSymmMatDouble::size(int d) const
{
   return(lower_data_.size(d));
}

inline int LaSymmMatDouble::inc(int d) const
{
   return(lower_data_.inc(d));
}

inline int LaSymmMatDouble::gdim(int d) const
{
   return(lower_data_.gdim(d));
}

#endif
