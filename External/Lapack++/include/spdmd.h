//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.


#ifndef _LA_SPD_MAT_DOUBLE_H_
#define _LA_SPD_MAT_DOUBLE_H_

//#include LA_LOWER_TRIANG_MAT_DOUBLE_H // changed of VC++
#include "ltgmd.h"


class LaSpdMatDouble 
{
  LaLowerTriangMatDouble lower_data_;
  static int debug_;        // print out debug info.
  static int* info_;        // print out information on matrix.

public:

        /*::::::::::::::::::::::::::*/

        /* Constructors/Destructors */

        /*::::::::::::::::::::::::::*/


  inline LaSpdMatDouble();
  inline LaSpdMatDouble(int,int);
  inline LaSpdMatDouble(double*,int,int);
  inline LaSpdMatDouble(const LaSpdMatDouble &);
  inline ~LaSpdMatDouble();

  // operators

  inline LaSpdMatDouble& ref(LaSpdMatDouble &);
  LaSpdMatDouble& copy(const LaSpdMatDouble &);
  LaSpdMatDouble& operator=(const double &);
  inline LaSpdMatDouble& operator=(LaSpdMatDouble &);
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

  inline LaSpdMatDouble& resize(const LaSpdMatDouble&);

  inline const LaSpdMatDouble& info() const {
        int *t = info_;
        *t = 1;
        return *this;}

  //* I/O *//

  friend std::ostream& operator<<(std::ostream&, const LaSpdMatDouble&);

};


  // constructor functions

inline LaSpdMatDouble::LaSpdMatDouble() : lower_data_()
{
}

inline LaSpdMatDouble::LaSpdMatDouble(int i,int j) : lower_data_(i,j)
{
}

inline LaSpdMatDouble::LaSpdMatDouble(double *d,int i,int j):lower_data_(d,i,j)
{
}

inline LaSpdMatDouble::LaSpdMatDouble(const LaSpdMatDouble &S)
{
  lower_data_.copy(S.lower_data_);
}

  // destructor function

inline LaSpdMatDouble::~LaSpdMatDouble()
{
        // automatically calls the destructor for LaLowerTriangMatDouble
}

  // operator functions

inline double& LaSpdMatDouble::operator()(int i, int j)
{

  if (j>i)
     return (lower_data_(j,i));
  else
     return (lower_data_(i,j));

}

inline double& LaSpdMatDouble::operator()(int i, int j) const
{

  if (j>i)
     return (lower_data_(j,i));
  else
     return (lower_data_(i,j));

}

inline LaSpdMatDouble& LaSpdMatDouble::ref(LaSpdMatDouble &S)
{
  lower_data_.ref(S.lower_data_);

  return *this;
}

inline LaSpdMatDouble& LaSpdMatDouble::resize(const LaSpdMatDouble &S)
{
    lower_data_.resize(S.lower_data_);

    return *this;
}

inline LaSpdMatDouble& LaSpdMatDouble::operator=(LaSpdMatDouble &S)
{
        return copy(S);
}

inline int LaSpdMatDouble::size(int d) const
{
   return(lower_data_.size(d));
}

inline int LaSpdMatDouble::inc(int d) const
{
   return(lower_data_.inc(d));
}

inline int LaSpdMatDouble::gdim(int d) const
{
   return(lower_data_.gdim(d));
}

#endif 
// _LA_SPD_MAT_DOUBLE_H_
