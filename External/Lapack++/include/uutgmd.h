//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.


#ifndef _LA_UNIT_UPPER_TRIANG_MAT_DOUBLE_H_
#define _LA_UNIT_UPPER_TRIANG_MAT_DOUBLE_H_

#include "lafnames.h"
//#include LA_GEN_MAT_DOUBLE_H // changed of VC++
#include "gmd.h"


//#define UNIT_UPPER_INDEX_CHK

class LaUnitUpperTriangMatDouble
{
  LaGenMatDouble data_;
  static double outofbounds_;
  static int debug_;         // print debug info. 
  static int *info_;         // print matrix info only, not values
                             //   originally 0, set to 1, and then
                             //   reset to 0 after use.
public:

  // constructors

  inline LaUnitUpperTriangMatDouble();
  inline LaUnitUpperTriangMatDouble(int,int);
  inline LaUnitUpperTriangMatDouble(double*,int,int);
  inline LaUnitUpperTriangMatDouble(const LaUnitUpperTriangMatDouble &);

  // operators

  inline LaUnitUpperTriangMatDouble& ref(LaUnitUpperTriangMatDouble &);
  inline LaUnitUpperTriangMatDouble& ref(LaGenMatDouble &);
  LaUnitUpperTriangMatDouble& copy(LaUnitUpperTriangMatDouble &);
        LaUnitUpperTriangMatDouble& operator=(const double &);
  inline LaUnitUpperTriangMatDouble& operator=(const LaUnitUpperTriangMatDouble &);
  double& operator()(int,int);
  double& operator()(int,int) const;

  inline operator LaGenMatDouble();

  inline int size(int) const;           // submatrix size
  inline int inc(int d) const;          // explicit increment
  inline int gdim(int d) const;         // global dimensions
  inline double* addr() const {        // return address of matrix.
        return data_.addr();}
  inline int ref_count() const {        // return ref_count of matrix.
        return data_.ref_count();}
  inline LaIndex index(int d) const {     // return indices of matrix.
        return data_.index(d);}
  inline int shallow() const {      // return indices of matrix.
        return data_.shallow();}
  inline int debug() const {    // return debug flag.
        return debug_;}
  inline int debug(int d) { // set debug flag.
        return debug_ = d;}

  inline LaUnitUpperTriangMatDouble& resize(const LaUnitUpperTriangMatDouble&);

  inline const LaUnitUpperTriangMatDouble& info() const {
        int *t = info_;
        *t = 1;
        return *this;};



  friend std:: ostream &operator<<(std::ostream &, const LaUnitUpperTriangMatDouble &);

  // destructor

  inline ~LaUnitUpperTriangMatDouble();
};

  // constructor functions

inline LaUnitUpperTriangMatDouble::LaUnitUpperTriangMatDouble() : data_()
{
    *info_ = 0;
}

inline LaUnitUpperTriangMatDouble::LaUnitUpperTriangMatDouble(int i,int j):
    data_(i,j)
{
    *info_ = 0;
}

inline LaUnitUpperTriangMatDouble::LaUnitUpperTriangMatDouble(double *d,int i,int j):data_(d,i,j)
{
    *info_ = 0;
}

inline LaUnitUpperTriangMatDouble::LaUnitUpperTriangMatDouble(const LaUnitUpperTriangMatDouble &A)
{

  data_.copy(A.data_);
}

  
  // operator functions

inline LaUnitUpperTriangMatDouble& LaUnitUpperTriangMatDouble::ref(LaUnitUpperTriangMatDouble &ob)
{

  data_.ref(ob.data_);

  return *this;
}
  
inline LaUnitUpperTriangMatDouble& LaUnitUpperTriangMatDouble::ref(LaGenMatDouble &ob)
{

  data_.ref(ob);

  return *this;
}
 
inline LaUnitUpperTriangMatDouble& LaUnitUpperTriangMatDouble::resize(const LaUnitUpperTriangMatDouble &ob)
{

  data_.resize(ob.data_);

  return *this;
}


     
inline LaUnitUpperTriangMatDouble& LaUnitUpperTriangMatDouble::operator=(const LaUnitUpperTriangMatDouble &U)
{

    data_ = U.data_;

    return *this;
}



  // destructor function

inline LaUnitUpperTriangMatDouble::~LaUnitUpperTriangMatDouble()
{
}

inline int LaUnitUpperTriangMatDouble::size(int d) const
{
   return(data_.size(d));
}

inline int LaUnitUpperTriangMatDouble::inc(int d) const
{
   return(data_.inc(d));
}

inline int LaUnitUpperTriangMatDouble::gdim(int d) const
{
   return(data_.gdim(d));
}


// type conversions between LaGenMat and LaUnitUpTriMat

inline LaUnitUpperTriangMatDouble::operator LaGenMatDouble()
{
  LaGenMatDouble G;

  G.ref((*this).data_);

  return G;
}

#endif 
// _LA_UNIT_UPPER_TRIANG_MAT_DOUBLE_H_
