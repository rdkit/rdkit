//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.


#ifndef _LA_SYMM_BAND_MAT_DOUBLE_H_
#define _LA_SYMM_BAND_MAT_DOUBLE_H_

#include "lafnames.h"
//#include LA_GEN_MAT_DOUBLE_H // changed of VC++
#include "gmd.h"
//#include LA_VECTOR_DOUBLE_H // changed of VC++
#include "lavd.h"


class LaSymmBandMatDouble
{
  LaGenMatDouble data_;  // internal storage.

  int N_;       // N_ is (NxN)
  int kl_;      // kl_ = # subdiags
  static double outofbounds_; // out of range value returned.
  static int debug_;         // print debug info.
  static int *info_;         // print matrix info only, not values
                             //   originally 0, set to 1, and then
                             //   reset to 0 after use.


public:

  // constructors

  inline LaSymmBandMatDouble();
  inline LaSymmBandMatDouble(int,int);
  inline LaSymmBandMatDouble(const LaSymmBandMatDouble &);

  // operators

  LaSymmBandMatDouble& operator=(double);
  inline LaSymmBandMatDouble& operator=(const LaSymmBandMatDouble&);
  double& operator()(int,int);
  double& operator()(int,int) const;
  friend std::ostream& operator<<(std::ostream &, const LaSymmBandMatDouble &);


  // member functions

  inline int size(int) const;           // submatrix size
  inline int inc(int d) const;          // explicit increment
  inline int gdim(int d) const;         // global dimensions

  inline LaSymmBandMatDouble& ref(LaSymmBandMatDouble &);
  inline LaSymmBandMatDouble& copy(const LaSymmBandMatDouble &);
  inline double* addr() const {        // return address of matrix.
        return data_.addr();}
  inline int ref_count() const {        // return ref_count of matrix.
        return data_.ref_count();}
  inline LaIndex index(int d) const {     // return indices of matrix.
        return data_.index(d);}
  inline int subdiags() {     // return # of subdiags of matrix.
        return (kl_);}
  inline int subdiags() const {     // return # of subdiags of const matrix.
        return (kl_);}
  inline int shallow() const {      // return shallow flag.
        return data_.shallow();}
  inline int debug() const {    // return debug flag.
        return debug_;}
  inline int debug(int d) { // set debug flag for lagenmat.
        return debug_ = d;}

  inline LaSymmBandMatDouble& resize(const LaSymmBandMatDouble&);

  inline const LaSymmBandMatDouble& info() const {
        int *t = info_;
        *t = 1;
        return *this;};

  inline LaSymmBandMatDouble print_data() const 
    { std::cout << data_; return *this;}

  // destructor

  inline ~LaSymmBandMatDouble();
};

  // constructors 

inline LaSymmBandMatDouble::LaSymmBandMatDouble()
    :data_()
{

  N_ = kl_ = 0;
}

inline LaSymmBandMatDouble::LaSymmBandMatDouble(int n,int l)
    :data_(n,2*l+1)
{

  N_ = n;
  kl_ = l;
}

inline LaSymmBandMatDouble::LaSymmBandMatDouble(const LaSymmBandMatDouble &A)
{

  data_.copy(A.data_);
  N_ = A.N_;
  kl_ = A.kl_;
}

  // destructor 

inline LaSymmBandMatDouble::~LaSymmBandMatDouble()
{
}

  
  // member functions and operators

inline LaSymmBandMatDouble& LaSymmBandMatDouble::ref(LaSymmBandMatDouble &ob)
{

  data_.ref(ob.data_);
  N_ = ob.N_;
  kl_ = ob.kl_;

  return *this;
}

inline LaSymmBandMatDouble& LaSymmBandMatDouble::resize(const LaSymmBandMatDouble &ob)
{

  data_.resize(ob.data_);

  return *this;
}


inline LaSymmBandMatDouble& LaSymmBandMatDouble::operator=(const LaSymmBandMatDouble &B)
{
    data_ = B.data_;
    N_ = B.N_;
    kl_ = B.kl_;

    return *this;
}

inline int LaSymmBandMatDouble::size(int d) const
{
   return(data_.size(d));
}

inline int LaSymmBandMatDouble::inc(int d) const
{
   return(data_.inc(d));
}

inline int LaSymmBandMatDouble::gdim(int d) const
{
   return(data_.gdim(d));
}

#endif 
// _LA_SYMM_BAND_MAT_DOUBLE_H_

