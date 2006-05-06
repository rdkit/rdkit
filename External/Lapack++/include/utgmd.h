//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.


#ifndef _LA_UPPER_TRIANG_MAT_DOUBLE_H_
#define _LA_UPPER_TRIANG_MAT_DOUBLE_H_

#include "lafnames.h"
//#include LA_GEN_MAT_DOUBLE_H // changed of VC++
#include "gmd.h"

//#define UPPER_INDEX_CHK

class LaUpperTriangMatDouble
{
  LaGenMatDouble data_;
  static double outofbounds_;
  static int debug_;         // print debug info. 
  static int *info_;         // print matrix info only, not values
                             //   originally 0, set to 1, and then
                             //   reset to 0 after use.
public:

  // constructors

  inline LaUpperTriangMatDouble();
  inline LaUpperTriangMatDouble(int,int);
  inline LaUpperTriangMatDouble(double*,int,int);
  inline LaUpperTriangMatDouble(const LaUpperTriangMatDouble &);

  // operators

  inline LaUpperTriangMatDouble& ref(LaUpperTriangMatDouble &);
  inline LaUpperTriangMatDouble& ref(LaGenMatDouble &);
  LaUpperTriangMatDouble& copy(LaUpperTriangMatDouble &);
        LaUpperTriangMatDouble& operator=(const double &);
  inline LaUpperTriangMatDouble& operator=(const LaUpperTriangMatDouble &);
  inline double& operator()(int,int);
  inline double& operator()(int,int) const;

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

  inline LaUpperTriangMatDouble& resize(const LaUpperTriangMatDouble&);

  inline const LaUpperTriangMatDouble& info() const {
        int *t = info_;
        *t = 1;
        return *this;};



  friend std::ostream &operator<<(std::ostream &, const LaUpperTriangMatDouble &);

  // destructor

  inline ~LaUpperTriangMatDouble();
};

  // constructor functions

inline LaUpperTriangMatDouble::LaUpperTriangMatDouble() : data_()
{
    *info_ = 0;
}

inline LaUpperTriangMatDouble::LaUpperTriangMatDouble(int i,int j):data_(i,j)
{
    *info_ = 0;
}

inline LaUpperTriangMatDouble::LaUpperTriangMatDouble(double *d,int i,int j):
    data_(d,i,j)
{
    *info_ = 0;
}

inline LaUpperTriangMatDouble::LaUpperTriangMatDouble(const LaUpperTriangMatDouble &A)
{

  data_.copy(A.data_);
}

  
  // operator functions

inline LaUpperTriangMatDouble& LaUpperTriangMatDouble::ref(LaUpperTriangMatDouble &ob)
{

  data_.ref(ob.data_);

  return *this;
}

inline LaUpperTriangMatDouble& LaUpperTriangMatDouble::ref(LaGenMatDouble &ob)
{

  data_.ref(ob);

  return *this;
}


inline LaUpperTriangMatDouble& LaUpperTriangMatDouble::resize(const LaUpperTriangMatDouble &ob)
{

    data_.resize(ob.data_);

    return *this;
}

  
     
inline LaUpperTriangMatDouble& LaUpperTriangMatDouble::operator=(const LaUpperTriangMatDouble &U)
{

    data_ = U.data_;
    
    return *this;
}


inline double& LaUpperTriangMatDouble::operator()(int i, int j)
{

#ifdef UPPER_INDEX_CHK
  if (j<i)
    {
     cout << "Warning, index to Upper Triular matrix out of range!\n";
     cout << " i = " << i << " " <<" j = " << j << endl;
    }
#endif

  if (j<i)
    return outofbounds_;
  else
    return data_(i,j);
}


inline double& LaUpperTriangMatDouble::operator()(int i, int j) const
{

#ifdef UPPER_INDEX_CHK
  if (j<i)
    {
     cout << "Warning, index to Upper Triular matrix out of range!\n";
     cout << " i = " << i << " " <<" j = " << j << endl;
    }
#endif

  if (j<i)
    return outofbounds_;
  else
    return data_(i,j);
}


  // destructor function

inline LaUpperTriangMatDouble::~LaUpperTriangMatDouble()
{
}

inline int LaUpperTriangMatDouble::size(int d) const
{
   return(data_.size(d));
}

inline int LaUpperTriangMatDouble::inc(int d) const
{
   return(data_.inc(d));
}

inline int LaUpperTriangMatDouble::gdim(int d) const
{
   return(data_.gdim(d));
}


// type conversions between LaGenMat and LaUpTriangMat

inline LaUpperTriangMatDouble::operator LaGenMatDouble()
{
  LaGenMatDouble G;

  G.ref((*this).data_);

  return G;
}

#endif 
// _LA_UPPER_TRIANG_MAT_DOUBLE_H_
