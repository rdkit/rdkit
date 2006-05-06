//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.


#ifndef _LA_LOWER_TRIANG_MAT_DOUBLE_H_
#define _LA_LOWER_TRIANG_MAT_DOUBLE_H_

#ifndef _LA_GEN_MAT_DOUBLE_H_
//#include LA_GEN_MAT_DOUBLE_H // changed of VC++
#include "gmd.h"
#endif

//#define LOWER_INDEX_CHK

class LaLowerTriangMatDouble
{
  LaGenMatDouble data_;
  static double outofbounds_;
  static int debug_;         // print debug info. 
  static int *info_;         // print matrix info only, not values
                             //   originally 0, set to 1, and then
                             //   reset to 0 after use.
public:

  // constructors

  inline LaLowerTriangMatDouble();
  inline LaLowerTriangMatDouble(int,int);
  inline LaLowerTriangMatDouble(double*,int,int);
  inline LaLowerTriangMatDouble(const LaLowerTriangMatDouble &);

  // operators

  inline LaLowerTriangMatDouble& ref(LaLowerTriangMatDouble &);
  inline LaLowerTriangMatDouble& ref(LaGenMatDouble &);
  LaLowerTriangMatDouble& copy(const LaLowerTriangMatDouble &);
         LaLowerTriangMatDouble& operator=(double);
  inline LaLowerTriangMatDouble& operator=(const LaLowerTriangMatDouble &);
  inline double& operator()(int,int);
  inline double& operator()(int,int) const;

//  inline operator LaGenMatDouble();

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
  inline int debug(int d) { // set debug flag for lagenmat.
        return debug_ = d;}

  inline LaLowerTriangMatDouble& resize(const LaLowerTriangMatDouble&); 

  inline const LaLowerTriangMatDouble& info() const {
        int *t = info_;
        *t = 1;
        return *this;};


  friend std::ostream &operator<<(std::ostream &, const LaLowerTriangMatDouble &);

  // destructor

  inline ~LaLowerTriangMatDouble();
};

  // constructor functions

inline LaLowerTriangMatDouble::LaLowerTriangMatDouble() : data_()
{
    *info_ = 0;
}

inline LaLowerTriangMatDouble::LaLowerTriangMatDouble(int i,int j):data_(i,j)
{
    *info_ = 0;
}

inline LaLowerTriangMatDouble::LaLowerTriangMatDouble(double *d,int i,int j):
    data_(d,i,j)
{
    *info_ = 0;
}

inline LaLowerTriangMatDouble::LaLowerTriangMatDouble(const LaLowerTriangMatDouble &A)
{

  data_.copy(A.data_);
}

  
  // operator functions

inline LaLowerTriangMatDouble& LaLowerTriangMatDouble::ref(LaLowerTriangMatDouble &ob)
{

  data_.ref(ob.data_);

  return *this;
}
  
 
inline LaLowerTriangMatDouble& LaLowerTriangMatDouble::ref(LaGenMatDouble &ob)
{

  data_.ref(ob);

  return *this;
}


inline LaLowerTriangMatDouble& LaLowerTriangMatDouble::resize(const LaLowerTriangMatDouble &ob)
{

    data_.resize(ob.data_);

    return *this;
}


     
inline LaLowerTriangMatDouble& LaLowerTriangMatDouble::operator=(const LaLowerTriangMatDouble &L)
{

    data_ = L.data_;

    return *this;
}
    

inline double& LaLowerTriangMatDouble::operator()(int i, int j)
{


#ifdef LOWER_INDEX_CHK
    if (i<j)
    {
     std::cout << "Warning, index to Lower Triular matrix out of range!\n";
     std::cout << " i = " << i << " " <<" j = " << j << std::endl;
    }
#endif

    if (i<j)
        return outofbounds_;
    else
        return data_(i,j);
}


inline double& LaLowerTriangMatDouble::operator()(int i, int j) const
{

#ifdef LOWER_INDEX_CHK
    if (i<j)
    {
     std::cout << "Warning, index to Lower Triular matrix out of range!\n";
     std::cout << " i = " << i << " " <<" j = " << j << std::endl;
    }
#endif

    if (i<j)
        return outofbounds_;
    else
        return data_(i,j);
}


  // destructor function

inline LaLowerTriangMatDouble::~LaLowerTriangMatDouble()
{
}

inline int LaLowerTriangMatDouble::size(int d) const
{
   return(data_.size(d));
}

inline int LaLowerTriangMatDouble::inc(int d) const
{
   return(data_.inc(d));
}

inline int LaLowerTriangMatDouble::gdim(int d) const
{
   return(data_.gdim(d));
}


#if 0
// type conversions between LaGenMat and LaUpTriMat

inline LaLowerTriangMatDouble::operator LaGenMatDouble()
{
  LaGenMatDouble G;

  G.ref((*this).data_);

  return G;
}
#endif


#endif 
// _LA_LOWER_TRIANG_MAT_DOUBLE_H_
