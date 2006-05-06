//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.


#ifndef _LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H_
#define _LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H_

#include "lafnames.h"
//#include LA_GEN_MAT_DOUBLE_H // changed of VC++
#include "gmd.h"

//#define UNIT_LOWER_INDEX_CHK

class LaUnitLowerTriangMatDouble
{
  LaGenMatDouble data_;
  static double outofbounds_;
  static int debug_;         // print debug info. 
  static int *info_;         // print matrix info only, not values
                             //   originally 0, set to 1, and then
                             //   reset to 0 after use.
public:

  // constructors

  inline LaUnitLowerTriangMatDouble();
  inline LaUnitLowerTriangMatDouble(int,int);
  inline LaUnitLowerTriangMatDouble(double*,int,int);
  inline LaUnitLowerTriangMatDouble(const LaUnitLowerTriangMatDouble &);

  // operators

  inline LaUnitLowerTriangMatDouble& ref(LaUnitLowerTriangMatDouble &);
  inline LaUnitLowerTriangMatDouble& ref(LaGenMatDouble &);
  LaUnitLowerTriangMatDouble& copy(LaUnitLowerTriangMatDouble &);
         LaUnitLowerTriangMatDouble& operator=(double );
  inline LaUnitLowerTriangMatDouble& operator=(const LaUnitLowerTriangMatDouble &);
  inline double& operator()(int,int);
  inline double& operator()(int,int) const;

  inline operator LaGenMatDouble(void);

  inline int size(int) const;           // submatrix size
  inline int inc(int d) const;          // explicit increment
  inline int gdim(int d) const;         // global dimensions
  inline double* addr() const {        // return address of data.
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
  
  inline LaUnitLowerTriangMatDouble& resize(const LaUnitLowerTriangMatDouble&);

  inline const LaUnitLowerTriangMatDouble& info() const {
        int *t = info_;
        *t = 1;
        return *this;};



		friend std::ostream &operator<<(std::ostream &, const LaUnitLowerTriangMatDouble &);

  // destructor

  inline ~LaUnitLowerTriangMatDouble();
};

  // constructor functions

inline LaUnitLowerTriangMatDouble::LaUnitLowerTriangMatDouble() : data_()
{
    *info_ = 0;
}

inline LaUnitLowerTriangMatDouble::LaUnitLowerTriangMatDouble(int i,int j):
    data_(i,j)
{
    *info_ = 0;
}

inline LaUnitLowerTriangMatDouble::LaUnitLowerTriangMatDouble(double *d,int i,int j):data_(d,i,j)
{
    *info_ = 0;
}

inline LaUnitLowerTriangMatDouble::LaUnitLowerTriangMatDouble(const LaUnitLowerTriangMatDouble &A)
{

  data_.copy(A.data_);
}

  
  // operator functions

inline LaUnitLowerTriangMatDouble& LaUnitLowerTriangMatDouble::ref(LaUnitLowerTriangMatDouble &ob)
{

  data_.ref(ob.data_);

  return *this;
}
  
inline LaUnitLowerTriangMatDouble& LaUnitLowerTriangMatDouble::ref(LaGenMatDouble &ob)
{

  data_.ref(ob);

  return *this;
}

 
inline LaUnitLowerTriangMatDouble& LaUnitLowerTriangMatDouble::resize(const LaUnitLowerTriangMatDouble &ob)
{

  data_.resize(ob.data_);

  return *this;
}

     
inline LaUnitLowerTriangMatDouble& LaUnitLowerTriangMatDouble::operator=(const LaUnitLowerTriangMatDouble &L)
{

    data_ = L.data_;

    return *this;
}
     

inline double& LaUnitLowerTriangMatDouble::operator()(int i, int j)
{

#ifdef UNIT_LOWER_INDEX_CHK
  if (j>=i)
   { 
     cout << "Warning, index to Lower Triular matrix out of range!\n";
     cout << " i = " << i << " " <<" j = " << j << endl;
   }
#endif

  if ((j==0)&&(i==0)) // special case, allows us to access beginning of
    return data_(0,0);// matrix without getting default outofbounds_.
  else if (j>=i)
    return outofbounds_;
  else
    return data_(i,j);
}


inline double& LaUnitLowerTriangMatDouble::operator()(int i, int j) const
{

#ifdef UNIT_LOWER_INDEX_CHK
  if (j>=i)
   {
     cout << "Warning, index to Lower Triangular matrix out of range!\n";
     cout << " i = " << i << " " <<" j = " << j << endl;
   }
#endif

  if ((j==0)&&(i==0)) // special case, allows us to access beginning of
    return data_(0,0);// matrix without getting default outofbounds_.
  else if (j>=i)
    return outofbounds_;
  else
    return data_(i,j);
}


  // destructor function

inline LaUnitLowerTriangMatDouble::~LaUnitLowerTriangMatDouble()
{
}

inline int LaUnitLowerTriangMatDouble::size(int d) const
{
   return(data_.size(d));
}

inline int LaUnitLowerTriangMatDouble::inc(int d) const
{
   return(data_.inc(d));
}

inline int LaUnitLowerTriangMatDouble::gdim(int d) const
{
   return(data_.gdim(d));
}


// type conversions between LaGenMat and LaUpTriMat

inline LaUnitLowerTriangMatDouble::operator LaGenMatDouble()
{
  LaGenMatDouble G;

  G.ref((*this).data_);

  return G;
}

#endif 
// _LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H_
