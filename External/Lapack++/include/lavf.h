//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.


#ifndef _LA_VECTOR_DOUBLE_H_
#define _LA_VECTOR_DOUBLE_H_

#include "conversion.h"

#ifndef _LA_GEN_MAT_DOUBLE_H_
//#include LA_GEN_MAT_DOUBLE_H // changed of VC++
#include "gmd.h"
#endif


// a vector is simply an nx1 or 1xn, matrix, only that it can
// be constructed and accessed by a single dimension

class LaVectorDouble: public LaGenMatDouble
{
public:

    inline LaVectorDouble();
    inline LaVectorDouble(int);
    inline LaVectorDouble(int, int);  
    inline LaVectorDouble(double*, int);
    inline LaVectorDouble(double*, int, int);
    inline LaVectorDouble(const LaGenMatDouble&);

    inline int size() const;
    inline inc() const;
    inline LaIndex index() const;
    inline int start() const;
    inline int end() const;

    inline LaVectorDouble& ref(const LaGenMatDouble &);
    inline LaVectorDouble& inject(const LaGenMatDouble &);
    inline LaVectorDouble& copy(const LaGenMatDouble &);

    inline double& operator()(int i);
    inline const double& operator()(int i) const ;
    inline LaVectorDouble operator()(const LaIndex&);

    
    inline LaVectorDouble& operator=(double);
    inline LaVectorDouble& operator=(const LaGenMatDouble&);

};

// NOTE: we default to column vectors, since matrices are column
//  oriented.

inline LaVectorDouble::LaVectorDouble() : LaGenMatDouble(0,1) {}
inline LaVectorDouble::LaVectorDouble(int i) : LaGenMatDouble(i,1) {}

// NOTE: one shouldn't be using this method to initalize, but
// it is here so that the constructor can be overloaded with 
// a runtime test.
//
inline LaVectorDouble::LaVectorDouble(int m, int n) : LaGenMatDouble(m,n)
{
    assert(n==1 || m==1);
}

inline LaVectorDouble::LaVectorDouble(double *d, int m) : 
    LaGenMatDouble(d,m,1) {}

inline LaVectorDouble::LaVectorDouble(double *d, int m, int n) : 
    LaGenMatDouble(d,m,n) {}

inline LaVectorDouble::LaVectorDouble(const LaGenMatDouble& G)
{
        assert(G.size(0)==1 || G.size(1)==1);

        (*this).ref(G);
}
        
//note that vectors can be either stored columnwise, or row-wise

// this will handle the 0x0 case as well.

inline int LaVectorDouble::size() const 
{ return LaGenMatDouble::size(0)*LaGenMatDouble::size(1); }

inline double& LaVectorDouble::operator()(int i)
{ if (LaGenMatDouble::size(0)==1 )
    return LaGenMatDouble::operator()(0,i);
  else
    return LaGenMatDouble::operator()(i,0);
}

inline const double& LaVectorDouble::operator()(int i) const
{ if (LaGenMatDouble::size(0)==1 )
    return LaGenMatDouble::operator()(0,i);
  else
    return LaGenMatDouble::operator()(i,0);
}

inline LaVectorDouble LaVectorDouble::operator()(const LaIndex& I)
{ if (LaGenMatDouble::size(0)==1)
    return LaGenMatDouble::operator()(LaIndex(0,0),I); 
  else
    return LaGenMatDouble::operator()(I,LaIndex(0,0)); 
}


inline LaVectorDouble& LaVectorDouble::copy(const LaGenMatDouble &A)
{
    assert(A.size(0) == 1 || A.size(1) == 1);   //make sure rhs is a
                                                // a vector.
    LaGenMatDouble::copy(A);
    return *this;
}


inline LaVectorDouble& LaVectorDouble::operator=(const LaGenMatDouble &A)
{
    return inject(A);
}

inline LaVectorDouble& LaVectorDouble::ref(const LaGenMatDouble &A)
{
    assert(A.size(0) == 1 || A.size(1) == 1);
    LaGenMatDouble::ref(A);
    return *this;
}

inline LaVectorDouble& LaVectorDouble::operator=(double d)
{
    LaGenMatDouble::operator=(d);
    return *this;
}

inline LaVectorDouble& LaVectorDouble::inject(const LaGenMatDouble &A)
{
    assert(A.size(0) == 1 || A.size(1) == 1);
    LaGenMatDouble::inject(A);
    return *this;
}
    
inline int LaVectorDouble::inc() const
{
  if (LaGenMatDouble::size(1)==1 )
    return LaGenMatDouble::inc(0);
  else
    return LaGenMatDouble::inc(1);
}

inline LaIndex LaVectorDouble::index() const
{
  if (LaGenMatDouble::size(1)==1 )
    return LaGenMatDouble::index(0);
  else
    return LaGenMatDouble::index(1);
}

inline int LaVectorDouble::start() const
{
  if (LaGenMatDouble::size(1)==1 )
    return LaGenMatDouble::start(0);
  else
    return LaGenMatDouble::start(1);
}

inline int LaVectorDouble::end() const
{
  if (LaGenMatDouble::size(1)==1 )
    return LaGenMatDouble::end(0);
  else
    return LaGenMatDouble::end(1);
}

#endif 
// _LA_VECTOR_DOUBLE_H_
