//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.


#ifndef _LA_VECTOR_INT_H_
#define _LA_VECTOR_INT_H_

#include "lafnames.h"

#ifndef _LA_GEN_MAT_INT_H_
//#include LA_GEN_MAT_INT_H // changed of VC++
#include "gmi.h"
#endif


// a vector is simply an nx1 or 1xn, matrix, only that it can
// be constructed and accessed by a single dimension

class LaVectorInt: public LaGenMatInt
{
public:

    inline LaVectorInt();
    inline LaVectorInt(int);
    inline LaVectorInt(int, int);  
    inline LaVectorInt(int*, int);
    inline LaVectorInt(int*, int, int);
    inline LaVectorInt(const LaGenMatInt&);

    inline int size() const;
    inline int inc() const;
    inline LaIndex index() const;
    inline int start() const;
    inline int end() const;

    inline LaVectorInt& ref(const LaGenMatInt &);
    inline LaVectorInt& inject(const LaGenMatInt &);
    inline LaVectorInt& copy(const LaGenMatInt &);

    inline int& operator()(int i);
    inline int& operator()(int i) const ;
    inline LaVectorInt operator()(const LaIndex&);

    
    inline LaVectorInt& operator=(int);
    inline LaVectorInt& operator=(const LaGenMatInt&);

};

// NOTE: we default to column vectors, since matrices are column
//  oriented.

inline LaVectorInt::LaVectorInt() : LaGenMatInt(0,1) {}
inline LaVectorInt::LaVectorInt(int i) : LaGenMatInt(i,1) {}

// NOTE: one shouldn't be using this method to initalize, but
// it is here so that the constructor can be overloaded with 
// a runtime test.
//
inline LaVectorInt::LaVectorInt(int m, int n) : LaGenMatInt(m,n)
{
    assert(n==1 || m==1);
}

inline LaVectorInt::LaVectorInt(int *d, int n) : 
    LaGenMatInt(d,n,1) {}

inline LaVectorInt::LaVectorInt(int *d, int n, int m) : 
    LaGenMatInt(d,n,m) {}

inline LaVectorInt::LaVectorInt(const LaGenMatInt &G) 
{
        assert(G.size(0)==1 || G.size(1)==1);

        (*this).ref(G);
}


//note that vectors can be either stored columnwise, or row-wise

// this will handle the 0x0 case as well.

inline int LaVectorInt::size() const 
{ return LaGenMatInt::size(0)*LaGenMatInt::size(1); }

inline int& LaVectorInt::operator()(int i)
{ if (LaGenMatInt::size(0)==1 )
    return LaGenMatInt::operator()(0,i);
  else
    return LaGenMatInt::operator()(i,0);
}

inline int& LaVectorInt::operator()(int i) const
{ if (LaGenMatInt::size(0)==1 )
    return LaGenMatInt::operator()(0,i);
  else
    return LaGenMatInt::operator()(i,0);
}

inline LaVectorInt LaVectorInt::operator()(const LaIndex& I)
{ if (LaGenMatInt::size(0)==1)
    return LaGenMatInt::operator()(LaIndex(0,0),I); 
  else
    return LaGenMatInt::operator()(I,LaIndex(0,0)); 
}


inline LaVectorInt& LaVectorInt::copy(const LaGenMatInt &A)
{
    assert(A.size(0) == 1 || A.size(1) == 1);   //make sure rhs is a
                                                // a vector.
    LaGenMatInt::copy(A);
    return *this;
}

inline LaVectorInt& LaVectorInt::operator=(const  LaGenMatInt &A)
{
    return inject(A);
}

inline LaVectorInt& LaVectorInt::ref(const LaGenMatInt &A)
{
    assert(A.size(0) == 1 || A.size(1) == 1);
    LaGenMatInt::ref(A);
    return *this;
}

inline LaVectorInt& LaVectorInt::operator=(int d)
{
    LaGenMatInt::operator=(d);
    return *this;
}

inline LaVectorInt& LaVectorInt::inject(const LaGenMatInt &A)
{
    assert(A.size(0) == 1 || A.size(1) == 1);
    LaGenMatInt::inject(A);
    return *this;
}
    
inline int LaVectorInt::inc() const
{
  if (LaGenMatInt::size(1)==1 )
    return LaGenMatInt::inc(0);
  else
    return LaGenMatInt::inc(1);
}

inline LaIndex LaVectorInt::index() const
{
  if (LaGenMatInt::size(1)==1 )
    return LaGenMatInt::index(0);
  else
    return LaGenMatInt::index(1);
}

inline int LaVectorInt::start() const
{
  if (LaGenMatInt::size(1)==1 )
    return LaGenMatInt::start(0);
  else
    return LaGenMatInt::start(1);
}

inline int LaVectorInt::end() const
{
  if (LaGenMatInt::size(1)==1 )
    return LaGenMatInt::end(0);
  else
    return LaGenMatInt::end(1);
}

#endif 
// _LA_VECTOR_INT_H_
