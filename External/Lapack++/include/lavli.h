//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.


#ifndef _LA_VECTOR_LONG_INT_H_
#define _LA_VECTOR_LONG_INT_H_

#include "lafnames.h"

#ifndef _LA_GEN_MAT_LONG_INT_H_
//#include LA_GEN_MAT_LONG_INT_H // changed of VC++
#include "gmli.h"
#endif


// a vector is simply an nx1 or 1xn, matrix, only that it can
// be constructed and accessed by a single dimension

class LaVectorLongInt: public LaGenMatLongInt
{
public:

    inline LaVectorLongInt();
    inline LaVectorLongInt(int);
    inline LaVectorLongInt(int, int);  
    inline LaVectorLongInt(long int*, int);
    inline LaVectorLongInt(long int*, int, int);
    inline LaVectorLongInt(const LaGenMatLongInt&);

    inline int size() const;
    inline int inc() const;
    inline LaIndex index() const;
    inline int start() const;
    inline int end() const;

    inline LaVectorLongInt& ref(const LaGenMatLongInt &);
    inline LaVectorLongInt& inject(const LaGenMatLongInt &);
    inline LaVectorLongInt& copy(const LaGenMatLongInt &);

    inline long int& operator()(int i);
    inline const long int& operator()(int i) const ;
    inline LaVectorLongInt operator()(const LaIndex&);

    
    inline LaVectorLongInt& operator=(long int);
    inline LaVectorLongInt& operator=(const LaGenMatLongInt&);

};

// NOTE: we default to column vectors, since matrices are column
//  oriented.

inline LaVectorLongInt::LaVectorLongInt() : LaGenMatLongInt(0,1) {}
inline LaVectorLongInt::LaVectorLongInt(int i) : LaGenMatLongInt(i,1) {}

// NOTE: one shouldn't be using this method to initalize, but
// it is here so that the constructor can be overloaded with 
// a runtime test.
//
inline LaVectorLongInt::LaVectorLongInt(int m, int n) : LaGenMatLongInt(m,n)
{
    assert(n==1 || m==1);
}

inline LaVectorLongInt::LaVectorLongInt(long int *d, int m) : 
    LaGenMatLongInt(d,m,1) {}

#if 0
inline LaVectorLongInt::LaVectorLongInt(long int *d, int m, int n) : 
    LaGenMatLongInt(d,m,n) {}
#endif

inline LaVectorLongInt::LaVectorLongInt(const LaGenMatLongInt& G) : 
        LaGenMatLongInt(G)
{
        assert(G.size(0)==1 || G.size(1)==1);

}
        
//note that vectors can be either stored columnwise, or row-wise

// this will handle the 0x0 case as well.

inline int LaVectorLongInt::size() const 
{ return LaGenMatLongInt::size(0)*LaGenMatLongInt::size(1); }

inline long int& LaVectorLongInt::operator()(int i)
{ if (LaGenMatLongInt::size(0)==1 )
    return LaGenMatLongInt::operator()(0,i);
  else
    return LaGenMatLongInt::operator()(i,0);
}

inline const long int& LaVectorLongInt::operator()(int i) const
{ if (LaGenMatLongInt::size(0)==1 )
    return LaGenMatLongInt::operator()(0,i);
  else
    return LaGenMatLongInt::operator()(i,0);
}

inline LaVectorLongInt LaVectorLongInt::operator()(const LaIndex& I)
{ if (LaGenMatLongInt::size(0)==1)
    return LaGenMatLongInt::operator()(LaIndex(0,0),I); 
  else
    return LaGenMatLongInt::operator()(I,LaIndex(0,0)); 
}


inline LaVectorLongInt& LaVectorLongInt::copy(const LaGenMatLongInt &A)
{
    assert(A.size(0) == 1 || A.size(1) == 1);   //make sure rhs is a
                                                // a vector.
    LaGenMatLongInt::copy(A);
    return *this;
}

inline LaVectorLongInt& LaVectorLongInt::operator=(const LaGenMatLongInt &A)
{
    return copy(A);
}

inline LaVectorLongInt& LaVectorLongInt::ref(const LaGenMatLongInt &A)
{
    assert(A.size(0) == 1 || A.size(1) == 1);
    LaGenMatLongInt::ref(A);
    return *this;
}

inline LaVectorLongInt& LaVectorLongInt::operator=(long int d)
{
    LaGenMatLongInt::operator=(d);
    return *this;
}

inline LaVectorLongInt& LaVectorLongInt::inject(const LaGenMatLongInt &A)
{
    assert(A.size(0) == 1 || A.size(1) == 1);
    LaGenMatLongInt::inject(A);
    return *this;
}
    
inline int LaVectorLongInt::inc() const
{
  if (LaGenMatLongInt::size(1)==1 )
    return LaGenMatLongInt::inc(0);
  else
    return LaGenMatLongInt::inc(1);
}

inline LaIndex LaVectorLongInt::index() const
{
  if (LaGenMatLongInt::size(1)==1 )
    return LaGenMatLongInt::index(0);
  else
    return LaGenMatLongInt::index(1);
}

inline int LaVectorLongInt::start() const
{
  if (LaGenMatLongInt::size(1)==1 )
    return LaGenMatLongInt::start(0);
  else
    return LaGenMatLongInt::start(1);
}

inline int LaVectorLongInt::end() const
{
  if (LaGenMatLongInt::size(1)==1 )
    return LaGenMatLongInt::end(0);
  else
    return LaGenMatLongInt::end(1);
}

#endif 
    // _LA_VECTOR_LONG_INT_H_
