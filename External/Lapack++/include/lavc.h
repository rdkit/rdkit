//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.


#ifndef _LA_VECTOR_COMPLEX_H_
#define _LA_VECTOR_COMPLEX_H_

#include "lafnames.h"

#ifndef _LA_GEN_MAT_COMPLEX_H_
//#include LA_GEN_MAT_COMPLEX_H // changed of VC++
#include "gmc.h"
#endif


// a vector is simply an nx1 or 1xn, matrix, only that it can
// be constructed and accessed by a single dimension

class LaVectorComplex: public LaGenMatComplex
{
public:

    inline LaVectorComplex();
    inline LaVectorComplex(int);
    inline LaVectorComplex(int, int);  
    inline LaVectorComplex(COMPLEX*, int);
    inline LaVectorComplex(COMPLEX*, int, int);
    inline LaVectorComplex(const LaGenMatComplex&);

    inline int size() const;
    inline inc() const;
    inline LaIndex index() const;
    inline int start() const;
    inline int end() const;

    inline LaVectorComplex& ref(const LaGenMatComplex &);
    inline LaVectorComplex& inject(const LaGenMatComplex &);
    inline LaVectorComplex& copy(const LaGenMatComplex &);

    inline COMPLEX& operator()(int i);
    inline const COMPLEX& operator()(int i) const ;
    inline LaVectorComplex operator()(const LaIndex&);

    
    inline LaVectorComplex& operator=(COMPLEX);
    inline LaVectorComplex& operator=(const LaGenMatComplex&);

};

// NOTE: we default to column vectors, since matrices are column
//  oriented.

inline LaVectorComplex::LaVectorComplex() : LaGenMatComplex(0,1) {}
inline LaVectorComplex::LaVectorComplex(int i) : LaGenMatComplex(i,1) {}

// NOTE: one shouldn't be using this method to initalize, but
// it is here so that the constructor can be overloaded with 
// a runtime test.
//
inline LaVectorComplex::LaVectorComplex(int m, int n) : LaGenMatComplex(m,n)
{
    assert(n==1 || m==1);
}

inline LaVectorComplex::LaVectorComplex(COMPLEX *d, int m) : 
    LaGenMatComplex(d,m,1) {}

inline LaVectorComplex::LaVectorComplex(COMPLEX *d, int m, int n) : 
    LaGenMatComplex(d,m,n) {}

inline LaVectorComplex::LaVectorComplex(const LaGenMatComplex& G)
{
        assert(G.size(0)==1 || G.size(1)==1);

        (*this).ref(G);
}
        
//note that vectors can be either stored columnwise, or row-wise

// this will handle the 0x0 case as well.

inline int LaVectorComplex::size() const 
{ return LaGenMatComplex::size(0)*LaGenMatComplex::size(1); }

inline COMPLEX& LaVectorComplex::operator()(int i)
{ if (LaGenMatComplex::size(0)==1 )
    return LaGenMatComplex::operator()(0,i);
  else
    return LaGenMatComplex::operator()(i,0);
}

inline const COMPLEX& LaVectorComplex::operator()(int i) const
{ if (LaGenMatComplex::size(0)==1 )
    return LaGenMatComplex::operator()(0,i);
  else
    return LaGenMatComplex::operator()(i,0);
}

inline LaVectorComplex LaVectorComplex::operator()(const LaIndex& I)
{ if (LaGenMatComplex::size(0)==1)
    return LaGenMatComplex::operator()(LaIndex(0,0),I); 
  else
    return LaGenMatComplex::operator()(I,LaIndex(0,0)); 
}


inline LaVectorComplex& LaVectorComplex::copy(const LaGenMatComplex &A)
{
    assert(A.size(0) == 1 || A.size(1) == 1);   //make sure rhs is a
                                                // a vector.
    LaGenMatComplex::copy(A);
    return *this;
}


inline LaVectorComplex& LaVectorComplex::operator=(const LaGenMatComplex &A)
{
    return inject(A);
}

inline LaVectorComplex& LaVectorComplex::ref(const LaGenMatComplex &A)
{
    assert(A.size(0) == 1 || A.size(1) == 1);
    LaGenMatComplex::ref(A);
    return *this;
}

inline LaVectorComplex& LaVectorComplex::operator=(COMPLEX d)
{
    LaGenMatComplex::operator=(d);
    return *this;
}

inline LaVectorComplex& LaVectorComplex::inject(const LaGenMatComplex &A)
{
    assert(A.size(0) == 1 || A.size(1) == 1);
    LaGenMatComplex::inject(A);
    return *this;
}
    
inline int LaVectorComplex::inc() const
{
  if (LaGenMatComplex::size(1)==1 )
    return LaGenMatComplex::inc(0);
  else
    return LaGenMatComplex::inc(1);
}

inline LaIndex LaVectorComplex::index() const
{
  if (LaGenMatComplex::size(1)==1 )
    return LaGenMatComplex::index(0);
  else
    return LaGenMatComplex::index(1);
}

inline int LaVectorComplex::start() const
{
  if (LaGenMatComplex::size(1)==1 )
    return LaGenMatComplex::start(0);
  else
    return LaGenMatComplex::start(1);
}

inline int LaVectorComplex::end() const
{
  if (LaGenMatComplex::size(1)==1 )
    return LaGenMatComplex::end(0);
  else
    return LaGenMatComplex::end(1);
}

#endif 
    // _LA_VECTOR_COMPLEX_H_
