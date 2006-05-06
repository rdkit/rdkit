//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.
//
//      Lapack++ "Shared" Vector Complex Class
//
//      A lightweight vector class with minimal overhead.
//
//      shallow assignment
//      unit stride
//      inlined access A(i)
//      optional (compile-time) array bounds checking through 
//              VECTOR_COMPLEX_BOUNDS_CHECK
//      A(i) is the same as A[i]
//      auto conversion to complex*
//      a null vector has size of 0, but has the ref_count structure
//              has been initalized
//

#ifndef _VECTOR_COMPLEX_H_
#define _VECTOR_COMPLEX_H_    

#include <iostream>       // for formatted printing of matrices

#include <complex>
//#include <math.h>

#ifndef __ASSERT_H
#include <assert.h>     // cheap "error" protection used in checking
#endif                  // checking array bounds.

#ifdef WIN32
typedef std::complex<double> COMPLEX;
#endif

typedef  struct {
    int        sz;                                        
    COMPLEX *data;                                       
    int        ref_count;
} vrefComplex;
                        


class VectorComplex
{                                                                      
    private:                                                           
           vrefComplex *p;
           COMPLEX *data;            // performance hack, avoid COMPLEX
                                    // indirection to data.
    public:                                                            
                                                                       
        /*::::::::::::::::::::::::::*/                                 
        /* Constructors/Destructors */                                 
        /*::::::::::::::::::::::::::*/                                 
                                                                       
    //inline VectorComplex();     // this should behave as VectorComplex(0)
    VectorComplex(int);                             
    VectorComplex(int, COMPLEX);   // can't be inlined because of 'for'
                                       // statement.
    VectorComplex(COMPLEX*, int);
    VectorComplex(const VectorComplex&); 
    ~VectorComplex() ;                              
                                                                       
        /*::::::::::::::::::::::::::::::::*/                           
        /*  Indices and access operations */                           
        /*::::::::::::::::::::::::::::::::*/                           
                                                                       
    inline COMPLEX&     operator[](int); 
    inline COMPLEX&     operator[](int) const;  // read only
    inline COMPLEX&     operator()(int); 
    inline COMPLEX&     operator()(int) const; // read only
    inline              operator    COMPLEX*(); 
    inline int          size() const;
    inline int          null() const;
           int          resize(int d);
    inline int          ref_count() const;  // return the number of ref counts
    inline COMPLEX*     addr() const;
                                                                       
        /*::::::::::::::*/                                             
        /*  Assignment  */                                             
        /*::::::::::::::*/                                             
                                                                       
    inline  VectorComplex& operator=(const VectorComplex&);
            VectorComplex& operator=(COMPLEX);
    inline  VectorComplex& ref(const VectorComplex &);
            VectorComplex& inject(VectorComplex&);
            VectorComplex& copy(const VectorComplex&);

    /* I/O */                                                      
    friend std::ostream&   operator<<(std::ostream&, const VectorComplex&);       

};                                                                     


    // operators and member functions

inline int VectorComplex::null()    const
{
    return (size() == 0) ;
}

inline int VectorComplex::size() const
{
    return   p-> sz;
}


inline int VectorComplex::ref_count() const
{
    return p->ref_count;
}

inline COMPLEX* VectorComplex::addr() const
{
    return data;
}

inline VectorComplex::operator COMPLEX*() 
{
    return data;
}


inline COMPLEX& VectorComplex::operator()(int i)
{
#ifdef VECTOR_COMPLEX_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif 
    return data[i];
}

inline COMPLEX& VectorComplex::operator()(int i) const
{
#ifdef VECTOR_COMPLEX_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif
    return data[i];
}

//  [] *always* performs bounds-check 
//  *CHANGE*  [] is the same as ()
inline COMPLEX& VectorComplex::operator[](int i)
{
#ifdef VECTOR_COMPLEX_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif  
    return data[i];
}

//  [] *always* performs bounds-check 
//  *CHANGE*  [] is the same as ()
inline COMPLEX& VectorComplex::operator[](int i) const
{
#ifdef VECTOR_COMPLEX_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif  
    return data[i];
}

inline VectorComplex& VectorComplex::ref(const VectorComplex& m)
{
    // always check that the p field has been initialized.
    // Either the lhs or rhs could be a NULL VectorComplex...
    
        if (&m != this)         // not really necessary...
        {
            m.p->ref_count++;
            if (--(p->ref_count) == 0)              // perform garbage col.
            {
                delete [] ( p->data);
                delete p;
            }
            p = m.p;
            data = m.data;
        }
        return *this;
}


inline VectorComplex& VectorComplex::operator=(const VectorComplex& m)
{

    return  ref(m);
}




#endif 
// _VECTOR_COMPLEX_H_

