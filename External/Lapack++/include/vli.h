//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.
//
//      Lapack++ "Shared" Vector LongInt Class
//
//      A lightweight vector class with minimal overhead.
//
//      shallow assignment
//      unit stride
//      inlined access A(i)
//      optional (compile-time) array bounds checking through 
//              VECTOR_LONG_INT_BOUNDS_CHECK
//      A(i) is the same as A[i]
//      auto conversion to long int*
//      a null vector has size of 0, but has the ref_count structure
//              has been initalized
//

#ifndef _VECTOR_LONG_INT_H_
#define _VECTOR_LONG_INT_H_    

#include <iostream>       // for formatted printing of matrices
#include <fstream>

#ifndef __ASSERT_H
#include <assert.h>     // cheap "error" protection used in checking
#endif                  // checking array bounds.


typedef  struct {
    int        sz;                                        
    long int*    data;                                       
    int        ref_count;
} vrefLongInt;
                        


class VectorLongInt
{                                                                      
    private:                                                           
           vrefLongInt *p;
           long int *data;            // performance hack, avoid long int
                                    // indirection to data.
    public:                                                            
                                                                       
        /*::::::::::::::::::::::::::*/                                 
        /* Constructors/Destructors */                                 
        /*::::::::::::::::::::::::::*/                                 
                                                                       
    //inline VectorLongInt();     // this should behave as VectorLongInt(0)
    VectorLongInt(int);                             
    VectorLongInt(int, long int);   // can't be inlined because of 'for'
                                       // statement.
    VectorLongInt(long int*, int);
    VectorLongInt(const VectorLongInt&); 
    ~VectorLongInt() ;                              
                                                                       
        /*::::::::::::::::::::::::::::::::*/                           
        /*  Indices and access operations */                           
        /*::::::::::::::::::::::::::::::::*/                           
                                                                       
    inline long int&        operator[](int); 
    inline long int&        operator[](int) const;  // read only
    inline long int&        operator()(int); 
    inline long int&        operator()(int) const; // read only
    inline              operator    long int*(void); 
    inline int          size() const;
    inline int          null() const;
           int          resize(int d);
    inline int          ref_count() const;  // return the number of ref counts
    inline long int*        addr() const;
                                                                       
        /*::::::::::::::*/                                             
        /*  Assignment  */                                             
        /*::::::::::::::*/                                             
                                                                       
    inline  VectorLongInt& operator=(const VectorLongInt&);
            VectorLongInt& operator=(long int);
    inline  VectorLongInt& ref(const VectorLongInt &);
            VectorLongInt& inject(VectorLongInt&);
            VectorLongInt& copy(const VectorLongInt&);

    /* I/O */                                                      
    friend std::ostream&   operator<<(std::ostream&, const VectorLongInt&);       

};                                                                     


    // operators and member functions

inline int VectorLongInt::null()    const
{
    return (size() == 0) ;
}

inline int VectorLongInt::size() const
{
    return   p-> sz;
}


inline int VectorLongInt::ref_count() const
{
    return p->ref_count;
}

inline long int* VectorLongInt::addr() const
{
    return data;
}

inline VectorLongInt::operator long int*(void) 
{
    return data;
}


inline long int& VectorLongInt::operator()(int i)
{
#ifdef VECTOR_LONG_INT_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif 
    return data[i];
}

inline long int& VectorLongInt::operator()(int i) const
{
#ifdef VECTOR_LONG_INT_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif
    return data[i];
}

//  *CHANGE*  [] is the same as ()
inline long int& VectorLongInt::operator[](int i)
{
#ifdef VECTOR_LONG_INT_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif  
    return data[i];
}

//  *CHANGE*  [] is the same as ()
inline long int& VectorLongInt::operator[](int i) const
{
#ifdef VECTOR_LONG_INT_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif  
    return data[i];
}

inline VectorLongInt& VectorLongInt::ref(const VectorLongInt& m)
{
    // always check that the p field has been initialized.
    // Either the lhs or rhs could be a NULL VectorLongInt...
    
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


inline VectorLongInt& VectorLongInt::operator=(const VectorLongInt& m)
{

    return  ref(m);
}




#endif 
// _VECTOR_LONG_INT_H_

