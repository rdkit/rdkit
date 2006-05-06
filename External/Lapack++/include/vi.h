//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.
//
//      Lapack++ "Shared" Vector Int Class
//
//      A lightweight vector class with minimal overhead.
//
//      shallow assignment
//      unit stride
//      inlined access A(i)
//      optional (compile-time) array bounds checking through 
//              VECTOR_INT_BOUNDS_CHECK
//      A(i) is the same as A[i]
//      auto conversion to int*
//      a null vector has size of 0, but has the ref_count structure
//              has been initalized
//

#ifndef _VECTOR_INT_H_
#define _VECTOR_INT_H_    

#include <iostream>       // for formatted printing of matrices

#ifndef __ASSERT_H
#include <assert.h>     // cheap "error" protection used in checking
#endif                  // checking array bounds.


typedef  struct {
    int        sz;                                        
    int*    data;                                       
    int        ref_count;
} vrefInt;
                        


class VectorInt
{                                                                      
    private:                                                           
           vrefInt *p;
           int *data;            // performance hack, avoid int
                                    // indirection to data.
    public:                                                            
                                                                       
        /*::::::::::::::::::::::::::*/                                 
        /* Constructors/Destructors */                                 
        /*::::::::::::::::::::::::::*/                                 
                                                                       
    //inline VectorInt();     // this should behave as VectorInt(0)
    VectorInt(int);                             
    VectorInt(int, int);   // can't be inlined because of 'for'
                                       // statement.
    VectorInt(int*, int);
    VectorInt(const VectorInt&); 
    ~VectorInt() ;                              
                                                                       
        /*::::::::::::::::::::::::::::::::*/                           
        /*  Indices and access operations */                           
        /*::::::::::::::::::::::::::::::::*/                           
                                                                       
    inline int&     operator[](int); 
    inline int&     operator[](int) const; // read only
    inline int&     operator()(int); 
    inline int&     operator()(int) const; // read only
    inline              operator    int*(); 
    inline int          size() const;
    inline int          null() const;
           int          resize(int d);
    inline int          ref_count() const;  // return the number of ref counts
    inline int*     addr() const;
                                                                       
        /*::::::::::::::*/                                             
        /*  Assignment  */                                             
        /*::::::::::::::*/                                             
                                                                       
    inline  VectorInt& operator=(const VectorInt&);
            VectorInt& operator=(int);
    inline  VectorInt& ref(const VectorInt &);
            VectorInt& inject(VectorInt&);
            VectorInt& copy(const VectorInt&);

    /* I/O */                                                      
    friend std::ostream&   operator<<(std::ostream&, const VectorInt&);       

};                                                                     


    // operators and member functions

inline int VectorInt::null()    const
{
    return (size() == 0) ;
}

inline int VectorInt::size() const
{
    return   p-> sz;
}


inline int VectorInt::ref_count() const
{
    return p->ref_count;
}

inline int* VectorInt::addr() const
{
    return data;
}

inline VectorInt::operator int*() 
{
    return data;
}


inline int& VectorInt::operator()(int i)
{
#ifdef VECTOR_INT_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif 
    return data[i];
}

inline int& VectorInt::operator()(int i) const
{
#ifdef VECTOR_INT_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif
    return data[i];
}

//  *CHANGE*  [] is the same as ()
inline int& VectorInt::operator[](int i)
{
#ifdef VECTOR_INT_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif  
    return data[i];
}

//  *CHANGE*  [] is the same as ()
inline int& VectorInt::operator[](int i) const
{
#ifdef VECTOR_INT_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif  
    return data[i];
}

inline VectorInt& VectorInt::ref(const VectorInt& m)
{
    // always check that the p field has been initialized.
    // Either the lhs or rhs could be a NULL VectorInt...
    
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


inline VectorInt& VectorInt::operator=(const VectorInt& m)
{

    return  ref(m);
}




#endif 
// _VECTOR_INT_H_

