//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.
//
//      Lapack++ "Shared" Vector Float Class
//
//      A lightweight vector class with minimal overhead.
//
//      shallow assignment
//      unit stride
//      inlined access A(i)
//      optional (compile-time) array bounds checking through 
//              VECTOR_FLOAT_BOUNDS_CHECK
//      A(i) is the same as A[i]
//      auto conversion to float*
//      a null vector has size of 0, but has the ref_count structure
//              has been initalized
//

#ifndef _VECTOR_FLOAT_H_
#define _VECTOR_FLOAT_H_    

#include <iostream>       // for formatted printing of matrices

#ifndef __ASSERT_H
#include <assert.h>     // cheap "error" protection used in checking
#endif                  // checking array bounds.


typedef  struct {
    int        sz;                                        
    float*    data;                                       
    int        ref_count;
} vrefFloat;
                        


class VectorFloat
{                                                                      
    private:                                                           
           vrefFloat *p;
           float *data;            // performance hack, avoid float
                                    // indirection to data.
    public:                                                            
                                                                       
        /*::::::::::::::::::::::::::*/                                 
        /* Constructors/Destructors */                                 
        /*::::::::::::::::::::::::::*/                                 
                                                                       
    //inline VectorFloat();     // this should behave as VectorFloat(0)
    VectorFloat(int);                             
    VectorFloat(int, float);   // can't be inlined because of 'for'
                                       // statement.
    VectorFloat(float*, int);
    VectorFloat(const VectorFloat&); 
    ~VectorFloat() ;                              
                                                                       
        /*::::::::::::::::::::::::::::::::*/                           
        /*  Indices and access operations */                           
        /*::::::::::::::::::::::::::::::::*/                           
                                                                       
    inline float&       operator[](int); 
    inline float&       operator[](int) const;  //read only
    inline float&       operator()(int); 
    inline float&       operator()(int) const; // read only
    inline              operator    float*(); 
    inline int          size() const;
    inline int          null() const;
           int          resize(int d);
    inline int          ref_count() const;  // return the number of ref counts
    inline float*       addr() const;
                                                                       
        /*::::::::::::::*/                                             
        /*  Assignment  */                                             
        /*::::::::::::::*/                                             
                                                                       
    inline  VectorFloat& operator=(const VectorFloat&);
            VectorFloat& operator=(float);
    inline  VectorFloat& ref(const VectorFloat &);
            VectorFloat& inject(VectorFloat&);
            VectorFloat& copy(const VectorFloat&);

    /* I/O */                                                      
    friend std::ostream&   operator<<(std::ostream&, const VectorFloat&);       

};                                                                     


    // operators and member functions

inline int VectorFloat::null()  const
{
    return (size() == 0) ;
}

inline int VectorFloat::size() const
{
    return   p-> sz;
}


inline int VectorFloat::ref_count() const
{
    return p->ref_count;
}

inline float* VectorFloat::addr() const
{
    return data;
}

inline VectorFloat::operator float*() 
{
    return data;
}


inline float& VectorFloat::operator()(int i)
{
#ifdef VECTOR_FLOAT_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif 
    return data[i];
}

inline float& VectorFloat::operator()(int i) const
{
#ifdef VECTOR_FLOAT_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif
    return data[i];
}

//  [] *always* performs bounds-check 
//  *CHANGE*  [] is the same as ()
inline float& VectorFloat::operator[](int i)
{
#ifdef VECTOR_FLOAT_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif  
    return data[i];
}

//  [] *always* performs bounds-check 
//  *CHANGE*  [] is the same as ()
inline float& VectorFloat::operator[](int i) const
{
#ifdef VECTOR_FLOAT_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif  
    return data[i];
}

inline VectorFloat& VectorFloat::ref(const VectorFloat& m)
{
    // always check that the p field has been initialized.
    // Either the lhs or rhs could be a NULL VectorFloat...
    
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


inline VectorFloat& VectorFloat::operator=(const VectorFloat& m)
{

    return  ref(m);
}




#endif 
// _VECTOR_FLOAT_H_

