//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.
//
//      Lapack++ "Shared" Vector Double Class
//
//      A lightweight vector class with minimal overhead.
//
//      shallow assignment
//      unit stride
//      inlined access A(i)
//      optional (compile-time) array bounds checking through 
//              VECTOR_DOUBLE_BOUNDS_CHECK
//      A(i) is the same as A[i]
//      auto conversion to double*
//      a null vector has size of 0, but has the ref_count structure
//              has been initalized
//

#ifndef _VECTOR_DOUBLE_H_
#define _VECTOR_DOUBLE_H_    

#include <iostream>       // for formatted printing of matrices

#ifndef __ASSERT_H
#include <assert.h>     // cheap "error" protection used in checking
#endif                  // checking array bounds.




typedef  struct {
    int        sz;                                        
    double*    data;                                       
    int        ref_count;
    int        d_ref_count;
} vrefDouble;
                        


class VectorDouble
{                                                                      
    private:                                                           
           vrefDouble *p;
           double *data;            // performance hack, avoid double
                                    // indirection to data.
    public:                                                            
                                                                       
        /*::::::::::::::::::::::::::*/                                 
        /* Constructors/Destructors */                                 
        /*::::::::::::::::::::::::::*/                                 
                                                                       
    //inline VectorDouble();     // this should behave as VectorDouble(0)
    VectorDouble(int);                             
    VectorDouble(int, double);   // can't be inlined because of 'for'
                                       // statement.
    VectorDouble(double*, int);
    VectorDouble(const VectorDouble&); 
    ~VectorDouble() ;                              
                                                                       
        /*::::::::::::::::::::::::::::::::*/                           
        /*  Indices and access operations */                           
        /*::::::::::::::::::::::::::::::::*/                           
                                                                       
    inline double&      operator[](int); 
    inline double&      operator[](int) const;  // read only
    inline double&      operator()(int); 
    inline double&      operator()(int) const; // read only
    inline              operator    double*(); 
    inline int          size() const;
    inline int          null() const;
           int          resize(int d);
    inline int          ref_count() const;  // return the number of ref counts
    inline double*      addr() const;
                                                                       
        /*::::::::::::::*/                                             
        /*  Assignment  */                                             
        /*::::::::::::::*/                                             
                                                                       
    inline  VectorDouble& operator=(const VectorDouble&);
            VectorDouble& operator=(double);
    inline  VectorDouble& ref(const VectorDouble &);
            VectorDouble& inject(VectorDouble&);
            VectorDouble& copy(const VectorDouble&);

    /* I/O */                                                      
    friend std::ostream&   operator<<(std::ostream&, const VectorDouble&);       

};                                                                     


    // operators and member functions

inline int VectorDouble::null() const
{
    return (size() == 0) ;
}

inline int VectorDouble::size() const
{
    return   p-> sz;
}


inline int VectorDouble::ref_count() const
{
    return p->ref_count;
}

inline double* VectorDouble::addr() const
{
    return data;
}

inline VectorDouble::operator double*() 
{
    return data;
}


inline double& VectorDouble::operator()(int i)
{
#ifdef VECTOR_DOUBLE_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif 
    return data[i];
}

inline double& VectorDouble::operator()(int i) const
{
#ifdef VECTOR_DOUBLE_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif
    return data[i];
}

//  [] *always* performs bounds-check 
//  *CHANGE*  [] is the same as ()
inline double& VectorDouble::operator[](int i)
{
#ifdef VECTOR_DOUBLE_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif  
    return data[i];
}

inline double& VectorDouble::operator[](int i) const
{
#ifdef VECTOR_DOUBLE_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif  
    return data[i];
}

inline VectorDouble& VectorDouble::ref(const VectorDouble& m)
{
    // always check that the p field has been initialized.
    // Either the lhs or rhs could be a NULL VectorDouble...
    
        if (&m != this)         // not really necessary...
        {
            m.p->ref_count++;
            m.p->d_ref_count++;

	    p->ref_count--;
	    p->d_ref_count--;
	    if(!p->d_ref_count) delete [] p->data;
	    if(!p->ref_count) delete p;
            p = m.p;
            data = m.data;
        }
        return *this;
}


inline VectorDouble& VectorDouble::operator=(const VectorDouble& m)
{

    return  ref(m);
}




#endif 
// _VECTOR_DOUBLE_H_

