//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.

//  Lapack++ templated  vector class

#ifndef _VECTOR_TEMPLATE_H_
#define _VECTOR_TEMPLATE_H_    

#include <iostream>       // for formatted printing of matrices

//#include <stdlib.h>

#ifndef __ASSERT_H
#include <assert.h>     // cheap "error" protection used in checking
#endif                  // checking array bounds.


template<class type>
typedef  struct vrep{
    int        sz;                                        
    type*    data;                                       
    int        ref_count;
} vref<type>;
                        

template<class type>
class Vector
{                                                                      
    private:                                                           
           vref<type> *p;
           type *data;            // performance hack, avoid type
                                    // indirection to data.
    public:                                                            
                                                                       
        /*::::::::::::::::::::::::::*/                                 
        /* Constructors/Destructors */                                 
        /*::::::::::::::::::::::::::*/                                 
                                                                       
    //inline Vector<type>();     // this should behave as Vector<type>(0)
    inline Vector(int);                             
    Vector(int, const type);   // can't be inlined because of 'for'
                                       // statement.
    inline Vector(type*, int);
    Vector(const Vector<type>&); 
    inline ~Vector() ;                              
                                                                       
        /*::::::::::::::::::::::::::::::::*/                           
        /*  Indices and access operations */                           
        /*::::::::::::::::::::::::::::::::*/                           
                                                                       
    inline type&    operator[](int); 
    inline type&    operator()(int); 
    inline type&    operator()(int) const; // read only
    inline operator type*(); 
    inline int     size() const;
    inline int      null();
    int     resize(int d);
    inline int      ref_count() const;  // return the number of ref counts
    inline type*    addr() const;
                                                                       
        /*::::::::::::::*/                                             
        /*  Assignment  */                                             
        /*::::::::::::::*/                                             
                                                                       
    inline Vector<type>& operator=(const Vector<type>&);
    inline Vector<type>& operator=(type);
    inline Vector<type>& ref(const Vector<type> &);
    Vector<type>& inject(Vector<type>&);
    Vector<type>& copy(const Vector<type>&);

    /* I/O */                                                      
    friend std::ostream&   operator<<(std::ostream&, const Vector<type>&);       

};                                                                     

    // constructors

template<class type>
inline Vector<type>::Vector(int n=0)
{                                                                      
    assert(n>=0);
    p = new vref<type>;
    p->sz = n;                                                          
    p->ref_count = 1;                        
    p->data = data = new type[n]; 
}                                                                      


template<class type>
inline Vector<type>::Vector(type *d, int n)
{                                                                      
    p = new vref<type>;
    p->sz = n;                                                          
    p->ref_count = 2; 
    p-> data = data = d;                        
}                                                                      
                                                                       
                                                                       

template<class type>
inline Vector<type>::~Vector()
{

        if (--(p->ref_count) == 0)              // perform garbage col.
        {
           delete [] ( (type*) p->data);
           delete p;
        }
}


    // operators and member functions

template<class type>
inline int Vector<type>::null()
{
    return (size() == 0) ;
}

template<class type>
inline int Vector<type>::size() const
{
    return   p-> sz;
}


template<class type>
inline int Vector<type>::ref_count() const
{
    return p->ref_count;
}

template<class type>
inline type* Vector<type>::addr() const
{
    return data;
}

template<class type>
inline Vector<type>::operator type*() 
{
    return &((*this)(0));
}

template<class type>
inline type& Vector<type>::operator()(int i)
{
#if 0   // already done by LA_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif 
    return data[i];
}

template<class type>
inline type& Vector<type>::operator()(int i) const
{
#if 0       // already done by LA_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif
    return data[i];
}

//  [] *always* performs bounds-check 
//
template<class type>
inline type& Vector<type>::operator[](int i)
{
    assert(0<=i && i<size());
    return data[i];
}

template<class type>
inline Vector<type>& Vector<type>::ref(const Vector<type>& m)
{
    // always check that the p field has been initialized.
    // Either the lhs or rhs could be a NULL Vector<type>...
        m.p->ref_count++;
        if (--(p->ref_count) == 0)              // perform garbage col.
        {
           delete [] ( (type*) p->data);
           delete p;
        }
        p = m.p;
        data = m.data;
        return *this;
}


template<class type>
inline Vector<type>& Vector<type>::operator=(const Vector<type>& m)
{

    return  copy(m);
}


template<class type>
inline Vector<type>& Vector<type>::operator=(type scalar)
{
    for (int i=0; i<size(); i++)
            (*this)(i) = scalar;
    return *this;
}


#endif 
// _VECTOR_TEMPLATE_H_

