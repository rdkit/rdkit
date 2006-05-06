//
//              LAPACK++ 1.1 Linear Algebra Package 1.1
//               University of Tennessee, Knoxvilee, TN.
//            Oak Ridge National Laboratory, Oak Ridge, TN.
//        Authors: J. J. Dongarra, E. Greaser, R. Pozo, D. Walker
//                 (C) 1992-1996 All Rights Reserved
//
//                             NOTICE
//
// Permission to use, copy, modify, and distribute this software and
// its documentation for any purpose and without fee is hereby granted
// provided that the above copyright notice appear in all copies and
// that both the copyright notice and this permission notice appear in
// supporting documentation.
//
// Neither the Institutions (University of Tennessee, and Oak Ridge National
// Laboratory) nor the Authors make any representations about the suitability 
// of this software for any purpose.  This software is provided ``as is'' 
// without express or implied warranty.
//
// LAPACK++ was funded in part by the U.S. Department of Energy, the
// National Science Foundation and the State of Tennessee.


#include "lafnames.h"
//#include VECTOR_COMPLEX_H // changed of VC++
#include "vc.h"

 VectorComplex::VectorComplex(int n=0)
{                                                                      
    assert(n>=0);
    p = new vrefComplex;
    p->sz = n;                                                          
    p->data = data = new COMPLEX[n]; 
    p->ref_count = 1;                        
}                                                                      


VectorComplex::VectorComplex( const VectorComplex& m)
{
// shallow assignment semantics

    p = m.p;
    data = m.data;
    p->ref_count++;
}

 VectorComplex::VectorComplex(COMPLEX *d, int n)
{                                                                      
    p = new vrefComplex;
    p->sz = n;                                                          
    p->ref_count = 2; 
    p-> data = data = d;                        
}                                                                      
                                                                       
                                                                       

 VectorComplex::~VectorComplex()
{

        if (--(p->ref_count) == 0)              // perform garbage col.
        {
           delete [] p->data;
           delete p;
        }
}

VectorComplex::VectorComplex(int n, COMPLEX scalar)
{                                                                      
    p = new vrefComplex;
    p->sz = n;                                                          
    p->ref_count = 1;                        
    p->data = data = new COMPLEX[n]; 
    for (int i=0; i<n; i++)                                            
        data[i] = scalar;                                            
}                                                                      
                                                                       


// this actually frees memory first, then resizes it.  it reduces
// internal fragmentation of memory pool, and the resizing of
// matrices > 1/2 available memory.

int VectorComplex::resize(int d)
{

    if (d<0)                // do nothing if invalid size
    {
        return size();
    }
    else
    {
        ref(VectorComplex(0));  // possibilby free up destination
        if (d>0)
            ref(VectorComplex(d));
    }
    return size();
}


VectorComplex& VectorComplex::inject( VectorComplex& m)
{
    if (m.size() != size())
    {
		std::cerr << "VectorComplex::inject(): vector sizes do not match.\n";
      return *this;
    }
    int N = size();
    for (int i=0; i<N; i++)
        (*this)(i) = m(i);

    return *this;
}


VectorComplex& VectorComplex::copy(const VectorComplex &m)
{
#if 0
        if (null()) resize(m.size());

        if (size() != m.size())
           cerr << "VectorComplex::copy(VectorComplex &): incompatible vector \
                sizes : "<< size() << " vs. " << m.size() << ".\n";
        else
#endif // 0

    resize(0);                  // free up destination

    int N = m.size();
    VectorComplex tmp(N);

    for (int i=0; i<N; i++)     // should use memcpy() here...
       (*this)(i) = m(i);

    ref(tmp);
    return *this;
}


std::ostream& operator<<(std::ostream& s, const VectorComplex& m)
{
        if (m.p)
        {
                int n = m.size();
                for (int i=0; i<n; i++)
                  s << m(i).real() << " i" << m(i).imag() << "  ";
                s << "\n";
        }
        else s << "NULL VectorComplex.\n";
    return s;
}



 VectorComplex& VectorComplex::operator=(COMPLEX scalar)
{
    for (int i=0; i<size(); i++)
            data[i] = scalar;
    return *this;
}
