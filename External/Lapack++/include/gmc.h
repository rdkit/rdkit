//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.
//
//      Lapack++ Rectangular Matrix Class
//
//      Dense (nonsingular) matrix, assumes no special structure or properties.
//
//      ) allows 2-d indexing
//      ) non-unit strides
//      ) deep (copy) assignment
//      ) cout << A.info()  prints out internal states of A
//      ) indexing via A(i,j) where i,j are either integers or
//              LaIndex         

#ifndef _LA_GEN_MAT_COMPLEX_H_
#define _LA_GEN_MAT_COMPLEX_H_

#include "lafnames.h"
//#include VECTOR_COMPLEX_H // changed of VC++
#include "vc.h"
//#include LA_INDEX_H // changed of VC++
#include "laindex.h"


class LaGenMatComplex
{
    VectorComplex     v;
    LaIndex           ii[2];
    int             dim[2];  // size of original matrix, not submatrix
    int             sz[2];   // size of this submatrix
    static int  debug_; // trace all entry and exits into methods and 
                        // operators of this class.  This variable is
                        // explicitly initalized in lagenmatCOMPLEX.cc

    static int      *info_;   // print matrix info only, not values
                             //   originally 0, set to 1, and then
                             //   reset to 0 after use.
                // use as in
                //
                //    cout << B.info() << endl;
                //
                // this *info_ member is unique in that it really isn't
                // part of the matrix info, just a flag as to how
                // to print it.   We've included in this beta release
                // as part of our testing, but we do not expect it 
                // to be user accessable.
                // It has to be declared as global static
                // so that we may monitor expresssions like
                // X::(const &X) and still utilize without violating
                // the "const" condition.
                // Because this *info_ is used at most one at a time,
                // there is no harm in keeping only one copy of it,
                // also, we do not need to malloc free space every time
                // we call a matrix constructor.


    int shallow_; // set flag to '0' in order to return matrices
                    // by value from functions without unecessary
                    // copying.


    // users shouldn't be able to modify assignment semantics..
    //
    //LaGenMatComplex& shallow_assign();

public:



        /*::::::::::::::::::::::::::*/

        /* Constructors/Destructors */

        /*::::::::::::::::::::::::::*/


        LaGenMatComplex();
        LaGenMatComplex(int, int);
        LaGenMatComplex(COMPLEX*, int, int);
        LaGenMatComplex(const LaGenMatComplex&);
    virtual ~LaGenMatComplex();


        /*::::::::::::::::::::::::::::::::*/

        /*  Indices and access operations */

        /*::::::::::::::::::::::::::::::::*/

    inline int size(int d) const;   // submatrix size
    inline int inc(int d) const;    // explicit increment
    inline int gdim(int d) const;   // global dimensions
    inline int start(int d) const;  // return ii[d].start()
    inline int end(int d) const;    // return ii[d].end()
    inline LaIndex index(int d) const;// index
    inline int ref_count() const;
    inline LaGenMatComplex& shallow_assign();
    inline COMPLEX* addr() const;       // begining addr of data space
    
    inline COMPLEX& operator()(int i, int j);
    inline COMPLEX& operator()(int i, int j) const;
    LaGenMatComplex operator()(const LaIndex& I, const LaIndex& J) ;
    LaGenMatComplex operator()(const LaIndex& I, const LaIndex& J) const;

            LaGenMatComplex& operator=(COMPLEX s);
    inline  LaGenMatComplex& operator=(const LaGenMatComplex& s); //copy

    LaGenMatComplex& resize(int m, int n);
    LaGenMatComplex& resize(const LaGenMatComplex& s);
    LaGenMatComplex& ref(const LaGenMatComplex& s);
    LaGenMatComplex& inject(const LaGenMatComplex& s);
    LaGenMatComplex& copy(const LaGenMatComplex& s);

    inline int shallow() const      // read global shallow flag
        { return shallow_;}
    inline int debug() const;       // read global debug flag
    inline int debug(int d);        // set global debug flag
    inline const LaGenMatComplex& info() const { 
            int *t = info_; 
            *t = 1; 
            return *this;};

    //* I/O *//
    friend std::ostream& operator<<(std::ostream&, const LaGenMatComplex&);
    std::ostream& Info(std::ostream& s)
    {
        s << "Size: (" << size(0) << "x" << size(1) << ") " ;
        s << "Indeces: " << ii[0] << " " << ii[1];
        s << "#ref: " << ref_count() << "addr: " << (unsigned) addr() << std::endl;
        return s;
    };

};  //* End of LaGenMatComplex Class *//


        

    //* Member Functions *//


 
inline int LaGenMatComplex::size(int d) const
{
    return sz[d];
}

inline int LaGenMatComplex::inc(int d) const
{
    return ii[d].inc();
}

inline int LaGenMatComplex::gdim(int d) const
{
    return dim[d];
}

inline int LaGenMatComplex::start(int d) const
{
    return ii[d].start();
}

inline int LaGenMatComplex::end(int d) const
{
    return ii[d].end();
}

inline int LaGenMatComplex::ref_count() const
{
    return v.ref_count();
}


inline LaIndex LaGenMatComplex::index(int d)  const
{
    return ii[d];
}

inline COMPLEX* LaGenMatComplex::addr() const
{
    return  v.addr();
}

inline int LaGenMatComplex::debug() const
{
    return debug_;
}

inline int LaGenMatComplex::debug(int d)
{
    return debug_ = d;
}

inline COMPLEX& LaGenMatComplex::operator()(int i, int j)
{

#ifdef LA_BOUNDS_CHECK
    assert(i>=0);
    assert(i<size(0));
    assert(j>=0);
    assert(j<size(1));
#endif
    return v( dim[0]*(ii[1].start() + j*ii[1].inc()) + 
                ii[0].start() + i*ii[0].inc());
}

inline COMPLEX& LaGenMatComplex::operator()(int i, int j) const
{

#ifdef LA_BOUNDS_CHECK
    assert(i>=0);
    assert(i<size(0));
    assert(j>=0);
    assert(j<size(1));
#endif

    return v( dim[0]*(ii[1].start() + j*ii[1].inc()) + 
                ii[0].start() + i*ii[0].inc());
}




inline LaGenMatComplex& LaGenMatComplex::operator=(const LaGenMatComplex& s)
{

    return copy(s);
}


inline  LaGenMatComplex&  LaGenMatComplex::shallow_assign()
{
    shallow_ = 1;
    return *this;
}





#endif 
// _LA_GEN_MAT_H_
