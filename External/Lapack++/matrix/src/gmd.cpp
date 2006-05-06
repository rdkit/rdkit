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
//#include LA_GEN_MAT_DOUBLE_H // changed to statement below for VC++
#include "gmd.h"

int LaGenMatDouble::debug_ = 0;     // turn off global deubg flag initially.
                                // use A.debug(1) to turn on/off,
                                // and A.debug() to check current status.

int* LaGenMatDouble::info_= new int;        // turn off info print flag.

LaGenMatDouble::~LaGenMatDouble()
{}

 LaGenMatDouble::LaGenMatDouble() : v(0)
{

    dim[0] = dim[1] = 0;
    sz[0] = sz[1] = 0;
    *info_ = 0;
    shallow_= 0;
}



  LaGenMatDouble::LaGenMatDouble(int m, int n) :v(m*n)
{
    ii[0](0,m-1);
    ii[1](0,n-1);
    dim[0] = m;
    dim[1] = n;
    sz[0] = m;
    sz[1] = n;
    *info_ = 0;
    shallow_= 0;
}

  LaGenMatDouble::LaGenMatDouble(double *d, int m, int n) :v(d, m*n)
{
    if (debug())
    {
		std::cout << ">>> LaGenMatDouble::LaGenMatDouble(double *d, int m, int n) :v(d, m*n)\n";
    }

    ii[0](0,m-1);
    ii[1](0,n-1);
    dim[0] = m;
    dim[1] = n;
    sz[0] = m;
    sz[1] = n;
    *info_ = 0;
    shallow_= 0;  

    if (debug())
    {
		std::cout << "<<< LaGenMatDouble::LaGenMatDouble(double *d, int m, int n) :v(d, m*n)\n";
    }
}

 LaGenMatDouble& LaGenMatDouble::operator=(double s)
{
if (debug())
{
	std::cout << ">>> LaGenMatDouble& LaGenMatDouble::operator=(double s)\n";
	std::cout << "       s = " << s << std::endl;
}

    int i,j;

    for (j=0; j<size(1); j++)
    {
        for (i=0; i<size(0); i++)
        {
            (*this)(i,j) = s;
        }
    }

if (debug())
{
	std::cout << "  *this: " << this->info() << std::endl;
	std::cout << " >>> LaGenMatDouble& LaGenMatDouble::operator=(double s)\n";
}
    
    return *this;
}

LaGenMatDouble& LaGenMatDouble::ref(const LaGenMatDouble& s)
{

        // handle trivial M.ref(M) case
        if (this == &s) return *this;
        else
        {
        ii[0] = s.ii[0];
        ii[1] = s.ii[1];
        dim[0] = s.dim[0];
        dim[1] = s.dim[1];
        sz[0] = s.sz[0];
        sz[1] = s.sz[1];
        shallow_ = 0;

        v.ref(s.v);

        return *this;
        }
}

LaGenMatDouble&  LaGenMatDouble::resize(int m, int n)
{
    if (debug())
    {
		std::cout << ">>> resize("<< m << "," << n << ")\n";
    }

    // first, reference 0x0 matrix, potentially freeing memory
    // this allows one to resize a matrix > 1/2 of the available
    // memory

    LaGenMatDouble tmp1(0,0);
    ref(tmp1);

    // now, reference an MxN matrix
    LaGenMatDouble tmp(m,n);
    ref(tmp);

if (debug())
{
	std::cout << "<<< resize: *this: " << this->info() << std::endl;
}
    return *this;

}

 LaGenMatDouble&  LaGenMatDouble::resize(const LaGenMatDouble &A)
{
    if (debug())
    {
		std::cout << ">>> resize(const LaGenMatDouble &A)" << std::endl;
		std::cout << "      A: " << A.info() << std::endl;
    }

    resize(A.size(0), A.size(1));
    if (debug())
    {
		std::cout << "<<< resize(const LaGenMatDouble &A) " << std::endl;
    }
    return *this;
}


LaGenMatDouble::LaGenMatDouble(const LaGenMatDouble& X) : v(0)
{
    if (X.debug())
    {
		std::cout << ">>> LaGenMatDouble::LaGenMatDouble(const LaGenMatDouble& X)\n";
		std::cout << "X: " << X.info() << std::endl;

    }
        debug_ = X.debug_;

        shallow_ = 0;  // do not perpeturate shallow copies, otherwise
                       //  B = A(I,J) does not work properly...

        if (X.shallow_)
        {
            v.ref(X.v);
            dim[0] = X.dim[0];
            dim[1] = X.dim[1];
            sz[0] = X.sz[0];
            sz[1] = X.sz[1];
            ii[0] = X.ii[0];
            ii[1] = X.ii[1];
        }
        else
        {
            if (X.debug())
				std::cout << std::endl << "Data is being copied!\n" << std::endl;
            v.resize(X.size(0)*X.size(1)); 
            ii[0](0,X.size(0)-1);
            ii[1](0,X.size(1)-1);
            dim[0] = sz[0] = X.size(0);
            dim[1] = sz[1] = X.size(1);
            int i, j, M = X.size(0), N = X.size(1);
                for (j=0; j<N; j++)
                    for (i=0; i<M; i++)
                        (*this)(i,j) = X(i,j);
        }  

    if (debug())
    {
		std::cout << "*this: " << info() << std::endl;
		std::cout << "<<< LaGenMatDouble::LaGenMatDouble(const LaGenMatDouble& X)\n";
    }
}


LaGenMatDouble& LaGenMatDouble::copy(const LaGenMatDouble& X) 
{
    if (debug())
    {
	std::cout << ">>> LaGenMatDouble& LaGenMatDouble::copy\
(const LaGenMatDouble& X)\n";
	std::cout << " X: " << X.info() << std::endl;
    }

        // current scheme in copy() is to detach the left-hand-side
        // from whatever it was pointing to.
        //
        resize(X);

        int i,j, M = X.size(0),  N = X.size(1);
        for (i=0; i<M; i++)
                for (j=0; j<N; j++)
                    (*this)(i,j) = X(i,j);

    if ((debug()))
    {
		std::cout << " *this: " << this->info() << std::endl;
        std::cout << " <<< LaGenMatDouble& LaGenMatDouble::copy \
(const LaGenMatDouble& s)\n";
    }

        return *this;
}


LaGenMatDouble& LaGenMatDouble::inject(const LaGenMatDouble& s)
{
    assert(s.size(0) == size(0));
    assert(s.size(1) == size(1));

    int i, j,  M=size(0), N=size(1);

    for (j=0;j<N;j++)
        for (i=0;i<M; i++)
            operator()(i,j) = s(i,j);

    return *this;
}


LaGenMatDouble LaGenMatDouble::operator()(const LaIndex& II, const LaIndex& JJ) const
{
if (debug())
{
	std::cout << ">>> LaGenMatDouble::operator(const LaIndex& const LaIndex&)\n";
}
    LaIndex I, J;

    if (II.null())
        I(0,size(0)-1);
    else
        I = II;
    if (JJ.null())
        J(0,size(1)-1);
    else
        J = JJ;

    //(JJ.null()) ?  J(0,size(1)) : J = JJ;

    assert(I.inc() != 0);
    assert(J.inc() != 0);

    if (I.inc() > 0)
    {
        assert(I.start() >= 0);
        assert(I.start() <= I.end());
        assert(I.end() < size(0));
    }
    else // I.inc() < 0
    {
        assert(I.start() < size(0));
        assert(I.start() >= I.end());
        assert(I.end() >= 0);
    }

    if (J.inc() > 0)
    {
        assert(J.start() >= 0);
        assert(J.start() <= J.end());
        assert(J.end() < size(1));
    }
    else  // J.inc() < 0
    {
        assert(J.start() < size(1));
        assert(J.start() >= J.end());
        assert(J.end() >= 0);
    }

    LaGenMatDouble tmp;

    tmp.dim[0] = dim[0];
    tmp.dim[1] = dim[1];
    tmp.sz[0] = (I.end() - I.start())/I.inc() + 1;
    tmp.sz[1] = (J.end() - J.start())/J.inc() + 1;

    tmp.ii[0].start() =  ii[0].start() + I.start()*ii[0].inc();
    tmp.ii[0].inc() = ii[0].inc() * I.inc();
    tmp.ii[0].end() = (I.end() - I.start())/ I.inc() * tmp.ii[0].inc() 
            + tmp.ii[0].start();

    tmp.ii[1].start() =  ii[1].start() + J.start()*ii[1].inc();
    tmp.ii[1].inc() = ii[1].inc() * J.inc();
    tmp.ii[1].end() = (J.end() - J.start())/ J.inc() * tmp.ii[1].inc() 
            + tmp.ii[1].start();

    tmp.v.ref(v);
    tmp.shallow_assign();

if (debug())
{
	std::cout << "    return value: " << tmp.info() << std::endl;
	std::cout << "<<< LaGenMatDouble::operator(const LaIndex& const LaIndex&)\n";
}

    return tmp;

}

LaGenMatDouble LaGenMatDouble::operator()(const LaIndex& II, const LaIndex& JJ) 
{
if (debug())
{
	std::cout << ">>> LaGenMatDouble::operator(const LaIndex& const LaIndex&)\n";
}
    LaIndex I, J;

    if (II.null())
        I(0,size(0)-1);
    else
        I = II;
    if (JJ.null())
        J(0,size(1)-1);
    else
        J = JJ;

    //(JJ.null()) ?  J(0,size(1)) : J = JJ;

    assert(I.inc() != 0);
    assert(J.inc() != 0);

    if (I.inc() > 0)
    {
        assert(I.start() >= 0);
        assert(I.start() <= I.end());
        assert(I.end() < size(0));
    }
    else // I.inc() < 0
    {
        assert(I.start() < size(0));
        assert(I.start() >= I.end());
        assert(I.end() >= 0);
    }

    if (J.inc() > 0)
    {
        assert(J.start() >= 0);
        assert(J.start() <= J.end());
        assert(J.end() < size(1));
    }
    else  // J.inc() < 0
    {
        assert(J.start() < size(1));
        assert(J.start() >= J.end());
        assert(J.end() >= 0);
    }

    LaGenMatDouble tmp;

    tmp.dim[0] = dim[0];
    tmp.dim[1] = dim[1];
    tmp.sz[0] = (I.end() - I.start())/I.inc() + 1;
    tmp.sz[1] = (J.end() - J.start())/J.inc() + 1;

    tmp.ii[0].start() =  ii[0].start() + I.start()*ii[0].inc();
    tmp.ii[0].inc() = ii[0].inc() * I.inc();
    tmp.ii[0].end() = (I.end() - I.start())/ I.inc() * tmp.ii[0].inc() 
            + tmp.ii[0].start();

    tmp.ii[1].start() =  ii[1].start() + J.start()*ii[1].inc();
    tmp.ii[1].inc() = ii[1].inc() * J.inc();
    tmp.ii[1].end() = (J.end() - J.start())/ J.inc() * tmp.ii[1].inc() 
            + tmp.ii[1].start();

    tmp.v.ref(v);
    tmp.shallow_assign();

if (debug())
{
	std::cout << "    return value: " << tmp.info() << std::endl;
	std::cout << "<<< LaGenMatDouble::operator(const LaIndex& const LaIndex&)\n";
}

    return tmp;

}


std::ostream& operator<<(std::ostream& s, const LaGenMatDouble& G)
{
    if (*(G.info_))     // print out only matrix info, not actual values
    {
        *(G.info_) = 0; // reset the flag
        s << "(" << G.size(0) << "x" << G.size(1) << ") " ;
        s << "Indices: " << G.ii[0] << " " << G.ii[1];
        s << " #ref: " << G.ref_count();
        s << " shallow:" << G.shallow_  ;
    }
    else 
    {

        int i,j;

        for (i=0; i<G.size(0); i++)
        {
            for (j=0; j<G.size(1); j++)
                s << G(i,j) << "  ";
            s << "\n";
        }
    }
    return s;
}

