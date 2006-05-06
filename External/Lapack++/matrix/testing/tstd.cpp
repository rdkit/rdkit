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


//#define LA_VECTOR_DEBUG

#include <stdlib.h>
#include "lafnames.h"
//#include LA_SYMM_TRIDIAG_MAT_DOUBLE_H // changed of VC++
#include "sytrmd.h"


int main(int argc, char *argv[])
{
    int N;
    int i,j;

    if (argc < 2)
    {
        std::cerr << "Usage:  " << argv[0] << "N " << std::endl;
        exit(1);
    }

    N = atoi(argv[1]);

    // Test constructors
    //
    LaSymmTridiagMatDouble A;
    std::cout << std::endl << "null consturctor " << std::endl;
    std::cout << "A:\n" << A.info() << std::endl;

    std::cout << std::endl;
    LaSymmTridiagMatDouble C(N);
    std::cout << std::endl << "(int, int) constructor " << std::endl;
    std::cout << "C(N):\n" << C.info() << std::endl;

    std::cout << std::endl;
    std::cout << " &C(0,0): " << &C(0,0) << std::endl;

    std::cout << std::endl;
    LaSymmTridiagMatDouble D(C);        // D is also N,N
    std::cout << std::endl << "X(const &X) constructor " << std::endl;
    std::cout << "D(C):\n" << D.info() << std::endl;

    std::cout << std::endl;
    std::cout << "test A.ref(C)\n";
    A.ref(C);
    std::cout << "A:\n" << A.info() << std::endl;

    std::cout << "D(i,i) = 3.3" << std::endl;

    for (j=0;j<N; j++)
        for (i=0; i<N; i++)
        {
            D(i,i) = 3.3;
        }
    std::cout << std::endl;
    std::cout << "D:\n" << D << std::endl;

    std::cout << "D(3,2) = 9.0" << std::endl;
    D(3,2) = 9.0;

    std::cout << std::endl;
    std::cout << "test A.copy(D)\n";
    A.copy(D);
    std::cout << "A:\n" << A.info() << std::endl;
    std::cout << "A:\n" << A << std::endl;
    

    LaVectorDouble E;
    E.ref(D.diag(0));
    std::cout << std::endl;
    std::cout << "test E.ref(D.diag(0))\n";
    std::cout << "E:\n" << E << std::endl;
    std::cout << std::endl;
    std::cout << "D:\n" << D.info() << std::endl;

    std::cout << "test error message: E.ref(D.diag(2))\n";
    std::cout << std::endl;
    E.ref(D.diag(2));
    
}
