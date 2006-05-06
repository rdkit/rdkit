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


#include <stdlib.h>
#include "lafnames.h"
//#include LA_TRIDIAG_MAT_DOUBLE_H // changed of VC++
#include "trmd.h"

int main(int argc, char *argv[])
{
    int N;

    if (argc < 2)
    {
         std::cerr << "Usage:  " << argv[0] << " N " << std::endl;
         exit(1);
    }

    N = atoi(argv[1]);

    // Test constructors
    //
    LaTridiagMatDouble A;
    std::cout << std::endl << "null consturctor " << std::endl;
    std::cout << "A:\n" << A.info() << std::endl;

    std::cout << std::endl;
    LaTridiagMatDouble C(N);
    std::cout << std::endl << "(int, int) constructor " << std::endl;
    std::cout << "C(N):\n" << C.info() << std::endl;

    std::cout << std::endl;
    std::cout << " &C(0,0): " << (long) &C(0,0) << std::endl;

    std::cout << std::endl;
    LaTridiagMatDouble D(C);        // D is also N,N
    std::cout << std::endl << "X(const &X) constructor " << std::endl;
    std::cout << "D(C):\n" << D.info() << std::endl;

    std::cout << std::endl;
    std::cout << "test A.ref(C)\n";
    A.ref(C);
    std::cout << "A:\n" << A.info() << std::endl;

    std::cout << "D.diag(0) = 3.3" << std::endl;

    D.diag(0) = 3.3;
   
    std::cout << std::endl;
    std::cout << "D:\n" << D << std::endl;

    std::cout << std::endl;
    std::cout << "test A.copy(D)\n";
    A.copy(D);
    std::cout << "A:\n" << A.info() << std::endl;
    std::cout << "A:\n" << A << std::endl;
    

    LaVectorDouble tmp(3*N-2);
    tmp(LaIndex(0,N-2)) = 9.9;

    C.diag(-1)(LaIndex(0,N-2)) = 1.1;
    std::cout << "\nC:\n" << C << std::endl;
    
    C.diag(-1)(LaIndex(0,N-2)) = tmp(LaIndex(0,N-2));
    std::cout << std::endl;
    std::cout << "test C.diag(-1)(LaIndex(0,N-2)) = tmp(LaIndex(0,N-2))\n";
    std::cout << "\nC:\n" << C << std::endl;
    std::cout << std::endl;

    std::cout << "\ntest error message: C.diag(3))\n";
    C.diag(3) = 5.0;
    std::cout << std::endl;
    
}
