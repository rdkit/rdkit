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

#include <iostream>
#include <stdlib.h>

#include "lafnames.h"       /* macros for LAPACK++ filenames */
//#include LA_GEN_MAT_DOUBLE_H // changed of VC++
#include "gmd.h"
//#include LA_UTIL_H // changed of VC++
#include "lautil.h"

main(int argc, char **argv)
{

    std::cout.precision(4);
    std::cout.setf(std::ios::scientific, std::ios::floatfield);

    if (argc < 2)
    {
        std::cerr << "Usage " << argv[0] << " N " << std::endl;
        exit(1);
    }
    int N = atoi(argv[1]);

    char fname[100];

    LaGenMatDouble A(N,N);

    for(;;)
    {
        std::cout << "Name of routine : ";
        std::cin >> fname;
        std::cout << "entered fname: " << fname << std::endl;

    std::cout << "Testing block size:  A(" << N << "x" << N << ")  , "
        << LaEnvBlockSize(fname, A) << std::endl;

    }
}

