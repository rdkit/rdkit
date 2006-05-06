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


//#include "lapack++.h"

//#define VECTOR_DEBUG
#define LA_BOUNDS_CHECK

#include "lafnames.h"
//#include LA_BAND_MAT_DOUBLE_H // changed of VC++
#include "bmd.h"

#ifdef WIN32
#include <string>
#endif

void mult2(double* v, int len)
{
    for (int i=0; i<len; i++)
        v[i] *= 2;
}

LaBandMatDouble madd(const LaBandMatDouble &A, const LaBandMatDouble &B)
{
if (A.debug())
{
    std::cout << ">>> madd(A,B) \n";
    std::cout << "    A: " << A.info() << std::endl;
    std::cout << "    B: " << B.info() << std::endl;
}

    if (A.size(0) != B.size(0) || A.size(1) != B.size(1))
    {
    std::cerr << "LaBandMatDouble madd(LaBandMatDouble &A, LaBandMatDouble &B): \
            non-conformant arrays.\n";

        return LaBandMatDouble(0,0,-1);     // 0x0 matrix
    }

    int i, j;
    LaBandMatDouble result(A);


    for (j=0; j<A.size(1); j++)
        for (i=0; i<A.size(0); i++)
            result(i,j) = A(i,j) + B(i,j);

if (A.debug())
{
    std::cout << "   A+B: " << result.info() << std::endl;
    std::cout << "<<< madd(A,B)\n";
}
    return result;

}
            
int main(int argc, char *argv[])
{
    int N,l,u;
    int i,j;

    if (argc<4)
    {
        std::cerr << "Usage:  " << argv[0] << " N u l\n" << std::endl;
		exit(1);
    }
    N = atoi(argv[1]);
    u = atoi(argv[2]);
    l = atoi(argv[3]);
    
    double v[100]; // = {4.0};
    for (int k=0;k<100;k++) v[k] = 4.0;
    
    // Test constructors
    
    LaBandMatDouble A;
    std::cout << std::endl << "null consturctor " << std::endl;
    std::cout << "A(): " << A.info() << std::endl;

    LaBandMatDouble C(N,l,u);
    std::cout << std::endl << "(int,int,int) constructor " << std::endl;
    std::cout << "C(N,l,u): " << C.info() << std::endl;

    C = 0.0;
    std::cout << std::endl << "C = 0.0 " << std::endl;
    std::cout << "C:\n" << C << std::endl;

    LaBandMatDouble D(C);       // D is also N,N
    std::cout << std::endl << "X(const &X) constructor " << std::endl;
    std::cout << "D(C) :" << D.info() << std::endl;

    //LaBandMatDouble F(v, 10, 10);     //create a 10x10 column-major
    //std::cout << std::endl << "(double *, int , int) constructor" << std::endl;
    //std::cout << "F(v,10,10): " << F.info() << std::endl;

    LaBandMatDouble E(C);
    std::cout << std::endl << "X(const &X) constructor with submatrices " << std::endl;
    std::cout << "E(C): " << E.info() << std::endl;


    for (j=0;j<N; j++)
        for (i=0; i<N; i++)
          if(((i>=j)&&(i-j<=l))||((j>i)&&(j-i<=u)))
            C(i,j) = i + j/100.0;

    std::cout << std::endl;   
    std::cout << "test operator(int, int)"  << std::endl;
    std::cout << "Initalize C(i,j) = i + j/100.0 " << std::endl;
    std::cout << "C: " << C.info() << std::endl;
    std::cout << std::endl;   
    std::cout << "C: " << std::endl;
    std::cout <<  C << std::endl;

    std::cout << std::endl;   
    E = 0.0;
    std::cout << "E:\n" << E << std::endl;
    for (i=0;i<N; i++)
        for (j=0;j<N; j++)
            if ((i>=j)&&(i-j<=E.subdiags()))
                E(i,j) = 9.0;
    std::cout << "test operator(int, int) assignment"  << std::endl;
    std::cout << "E(i,i) = 9.0\n" << E.info() << std::endl;
    std::cout << "E:\n" << E << std::endl;
    std::cout << std::endl;   
    E.print_data();
    E(2,2) = 3;
    E(3,2) = 5;
    std::cout << std::endl;   
    E.print_data();
    std::cout << std::endl;   
    std::cout << "E:\n" << E << std::endl;
    std::cout <<"E(2,2) = 3\n";
    double hold = E(2,4);
    std::cout << "double hold = E(2,4) : "<< hold << std::endl;
    
    std::cout << std::endl;
    std::cout << "test debug and copy(): " << std::endl;
    D.copy(C);
    std::cout << std::endl;
    D.debug(1);
    std::cout << "D.copy(C): " << D.info() << std::endl;
    D.debug(0);
    std::cout << "D:\n" << D << std::endl;
    std::cout << "C:\n" << C << std::endl;

    E = 5.55;
    std::cout << std::endl;
    std::cout << "operator=(double) " << std::endl;
    std::cout << "E = 5.55 : " << E.info() << std::endl;
    std::cout << "E:\n" << E << std::endl;
    std::cout << std::endl;   

    D.ref(C);
    std::cout << std::endl;
    std::cout << "D.ref(C): " << D.info() << std::endl;
    std::cout << std::endl;   
    std::cout << D << std::endl;

    std::cout << std::endl;
    std::cout << "test automatic destructuion of temporaries: " << std::endl;
    LaBandMatDouble B;
    for (int c=0; c<10; c++)
    {
        B.ref(C);
        std::cout << "B.ref(C)): " << B.info() << std::endl;
    }

    C.ref(C);
    std::cout << std::endl;
    std::cout <<"test C.ref(C) case works correctly." << std::endl;
    std::cout << "C.ref(C) " << C.info() << std::endl;
    std::cout << std::endl;   

    return 0;
}
