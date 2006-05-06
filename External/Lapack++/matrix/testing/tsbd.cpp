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
#define LA_BOUNDS_CHECK

#include "lafnames.h"
//#include LA_SYMM_BAND_MAT_DOUBLE_H // changed of VC++
#include "sybmd.h"

void mult2(double* v, int len)
{
    for (int i=0; i<len; i++)
        v[i] *= 2;
}

LaSymmBandMatDouble madd(const LaSymmBandMatDouble &A, const LaSymmBandMatDouble &B)
{
if (A.debug())
{
    cout << ">>> madd(A,B) \n";
    cout << "    A: " << A.info() << endl;
    cout << "    B: " << B.info() << endl;
}

    if (A.size(0) != B.size(0) || A.size(1) != B.size(1))
    {
    cerr << "LaSymmBandMatDouble madd(LaSymmBandMatDouble &A, LaSymmBandMatDouble &B): \
            non-conformant arrays.\n";

        return LaSymmBandMatDouble(0,0);        // 0x0 matrix
    }

    int i, j;
    LaSymmBandMatDouble result(A);


    for (j=0; j<A.size(1); j++)
        for (i=0; i<A.size(0); i++)
            result(i,j) = A(i,j) + B(i,j);

if (A.debug())
{
    cout << "   A+B: " << result.info() << endl;
    cout << "<<< madd(A,B)\n";
}
    return result;

}
            
int main(int argc, char *argv[])
{
    int N,l;
    int i,j;

    if (argc < 2)
    {
        cerr << "Usage:  " << argv[0] << " N l " << endl;
        exit(1);
    }
    N = atoi(argv[1]);
    l = atoi(argv[2]);
    
    double v[100]; // = {4.0};
    for (int k=0;k<100;k++) v[k] = 4.0;
    
    // Test constructors
    
    LaSymmBandMatDouble A;
    cout << endl << "null consturctor " << endl;
    cout << "A(): " << A.info() << endl;

    LaSymmBandMatDouble C(N,l);
    cout << endl << "(int,int,int) constructor " << endl;
    cout << "C(N,l,u): " << C.info() << endl;

    C = 0.0;
    cout << endl << "C = 0.0 " << endl;
    cout << "C:\n" << C << endl;

    LaSymmBandMatDouble D(C);       // D is also N,N
    cout << endl << "X(const &X) constructor " << endl;
    cout << "D(C) :" << D.info() << endl;

    LaSymmBandMatDouble E(C);
    cout << endl << "X(const &X) constructor with submatrices " << endl;
    cout << "E(C): " << E.info() << endl;


    for (j=0;j<N; j++)
        for (i=0; i<N; i++)
          if((i>=j)&&(i-j<=l))
            C(i,j) = i + j/100.0;

    cout << endl;   
    cout << "test operator(int, int)"  << endl;
    cout << "Initalize C(i,j) = i + j/100.0 " << endl;
    cout << "C: " << C.info() << endl;
    cout << endl;   
    cout << "C: " << endl;
    cout <<  C << endl;

    cout << endl;   
    E = 0.0;
    cout << "E:\n" << E << endl;
    for (i=0;i<N; i++)
        for (j=0;j<N; j++)
            //if ((i>=j)&&(i-j<=E.subdiags()))
                E(i,j) = 9.0;
    cout << "test operator(int, int) assignment"  << endl;
    cout << "E(i,i) = 9.0\n" << E.info() << endl;
    cout << "E:\n" << E << endl;
    cout << endl;   
    E.print_data();
    E(2,2) = 3;
    E(3,2) = 5;
    cout << endl;   
    E.print_data();
    cout << endl;   
    cout << "E:\n" << E << endl;
    cout <<"E(2,2) = 3\n";
    double hold = E(1,4);
    cout << "hold = E(1,4) : "<< hold << endl;
    hold = E(2,4);
    cout << "hold = E(2,4) : "<< hold << endl;
    hold = E(3,4);
    cout << "hold = E(3,4) : "<< hold << endl;
    hold = E(4,4);
    cout << "hold = E(4,4) : "<< hold << endl;
    
    cout << endl;
    cout << "test debug and copy(): " << endl;
    D.copy(C);
    cout << endl;
    D.debug(1);
    cout << "D.copy(C): " << D.info() << endl;
    D.debug(0);
    cout << "D:\n" << D << endl;
    cout << "C:\n" << C << endl;

    E = 5.55;
    cout << endl;
    cout << "operator=(double) " << endl;
    cout << "E = 5.55 : " << E.info() << endl;
    cout << "E:\n" << E << endl;
    cout << endl;   

    D.ref(C);
    cout << endl;
    cout << "D.ref(C): " << D.info() << endl;
    cout << endl;   
    cout << D << endl;

    cout << endl;
    cout << "test automatic destructuion of temporaries: " << endl;
    LaSymmBandMatDouble B;
    for (int c=0; c<10; c++)
    {
        B.ref(C);
        cout << "B.ref(C)): " << B.info() << endl;
    }

    C.ref(C);
    cout << endl;
    cout <<"test C.ref(C) case works correctly." << endl;
    cout << "C.ref(C) " << C.info() << endl;
    cout << endl;   
}
