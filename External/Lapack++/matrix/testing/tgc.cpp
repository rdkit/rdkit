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

#define LA_BOUNDS_CHECK

#include <stdlib.h>
#include "lafnames.h"
//#include LA_GEN_MAT_COMPLEX_H // changed of VC++
#include "gmc.h"

            
int main(int argc, char *argv[])
{
    int M, N;
    int i,j;

    if (argc <3)
    {
        cout << "Usage " << argv[0] << " M N " << endl;
        exit(1);
    }

    M = atoi(argv[1]);
    N = atoi(argv[2]);
    
    COMPLEX v[100]; // = {4.0};
    for (int k=0;k<100;k++) v[k] = COMPLEX(4.0, 1.1);
    
    // Test constructors
    //

    LaGenMatComplex A;
    cout << endl << "null consturctor " << endl;
    cout << "A(): " << A.info() << endl;

    cout << endl << "(int, int) constructor " << endl;
    LaGenMatComplex C(M,N);
    cout << "C(M,N): " << C.info() << endl;

    // C.debug(1);
    cout << endl << "X(const &X) constructor " << endl;
    LaGenMatComplex D(C);        // D is also N,N
    cout << "D(C) :" << D.info() << endl;


    cout << endl; 
    LaGenMatComplex K, O(100,100);
    LaGenMatComplex L(O(LaIndex(2,4),LaIndex(2,8)));
    cout << "L(O(LaIndex(2,4),LaIndex(2,8)))\n";
    L = 0.0;
    cout << "L: " << L.info() << endl;

    cout << endl <<"K.copy(L) " << endl;
    K.copy(L);
    cout << "K: " << K.info() << endl;
    cout << "K:\n" << K << endl;

    
    LaIndex I(2,M-1), J(1,N-1);       // evens, odd
    cout << endl << "create indices  I=" << I << ", J=" << J << endl;

    LaGenMatComplex E(C(I, J));
    cout << endl << "X(const &X) constructor with submatrices " << endl;
    cout << "E(C(I,J)): " << E.info() << endl;


    for (j=0;j<N; j++)
        for (i=0; i<M; i++)
            C(i,j) = i + j/100.0;

    cout << endl;   
    cout << "test operator(int, int)"  << endl;
    cout << "Initalize C(i,j) = i + j/100.0 " << endl;
    cout << "C: " << C.info() << endl;
    cout <<  C << endl;

    cout << endl;
    cout <<  "operator(LaIndex, LaIndex)" << endl;  
    cout << "C(I,J) " << C(I,J).info() << endl;
    cout <<  C(I,J) << endl;

    cout << endl;
    cout << "test missing indices (default to whole row or column" << endl;
    cout << "C(LaIndex(),J) " << C(LaIndex(),J).info() << endl;
    cout << C(LaIndex(),J) << endl;
    cout << endl;
    cout << "C(I,LaIndex()) " << C(I,LaIndex()).info() << endl;
    cout << C(I,LaIndex()) << endl;

    cout << endl;
    LaGenMatComplex F;
    cout << "F.ref(C(I,J))\n";
    F.ref(C(I,J));
    cout << F.info() << endl;
    F = 4.44;
    cout <<"F:\n" << endl;
    cout << F << endl;

    E = F;
    cout << endl;
    cout << "operator=() " << endl;
    cout << "E = F : " << E.info() << endl;

    D = C;
    cout << endl;
    cout << "operator=(const Matrix&) "<< endl;
    cout << "D = C : " << D.info() << endl;
    cout << D << endl;


    cout << endl;
    cout << "test automatic destructuion of temporaries: " << endl;
    LaGenMatComplex B;
    for (int c=0; c<10; c++)
    {
        B.ref(C(I,J));
        cout << "B.ref(C(I,J)): " << B.info() << endl;
    }

    C.ref(C);
    cout << endl;
    cout <<"test C.ref(C) case works correctly." << endl;
    cout << "C.ref(C) " << C.info() << endl;

    return 0;
}
