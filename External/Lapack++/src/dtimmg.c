/* dtimmg.f -- translated by f2c (version of 20 August 1993  13:15:44).
   You must link the resulting object file with the libraries:
    -lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;

/* Subroutine */ int dtimmg_(iflag, m, n, a, lda, kl, ku)
integer *iflag, *m, *n;
doublereal *a;
integer *lda, *kl, *ku;
{
    /* Initialized data */

    static integer iseed[4] = { 0,0,0,1 };

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Builtin functions */
    double d_sign();

    /* Local variables */
    static integer i, j, k;
    extern /* Subroutine */ int dcopy_();
    static integer jj, jn, mj, mu;
    extern /* Subroutine */ int dlarnv_();


/*  -- LAPACK timing routine (version 1.1b) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     February 29, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DTIMMG generates a real test matrix whose type is given by IFLAG. */
/*  All the matrices are Toeplitz (constant along a diagonal), with */
/*  random elements on each diagonal. */

/*  Arguments */
/*  ========= */

/*  IFLAG   (input) INTEGER */
/*          The type of matrix to be generated. */
/*          = 0 or 1:   General matrix */
/*          = 2 or -2:  General banded matrix */
/*          = 3 or -3:  Symmetric positive definite matrix */
/*          = 4 or -4:  Symmetric positive definite packed */
/*          = 5 or -5:  Symmetric positive definite banded */
/*          = 6 or -6:  Symmetric indefinite matrix */
/*          = 7 or -7:  Symmetric indefinite packed */
/*          = 8 or -8:  Symmetric indefinite banded */
/*          = 9 or -9:  Triangular */
/*          = 10 or -10:  Triangular packed */
/*          = 11 or -11:  Triangular banded */
/*          For symmetric or triangular matrices, IFLAG > 0 indicates */
/*          upper triangular storage and IFLAG < 0 indicates lower */
/*          triangular storage. */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix to be generated. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix to be generated. */

/*  A       (output) DOUBLE PRECISION array, dimension (LDA,N) */
/*          The generated matrix. */

/*          If the absolute value of IFLAG is 1, 3, or 6, the leading */
/*          M x N (or N x N) subblock is used to store the matrix. */
/*          If the matrix is symmetric, only the upper or lower triangle 
*/
/*          of this block is referenced. */

/*          If the absolute value of IFLAG is 4 or 7, the matrix is */
/*          symmetric and packed storage is used for the upper or lower */
/*          triangle.  The triangular matrix is stored columnwise as a */
/*          inear array, and the array A is treated as a vector of */
/*          length LDA.  LDA must be set to at least N*(N+1)/2. */

/*          If the absolute value of IFLAG is 2 or 5, the matrix is */
/*          returned in band format.  The columns of the matrix are */
/*          specified in the columns of A and the diagonals of the */
/*          matrix are specified in the rows of A, with the leading */
/*          diagonal in row */
/*              KL + KU + 1,  if IFLAG = 2 */
/*              KU + 1,       if IFLAG = 5 or -2 */
/*              1,            if IFLAG = -5 */
/*          If IFLAG = 2, the first KL rows are not used to leave room */
/*          for pivoting in DGBTRF. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of A.  If the generated matrix is */
/*          packed, LDA >= N*(N+1)/2, otherwise LDA >= max(1,M). */

/*  KL      (input) INTEGER */
/*          The number of subdiagonals if IFLAG = 2, 5, or -5. */

/*  KU      (input) INTEGER */
/*          The number of superdiagonals if IFLAG = 2, 5, or -5. */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Data statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
/*     .. */
/*     .. Executable Statements .. */

    if (*m <= 0 || *n <= 0) {
    return 0;

    } else if (*iflag == 0 || *iflag == 1) {

/*        General matrix */

/*        Set first column and row to random values. */

    dlarnv_(&c__2, iseed, m, &a[a_dim1 + 1]);
    i__1 = *n;
    i__2 = *m;
    for (j = 2; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
        i__3 = *m, i__4 = *n - j + 1;
        mj = min(i__3,i__4);
        dlarnv_(&c__2, iseed, &mj, &a[j * a_dim1 + 1]);
        if (mj > 1) {
        i__3 = mj - 1;
        dcopy_(&i__3, &a[j * a_dim1 + 2], &c__1, &a[(j + 1) * a_dim1 
            + 1], lda);
        }
/* L10: */
    }

/*        Fill in the rest of the matrix. */

    i__2 = *n;
    for (j = 2; j <= i__2; ++j) {
        i__1 = *m;
        for (i = 2; i <= i__1; ++i) {
        a[i + j * a_dim1] = a[i - 1 + (j - 1) * a_dim1];
/* L20: */
        }
/* L30: */
    }

    } else if (*iflag == 2 || *iflag == -2) {

/*        General band matrix */

    if (*iflag == 2) {
        k = *kl + *ku + 1;
    } else {
        k = *ku + 1;
    }
/* Computing MIN */
    i__1 = *m, i__3 = *kl + 1;
    i__2 = min(i__1,i__3);
    dlarnv_(&c__2, iseed, &i__2, &a[k + a_dim1]);
/* Computing MIN */
    i__2 = *n - 1;
    mu = min(i__2,*ku);
    dlarnv_(&c__2, iseed, &mu, &a[k - mu + *n * a_dim1]);
    a[k + *n * a_dim1] = a[k + a_dim1];
    i__2 = *n - 1;
    for (j = 2; j <= i__2; ++j) {
/* Computing MIN */
        i__1 = j - 1;
        mu = min(i__1,*ku);
        dcopy_(&mu, &a[k - mu + *n * a_dim1], &c__1, &a[k - mu + j * 
            a_dim1], &c__1);
/* Computing MIN */
        i__3 = *m - j + 1, i__4 = *kl + 1;
        i__1 = min(i__3,i__4);
        dcopy_(&i__1, &a[k + a_dim1], &c__1, &a[k + j * a_dim1], &c__1);
/* L40: */
    }

    } else if (*iflag == 3) {

/*        Symmetric positive definite, upper triangle */

    i__2 = *n - 1;
    dlarnv_(&c__2, iseed, &i__2, &a[*n * a_dim1 + 1]);
    a[*n + *n * a_dim1] = (doublereal) (*n);
    for (j = *n - 1; j >= 1; --j) {
        dcopy_(&j, &a[*n - j + 1 + *n * a_dim1], &c__1, &a[j * a_dim1 + 1]
            , &c__1);
/* L50: */
    }

    } else if (*iflag == -3) {

/*        Symmetric positive definite, lower triangle */

    a[a_dim1 + 1] = (doublereal) (*n);
    if (*n > 1) {
        i__2 = *n - 1;
        dlarnv_(&c__2, iseed, &i__2, &a[a_dim1 + 2]);
    }
    i__2 = *n;
    for (j = 2; j <= i__2; ++j) {
        i__1 = *n - j + 1;
        dcopy_(&i__1, &a[a_dim1 + 1], &c__1, &a[j + j * a_dim1], &c__1);
/* L60: */
    }

    } else if (*iflag == 4) {

/*        Symmetric positive definite packed, upper triangle */

    jn = (*n - 1) * *n / 2 + 1;
    i__2 = *n - 1;
    dlarnv_(&c__2, iseed, &i__2, &a[jn + a_dim1]);
    a[jn + *n - 1 + a_dim1] = (doublereal) (*n);
    jj = jn;
    for (j = *n - 1; j >= 1; --j) {
        jj -= j;
        ++jn;
        dcopy_(&j, &a[jn + a_dim1], &c__1, &a[jj + a_dim1], &c__1);
/* L70: */
    }

    } else if (*iflag == -4) {

/*        Symmetric positive definite packed, lower triangle */

    a[a_dim1 + 1] = (doublereal) (*n);
    if (*n > 1) {
        i__2 = *n - 1;
        dlarnv_(&c__2, iseed, &i__2, &a[a_dim1 + 2]);
    }
    jj = *n + 1;
    i__2 = *n;
    for (j = 2; j <= i__2; ++j) {
        i__1 = *n - j + 1;
        dcopy_(&i__1, &a[a_dim1 + 1], &c__1, &a[jj + a_dim1], &c__1);
        jj = jj + *n - j + 1;
/* L80: */
    }

    } else if (*iflag == 5) {

/*        Symmetric positive definite banded, upper triangle */

    k = *kl;
/* Computing MIN */
    i__2 = *n - 1;
    mu = min(i__2,k);
    dlarnv_(&c__2, iseed, &mu, &a[k + 1 - mu + *n * a_dim1]);
    a[k + 1 + *n * a_dim1] = (doublereal) (*n);
    for (j = *n - 1; j >= 1; --j) {
/* Computing MIN */
        i__2 = j, i__1 = k + 1;
        mu = min(i__2,i__1);
        dcopy_(&mu, &a[k + 2 - mu + *n * a_dim1], &c__1, &a[k + 2 - mu + 
            j * a_dim1], &c__1);
/* L90: */
    }

    } else if (*iflag == -5) {

/*        Symmetric positive definite banded, lower triangle */

    k = *kl;
    a[a_dim1 + 1] = (doublereal) (*n);
/* Computing MIN */
    i__1 = *n - 1;
    i__2 = min(i__1,k);
    dlarnv_(&c__2, iseed, &i__2, &a[a_dim1 + 2]);
    i__2 = *n;
    for (j = 2; j <= i__2; ++j) {
/* Computing MIN */
        i__3 = *n - j + 1, i__4 = k + 1;
        i__1 = min(i__3,i__4);
        dcopy_(&i__1, &a[a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &c__1);
/* L1.1: */
    }

    } else if (*iflag == 6) {

/*        Symmetric indefinite, upper triangle */

    dlarnv_(&c__2, iseed, n, &a[*n * a_dim1 + 1]);
    for (j = *n - 1; j >= 1; --j) {
        dcopy_(&j, &a[*n - j + 1 + *n * a_dim1], &c__1, &a[j * a_dim1 + 1]
            , &c__1);
/* L1.1: */
    }

    } else if (*iflag == -6) {

/*        Symmetric indefinite, lower triangle */

    dlarnv_(&c__2, iseed, n, &a[a_dim1 + 1]);
    i__2 = *n;
    for (j = 2; j <= i__2; ++j) {
        i__1 = *n - j + 1;
        dcopy_(&i__1, &a[a_dim1 + 1], &c__1, &a[j + j * a_dim1], &c__1);
/* L1.1: */
    }

    } else if (*iflag == 7) {

/*        Symmetric indefinite packed, upper triangle */

    jn = (*n - 1) * *n / 2 + 1;
    dlarnv_(&c__2, iseed, n, &a[jn + a_dim1]);
    jj = jn;
    for (j = *n - 1; j >= 1; --j) {
        jj -= j;
        ++jn;
        dcopy_(&j, &a[jn + a_dim1], &c__1, &a[jj + a_dim1], &c__1);
/* L1.1: */
    }

    } else if (*iflag == -7) {

/*        Symmetric indefinite packed, lower triangle */

    dlarnv_(&c__2, iseed, n, &a[a_dim1 + 1]);
    jj = *n + 1;
    i__2 = *n;
    for (j = 2; j <= i__2; ++j) {
        i__1 = *n - j + 1;
        dcopy_(&i__1, &a[a_dim1 + 1], &c__1, &a[jj + a_dim1], &c__1);
        jj = jj + *n - j + 1;
/* L1.1: */
    }

    } else if (*iflag == 8) {

/*        Symmetric indefinite banded, upper triangle */

    k = *kl;
/* Computing MIN */
    i__2 = *n, i__1 = k + 1;
    mu = min(i__2,i__1);
    dlarnv_(&c__2, iseed, &mu, &a[k + 2 - mu + *n * a_dim1]);
    for (j = *n - 1; j >= 1; --j) {
/* Computing MIN */
        i__2 = j, i__1 = k + 1;
        mu = min(i__2,i__1);
        dcopy_(&mu, &a[k + 2 - mu + *n * a_dim1], &c__1, &a[k + 2 - mu + 
            j * a_dim1], &c__1);
/* L1.1: */
    }

    } else if (*iflag == -8) {

/*        Symmetric indefinite banded, lower triangle */

    k = *kl;
/* Computing MIN */
    i__1 = *n, i__3 = k + 1;
    i__2 = min(i__1,i__3);
    dlarnv_(&c__2, iseed, &i__2, &a[a_dim1 + 1]);
    i__2 = *n;
    for (j = 2; j <= i__2; ++j) {
/* Computing MIN */
        i__3 = *n - j + 1, i__4 = k + 1;
        i__1 = min(i__3,i__4);
        dcopy_(&i__1, &a[a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &c__1);
/* L1.1: */
    }

    } else if (*iflag == 9) {

/*        Upper triangular */

    dlarnv_(&c__2, iseed, n, &a[*n * a_dim1 + 1]);
    d__1 = (doublereal) (*n);
    a[*n + *n * a_dim1] = d_sign(&d__1, &a[*n + *n * a_dim1]);
    for (j = *n - 1; j >= 1; --j) {
        dcopy_(&j, &a[*n - j + 1 + *n * a_dim1], &c__1, &a[j * a_dim1 + 1]
            , &c__1);
/* L1.1: */
    }

    } else if (*iflag == -9) {

/*        Lower triangular */

    dlarnv_(&c__2, iseed, n, &a[a_dim1 + 1]);
    d__1 = (doublereal) (*n);
    a[a_dim1 + 1] = d_sign(&d__1, &a[a_dim1 + 1]);
    i__2 = *n;
    for (j = 2; j <= i__2; ++j) {
        i__1 = *n - j + 1;
        dcopy_(&i__1, &a[a_dim1 + 1], &c__1, &a[j + j * a_dim1], &c__1);
/* L1.1: */
    }

    } else if (*iflag == 10) {

/*        Upper triangular packed */

    jn = (*n - 1) * *n / 2 + 1;
    dlarnv_(&c__2, iseed, n, &a[jn + a_dim1]);
    d__1 = (doublereal) (*n);
    a[jn + *n - 1 + a_dim1] = d_sign(&d__1, &a[jn + *n - 1 + a_dim1]);
    jj = jn;
    for (j = *n - 1; j >= 1; --j) {
        jj -= j;
        ++jn;
        dcopy_(&j, &a[jn + a_dim1], &c__1, &a[jj + a_dim1], &c__1);
/* L1.1: */
    }

    } else if (*iflag == -10) {

/*        Lower triangular packed */

    dlarnv_(&c__2, iseed, n, &a[a_dim1 + 1]);
    d__1 = (doublereal) (*n);
    a[a_dim1 + 1] = d_sign(&d__1, &a[a_dim1 + 1]);
    jj = *n + 1;
    i__2 = *n;
    for (j = 2; j <= i__2; ++j) {
        i__1 = *n - j + 1;
        dcopy_(&i__1, &a[a_dim1 + 1], &c__1, &a[jj + a_dim1], &c__1);
        jj = jj + *n - j + 1;
/* L200: */
    }

    } else if (*iflag == 11) {

/*        Upper triangular banded */

    k = *kl;
/* Computing MIN */
    i__2 = *n, i__1 = k + 1;
    mu = min(i__2,i__1);
    dlarnv_(&c__2, iseed, &mu, &a[k + 2 - mu + *n * a_dim1]);
    d__1 = (doublereal) (k + 1);
    a[k + 1 + *n * a_dim1] = d_sign(&d__1, &a[k + 1 + *n * a_dim1]);
    for (j = *n - 1; j >= 1; --j) {
/* Computing MIN */
        i__2 = j, i__1 = k + 1;
        mu = min(i__2,i__1);
        dcopy_(&mu, &a[k + 2 - mu + *n * a_dim1], &c__1, &a[k + 2 - mu + 
            j * a_dim1], &c__1);
/* L210: */
    }

    } else if (*iflag == -11) {

/*        Lower triangular banded */

    k = *kl;
/* Computing MIN */
    i__1 = *n, i__3 = k + 1;
    i__2 = min(i__1,i__3);
    dlarnv_(&c__2, iseed, &i__2, &a[a_dim1 + 1]);
    d__1 = (doublereal) (k + 1);
    a[a_dim1 + 1] = d_sign(&d__1, &a[a_dim1 + 1]);
    i__2 = *n;
    for (j = 2; j <= i__2; ++j) {
/* Computing MIN */
        i__3 = *n - j + 1, i__4 = k + 1;
        i__1 = min(i__3,i__4);
        dcopy_(&i__1, &a[a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &c__1);
/* L220: */
    }
    }

    return 0;

/*     End of DTIMMG */

} /* dtimmg_ */

