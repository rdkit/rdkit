/* dopla.f -- translated by f2c (version of 20 August 1993  13:15:44).
   You must link the resulting object file with the libraries:
    -lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__2 = 2;
static integer c__3 = 3;

doublereal dopla_(subnam, m, n, kl, ku, nb, subnam_len)
char *subnam;
integer *m, *n, *kl, *ku, *nb;
ftnlen subnam_len;
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal ret_val;

    /* Builtin functions */
    /* Subroutine */ int s_copy();

    /* Local variables */
    static doublereal adds;
    static logical sord, corz;
    static integer i;
    extern logical lsame_();
    static char c1[1], c2[2], c3[3];
    static doublereal mults, addfac, ek, em, en, wl, mulfac, wu;
    extern logical lsamen_();
    static doublereal emn;


/*  -- LAPACK timing routine (version 1.1b) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     February 29, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DOPLA computes an approximation of the number of floating point */
/*  operations used by the subroutine SUBNAM with the given values */
/*  of the parameters M, N, KL, KU, and NB. */

/*  This version counts operations for the LAPACK subroutines. */

/*  Arguments */
/*  ========= */

/*  SUBNAM  (input) CHARACTER*6 */
/*          The name of the subroutine. */

/*  M       (input) INTEGER */
/*          The number of rows of the coefficient matrix.  M >= 0. */

/*  N       (input) INTEGER */
/*          The number of columns of the coefficient matrix. */
/*          For solve routine when the matrix is square, */
/*          N is the number of right hand sides.  N >= 0. */

/*  KL      (input) INTEGER */
/*          The lower band width of the coefficient matrix. */
/*          If needed, 0 <= KL <= M-1. */
/*          For xGEQRS, KL is the number of right hand sides. */

/*  KU      (input) INTEGER */
/*          The upper band width of the coefficient matrix. */
/*          If needed, 0 <= KU <= N-1. */

/*  NB      (input) INTEGER */
/*          The block size.  If needed, NB >= 1. */

/*  Notes */
/*  ===== */

/*  In the comments below, the association is given between arguments */
/*  in the requested subroutine and local arguments.  For example, */

/*  xGETRS:  N, NRHS  =>  M, N */

/*  means that arguments N and NRHS in DGETRS are passed to arguments */
/*  M and N in this procedure. */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     -------------------------------------------------------- */
/*     Initialize DOPLA to 0 and do a quick return if possible. */
/*     -------------------------------------------------------- */

    ret_val = 0.;
    mults = 0.;
    adds = 0.;
    *c1 = *subnam;
    s_copy(c2, subnam + 1, 2L, 2L);
    s_copy(c3, subnam + 3, 3L, 3L);
    sord = lsame_(c1, "S", 1L, 1L) || lsame_(c1, "D", 1L, 1L);
    corz = lsame_(c1, "C", 1L, 1L) || lsame_(c1, "Z", 1L, 1L);
    if (*m <= 0 || ! (sord || corz)) {
    return ret_val;
    }

/*     --------------------------------------------------------- */
/*     If the coefficient matrix is real, count each add as 1 */
/*     operation and each multiply as 1 operation. */
/*     If the coefficient matrix is complex, count each add as 2 */
/*     operations and each multiply as 6 operations. */
/*     --------------------------------------------------------- */

    if (lsame_(c1, "S", 1L, 1L) || lsame_(c1, "D", 1L, 1L)) {
    addfac = 1.;
    mulfac = 1.;
    } else {
    addfac = 2.;
    mulfac = 6.;
    }
    em = (doublereal) (*m);
    en = (doublereal) (*n);
    ek = (doublereal) (*kl);

/*     --------------------------------- */
/*     GE:  GEneral rectangular matrices */
/*     --------------------------------- */

    if (lsamen_(&c__2, c2, "GE", 2L, 2L)) {

/*        xGETRF:  M, N  =>  M, N */

    if (lsamen_(&c__3, c3, "TRF", 3L, 3L)) {
        emn = (doublereal) min(*m,*n);
        adds = emn * (em * en - (em + en) * (emn + 1.) / 2. + (emn + 1.) *
             (emn * 2. + 1.) / 6.);
        mults = adds + emn * (em - (emn + 1.) / 2.);

/*        xGETRS:  N, NRHS  =>  M, N */

    } else if (lsamen_(&c__3, c3, "TRS", 3L, 3L)) {
        mults = en * em * em;
        adds = en * (em * (em - 1.));

/*        xGETRI:  N  =>  M */

    } else if (lsamen_(&c__3, c3, "TRI", 3L, 3L)) {
        mults = em * (em * (em * .66666666666666663 + .5) + 
            .83333333333333337);
        adds = em * (em * (em * .66666666666666663 - 1.5) + 
            .83333333333333337);

/*        xGEQRF or xGEQLF:  M, N  =>  M, N */

    } else if (lsamen_(&c__3, c3, "QRF", 3L, 3L) || lsamen_(&c__3, c3, 
        "QR2", 3L, 3L) || lsamen_(&c__3, c3, "QLF", 3L, 3L) || 
        lsamen_(&c__3, c3, "QL2", 3L, 3L)) {
        if (*m >= *n) {
        mults = en * (em + 3.8333333333333335 + en / 2. + en * (em - 
            en / 3.));
        adds = en * (en * (em - en / 3. + .5) + .83333333333333337);
        } else {
        mults = em * (en * 2. + 3.8333333333333335 - em / 2. + em * (
            en - em / 3.));
        adds = em * (en + .83333333333333337 - em / 2. + em * (en - 
            em / 3.));
        }

/*        xGERQF or xGELQF:  M, N  =>  M, N */

    } else if (lsamen_(&c__3, c3, "RQF", 3L, 3L) || lsamen_(&c__3, c3, 
        "RQ2", 3L, 3L) || lsamen_(&c__3, c3, "LQF", 3L, 3L) || 
        lsamen_(&c__3, c3, "LQ2", 3L, 3L)) {
        if (*m >= *n) {
        mults = en * (em + 4.833333333333333 + en / 2. + en * (em - 
            en / 3.));
        adds = en * (em + .83333333333333337 + en * (em - en / 3. - 
            .5));
        } else {
        mults = em * (en * 2. + 4.833333333333333 - em / 2. + em * (
            en - em / 3.));
        adds = em * (em / 2. + .83333333333333337 + em * (en - em / 
            3.));
        }

/*        xGEQPF: M, N => M, N */

    } else if (lsamen_(&c__3, c3, "QPF", 3L, 3L)) {
        emn = (doublereal) min(*m,*n);
        mults = en * 2 * en + emn * (em * 3 + en * 5 + em * 2 * en - (emn 
            + 1) * (en + 4 + em - (emn * 2 + 1) / 3));
        adds = en * en + emn * (em * 2 + en + em * 2 * en - (emn + 1) * (
            en + 2 + em - (emn * 2 + 1) / 3));

/*        xGEQRS or xGERQS:  M, N, NRHS  =>  M, N, KL */

    } else if (lsamen_(&c__3, c3, "QRS", 3L, 3L) || lsamen_(&c__3, c3, 
        "RQS", 3L, 3L)) {
        mults = ek * (en * (2. - ek) + em * (en * 2. + (em + 1.) / 2.));
        adds = ek * (en * (1. - ek) + em * (en * 2. + (em - 1.) / 2.));

/*        xGELQS or xGEQLS:  M, N, NRHS  =>  M, N, KL */

    } else if (lsamen_(&c__3, c3, "LQS", 3L, 3L) || lsamen_(&c__3, c3, 
        "QLS", 3L, 3L)) {
        mults = ek * (em * (2. - ek) + en * (em * 2. + (en + 1.) / 2.));
        adds = ek * (em * (1. - ek) + en * (em * 2. + (en - 1.) / 2.));

/*        xGEBRD:  M, N  =>  M, N */

    } else if (lsamen_(&c__3, c3, "BRD", 3L, 3L)) {
        if (*m >= *n) {
        mults = en * (en * (em * 2. - en * .66666666666666663 + 2.) + 
            6.666666666666667);
        adds = en * (en - em + 1.6666666666666667 + en * (em * 2. - 
            en * .66666666666666663));
        } else {
        mults = em * (em * (en * 2. - em * .66666666666666663 + 2.) + 
            6.666666666666667);
        adds = em * (em - en + 1.6666666666666667 + em * (en * 2. - 
            em * .66666666666666663));
        }

/*        xGEHRD:  N  =>  M */

    } else if (lsamen_(&c__3, c3, "HRD", 3L, 3L)) {
        if (*m == 1) {
        mults = 0.;
        adds = 0.;
        } else {
        mults = em * (em * (em * 1.6666666666666667 + .5) - 
            1.1666666666666667) - 13.;
        adds = em * (em * (em * 1.6666666666666667 - 1.) - 
            .66666666666666663) - 8.;
        }

    }

/*     ---------------------------- */
/*     GB:  General Banded matrices */
/*     ---------------------------- */
/*        Note:  The operation count is overestimated because */
/*        it is assumed that the factor U fills in to the maximum */
/*        extent, i.e., that its bandwidth goes from KU to KL + KU. */

    } else if (lsamen_(&c__2, c2, "GB", 2L, 2L)) {

/*        xGBTRF:  M, N, KL, KU  =>  M, N, KL, KU */

    if (lsamen_(&c__3, c3, "TRF", 3L, 3L)) {
        for (i = min(*m,*n); i >= 1; --i) {
/* Computing MAX */
/* Computing MIN */
        i__3 = *kl, i__4 = *m - i;
        i__1 = 0, i__2 = min(i__3,i__4);
        wl = (doublereal) max(i__1,i__2);
/* Computing MAX */
/* Computing MIN */
        i__3 = *kl + *ku, i__4 = *n - i;
        i__1 = 0, i__2 = min(i__3,i__4);
        wu = (doublereal) max(i__1,i__2);
        mults += wl * (wu + 1.);
        adds += wl * wu;
/* L10: */
        }

/*        xGBTRS:  N, NRHS, KL, KU  =>  M, N, KL, KU */

    } else if (lsamen_(&c__3, c3, "TRS", 3L, 3L)) {
/* Computing MAX */
/* Computing MIN */
        i__3 = *kl, i__4 = *m - 1;
        i__1 = 0, i__2 = min(i__3,i__4);
        wl = (doublereal) max(i__1,i__2);
/* Computing MAX */
/* Computing MIN */
        i__3 = *kl + *ku, i__4 = *m - 1;
        i__1 = 0, i__2 = min(i__3,i__4);
        wu = (doublereal) max(i__1,i__2);
        mults = en * (em * (wl + 1. + wu) - (wl * (wl + 1.) + wu * (wu + 
            1.)) * .5);
        adds = en * (em * (wl + wu) - (wl * (wl + 1.) + wu * (wu + 1.)) * 
            .5);

    }

/*     -------------------------------------- */
/*     PO:  POsitive definite matrices */
/*     PP:  Positive definite Packed matrices */
/*     -------------------------------------- */

    } else if (lsamen_(&c__2, c2, "PO", 2L, 2L) || lsamen_(&c__2, c2, "PP", 
        2L, 2L)) {

/*        xPOTRF:  N  =>  M */

    if (lsamen_(&c__3, c3, "TRF", 3L, 3L)) {
        mults = em * (em * (em * .16666666666666666 + .5) + 
            .33333333333333331);
        adds = em * .16666666666666666 * (em * em - 1.);

/*        xPOTRS:  N, NRHS  =>  M, N */

    } else if (lsamen_(&c__3, c3, "TRS", 3L, 3L)) {
        mults = en * (em * (em + 1.));
        adds = en * (em * (em - 1.));

/*        xPOTRI:  N  =>  M */

    } else if (lsamen_(&c__3, c3, "TRI", 3L, 3L)) {
        mults = em * (em * (em * .33333333333333331 + 1.) + 
            .66666666666666663);
        adds = em * (em * (em * .33333333333333331 - .5) + 
            .16666666666666666);

    }

/*     ------------------------------------ */
/*     PB:  Positive definite Band matrices */
/*     ------------------------------------ */

    } else if (lsamen_(&c__2, c2, "PB", 2L, 2L)) {

/*        xPBTRF:  N, K  =>  M, KL */

    if (lsamen_(&c__3, c3, "TRF", 3L, 3L)) {
        mults = ek * (ek * (ek * -.33333333333333331 - 1.) - 
            .66666666666666663) + em * (ek * (ek * .5 + 1.5) + 1.);
        adds = ek * (ek * (ek * -.33333333333333331 - .5) - 
            .16666666666666666) + em * (ek / 2. * (ek + 1.));

/*        xPBTRS:  N, NRHS, K  =>  M, N, KL */

    } else if (lsamen_(&c__3, c3, "TRS", 3L, 3L)) {
        mults = en * ((em * 2 - ek) * (ek + 1.));
        adds = en * (ek * (em * 2 - (ek + 1.)));

    }

/*     -------------------------------------------------------- */
/*     SY:  SYmmetric indefinite matrices */
/*     SP:  Symmetric indefinite Packed matrices */
/*     HE:  HErmitian indefinite matrices (complex only) */
/*     HP:  Hermitian indefinite Packed matrices (complex only) */
/*     -------------------------------------------------------- */

    } else if (lsamen_(&c__2, c2, "SY", 2L, 2L) || lsamen_(&c__2, c2, "SP", 
        2L, 2L) || lsamen_(&c__3, subnam, "ZHE", 6L, 3L) || lsamen_(&c__3,
         subnam, "ZHE", 6L, 3L) || lsamen_(&c__3, subnam, "ZHP", 6L, 3L) 
        || lsamen_(&c__3, subnam, "ZHP", 6L, 3L)) {

/*        xSYTRF:  N  =>  M */

    if (lsamen_(&c__3, c3, "TRF", 3L, 3L)) {
        mults = em * (em * (em * .16666666666666666 + .5) + 
            3.3333333333333335);
        adds = em / 6. * (em * em - 1.);

/*        xSYTRS:  N, NRHS  =>  M, N */

    } else if (lsamen_(&c__3, c3, "TRS", 3L, 3L)) {
        mults = en * em * em;
        adds = en * (em * (em - 1.));

/*        xSYTRI:  N  =>  M */

    } else if (lsamen_(&c__3, c3, "TRI", 3L, 3L)) {
        mults = em * (em * em * .33333333333333331 + .66666666666666663);
        adds = em * (em * em * .33333333333333331 - .33333333333333331);

/*        xSYTRD, xSYTD2:  N  =>  M */

    } else if (lsamen_(&c__3, c3, "TRD", 3L, 3L) || lsamen_(&c__3, c3, 
        "TD2", 3L, 3L)) {
        if (*m == 1) {
        mults = 0.;
        adds = 0.;
        } else {
        mults = em * (em * (em * .66666666666666663 + 2.5) - 
            .16666666666666666) - 15.;
        adds = em * (em * (em * .66666666666666663 + 1.) - 
            2.6666666666666665) - 4.;
        }
    }

/*     ------------------- */
/*     Triangular matrices */
/*     ------------------- */

    } else if (lsamen_(&c__2, c2, "TR", 2L, 2L) || lsamen_(&c__2, c2, "TP", 
        2L, 2L)) {

/*        xTRTRS:  N, NRHS  =>  M, N */

    if (lsamen_(&c__3, c3, "TRS", 3L, 3L)) {
        mults = en * em * (em + 1.) / 2.;
        adds = en * em * (em - 1.) / 2.;

/*        xTRTRI:  N  =>  M */

    } else if (lsamen_(&c__3, c3, "TRI", 3L, 3L)) {
        mults = em * (em * (em * .16666666666666666 + .5) + 
            .33333333333333331);
        adds = em * (em * (em * .16666666666666666 - .5) + 
            .33333333333333331);

    }

    } else if (lsamen_(&c__2, c2, "TB", 2L, 2L)) {

/*        xTBTRS:  N, NRHS, K  =>  M, N, KL */

    if (lsamen_(&c__3, c3, "TRS", 3L, 3L)) {
        mults = en * (em * (em + 1.) / 2. - (em - ek - 1.) * (em - ek) / 
            2.);
        adds = en * (em * (em - 1.) / 2. - (em - ek - 1.) * (em - ek) / 
            2.);
    }

/*     -------------------- */
/*     Trapezoidal matrices */
/*     -------------------- */

    } else if (lsamen_(&c__2, c2, "TZ", 2L, 2L)) {

/*        xTZRQF:  M, N => M, N */

    if (lsamen_(&c__3, c3, "RQF", 3L, 3L)) {
        emn = (doublereal) min(*m,*n);
        mults = em * 3 * (en - em + 1) + (en * 2 - em * 2 + 3) * (em * em 
            - emn * (emn + 1) / 2);
        adds = (en - em + 1) * (em + em * 2 * em - emn * (emn + 1));
    }

/*     ------------------- */
/*     Orthogonal matrices */
/*     ------------------- */

    } else if (sord && lsamen_(&c__2, c2, "OR", 2L, 2L) || corz && lsamen_(&
        c__2, c2, "UN", 2L, 2L)) {

/*        -MQR, -MLQ, -MQL, or -MRQ:  M, N, K, SIDE  =>  M, N, KL, KU 
*/
/*           where KU<= 0 indicates SIDE = 'L' */
/*           and   KU> 0  indicates SIDE = 'R' */

    if (lsamen_(&c__3, c3, "MQR", 3L, 3L) || lsamen_(&c__3, c3, "MLQ", 3L,
         3L) || lsamen_(&c__3, c3, "MQL", 3L, 3L) || lsamen_(&c__3, 
        c3, "MRQ", 3L, 3L)) {
        if (*ku <= 0) {
        mults = ek * en * (em * 2. + 2. - ek);
        adds = ek * en * (em * 2. + 1. - ek);
        } else {
        mults = ek * (em * (en * 2. - ek) + (em + en + (1. - ek) / 2.)
            );
        adds = ek * em * (en * 2. + 1. - ek);
        }

/*        -GQR or -GQL:  M, N, K  =>  M, N, KL */

    } else if (lsamen_(&c__3, c3, "GQR", 3L, 3L) || lsamen_(&c__3, c3, 
        "GQL", 3L, 3L)) {
        mults = ek * (en * 2. - ek - 1.6666666666666667 + (em * 2. * en + 
            ek * (ek * .66666666666666663 - em - en)));
        adds = ek * (en - em + .33333333333333331 + (em * 2. * en + ek * (
            ek * .66666666666666663 - em - en)));

/*        -GLQ or -GRQ:  M, N, K  =>  M, N, KL */

    } else if (lsamen_(&c__3, c3, "GLQ", 3L, 3L) || lsamen_(&c__3, c3, 
        "GRQ", 3L, 3L)) {
        mults = ek * (em + en - ek - .66666666666666663 + (em * 2. * en + 
            ek * (ek * .66666666666666663 - em - en)));
        adds = ek * (em - en + .33333333333333331 + (em * 2. * en + ek * (
            ek * .66666666666666663 - em - en)));

    }

    }

    ret_val = mulfac * mults + addfac * adds;

    return ret_val;

/*     End of DOPLA */

} /* dopla_ */

