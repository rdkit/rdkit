/* hc.f -- translated by f2c (version 20020208).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C */
/*                                                            C */
/*  HIERARCHICAL CLUSTERING using (user-specified) criterion. C */
/*                                                            C */
/*  Parameters:                                               C */
/*                                                            C */
/*  DISS(LEN)         dissimilarities in lower half diagonal  C */
/*                    storage; LEN = N.N-1/2,                 C */
/*  IOPT              clustering criterion to be used,        C */
/*  IA, IB, CRIT      history of agglomerations; dimensions   C */
/*                    N, first N-1 locations only used,       C */
/*  MEMBR, NN, DISNN  vectors of length N, used to store      C */
/*                    cluster cardinalities, current nearest  C */
/*                    neighbour, and the dissimilarity assoc. C */
/*                    with the latter.                        C */
/*  FLAG              boolean indicator of agglomerable obj./ C */
/*                    clusters.                               C */
/*                                                            C */
/*  F. Murtagh, ESA/ESO/STECF, Garching, February 1986.       C */
/*                                                            C */
/* ------------------------------------------------------------C */
/* Subroutine */ int hc_(n, len, iopt, ia, ib, crit, membr, nn, disnn, flag__,
	 diss)
integer *n, *len, *iopt, *ia, *ib;
doublereal *crit, *membr;
integer *nn;
doublereal *disnn;
logical *flag__;
doublereal *diss;
{
    /* Initialized data */

    static doublereal inf = 1e20;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal dmin__;
    static integer i__, j, k;
    static doublereal x;
    static integer i2, j2, jj, im, jm;
    static doublereal xx;
    static integer ind, ncl;
    extern integer ioffset_();
    static integer ind1, ind2, ind3;

    /* Parameter adjustments */
    --flag__;
    --disnn;
    --nn;
    --membr;
    --crit;
    --ib;
    --ia;
    --diss;

    /* Function Body */

/*  Initializations */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	membr[i__] = (float)1.;
	flag__[i__] = TRUE_;
    }
    ncl = *n;
    if (*iopt == 1) {
	i__1 = *n * (*n - 1) / 2;
	for (ind = 1; ind <= i__1; ++ind) {
	    diss[ind] /= (float)2.;
	}
    }
/*     (Above is done for the case of the min. var. method */
/*     where merging criteria are defined in terms of variances */
/*     rather than distances.) */

/*  Carry out an agglomeration - first create list of NNs */

    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dmin__ = inf;
	i__2 = *n;
	for (j = i__ + 1; j <= i__2; ++j) {
	    ind = ioffset_(n, &i__, &j);
	    if (diss[ind] >= dmin__) {
		goto L500;
	    }
	    dmin__ = diss[ind];
	    jm = j;
L500:
	    ;
	}
	nn[i__] = jm;
	disnn[i__] = dmin__;
    }

L400:
/*     Next, determine least diss. using list of NNs */
    dmin__ = inf;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (! flag__[i__]) {
	    goto L600;
	}
	if (disnn[i__] >= dmin__) {
	    goto L600;
	}
	dmin__ = disnn[i__];
	im = i__;
	jm = nn[i__];
L600:
	;
    }
    --ncl;

/*  This allows an agglomeration to be carried out. */

    i2 = min(im,jm);
    j2 = max(im,jm);
    ia[*n - ncl] = i2;
    ib[*n - ncl] = j2;
    crit[*n - ncl] = dmin__;
    ind1 = ioffset_(n, &i2, &j2);
/*      write(6,*) "agglom: ",i2,j2,dmin,diss(ind1) */

/*  Update dissimilarities from new cluster. */

    flag__[j2] = FALSE_;
    dmin__ = inf;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (! flag__[k]) {
	    goto L800;
	}
	if (k == i2) {
	    goto L800;
	}
	x = membr[i2] + membr[j2] + membr[k];
	if (i2 < k) {
	    ind1 = ioffset_(n, &i2, &k);
	} else {
	    ind1 = ioffset_(n, &k, &i2);
	}
	if (j2 < k) {
	    ind2 = ioffset_(n, &j2, &k);
	} else {
	    ind2 = ioffset_(n, &k, &j2);
	}
	ind3 = ioffset_(n, &i2, &j2);
	xx = diss[ind3];

/*  WARD'S MINIMUM VARIANCE METHOD - IOPT=1. */

	if (*iopt == 1) {
	    diss[ind1] = (membr[i2] + membr[k]) * diss[ind1] + (membr[j2] + 
		    membr[k]) * diss[ind2] - membr[k] * xx;
	    diss[ind1] /= x;
	}

/*  SINGLE LINK METHOD - IOPT=2. */

	if (*iopt == 2) {
/* Computing MIN */
	    d__1 = diss[ind1], d__2 = diss[ind2];
	    diss[ind1] = min(d__1,d__2);
	}

/*  COMPLETE LINK METHOD - IOPT=3. */

	if (*iopt == 3) {
/* Computing MAX */
	    d__1 = diss[ind1], d__2 = diss[ind2];
	    diss[ind1] = max(d__1,d__2);
	}

/*  AVERAGE LINK (OR GROUP AVERAGE) METHOD - IOPT=4. */

	if (*iopt == 4) {
	    diss[ind1] = (membr[i2] * diss[ind1] + membr[j2] * diss[ind2]) / (
		    membr[i2] + membr[j2]);
	}

/*  MCQUITTY'S METHOD - IOPT=5. */

	if (*iopt == 5) {
	    diss[ind1] = diss[ind1] * (float).5 + diss[ind2] * (float).5;
	}

/*  MEDIAN (GOWER'S) METHOD - IOPT=6. */

	if (*iopt == 6) {
	    diss[ind1] = diss[ind1] * (float).5 + diss[ind2] * (float).5 - xx 
		    * (float).25;
	}

/*  CENTROID METHOD - IOPT=7. */

	if (*iopt == 7) {
	    diss[ind1] = (membr[i2] * diss[ind1] + membr[j2] * diss[ind2] - 
		    membr[i2] * membr[j2] * xx / (membr[i2] + membr[j2])) / (
		    membr[i2] + membr[j2]);
	}

	if (i2 > k) {
	    goto L800;
	}
	if (diss[ind1] >= dmin__) {
	    goto L800;
	}
	dmin__ = diss[ind1];
	jj = k;
L800:
	;
    }
    membr[i2] += membr[j2];
    disnn[i2] = dmin__;
    nn[i2] = jj;

/*  Update list of NNs insofar as this is required. */

    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (! flag__[i__]) {
	    goto L900;
	}
	if (nn[i__] == i2) {
	    goto L850;
	}
	if (nn[i__] == j2) {
	    goto L850;
	}
	goto L900;
L850:
/*        (Redetermine NN of I:) */
	dmin__ = inf;
	i__2 = *n;
	for (j = i__ + 1; j <= i__2; ++j) {
	    ind = ioffset_(n, &i__, &j);
	    if (! flag__[j]) {
		goto L870;
	    }
	    if (i__ == j) {
		goto L870;
	    }
	    if (diss[ind] >= dmin__) {
		goto L870;
	    }
	    dmin__ = diss[ind];
	    jj = j;
L870:
	    ;
	}
	nn[i__] = jj;
	disnn[i__] = dmin__;
L900:
	;
    }

/*  Repeat previous steps until N-1 agglomerations carried out. */

    if (ncl > 1) {
	goto L400;
    }


    return 0;
} /* hc_ */



integer ioffset_(n, i__, j)
integer *n, *i__, *j;
{
    /* System generated locals */
    integer ret_val;

/*  Map row I and column J of upper half diagonal symmetric matrix */
/*  onto vector. */
    /*ret_val = *j + (*i__ - 1) * *n - *i__ * (*i__ + 1) / 2;*/
    if(*j > *i__) ret_val = (*j-1)*(*j-2)/2+*i__;
    else ret_val = (*i__-1)*(*i__-2)/2+*j;
    return ret_val;
} /* ioffset_ */

