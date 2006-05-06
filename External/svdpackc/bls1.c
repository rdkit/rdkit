/*************************************************************************
                           (c) Copyright 1993
                        University of Tennessee
                          All Rights Reserved                          
 *************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <fcntl.h>
#include "bls1.h"

#define  ZERO   0.0

/* #define  UNIX_CREAT */

#ifdef UNIX_CREAT
#define PERMS 0664
#endif


double ddot(long, double *,long, double *, long);
void   opa(long, long, double *, double *);
void   opat(long, double *, double *);
void   daxpy(long, double, double *,long, double *, long);
long   validate(FILE *, long, long, long, long, long, long, double);
float  timer(void);
long   blklan1(FILE *, long, long, long, long, double *, double *, 
	     double *, long, long, double, double *, long , long *,
	     long *, long *, long *);

/***********************************************************************
 *                                                                     *
 *                        main()                                       *
 *     Sparse SVD Via Hybrid Block Lanczos Procedure for Equivalent    *
 *               2-Cyclic Eigensystems  (double precision)             *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   This is a sample driver that invokes BLKLAN1 to compute singular 
   triplets of a large sparse matrix A.  In this test program, bls1,
   the Harwell-Boeing sparse matrix format is employed for accessing
   elements of the m by n sparse matrix A and its transpose.  Other 
   sparse matrix formats can be used, of course.  Approximate singular
   values of A and their residuals can be found in formatted output
   file blo1 and corresponding right and left singular vectors can be
   found in unformatted output file blv1.

   User supplied routines:  opa, opat, timer
   
   opa(m, n, x, y) takes an n-vector x and returns the m-vector y = A * x.  

   opat(n, m, x, y) takes an m-vector x and returns the n-vector y = A' * x,
   where A' denotes the transpose of the matrix A.

   User should edit timer() with an appropriate call to an intrinsic
   timing routine that returns elapsed user cpu time.


   External parameters
   -------------------

   Defined and documented in bls1.h


   Local parameters
   ----------------

  (input)
   maxit     maximum number of iterations
   nc        maximum subspace size
   nb        initial block size
   nums      number of singular values desired 
   tol       residual tolerance
   vtf       TRUE or FALSE to indicate whether both singular values and
		vectors are wanted or only singular values are wanted.
		If vtf is set to "TRUE", the unformatted output file blv1
		will contain the approximate singular vectors written in 
		the order u[0], v[0], u[1], v[1], ..., u[iko-1], v[iko-1].
		u[i], v[i] denote the left and right singular vectors,
		respectively, corresponding to the i-th approximate singular
		value sing[i], and iko is the number of singular values 
		approximated.

   (output)
   sing      array of singular values
   res       array of residual norms of the singular values calculated
   memory    total memory allocated in bytes


   Functions called
   --------------

   BLAS         ddot, daxpy
   USER         opa, opat, timer
   MISC         validate
   BLS1         blklan1


   BLS1 development
   ----------------

   BLS1 is a C translation of the Fortran-77 BLS1 from the SVDPACK
   library written by Michael W. Berry, University of Tennessee,
   Dept. of Computer Science, 107 Ayres Hall, Knoxville, TN, 37996-1301.
   It is a modified version of the block Lanczos algorithm first published
   by Golub, Luk, and Overton (ACM TOMS 7(2):149-169, 1981).
   This particular implementation is discussed in "Multiprocessor Sparse
   SVD Algorithms and Applications", Ph.D. thesis by M. Berry, University
   of Illinois at Urbana-Champaign, October 1990.

   Date written:  05 May 1992

   Theresa H. Do
   University of Tennessee
   Dept. of Computer Science
   107 Ayres Hall
   Knoxville, TN, 37996-1301
   internet: tdo@cs.utk.edu

 ***********************************************************************/

  main() 
{
   float exetime;
   long i, nrow, ncol, nnzero, nc, nb, nums, maxit, nn;
   long ncold, nbold, nk, size1, size2, indexu, indexv;
   long memory, vectors;
   double tol, *sing, *v, *u, *res, *tmpu, *tmpv, *tptr1;
   double tmp1, tmp2, tmp3, xnorm;
   char title[73], name[41], vtf[6];
   FILE *fp_in1, *fp_in2;
   FILE *fp_out1 = NULL;
   long  fp_out2;

   char *in1, *in2, *out1, *out2;
 
   in1 = "blp1";
   in2 = "matrix";
   out1 = "blo1";
   out2 = "blv1";

   /* open files for input/output */
    if (!(fp_in1 = fopen(in1, "r"))) { 
       printf("cannot open file %s for reading\n", in1);
       exit(-1);
    }
    if (!(fp_in2 = fopen(in2, "r"))) {
       printf("cannot open file %s for reading\n", in2);
       exit(-1);
    }
    if (!(fp_out1 = fopen(out1, "w"))) { 
       printf("cannot open output file %s \n", out1);
       exit(-1);
    }

   /* read data */
    fscanf (fp_in2,"%72c%*s%*s%*s%ld%ld%ld%*d",
	    title, &nrow, &ncol, &nnzero);
    nn=nrow+ncol;
    title[73] = '\0';
    fscanf (fp_in1,"%s%ld%ld%ld%ld%lf%s", name, &maxit, &nc, &nb, 
	    &nums, &tol, vtf);

    if (!(strcmp(vtf, "TRUE"))) {
       vectors = 1;

#if  !defined UNIX_CREAT
       if ((fp_out2 = open(out2, O_CREAT | O_RDWR)) == -1) {
          printf("cannot open output file %s \n", out2);
          exit(-1);
       }
#else
       if ((fp_out2 = creat(out2, PERMS )) == -1) {
          printf("cannot open output file %s \n", out2);
          exit(-1);
       }
#endif
    }
    else vectors = 0;

    /* even though the validity of the parameters will be checked in the
     * SVD code, some parameter checking should also be done before
     * allocating memory to ensure that they are nonnegative */
     if (validate(fp_out1, nnzero, nrow, ncol, nc, nb, nums, tol)) {
	fclose(fp_in1);
	fclose(fp_in2);
	fclose(fp_out1);
	if (vectors) close(fp_out2);
	exit(-1);
     }
    /*******************************************************************
     * allocate memory                                                 *
     * pointr - column start array of harwell-boeing sparse matrix     *
     *          format                                    (ncol + 1)   *
     * rowind - row indices array of harwell-boeing format  (nnzero)   *
     * value  - nonzero values array of harwell-boeing sparse matrix   *
     *          format                                      (nnzero)   *
     * sing   -                                               (nums)   *
     * res    -                                               (nums)   *
     * u      -                                        (nrow * nums)   *
     * v      -                                        (ncol * nums)   *
     *******************************************************************/

     size1 = sizeof(double) * (nnzero + nums * (nrow + ncol + 2));
     size2 = sizeof(long) * (ncol + nnzero + 1);
     if (!(rowind = (long *)   malloc(size2))   ||
         !(value  = (double *) malloc(size1))){
	 perror("MALLOC FAILED in MAIN()");
	 exit(errno);
     }
     tptr1 = value;

     /* calculate memory allocated for problem */
     memory = size1 + size2;

     pointr = rowind + nnzero;
     tptr1 += nnzero;
     u      = tptr1;
     tptr1 += nrow * nums;
     v      = tptr1;
     tptr1 += ncol * nums;
     sing   = tptr1;
     tptr1 += nums;
     res    = tptr1;

    /* skip data format line */
    fscanf(fp_in2,"%*s %*s %*s %*s");

    /* read data */
    for (i = 0; i <= ncol; i++) fscanf(fp_in2, "%ld", &pointr[i]);
    for (i = 0; i < ncol; i++) pointr[i] -= 1;

    /* define last element of pointr in case it is not */ 
    pointr[i] = nnzero;

    for (i = 0; i < nnzero; i++) fscanf(fp_in2, "%ld", &rowind[i]);
    for (i = 0; i < nnzero; i++) rowind[i] -= 1;
    for (i = 0; i < nnzero; i++) fscanf(fp_in2, "%lf", &value[i]);

    ncold = nc;
    nbold = nb;
    nk    = nums;
    mxvcount = 0;
    mtxvcount = 0;
    exetime = timer();

    /* make a lanczos run; exit upon error */
    if (blklan1(fp_out1, nnzero, nrow, ncol, nums, v, u, sing, ncold, 
               nbold, tol, res, maxit, &nk, &nc, &nb, &memory)) {

       free(value);
       free(rowind);
       fclose(fp_in1);
       fclose(fp_in2);
       fclose(fp_out1);
       if (vectors) close(fp_out2);
       exit(-1);
    }

    exetime = timer() - exetime;

    if (vectors) {
       size1 = sizeof(double) * nrow;
       size2 = sizeof(double) * ncol;
       memory += size1 + size2;

       /* extra memory required if singular vectors are desired */
       if (!(tmpu = (double *) malloc(size1 + size2))) {
	   perror("SECOND MALLOC FAILED in MAIN()");
	   exit(errno);
       }
       tmpv = tmpu + nrow;
       indexu = 0;
       indexv = 0;
       for (i = 0; i < nk; i++) {
	  tmp1 = ddot(ncol, &v[indexv], 1, &v[indexv], 1);
	  tmp2 = ddot(nrow, &u[indexu], 1, &u[indexu], 1);
	  tmp3 = sqrt(tmp1 + tmp2);
	  opa(nrow, ncol, &v[indexv], tmpu);
	  opat(ncol, &u[indexu], tmpv);
	  daxpy(nrow, -sing[i], &u[indexu], 1, tmpu, 1);
	  daxpy(ncol, -sing[i], &v[indexv], 1, tmpv, 1);
	  xnorm = ddot(ncol, tmpv, 1, tmpv, 1) +
		  ddot(nrow, tmpu, 1, tmpu, 1);
          xnorm = sqrt(xnorm);
	  res[i] = xnorm / tmp3;

	  /* write vectors to binary output file */
	  write (fp_out2, (char *)&u[indexu], size1);
	  write (fp_out2, (char *)&v[indexv], size2);
	  indexu += nrow;
	  indexv += ncol;
       }
    }

    /* write results to output file */
    fprintf(fp_out1,"\n");
    fprintf(fp_out1, " ... \n");
    fprintf(fp_out1, " ... HYBRID BLOCK LANCZOS (CYCLIC)\n");
    fprintf(fp_out1, " ... NO. OF EQUATIONS          =%10ld\n", nn);
    fprintf(fp_out1, " ... MAX. NO. OF ITERATIONS    =%10ld\n", maxit);
    fprintf(fp_out1, " ... NO. OF ITERATIONS TAKEN   =%10ld\n", iter);
    fprintf(fp_out1, " ... NO. OF TRIPLETS SOUGHT    =%10ld\n", nums);
    fprintf(fp_out1, " ... NO. OF TRIPLETS FOUND     =%10ld\n", nk);
    fprintf(fp_out1, " ... INITIAL BLOCKSIZE         =%10ld\n", nbold);
    fprintf(fp_out1, " ... FINAL   BLOCKSIZE         =%10ld\n", nb);
    fprintf(fp_out1, " ... MAXIMUM SUBSPACE BOUND    =%10ld\n", ncold);
    fprintf(fp_out1, " ... FINAL   SUBSPACE BOUND    =%10ld\n", nc);
    fprintf(fp_out1, " ... NO. MULTIPLICATIONS BY A  =%10ld\n", mxvcount);
    fprintf(fp_out1, " ... NO. MULT. BY TRANSPOSE(A) =%10ld\n", mtxvcount);
    fprintf(fp_out1, " ... TOTAL SPARSE MAT-VEC MULT.=%10ld\n", 
            mxvcount + mtxvcount);
    fprintf(fp_out1, " ... MEMORY NEEDED IN BYTES    =%10ld\n", memory);
    if (vectors) 
       fprintf(fp_out1, " ... WANT S-VECTORS ON OUTPUT?           T\n");
    else
       fprintf(fp_out1, " ... WANT S-VECTORS ON OUTPUT?           F\n");
    fprintf(fp_out1, " ... TOLERANCE                 =%10.2E\n\n", tol);
    fprintf(fp_out1, " %s\n", title);
    fprintf(fp_out1, "           %s\n", name);
    fprintf(fp_out1, " ... NO. OF TERMS (ROWS OF A)   = %8ld\n", nrow);
    fprintf(fp_out1, " ... NO. OF DOCS  (COLS OF A)   = %8ld\n", ncol);
    fprintf(fp_out1, " ... \n\n");
    fprintf(fp_out1, " ...... BLSVD EXECUTION TIME=%10.2E\n",exetime);
    fprintf(fp_out1, " ...... \n");
    fprintf(fp_out1, " ......         COMPUTED S-VALUES     (RES. NORMS)\n");
    fprintf(fp_out1, " ...... \n");
    for (i = 0; i < nk; i++)
	   fprintf(fp_out1," ...... %3ld   %22.14E  (%11.2E)\n",i+1,sing[i],res[i]);

    free(value);
    free(rowind);
    fclose(fp_in1);
    fclose(fp_in2);
    fclose(fp_out1);
    if (vectors) {
       free(tmpu);
       close(fp_out2);
    }
    exit(0);
}
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <fcntl.h>

extern double *tres, *y, *temp, *uu, *vv, *u0, *v0, *uvtmp, *pp, *qq;
extern double *alpha, *beta, *p, *q, *t, *z;
extern double **uup, **yp, **vvp, **uvtmpp, **ppp, **qqp;
extern long iconv, nn, iter;

#define   SEED      91827211
#define   MINBLKS   2
#define   TRANSP    1
#define   NTRANSP   0
#define   TRUE      1
#define   FALSE     0
#define   ZERO      0.0
#define   ONE       1.0

void   qriter2(long, double *, double *, double **, double **);
void   dgemv(long, long, long, double, double **, double *, double, double *);
void   dgemm(long, long, long, long, long,
	     double, double **, double *, double, double *);
long   imin(long, long);
double random(long *);
long   validate(FILE *, long, long, long, long, long, long, double);
long   point1(long *, long, long, double **, double *, double *, double);
void   block1(double *, double **, double **, double **, double **,
	    long, long, long, long, long *);
void   orthg(long, long, long, double **, double **, double *);

/***********************************************************************
 *                                                                     *
 *                        blklan1()                                    *
 *                  Block Lanczos SVD Alogrithm                        *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   This hybrid block Lanczos procedure consists of five phases:

   Phase 1:  Block Lanczos outer iteration to yield a block upper
             bidiagonal matrix S whose singular values approximate
	     those of the original sparse matrix A.  Total (or complete)
	     re- orthogonalization is used here.

   Phase 2:  Lanczos method for bi-diagonalizing the S matrix from
	     Phase 1 to yield the bi-diagonal matrix B whose singular
	     values also approximate those of A.  Total (or complete)
	     re-orthogonalization is used for this Lanczos recursion.  
	     a point Lanczos method (single vector) is used if a blocksize
	     of 1 is encountered via normal deflation.

   Phase 3:  Apply an appropriate QR iteration to diagonalize B and
	     hence produce approximate singular values (array alpha)
	     of the original matrix A.

   Phase 4:  Convergence test using a user-supplied residual tolerance.

   Phase 5:  Iteration restart with orthogonalization with respect
	     to any (all) converged singular vectors.

   This version of blklan1 is designed to approximate the ik-largest
   singular triplets of A.  Users interested in the ik-smallest
   singular triplets need only sort the alpha array in increasing
   (as opposed to the default ascending order) following the call to
   qriter2 in Phase 3.  Also, the rows of the two-dimensional arrays
   ppp and qqp must be reordered to reflect a one-to-one correspondence
   with the newly sorted elements of alpha (which are approximate
   singular values of the matrix A).


   Arguments
   ---------

   (input)
   fp       pointer to output file
   nnzero   number of nonzeros in matrix A
   m        row dimension of the sparse matrix A whose SVD is sought
   n        column dimension of the sparse matrix A whose SVD is sought
   ik       number of singular triplets desired
   ib       initial block size for outer iteration
   ic       upper bound for dimension of Krylov subspace generated via
	      outer iteration.  It is the maximum dimension for the block
	      upper-bidiagonal matrix S generated in Phase 1 above.
   tol     user-specified tolerance for approximate singular triplets
   maxit   maximum number of outer iterations allowed

   (output)
   iko     number of singular triplets approximated
   ico     last bound for dimension of Krylov subspace used within outer 
	      iteration
   ibo     final block size used in outer iteration
   sing    linear array containing the iko approximate singular values
   v       two-dimensional array containing the iko approximate right 
	      singular vectors corresponding to the approximate singular 
	      values in array sing
   u       two-dimensional array containing the iko approximate left
	   singular vectors corresponding to the approximate singular values 
	   in array sing
   res     linear array containing the iko residuals of the approximate 
	   singular triplets
   memory  memory storage needed in bytes


   External parameters
   -------------------

   Defined and documented in bls1.h


   Local parameters
   ----------------

    k        current # of desired triplets (exit when reduced to 0)
    k0       count of triplets found in current iteration
    nb       current block size
    nc       size of current subspace
    ns       number of blocks in current iteration


   Functions called
   --------------

   BLAS         dgemv, dgemm, qriter2, orthg
   MISC         validate, random, imin
   BLS1         point1, block1


   NOTE:  Unlike Fortran, C stores arrays in row-major order.  For the
	  sake of efficiency, most matrices and vectors in the algorithm
	  are transposed, i.e., each vector and each column of a
	  matrix are stored row-wise.

 ***********************************************************************/

long   blklan1(FILE *fp, long nnzero, long m, long n, long ik, double *v,
	       double *u, double *sing, long ic, long ib, double tol,
               double *res, long maxit, long *iko, long *ico, long *ibo, 
	       long *memory)

{
   long nbold, cont, irand, ns, nb, nc, k, k0, length, memsize;
   long jj, ii, kk, ll, final, flag, i, j, *index;
   double **wp, **sp, **rp, **bigsp, **workptr2, **tempptr2;
   double *w, *s, *r, *bigs, *workptr1, *tempptr1;
   double *ptr1, *ptr2, *ptr3, *ptr4;

   nb = ib;
   nc = ic;
   k  = ik;
   k0 = 0;
   ns = nc / nb;
   iconv = 0;

   irand = SEED;

   if (ns < MINBLKS) {
      nb = nc / MINBLKS;
      ns = nc / nb;
   }

   /* return upon error on input */
   if (validate(fp, nnzero, m, n, *ico, *ibo, ik, tol))
      return(-1);

   memsize = sizeof(double) * 
	     (nb * (m + nc + nc - nb + m + n + 1) +
              nc * (nb + m + n + m + 3 + nc + nc) +
	      ik * (m + m + n) + n + n + m + m);                   

   *memory += memsize;

 /***********************************************************************
  *          Allocate work area and initialize pointers                 *
  *     pointer              size                                       *
  *    w     (wp)          nb      by m                                 *
  *    s     (sp)          nc      by nb                                *
  *    r     (rp)          (nc-nb) by nb                                *
  *    bigs  (bigsp)       nc      by (nb+1)                            *
  *    temp                m       by nb                                *
  *    y     (yp)          nb      by n                                 *
  *    uu    (uup)         nc      by m                                 *
  *    vv    (vvp)         nc      by n                                 *
  *    u0                  ik      by m                                 *
  *    v0                  ik      by n                                 *
  *    tres                nb                                           *
  *    uvtmp (uvtmpp)      (nc+ik) by m  (assume max (m,n) = m)         *
  *    z                   n                                            *
  *    p                   n                                            *
  *    q                   m                                            *
  *    t                   m                                            *
  *    alpha               nc                                           *
  *    beta                nc                                           *
  *    pp    (ppp)         nc      by nc                                *
  *    qq    (qqp)         nc      by nc                                *
  *    u                   ik      by m (allocated in bls1.c)           *
  *    v                   ik      by n (allocated in bls1.c)           *
  ***********************************************************************/

   if (!(workptr1 = (double *)malloc(memsize)) ||
       !(index   = (long   *)malloc(sizeof(long) * ib))){
      perror("FIRST MALLOC FAILED in BLKLAN1()");
      exit(errno);
   }

   *memory += sizeof(long) * ib;
   memsize = sizeof(double *) * (8 * nc + nb + ik);
   *memory += memsize;
   if (!(workptr2 = (double **)malloc(memsize))){
      perror("SECOND MALLOC FAILED in BLKLAN1()");
      exit(errno);
   }
   tempptr1  = workptr1;
   tempptr2  = workptr2;

   length    = m * nb;
   w         = tempptr1;
   tempptr1 += length;
   wp        = tempptr2;
   tempptr2 += nb;
   j = 0;
   for (i = 0; i < length; i += m) wp[j++] = &w[i];

   length    = nc * nb;
   s         = tempptr1;
   tempptr1 += length;
   sp        = tempptr2;
   tempptr2 += nc;
   j = 0;
   for (i = 0; i < length; i += nb) sp[j++] = &s[i];

   length = (nc - nb) * nb;
   r         = tempptr1; 
   tempptr1 += length;
   rp        = tempptr2;
   tempptr2 += (nc - nb);
   j = 0;
   for (i = 0; i < length; i += nb) rp[j++] = &r[i];

   length    = nc * (nb + 1);
   bigs      = tempptr1;
   tempptr1 += length;
   bigsp     = tempptr2;
   tempptr2 += nc;
   j = 0;
   for (i = 0; i < length; i += nb + 1) bigsp[j++] = &bigs[i];

   temp     = tempptr1;
   tempptr1 += m * nb;

   length    = n * nb;
   y         = tempptr1;
   tempptr1 += length;
   yp        = tempptr2;
   tempptr2 += nb;
   j = 0;
   for (i = 0; i < length; i += n) yp[j++] = &y[i];

   length    = m * nc;
   uu        = tempptr1;
   tempptr1 += length; 
   uup       = tempptr2;
   tempptr2 += nc;
   j = 0;
   for (i = 0; i < length; i += m) uup[j++] = &uu[i];

   length    = n * nc;
   vv        = tempptr1;
   tempptr1 += length;
   vvp       = tempptr2;
   tempptr2 += nc;
   j = 0;
   for (i = 0; i < length; i += n) vvp[j++] = &vv[i];

   u0        = tempptr1;
   tempptr1 += m * ik;

   v0        = tempptr1;
   tempptr1 += n * ik;

   tres      = tempptr1;
   tempptr1 += nb;

   length    = m * (nc + ik);
   uvtmp     = tempptr1;
   tempptr1 += length;
   uvtmpp    = tempptr2;
   tempptr2 += nc + ik;
   j = 0;
   for (i = 0; i < length; i += m) uvtmpp[j++] = &uvtmp[i];

   z         = tempptr1;
   tempptr1 += n;

   p         = tempptr1;
   tempptr1 += n;

   q         = tempptr1;
   tempptr1 += m;

   t         = tempptr1;
   tempptr1 += m;

   alpha     = tempptr1;
   tempptr1 += nc;

   beta      = tempptr1;
   tempptr1 += nc;

   length    = nc * nc;
   pp        = tempptr1;
   tempptr1 += length;
   qq        = tempptr1;
   tempptr1 += length;
   ppp       = tempptr2;
   tempptr2 += nc;
   qqp       = tempptr2;
   tempptr2 += nc;
   j = 0;
   for (i = 0; i < length; i += nc) {
      ppp[j]   = &pp[i];
      qqp[j++] = &qq[i];
   }

   /* choose an initial random V[1] matrix */
   length = n * nb;
   for (i = 0; i < length; i++) vv[i] = random(&irand);
   orthg(nb, 0, n, yp, vvp, temp); 

   /* initialize iteration counter */
   iter = 0;

   while (iter < maxit) {
      nn = nb * ns;
      cont = TRUE;
      iter += 1;

      /*------------------------------------------------------------------*
       *              PHASE 1 and PHASE 2 (block algorithm)               *
       *------------------------------------------------------------------*/
      if (nb > 1) 
	 block1(w, wp, sp, rp, bigsp, m, n, nb, ns, &irand);
      else {
	 if (nbold != 1) k0 = 0;
         /*------------------------------------------------------------------*
          *                   PHASE 2A (point algorithm)                     *
          *------------------------------------------------------------------*/
	 cont = point1(&k, m, n, wp, res, sing, tol);
         if (cont) {
	    for (i = 0; i < nn; i++) {
	       for (j = 0; j < nn; j++) {
	          ppp[i][j] = ZERO;
	          qqp[i][j] = ZERO;
               }
            }
	    for (i = 0; i < nn; i++) {
	       ppp[i][i] = ONE;
	       qqp[i][i] = ONE;
            }
	 }
      }

      /*------------------------------------------------------------------*
       *                        PHASE 3                                   *
       *------------------------------------------------------------------*/

      if (!cont) break;
      qriter2(nn, alpha, beta, qqp, ppp);

      /*------------------------------------------------------------------*
       *                        PHASE 4                                   *
       *------------------------------------------------------------------*/

      nbold = nb;
      if (iter > 1 && nb > 1) {
	 k0 = 0;
	 final = imin(nb, k);
	 for (i = 0; i < final; i++) {
	    if (fabs(tres[i]) <= tol) {
	      index[k0] = i;
	      sing[iconv + k0] = alpha[i];
	      res[iconv + k0] = tres[i];
	      k0 += 1;
            }
         }
	 if (nb >= k) nb -= k0;
	 else nb = imin(nb, k - k0);
	 nc -= k0;
	 k -= k0;
	 if (k) ns = nc / nb;
         if (ns < MINBLKS) {
            nb = nc / MINBLKS;
            ns = nc / nb;
         }
      }

      /*------------------------------------------------------------------*
       *                  PHASE 5 (back transformation)                   *
       *------------------------------------------------------------------*/

      if (nbold > 1) {

         /* uu[nn x m],  qq[(nb+k0) x nn] */

         dgemm(NTRANSP, NTRANSP, nb + k0, m, nn, ONE, qqp, uu, ZERO, u);
         dgemm(NTRANSP, NTRANSP, nb + k0, n, nn, ONE, ppp, vv, ZERO, v);
      }
      else {
         dgemv(TRANSP, nn, m, ONE, uup, qq, ZERO, u);
         dgemv(TRANSP, nn, n, ONE, vvp, pp, ZERO, v);
      }
      if (k0) {
         final = iconv + k0;
         ii = 0;
         for (jj = iconv; jj < final; jj++) {
	    ptr1 =  &u[m * index[ii]];
	    ptr2 = &u0[m * jj];
            for (i = 0; i < m; i++) *ptr2++ = *ptr1++;
	    ptr1 =  &v[n * index[ii++]];
	    ptr2 = &v0[n * jj];
            for (i = 0; i < n; i++) *ptr2++ = *ptr1++;
         }
         iconv = final;
         if (k <= 0) {
	    iter -= 1;
	    cont = FALSE;
         }
         /* reload unconverged right S-vector */
         if (cont) {
            kk = 0;
            final = nb + k0;
            for (jj = 0; jj < final; jj++) {
	       flag = FALSE;
	       ll = 0;
	       while (ll < k0 && !flag) {
	          if (index[ll] != jj) ll++;
	          else flag = TRUE;
	       }
	       if (!flag) {
	          ptr1 = &vv[n * kk];
	          ptr2 =  &v[n * jj];
                  for (i = 0; i < n; i++) *ptr1++ = *ptr2++;
	          kk++;
	       }
            }
	 }
      }
      else { 
         ptr1 =  v;
         for (jj = 0; jj < nb; jj++)
            for (i = 0; i < n; i++) vvp[jj][i] = *ptr1++; 
      }
      if (!cont) break;
   }

   /* return all converged S-vectors (right and left) */
   ptr1 =  u;
   ptr2 = u0;
   ptr3 =  v;
   ptr4 = v0;
   for (j = 0; j < iconv; j++) {
      for (i = 0; i < m; i++) *ptr1++ = *ptr2++;
      for (i = 0; i < n; i++) *ptr3++ = *ptr4++;
   }
   *iko = ik - k;
   *ibo  = nb;
   *ico = nc;

   free(workptr1);
   free(workptr2);
   free(index);
   return(0);
}
#define       ZERO       0.0
#define       ONE        1.0
#define       TRANSP     1
#define       NTRANSP    0

extern double *tres, *temp, *uu, *vv, *u0, *v0, *y, *uvtmp, *pp, *qq;
extern double **uup, **vvp, **uvtmpp, **yp, **ppp, **qqp;
extern double *p, *q, *t, *z, *alpha, *beta;
extern long iconv, nn, iter;

double ddot(long, double *,long, double *, long);
double enorm(long, double *);
void   daxpy(long, double, double *,long, double *, long);
double random(long *);
void   dgemm(long, long, long, long, long, 
 	     double, double **, double *, double, double *);
void   dtbmv(long trans, long n, long k, double **a, double *x);
void   orthg(long, long, long, double **, double **, double *);
void   formbigs(long, long, double **, double **, double **);
void   opat(long, double *, double *);


/***********************************************************************
 *                                                                     *
 *                        block1()                                     *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   This function implements the first two phases of the hybrid block
   Lanczos procedure.  In the first phase, which is also known as the
   block Lanczos outer iteration, a block upper-bidiagonal matrix S
   is formed.  The singular values of the matrix S approximate those 
   of the original sparse matrix A.  Total (or complete) re-
   orthogonalization is used.

   In the second phase, single-vector Lanczos bidiagonalization is used
   to reduce (preserving singular values of S) the block matrix S to a 
   bi-diagonal matrix B. 


   Arguments
   ---------

   (input)
   w, wp     work space
   sp        diagonal blocks (upper-triangular submatrices) of the block 
		upper-bidiagonal matrix S 
   rp        super-diagonal blocks (upper-triangular submatrices) of the 
		block upper-bidiagonal matrix S 
   bigsp     block upper-bidiagonal matrix S
   m         row dimension of sparse matrix A
   n         column dimension of sparse matrix A
   nb        current block size
   ns        number of blocks in current iteration
   irand     seed for random number generator
   

   (output - globally defined)
   alpha     diagonal elements of bidiagonal matrix B (reduced from 
	        matrix S)
   beta      off-diagonal elements of bidiagonal matrix B (reduced from 
		matrix S)
   tres      residuals of approximate singular triplets determined from
	     a previous set of block Lanczos outer iterations.
   ppp       array of left singular vectors of S
   qqp       array of right singular vectors of S

   External parameters
   -------------------

   Defined and documented in bls1.h


   Functions called
   --------------

   BLAS         ddot, daxpy, enorm, dgemm, orthg, dtbmv
   USER         opat
   MISC         random
   BLS1         formbigs

 ***********************************************************************/
void block1(double *w, double **wp, double **sp, double **rp, double **bigsp,
      long m, long n, long nb, long ns, long *irand)

{
   long jinc, nk, nj, i, j, k, blkptr;
   double *ptr, dum, pnm;

   /*------------------------------------------------------------------*
    *                        PHASE 1                                   *
    *------------------------------------------------------------------*/
   /* ns (number of blocks) is assumed to be at least 2 */
   nk = nn - nb;

   /* initialize S and R arrays */
   for (i = 0; i < nn; i++)
      for (j = 0; j < nb; j++) sp[i][j] = ZERO;

   for (i = 0; i < nk; i++)
      for (j = 0; j < nb; j++) rp[i][j] = ZERO;

   for (i = 0; i < nb; i++) 
      opa(m, n, vvp[i], wp[i]);
   
   ptr = uvtmp;
   nk = iconv * m;
   if (iconv) 
      for (i = 0; i < nk; i++) *ptr++ = u0[i]; 

   nj = nb * m;
   for (i = 0; i < nj; i++) {
      *ptr++ = w[i];
      w[i] = ZERO;
   }

   orthg(nb, iconv, m, wp, uvtmpp, temp);

   for (i = 0; i < nb; i++) 
      for (j = i; j < nb; j++) sp[i][j] = wp[j][i];

   ptr = uvtmp + nk;
   for (i = 0; i < nj; i++) uu[i] = *ptr++;

   blkptr = 0;

   /* iterate for j = 1 to (ns-1) */ 
   for (j = 1; j < ns; j++) {

      /* compute Y[j] = A^T * U[j] - V[j] * S[j]^T */
      for (i = 0; i < nb; i++) 
         opat(n, uup[blkptr + i], yp[i]);

      /* compute -S * VVP + Y, Y is [nb,n] */
      dgemm(NTRANSP, NTRANSP, nb, n, nb, -ONE, &sp[blkptr],
            vvp[blkptr], ONE, y);

      /* store residuals in Y for convergence test in Phase 4 */
      if (j == 1 && iter > 1)
         for (i = 0; i < nb; i++) tres[i] = enorm(n, yp[i]);

      if (iconv) {
         /* orthogonalize Y[j] w.r.t. [V[0], ... V[j] | V0] */ 
         nk = nb * j;
         nj = nk + iconv;

         ptr = vv;
         for (k = 0; k < nk; k++) 
            for (i = 0; i < n; i++) uvtmpp[k][i] = *ptr++; 
         ptr = v0;
         for (k = nk; k < nj; k++) 
            for (i = 0; i < n; i++) uvtmpp[k][i] = *ptr++; 
      }
      else {
         /* orthogonalize Y[j] w.r.t. [V[0], ... V[j]] */ 
         nj = nb * j;
         ptr = vv;
         for (k = 0; k < nj; k++) 
            for (i = 0; i < n; i++) uvtmpp[k][i] = *ptr++;
      }
      ptr = y;
      for (k = 0; k < nb; k++) 
         for (i = 0; i < n; i++) {
	    uvtmpp[nj + k][i] = *ptr;
	    *ptr++ = ZERO;
         }
      orthg(nb, nj, n, yp, uvtmpp, temp);

      /* get QR factorization: Y[j] = V[j] * R[j], [R[j] := upper tri. ] */
      /* load R[j] submatrix */
      for (k = 0; k < nb; k++) 
         for (i = k; i < nb; i++) rp[blkptr + k][i] = yp[i][k];

      /* increment pointer for V[j] */
      jinc = blkptr + nb;

      /* load V[j] submatrix */
      ptr = vvp[jinc];
      for (k = nj; k < nj + nb; k++) 
         for (i = 0; i < n; i++) *ptr++ = uvtmpp[k][i]; 
  
      /* compute W[j] = A * V[j] - U[j] * R[j]^T */
      for (k = 0; k < nb; k++)
         opa(m, n, vvp[jinc + k], wp[k]);

      dgemm(NTRANSP, NTRANSP, nb, m, nb, -ONE, &rp[blkptr],
            uup[blkptr], ONE, w);

      if (iconv) {
         /* orthogonalize W[j] w.r.t. [U[0], ... U[j] | U0] */ 
	 nk = nb * j;
	 nj = nk + iconv;
	 ptr = uu;

         for (k = 0; k < nk; k++) 
            for (i = 0; i < m; i++) uvtmpp[k][i] = *ptr++; 
         ptr = u0;
         for (k = nk; k < nj; k++) 
            for (i = 0; i < m; i++) uvtmpp[k][i] = *ptr++; 
      }
      else {
         /* orthogonalize W[j] w.r.t. [U[0], ... U[j]] */ 
         nj = nb * j;
         ptr = uu;
         for (k = 0; k < nj; k++) 
            for (i = 0; i < m; i++) uvtmpp[k][i] = *ptr++;
      }
      ptr = w;
      for (k = 0; k < nb; k++) 
         for (i = 0; i < m; i++) {
	    uvtmpp[nj + k][i] = *ptr;
	    *ptr++ = ZERO;
         }
      orthg(nb, nj, m, wp, uvtmpp, temp);

      /* get QR factorization: W[j] = U[j] * S[j], [S[j] := upper tri. ] */
      /* load S[j] submatrix */
      for (k = 0; k < nb; k++) 
         for (i = k; i < nb; i++) sp[k + jinc][i] = wp[i][k];

      ptr = uup[jinc];
      for (k = 0; k < nb; k++) 
         for (i = 0; i < m; i++) *ptr++ = uvtmpp[nj + k][i];

      /* compute pointer to next block */
      blkptr += nb;
   }
   /* form the block upper-triangular S matrix in BIGS array;
    * array ppp will store left S-vectors of BIGS
    * array qqp will store right S-vectors of BIGS            */

   formbigs(nn, nb, rp, sp, bigsp);

   /*------------------------------------------------------------------*
    *                        PHASE 2 (block size > 1)                  *
    *------------------------------------------------------------------*/

   /* choose starting vector P[0] having unity 2-norm */
   for (i = 0; i < nn; i++) p[i] = random(irand);
   pnm = enorm(nn,p);
   ptr = pp;
   for (i = 0; i < nn; i++) {
      p[i] /= pnm;
      t[i] = p[i];
      *ptr++ = p[i];
   }
   /* compute T[0] := S * P[0] */
   dtbmv(NTRANSP, nn, nb, bigsp, t);

   /* compute alpha[0] */
   alpha[0] = enorm(nn,t);

   for (j = 0; j < nn; j++) {
      if (alpha[j] != ZERO) {
	 /* compute Q[j] := T[j] / alpha[j] */
         ptr = qqp[j]; 
         for (i = 0; i < nn; i++) {
	    q[i] = t[i] / alpha[j];
	    z[i] = q[i];
	    *ptr++ = q[i];
	 }
	 if (j == nn - 1) return;

	 /* compute Z[j] := S^T * Q[J] - alpha[j] * P[j] */
	 dtbmv(TRANSP, nn, nb, bigsp, z);
	 daxpy(nn, -alpha[j], p, 1, z, 1);

	 /* orthogonalize Z w.r.t. previous PP's */
         for (k = 0; k <= j; k++)
            for (i = 0; i < nn; i++) uvtmpp[k][i] = ppp[k][i];

         if (j) {
	    ptr = uvtmpp[j + 1];
            for (i = 0; i < nn; i++) *ptr++ = z[i];
	    orthg(1, j+1, nn, wp, uvtmpp, temp);
	    ptr = uvtmpp[j + 1];
            for (i = 0; i < nn; i++) z[i] = *ptr++;
	 }
	 else {
	    dum = -ddot(nn, uvtmp, 1, z, 1);
	    daxpy(nn, dum, uvtmp, 1, z, 1);
	 }
	 beta[j] = enorm(nn,z);
	 if (beta[j] != ZERO) {

	    /* compute p[j+1] := z[j] / beta[j] */
	    ptr = ppp[j + 1];
            for (i = 0; i < nn; i++) {
	       p[i] = z[i] / beta[j];
	       t[i] = p[i];
	       *ptr++ = p[i];
	    }
	    /* compute T[j+1] := S * P[j+1] - beta[j] * Q[j] */
	    dtbmv(NTRANSP, nn, nb, bigsp, t);
	    daxpy(nn, -beta[j], q, 1, t, 1);

	    /* orthogonalize Z w.r.t. previous QQ's */
            for (k = 0; k <= j; k++) 
               for (i = 0; i < nn; i++) uvtmpp[k][i] = qqp[k][i];
            if (j) {
	       ptr = uvtmpp[j + 1];
               for (i = 0; i < nn; i++) *ptr++ = t[i];
	       orthg(1, j+1, nn, wp, uvtmpp, temp);
	       ptr = uvtmpp[j + 1];
               for (i = 0; i < nn; i++) t[i] = *ptr++;
	    }
	    else {
	       dum = -ddot(nn, uvtmp, 1, t, 1);
	       daxpy(nn, dum, uvtmp, 1, t, 1);
	    }
	    /* compute alpha[j+1] */
	    alpha[j + 1] = enorm(nn, t);
	 }
	 else return;
      }
      else return;
   }
   return;
}
#define       ZERO       0.0
#define       ONE        1.0
#define       TRANSP     1
#define       NTRANSP    0
#define       CONTINUE   1
#define       DONE       0

extern double *u0, *v0, *temp, *uvtmp, *uu, *vv, *alpha, *beta;
extern double *p, *q, *t, *z;
extern double **uvtmpp, **uup, **vvp;
extern long iconv, nn, iter;

double ddot(long, double *,long, double *, long);
double enorm(long, double *);
void   daxpy(long, double, double *,long, double *, long);
void   opa(long, long, double *, double *);
void   opat(long, double *, double *);
void   orthg(long, long, long, double **, double **, double *);

/***********************************************************************
 *                                                                     *
 *                        point1()                                     *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   This function is a single-vector Lanczos bidiagonalization procedure
   (degenerate case of block size = 1 for block1.c) which is used when
   normal deflation reduces the current block size to 1.

   The function returns DONE when the number of remaining triplets to be 
   approximated is less than or equal to zero.  Otherwise, it returns
   status CONTINUE.
   

   Arguments
   ---------

   (input)
   k        current # of desired triplets
   m        row dimension of the sparse matrix A whose SVD is sought
   n        column dimension of the sparse matrix A whose SVD is sought
   tol     user-specified tolerance for approximate singular triplets
   wp      work space

   (output)
   sing    linear array containing the iko approximate singular values
   res     linear array containing the iko residuals of the approximate 
	   singular triplets
   alpha   diagonals of bidiagonal matrix from the inner recursion
   beta    super-diagonals of bidiagonal matrix from the inner recursion

   External parameters
   -------------------

   Defined and documented in bls1.h


   Functions called
   --------------

   BLAS         ddot, enorm, daxpy, orthg
   USER         opa, opat

 ***********************************************************************/

long point1(long *k, long m, long n, 
           double **wp, double *res, double *sing, double tol)

{
   double *ptr, dum, znm;
   long i, j, jj, convpj;

   for (j = 0; j < n; j++) p[j] = vv[j];
   opa(m, n, p, t);

   /* orthogonalize T w.r.t. converged left S-vectors (U0 array) */
   if (iconv) {
      ptr = u0;
      for (j = 0; j < iconv; j++)
         for (i = 0; i < m; i++) uvtmpp[j][i] = *ptr++;
      if (iconv > 1) {
	 ptr = uvtmpp[iconv];
         for (i = 0; i < m; i++) *ptr++ = t[i];
	 orthg(1, iconv, m, wp, uvtmpp, temp);
	 ptr = uvtmpp[iconv];
         for (i = 0; i < m; i++) t[i] = *ptr++;
      }
      else {
	 dum = -ddot(nn, uvtmp, 1, t, 1);
	 daxpy(nn, dum, uvtmp, 1, t, 1);
      }
   }
   /* compute alpha[0] */
   alpha[0] = enorm(m, t);
   for (j = 0; j < nn; j++) {
      if (alpha[j] != ZERO) {
	 
	 /* compute Q[j] := T[j] / alpha[j] */
	 ptr = uup[j];
         for (i = 0; i < m; i++) {
	    q[i] = t[i] / alpha[j];
	    *ptr++ = q[i];
	 }
	 if (j == nn - 1) return(CONTINUE);

	 /* compute Z[j] := A^T * Q[J] - alpha[j] * P[j] */
	 opat(n, q, z);
	 daxpy(n, -alpha[j], p, 1, z, 1);
	 if (!j) {
	    if ((znm = enorm(n, z)) <= tol) {
	       ptr = &u0[iconv * m];
               for (i = 0; i < m; i++) *ptr++ = q[i];
	       ptr = &v0[iconv * n];
               for (i = 0; i < n; i++) *ptr++ = p[i];
	       sing[iconv] = alpha[0];
	       res[iconv] = znm;
	       iconv += 1;
               *k -= 1;
	       if (!*k) {
		  iter -= 1;
		  return(DONE);
	       }
            }
	 }
         convpj = iconv + j;
	 ptr = v0;

	 /* orthogonalize Z w.r.t. converged right S-vectors and 
	  * previous VV's */
	 for (jj = 0; jj < iconv; jj++)
	    for (i = 0; i < n; i++)
	       uvtmpp[jj][i] = *ptr++;
	 ptr = vv;
	 for (jj = iconv; jj <= convpj; jj++)
	    for (i = 0; i < n; i++)
	       uvtmpp[jj][i] = *ptr++;
         if (convpj) {
	    ptr = uvtmpp[convpj + 1];
	    for (i = 0; i < n; i++) *ptr++ = z[i];
	    orthg(1, convpj + 1, n, wp, uvtmpp, temp);
	    ptr = uvtmpp[convpj + 1];
	    for (i = 0; i < n; i++) z[i] = *ptr++;
	 }
	 else {
	    dum = -ddot(n, uvtmp, 1, z, 1);
	    daxpy(n, dum, uvtmp, 1, z, 1);
	 }
	 beta[j] = enorm(n,z);
	 if (beta[j] != ZERO) {

	    /* compute P[j+1] := Z[j] / beta[j] */
	    ptr = vvp[j + 1];
            for (i = 0; i < n; i++) {
	       p[i] = z[i] /beta[j];
	       *ptr++ = p[i];
	    }
	    /* compute T[j+1] := A * P[j+1] - beta[j] * Q[j] */
	    opa(m, n, p, t);
	    daxpy(m, -beta[j], q, 1, t, 1);

	    /* orthogonalize T w.r.t. converged left S-vectors and 
	     * previous UU's */

            convpj = iconv + j;
	    ptr = u0;
	    for (jj = 0; jj < iconv; jj++)
	       for (i = 0; i < m; i++)
	          uvtmpp[jj][i] = *ptr++;
	    ptr = uu;
	    for (jj = iconv; jj <= convpj; jj++)
	       for (i = 0; i < m; i++)
	          uvtmpp[jj][i] = *ptr++;
            if (convpj) {
	       ptr = uvtmpp[convpj + 1];
	       for (i = 0; i < m; i++) *ptr++ = t[i];
	       orthg(1, convpj + 1, m, wp, uvtmpp, temp);
	       ptr = uvtmpp[convpj + 1];
	       for (i = 0; i < m; i++) t[i] = *ptr++;
	    }
	    else {
	       dum = -ddot(m, uvtmp, 1, t, 1);
	       daxpy(m, dum, uvtmp, 1, t, 1);
	    }
	    alpha[j + 1] = enorm(m, t);
         }
	 else return(CONTINUE);
      }
      else return(CONTINUE);
   }
}
#include <stdio.h>
#include <fcntl.h>

#define  NCMAX  5200
#define  NZMAX  800000
#define  ZERO   0.0

/***********************************************************************
 *								       *
 *		              validate()			       *
 *								       *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------
   Function validates input parameters and returns error code (long)  

   Arguments 
   ---------
  (input)
  fp_out1   pointer to output file
  nnzero    number of nonzero elements in input sparse matrix A
  nrow      row dimension of A
  ncol      column dimension of A
  maxsubsp  maximum dimension of Krylov subspace allowed
  blksize   initial block size
  nums      number of singular values desired
  tol       user-specified tolerance for approximate singular triplets

 ***********************************************************************/

long validate(FILE *fp_out1, long nnzero, long nrow, long ncol,
	      long maxsubsp, long blksize, long nums, double tol)

{
   long error_index;
   char *error[8];

   error_index = 0;
   error[1] = " ***** SORRY, YOUR MATRIX IS TOO BIG *****";
   error[2] = " ***** NCOL MUST NOT BE GREATER THAN NROW *****";
   error[3] = " ***** TOLERANCE IS INVALID *****";
   error[4] = " ***** MAXIMUM SUBSPACE DIMENSION IS INVALID *****";
   error[5] = " ***** INITIAL BLOCK SIZE MUST BE GREATER THAN 1 *****";
   error[6] = " ***** NUMBER OF SINGULAR VALUES DESIRED IS INVALID *****";
   error[7] = " ***** INIT BLK SIZE MUST BE LESS THAN NO. OF S-VALUES DESIRED *****";

   if (ncol > NCMAX || nnzero > NZMAX)           error_index = 1;
   else if (ncol > nrow)                         error_index = 2;
   else if (tol < ZERO)                          error_index = 3;
   else if (maxsubsp > ncol || maxsubsp <= 0)    error_index = 4;
   else if (blksize <= 1 || blksize > maxsubsp)  error_index = 5;
   else if (nums > maxsubsp || nums <= 0)        error_index = 6;
   else if (blksize > nums)                      error_index = 7;

   if (error_index) fprintf(fp_out1, "%s\n", error[error_index]);
   return(error_index);
}
#define       ZERO    0.0
extern long mxvcount;
extern long *pointr, *rowind;
extern double *value;

/**************************************************************
 *                                                            *
 * multiplication of matrix A by vector x, where A is m by n  *
 * (m >> n) and is stored using the Harwell-Boeing compressed *
 * column sparse matrix format.  y stores product vector.     *
 *                                                            *
 **************************************************************/

void opa(long m, long n, double *x, double *y)

{
   long end,i,j;
   
   mxvcount += 1;
   for (i = 0; i < m; i++) y[i] = ZERO;

   for (i = 0; i < n; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++)
	 y[rowind[j]] += value[j] * x[i]; 
   }
   return;
}
#define       ZERO    0.0
extern long mtxvcount;
extern long *pointr, *rowind;
extern double *value;

/**************************************************************
 *                                                            *
 * multiplication of an n by m matrix A' by a vector X, store *
 * in Y.                                                      *
 *                                                            *
 **************************************************************/

void opat(long n, double *x, double *y)
{
   long end,i,j;
   
   mtxvcount += 1;
   for (i = 0; i < n; i++) y[i] = ZERO;

   for (i = 0; i < n; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++)
	 y[i] += value[j] * x[rowind[j]]; 
   }
   return;
}
double fsign(double a,double b)
/************************************************************** 
 * returns |a| if b is positive; else fsign returns -|a|      *
 **************************************************************/ 
{

   if ((a>=0.0 && b>=0.0) || (a<0.0 && b<0.0))return(a);
   if  ((a<0.0 && b>=0.0) || (a>=0.0 && b<0.0))return(-a);
}

double dmax(double a, double b)
/************************************************************** 
 * returns the larger of two double precision numbers         *
 **************************************************************/ 
{

   if (a > b) return(a);
   else return(b);
}

double dmin(double a, double b)
/************************************************************** 
 * returns the smaller of two double precision numbers        *
 **************************************************************/ 
{

   if (a < b) return(a);
   else return(b);
}

long imin(long a, long b)
/************************************************************** 
 * returns the smaller of two integers                        *
 **************************************************************/ 
{

   if (a < b) return(a);
   else return(b);
}

long imax(long a,long b)
/************************************************************** 
 * returns the larger of two integers                         *
 **************************************************************/ 
{

   if (a > b) return(a);
   else return(b);
}
#include <math.h>
#include <stdio.h>

#define		TRUE	 1
#define		FALSE  	 0
#define		TRANSP   1
#define		NTRANSP  0
#define		ZERO	 0.0
#define		ONE	 1.0
#define		CONST    100.0

void   dgemv(long, long, long, double, double **, double *, double, double *);
double ddot(long, double *,long, double *, long);
void   dscal(long, double, double *,long);
void   daxpy(long, double, double *,long, double *, long);
void   dcopy(long, double *, long, double *, long);

/***********************************************************************
 *                                                                     *
 *                        orthg()                                      *
 *         Gram-Schmidt orthogonalization procedure                    *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   The p by n matrix Z stored row-wise in rows f to (f+p-1) of
   array X is reorthogonalized w.r.t. the first f rows of array X.
   The resulting matrix Z is then factored into the product of a
   p by n orthonormal matrix (stored over matrix Z) and a p by p
   upper-triangular matrix (stored in the first p rows and columns 
   of array B).  (Based on orthog from Rutishauser) 


   Parameters
   ----------

   (input)
   p           number of consecutive vectors of array x (stored row-wise)
	       to be orthogonalized
   f           number of rows against which the next p rows are to be
	       orthogonalized
   n           column dimension of x
   x           2-dimensional array whose p rows are to be orthogonalized
	       against its first f rows
   temp        work array


   (output)
   x           output matrix whose f+p rows are orthonormalized
   b           p by p upper-triangular matrix


   Functions called
   --------------

   BLAS         dgemv, ddot, dscal, daxpy, dcopy

 ***********************************************************************/

void orthg(long p, long f, long n, double **b, double **x, double *temp)

{
   long fp, k, km1;
   long orig, small;
   double t, s;

   if (!p) return;
   if (f == 0 && p > n) {
      fprintf(stderr,"%s\n",
         "*** ON ENTRY TO ORTHG, MATRIX TO BE ORTHONORMALIZED IS SINGULAR");
      exit(-1);
   }
   fp = f + p;

   for (k = f; k < fp; k++) {
      km1 = k - 1;
      orig = TRUE;

      while(TRUE) {
         t = ZERO;

	 if (km1 >= 0) {
	    if (km1 > 0) {
	       dgemv(NTRANSP, k, n, ONE, x, x[k], ZERO, temp);
	       t += ddot(k, temp, 1, temp, 1);
	    }

	    else {
	       temp[0] = ddot(n, x[0], 1, x[k], 1);
	       t += temp[0] * temp[0];
	    }

	    if (orig && km1 >= f) 
               dcopy(k - f, &temp[f], 1, &b[k - f][0], 1); 

            if (km1 > 0) 
	       dgemv(TRANSP, k, n, -ONE, x, temp, ONE, &x[k][0]);
            else
	       daxpy(n, -temp[0], x[0], 1, x[k], 1);
         }

	 if (km1 < 0 || p != 1) {
	    s = ddot(n, x[k], 1, x[k], 1);
	    t += s;
	    if (s > t/CONST) {
	       small = FALSE;
	       s = sqrt(s);
               b[k - f][k - f] = s;
	       if (s != ZERO) s = ONE/s;
	       dscal(n, s, x[k], 1);
	    }
	    else {
	       small = TRUE;
	       orig  = FALSE;
	    }
	 }
	 if (small == FALSE || p == 1) break;
      }
   }
}
#include <stdio.h>

#define		TRANSP   1
#define		NTRANSP  0
#define		ZERO	 0.0
#define		ONE	 1.0

/***********************************************************************
 *                                                                     *
 *                         dgemv()                                     *
 * A C-translation of the level 2 BLAS routine DGEMV by Dongarra,      *
 * du Croz, and Hammarling, and Hanson (see LAPACK Users' Guide).      *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   dgemv() performs one of the matrix-vector operations

   y := alpha * A * x + beta * y  or  y := alpha * A' * x + beta * y

   where alpha and beta are scalars, X, Y are vectors and A is an
   m by n matrix.

void dgemv(long transa, long m, long n, 
           double alpha, double **a, double *x, double beta, double *y)

   Parameters
   ----------

   (input)
   transa   TRANSP indicates op(A) = A' is to be used in the multiplication
	    NTRANSP indicates op(A) = A is to be used in the multiplication

   m        on entry, m specifies the number of rows of the matrix A.
	    m must be at least zero.  Unchanged upon exit.

   n        on entry, n specifies the number of columns of the matrix A.
	    n must be at least zero.  Unchanged upon exit.

   alpha    a scalar multiplier.  Unchanged upon exit.

   a        matrix A as a 2-dimensional array.  Before entry, the leading
	    m by n part of the array a must contain the matrix A.

   x        linear array of dimension of at least n if transa = NTRANSP
	    and at least m otherwise.

   beta     a scalar multiplier.  When beta is supplied as zero then y
	    need not be set on input.  Unchanged upon exit.

   y        linear array of dimension of at least m if transa = NTRANSP
	    and at leat n otherwise.  Before entry with beta nonzero,
	    the array y must contain the vector y.  On exit, y is 
	    overwritten by the updated vector y.


 ***********************************************************************/

void dgemv(long transa, long m, long n, 
           double alpha, double **a, double *x, double beta, double *y)

{
   long info, leny, i, j;
   double temp, *ptrtemp;

   info = 0;
   if      ( transa != TRANSP && transa != NTRANSP ) info = 1;
   else if ( m < 0 ) 				     info = 2;
   else if ( n < 0 )				     info = 3;

   if (info) {
      fprintf(stderr, "%s %1ld %s\n",
      "*** ON ENTRY TO DGEMV, PARAMETER NUMBER",info,"HAD AN ILLEGAL VALUE");
      exit(info);
   }

   if (transa) leny = n;
   else        leny = m;

   if (!m || !n || (alpha == ZERO && beta == ONE))
      return;

   ptrtemp = y; 

   /* form Y := beta * Y */
   if (beta == ZERO) 
      for (i = 0; i < leny; i++) *ptrtemp++ = ZERO;
   else if (beta != ONE) 
      for (i = 0; i < leny; i++) *ptrtemp++ *= beta;

   if (alpha == ZERO) return;

   switch(transa) {

      /* form Y := alpha * A * X + Y */
      case NTRANSP:  for(i = 0; i < m; i++) {
                        ptrtemp = *a++;
		        temp = ZERO;
		        for(j = 0; j < n; j++) 
			   temp += *ptrtemp++ * x[j];
			y[i] += alpha * temp;
		     }
		     break;
		     
      /* form Y := alpha * A' * X + Y */
      case TRANSP:   for(i = 0; i < m; i++) { 
                        ptrtemp = *a++;
			if (x[i] != ZERO) {
			   temp = alpha * x[i];
			   for(j = 0; j < n; j++)
			      y[j] += temp * (*ptrtemp++);
			}
		     }
		     break;
   }
}
#include <stdio.h>

#define		TRANSP   1
#define		NTRANSP  0
#define		ZERO	 0.0
#define		ONE	 1.0

/***********************************************************************
 *                                                                     *
 *                         dgemm()                                     *
 *                                                                     *
 * A C-translation of the level 3 BLAS routine DGEMM by Dongarra,      *
 * Duff, du Croz, and Hammarling (see LAPACK Users' Guide).            *
 * In this version, two of the three arrays which store the matrices   *
 * used in this matrix-matrix multiplication are accessed as linear    *
 * arrays.                                                             *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   dgemm() performs one of the matrix-matrix operations

	      C := alpha * op(A) * op(B) + beta * C,

   where op(X) = X or op(X) = X', alpha and beta are scalars, and A, B
   and C are matrices, with op(A) an m by k matrix, op(B) a k by n
   matrix and C an m by n matrix.

   Note that the arrays storing matrices B and C are linear arrays while
   the array of A is two-dimensional.


   Parameters
   ----------

   (input)
   transa   TRANSP indicates op(A) = A' is to be used in the multiplication
	    NTRANSP indicates op(A) = A is to be used in the multiplication

   transb   TRANSP indicates op(B) = B' is to be used in the multiplication
	    NTRANSP indicates op(B) = B is to be used in the multiplication

   m        on entry, m specifies the number of rows of the matrix op(A)
	    and of the matrix C.  m must be at least zero.  Unchanged
	    upon exit.

   n        on entry, n specifies the number of columns of the matrix op(B)
	    and of the matrix C.  n must be at least zero.  Unchanged
	    upon exit.

   k        on entry, k specifies the number of columns of the matrix op(A)
	    and the number of rows of the matrix B.  k must be at least 
	    zero.  Unchanged upon exit.

   alpha    a scalar multiplier

   a        matrix A as a 2-dimensional array.  When transa = NTRANSP, the
	    first k columns of the first m rows must contain the matrix A.
	    Otherwise, the first m columns of the first k rows must contain
	    the matrix A.

   b        matrix B as a linear array.  The leading (k * n) elements of
	    b must contain the matrix B.

   beta     a scalar multiplier.  When beta is supplied as zero then C
	    need not be set on input.

   c        matrix C as a linear array.  Before entry, the leading (m * n)
	    elements of c must contain the matrix C except when beta = 0.
	    In this case, c need not be set on entry.
	    On exit, c is overwritten by the (m * n) elements of matrix
	    (alpha * op(A) * op(B) + beta * C).

 ***********************************************************************/

void dgemm(long transa, long transb, long m, long n, long k, 
           double alpha, double **a, double *b, double beta, double *c)

{
   long info;
   long i, j, l, nrowa, ncola, nrowb, ncolb, nc;
   double temp, *atemp, *btemp1, *ptrtemp, *ctemp;

   info = 0;
   if      ( transa != TRANSP && transa != NTRANSP ) info = 1;
   else if ( transb != TRANSP && transb != NTRANSP ) info = 2;
   else if ( m < 0 ) 				     info = 3;
   else if ( n < 0 )				     info = 4;
   else if ( k < 0 )        			     info = 5;

   if (info) {
      fprintf(stderr, "%s %1ld %s\n",
      "*** ON ENTRY TO DGEMM, PARAMETER NUMBER",info,"HAD AN ILLEGAL VALUE");
      exit(info);
   }

   if (transa) {
      nrowa = k;
      ncola = m;
   }
   else { 
      nrowa = m;
      ncola = k;
   }
   if (transb) {
      nrowb = n;
      ncolb = k;
   }
   else {
      nrowb = k;
      ncolb = n;
   }
   nc = m * n;

   if (!m || !n || ((alpha == ZERO || !k) && beta == ONE))
      return;

   ctemp = c; 
   if (alpha == ZERO) {
      if (beta == ZERO)
         for (i = 0; i < nc; i++) *ctemp++ = ZERO;
      else if (beta != ONE)
         for (i = 0; i < nc; i++) *ctemp++ *= beta;
      return;
   }

   if (beta == ZERO) 
      for (i = 0; i < nc; i++) *ctemp++ = ZERO;
   else if (beta != ONE) 
      for (i = 0; i < nc; i++) *ctemp++ *= beta;

   if (!transb) { 

      switch(transa) {

	 /* form C := alpha * A * B + beta * C */
	 case NTRANSP:  ptrtemp = c;
		        for(l = 0; l < nrowa; l++) {
                           atemp = *a++;
      	       		   btemp1 = b;
      		  	   for(j = 0; j < ncola; j++) {
	 	     	      temp = *atemp * alpha;
	             	      ctemp = ptrtemp;
         	     	      for(i = 0; i < ncolb; i++) 
	    			 (*ctemp++) += temp * (*btemp1++);
	 	     	      atemp++;
      		  	   }
   	       		   ptrtemp = ctemp;
   	    	        }
			break;

	 /* form C := alpha * A' * B + beta * C */
	 case TRANSP:   ptrtemp = b;
	                for(l = 0; l < nrowa; l++) {
                           atemp = *a++;
      	       		   ctemp = c;
      		  	   for(j = 0; j < ncola; j++) {
	 	     	      temp = *atemp * alpha;
	             	      btemp1 = ptrtemp;
         	     	      for(i = 0; i < ncolb; i++) 
	    			 (*ctemp++) += temp * (*btemp1++);
	 	     	      atemp++;
      		  	   }
   	       		   ptrtemp = btemp1;
   	    		}
			break;
      }
   }
   else { 
      ctemp = c;

      switch(transa) {

	 /* form C := alpha * A * B' + beta * C */
	 case NTRANSP: for(l = 0; l < nrowa; l++) {
      	       		   btemp1 = b;
      		  	   for(j = 0; j < nrowb; j++) {
	 	     	      atemp = *a;
         	     	      for(i = 0; i < ncolb; i++) 
	    			 *ctemp += (*atemp++) * alpha * (*btemp1++);
	 	     	      ctemp++;
      		  	   }
			   a++;
   	    		}
			break;

	 /* form C := alpha * A' * B' + beta * C */
	 case TRANSP:   for(i = 0; i < ncola; i++) {
			   btemp1 = b;
			   for (l = 0; l < nrowb; l++) {
      	       		      temp = ZERO;
	 		      for(j = 0; j < nrowa; j++) 
			         temp += a[j][i] * (*btemp1++);
	    		      *ctemp++ += alpha * temp;
			   }
   	    		}
			break;
      }
   }
}
#include <math.h>

#define       ZERO       0.0
#define       ONE        1.0
#define       RDWARF     3.834e-20
#define       RGIANT     1.304e19

/***********************************************************************
 *                                                                     *
 *                        enorm()                                      *
 *  a C translation of the Fortran-77 version by Burton, Garbow,       *
 *  Hillstrom and More of Argonne National Laboratory.                 *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   given an n-vector x, this function calculates the Euclidean norm of x.
   The Euclidean norm is computed by accumulating the sum of squares in 
   three different sums.  The sums of squares for the small and large 
   components are scaled so that no overflows occur.  Non-destructive 
   underflows are permitted.  Underflows and overflows do not occur in the 
   computation of the unscaled sum of squares for the intermediate components. 
   The definitions of small, intermediate and large components depend on two
   constants, rdwarf and rgiant.  The restrictions on these constants are 
   that rdwarf**2 not underflow and rgiant**2 not overflow.  The constants
   given here are suitable for every known computer.
   The function returns the Euclidean norm of vector x in double precision.


   Parameters
   ----------

   n         number of elements in vector x
   x         linear array of vector x whose Euclidean norm is to be 
		calculated

 ***********************************************************************/

double enorm(long n, double *x)

{
   double norm2, agiant, floatn, s1, s2, s3, xabs, x1max, x3max;
   long i;

   s1     = ZERO;
   s2     = ZERO;
   s3     = ZERO;
   x1max  = ZERO;
   x3max  = ZERO;
   floatn = (double)n;
   agiant = RGIANT / floatn;

   for (i = 0; i < n; i++) {
      xabs = fabs(x[i]);
      /* summing components of vector that need no scaling */
      if (xabs > RDWARF && xabs < agiant)
	 s2 += xabs * xabs;
      else {
	 /* underflow... */
	 if (xabs <= RDWARF) {
            if (xabs > x3max) {
	       s3 = ONE + s3 * (x3max/xabs) * (x3max/xabs);
	       x3max = xabs;
	    }
	    else if (xabs != 0)
	       s3 += (xabs/x3max) * (xabs/x3max);
	 }
	 /* overflow... */
	 else {
	    /* summing large components of vector */
	    if (xabs <= x1max)
	       s1 += (xabs/x1max) * (xabs/x1max);
            else {
	       s1 = ONE + s1 * (x1max/xabs) * (x1max/xabs);
	       x1max = xabs;
	    }
	 }
      }
   }
   if (s1 != ZERO)
      norm2 = x1max * sqrt(s1 + (s2/x1max) / x1max);
   else if (s2 != ZERO) {
      if (s2 >= x3max)
	 norm2 = sqrt(s2 * (ONE + (x3max/s2) * (x3max*s3)));
      else 
	 norm2 = sqrt(x3max * ((s2/x3max) + (x3max*s3)));
   }
   else
      norm2 = x3max * sqrt(s3);
   return(norm2);
}
#define         ZERO     0.0

/***********************************************************************
 *                                                                     *
 *                        formbigs()                                   *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   This function forms the block upper-bidiagonal or the symmetric block
   tridiagonal matrix S from the block Lanczos algorithm in Phase 1 of 
   blklan1.c or blklan2.c, respectively. 


   Arguments
   ---------

   (input)
   r, s    submatrices from which the bidiagonal block matrix S
	   (Phase 1 of blklan1.c) is formed.
           The following data structure is assumed for the submatrices
	   s[j] and r[j], where j = 0, 1, ..., p-1.  For blklan1.c,
	   s[j] and r[j] are both upper-triangular.  For blklan2.c,
	   s[j] is dense and symmetric and r[j] is upper-triangular.

   p       number of block partitions used in blklan1.c.

		    s = s[0]               r= r[0]
			----                  ----
			s[1]                  r[1]
			----                  ----
			  .                     .
                          .                     .
			  .                     .
			s[p-1]                r[p-2]

   n       dimension of bigs

   (output)
   bigs    The 2-dimensional array bigs will contain this matrix in the 
	   following band matrix format:

           EXAMPLE WITH 4 SUP. DIAGONALS:
  
             TRANSPOSE OF   [0 0 0 0--1ST SUP. DIAGONAL----]
             ------------   [0 0 0 ---2ND SUP. DIAGONAL----]
                            [0 0 -----3RD SUP. DIAGONAL----]
                            [0 -------4TH SUP. DIAGONAL----]
                            [-------- MAIN DIAGONAL--------]
  
    Note that the super-diagonals and main diagonal of S are stored
    in COLUMNS of bigs (bigs is n by n).  

 ***********************************************************************/

void formbigs(long n, long ndp, double **r, double **s, double **bigs)

{
   long p, i, j, k, kk, m, row, col;

   p = n / ndp;

   /* load main diagonal of bigs */
   j = 0;
   for (i = 0; i < p; i++) 
      for (k = 0; k < ndp; k++) { 
	 bigs[j][ndp] = s[j][k];
	 j++;
      }

   /* load super-diagonals of bigs (from top to bottom) */
   for (i = 0; i < ndp; i++) {

      /* pad zeros at start of a column */
      for (kk = 0; kk < ndp - i; kk++) 
	 bigs[kk][i] = ZERO;

      /* load first row of bigs with main diagonals of r[j] */
      if (i == 0) {
	 j = 0;
         for (m = 0; m < p - 1; m++) 
            for (k = 0; k < ndp; k++) 
	       bigs[kk++][0] = r[j++][k];
      }
      else {
	 m = 0;
         for (j = 0; j < p; j++) {
	    row = m;
	    col = ndp - i;
	    /* load elements form s[j] submatrices */
	    while (col < ndp)
	       bigs[kk++][i] = s[row++][col++];

	    /* load elements form r[j] submatrices */
            if (j < p - 1) {
	       col = i;
	       row = m;
	       while (col < ndp)
	         bigs[kk++][i] = r[row++][col++];
            }
            m += ndp;
         }
      }
   }
}
#include <stdio.h>

#define         TRANSP   1
#define         NTRANSP  0

long imax(long, long);

/***********************************************************************
 *                                                                     *
 *                          dtbmv()                                    *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   The function performs one of the matrix-vector operations

	x := A * x,  or  x := A' * x,

   where A is an upper-triangular matrix.

   Parameters
   ----------

   trans     if trans = TRANSP, A' is to be used in the multiplication
             if trans = NTRANSP, A is to be used in the multiplication

   n         number of rows of matrix A; n must be at least 0.  Unchanged
	     upon exit.

   k         number of super-diagonals of matrix A

   a         2-dimensional array whose leading n by (k + 1) part must 
	     contain the upper triangular band part of the matrix of 
	     coefficients, supplied row by row, with the leading diagonal
	     of the matrix in column (k + 1) of the array, the first super-
	     diagonal starting at position 2 in column k, and so on.
	     The top left k by k triangle of the array A is not referenced.

   x         linear array of dimension of at least n.  Before entry,
	     x must contain the n elements of vector x.  On exit, x is
	     overwritten with the transformed vector x.


   Functions called
   --------------

   MISC      imax  

 ***********************************************************************/

void dtbmv(long trans, long n, long k, double **a, double *x)

{
   long info, j, i, l, end;
   double temp;

   info = 0;
   if      ( trans != TRANSP && trans != NTRANSP )   info = 1;
   else if ( n < 0 )                                 info = 2;
   else if ( k < 0 )                                 info = 3;

   if (info) {
      fprintf(stderr, "%s %1ld %s\n",
      "*** ON ENTRY TO DTBMV, PARAMETER NUMBER",info,"HAD AN ILLEGAL VALUE");
      exit(info);
   }

   switch(trans) {
      case NTRANSP:  for (j = 0; j < n; j++) {
                        temp = x[j];
                        l = k - j;
                        for (i = imax(0, j - k); i < j; i++) 
                           x[i] += temp * a[j][l+i];
                        x[j] *= a[j][k];
                     }
		     break;

      case TRANSP:   for (j = n - 1; j >= 0; j--) {
			temp = x[j] * a[j][k];
			l = k - j;
			end = imax(0, j - k);
			for (i = j - 1; i >= end; i--)
			   temp += x[i] * a[j][l+i];
                        x[j] = temp;
		     }
		     break;
   }
}
/************************************************************** 
 * Function interchanges two vectors		     	      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 

void dswap(long n,double *dx,long incx,double *dy,long incy)

{
   long i;
   double dtemp;

   if (n <= 0 || incx == 0 || incy == 0) return;
   if (incx == 1 && incy == 1) {
      for (i=0; i < n; i++) {
	 dtemp = *dy;
	 *dy++ = *dx;
	 *dx++ = dtemp;
      }	
   }
   else {
      if (incx < 0) dx += (-n+1) * incx;
      if (incy < 0) dy += (-n+1) * incy;
      for (i=0; i < n; i++) {
         dtemp = *dy;
         *dy = *dx;
         *dx = dtemp;
         dx += incx;
         dy += incy;
      }
   }
}
#include <math.h>

#define               TRUE      1
#define               FALSE     0
#define               MAXIT     30
#define               ONE       1.0
#define               ZERO      0.0
#define               CASE1     1
#define               CASE2     2
#define               CASE3     3
#define               CONVERGE  4

double dmax(double, double);
long   imin(long, long);
void   dscal(long, double, double *,long);
void   dswap(long, double *, long, double *, long);
void   drot(long, double *, double *, double, double);
void   drotg(double *, double *, double *, double *);

/***********************************************************************
 *                                                                     *
 *                        qriter2()                                    *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   This function reduces an upper bidiagonal matrix B to diagonal form.
   It is a C translation of a portion of DSVDC from Linpack.  In this 
   version, vectors are accumulated and B is assumed to be a square matrix.


   Parameters
   ---------

   (input)
   n           order of B
   s           linear array containing the diagonal elements of B
   e           linear array containing the off-diagonal elements of B

   (output)
   s           contains the singular values of the original matrix B
   up          2-dimensional array containing left singular vectors of B
   vp          2-dimensional array containing right singular vectors of B


   Functions called
   --------------

   BLAS         dscal, dswap, drot, drotg 
   MISC         dmax

 ***********************************************************************/
void qriter2(long n, double *s, double *e, double **up, double **vp)

{
   long negligible, iter, m, mm1, k, l, qrcase;
   double ztest, test, sl, g, t, smm1;
   double f, t1, cs, sn, scale, sm, el, emm1, b, c, shift;

   m     = n - 1;
   iter  = 0;
   while (m >= 0) {
      if (iter >= MAXIT) return;
      negligible = FALSE;

      /* this portion of the code inspects for negligible elements
       * in the s and e arrays.  On completion the variable qrcase
       * is set as follows:
       * qrcase = CASE1 if s[m] and e[l-1] are negligible and l < m
       * qrcase = CASE2 if s[l] is negligible and l < m
       * qrcase = CASE3 if s[l], ..., s[m] are not negligible (QR step),
       *                e[l-1] is negligible and l < m
       * qrcase = CONVERGENCE if e[m-1] is negligible */

      for (l = m - 1; l >= 0; l--) {
	 test = fabs(s[l]) + fabs(s[l + 1]);
	 ztest = test + fabs(e[l]);
	 if (ztest == test) {
	    e[l] = ZERO;
	    negligible = TRUE;
	 }
	 if (negligible) break;
      }

      if (l == m - 1) qrcase = CONVERGE;
      else {
         negligible = FALSE;
         for (k = m; k > l; k--) {
	    test = ZERO;
	    if (k != m) test += fabs(e[k]);
	    if (k != l + 1) test += fabs(e[k-1]);
	    ztest = test + fabs(s[k]);
	    if (ztest == test) {
	       s[k] = ZERO;
	       negligible = TRUE;
	    }
	    if (negligible) break;
	 }
	 if (k == l) qrcase = CASE3;
	 else if (k == m) qrcase = CASE1;
	 else {
	    qrcase = CASE2;
	    l = k;
	 }
      }
      l += 1;

      switch(qrcase) {

         /* deflate negligible s[m] */
	 case CASE1:    mm1 = m - 1;
			f = e[mm1];
			e[mm1] = ZERO;
			for (k = mm1; k >= l; k--) {
			   t1 = s[k];
			   drotg(&t1, &f, &cs, &sn);
			   s[k] = t1;
			   if (k != l) {
			      f = -sn * e[k - 1];
			      e[k - 1] *= cs;
			   }
			   drot(n, vp[k], vp[m], cs, sn);
                        }
			break;

         /* split at negligible s[l] */
	 case CASE2:    f = e[l - 1];
			e[l - 1] = ZERO;
			for (k = l; k <= m; k++) {
			   t1 = s[k];
			   drotg(&t1, &f, &cs, &sn);
			   s[k] = t1;
			   f = -sn * e[k];
			   e[k] *= cs;
			   drot(n, up[k], up[l - 1], cs, sn);
                        }
			break;

         /* perform one QR step */
	 case CASE3:    f = e[l - 1];
			/* calculate the shift */
			scale = dmax(fabs(s[m]), fabs(s[m - 1]));
			if (scale < fabs(e[m - 1])) scale = fabs(e[m - 1]);
			if (scale < fabs(s[l])) scale = fabs(s[l]);
			if (scale < fabs(e[l])) scale = fabs(e[l]);
			sm = s[m] / scale;
			smm1 = s[m - 1] / scale;
			emm1 = e[m - 1] / scale;
			sl = s[l] / scale;
			el = e[l] / scale;
			b = ((smm1 + sm) * (smm1 - sm) + emm1 * emm1) / 2.0;
			c = (sm * emm1);
			c *= c;
			shift = ZERO;
			if (b != ZERO || c !=ZERO) {
			   shift = sqrt(b * b + c);
			   if (b < ZERO) shift = -shift;
			   shift = c / (b + shift);
			}
			f = (sl + sm) * (sl - sm) + shift;
			g = sl * el;

			/* chase zeros */
			mm1 = m - 1;
			for (k = l; k <= mm1; k++) {
			   drotg(&f, &g, &cs, &sn);
			   if (k != l) e[k - 1] = f;
			   f = cs * s[k] + sn * e[k];
			   e[k] = cs * e[k] - sn * s[k];
			   g = sn * s[k + 1];
			   s[k + 1] = cs * s[k + 1];
			   drot(n, vp[k], vp[k + 1], cs, sn);
			   drotg(&f, &g, &cs, &sn);
			   s[k] = f;
			   f = cs * e[k] + sn * s[k + 1];
			   s[k + 1] = -sn * e[k] + cs * s[k + 1];
			   g = sn * e[k + 1];
			   e[k + 1] = cs * e[k + 1];
			   if (k < n - 1)
			      drot(n, up[k], up[k + 1], cs, sn);
			}
			e[mm1] = f;
			iter += 1;
			break;

         /* convergence */
	 case CONVERGE: if (s[l] < ZERO) {
			   /* make singular value positive */
	                   s[l] = -s[l];
			   dscal (n, -ONE, vp[l], 1); 
                        }
			/* order the singular value */
			while (l < n - 1) {
			   if (s[l] < s[l + 1]) {
			      t = s[l];
			      s[l] = s[l + 1];
			      s[l + 1] = t;
			      if (l < n - 1) 
				 dswap(n, vp[l], 1, vp[l + 1], 1);
			      if (l < n - 1) 
				 dswap(n, up[l], 1, up[l + 1], 1);
                              l += 1;
			   }
			   else break;
			}
			iter = 0;
			m -= 1;
			break;
      }
   }
}
/***************************************************************** 
 * applies a plane rotation;  assume a stride of 1 for dx and dy *
 * based on FORTRAN 77 routine from Linpack by J. Dongarra       *
 *****************************************************************/ 

void drot(long n, double *dx, double *dy, double c, double s)

{
   long i;
   double temp;

   if (n <= 0) return;

   for (i = 0; i < n; i++) {
      temp = c * (*dx) + s * (*dy);
      *dy = c * (*dy) - s * (*dx);
      dy++;
      *dx++ = temp;
   }
   return;
}
#include <math.h>

#define       ZERO       0.0
#define       ONE        1.0

double fsign(double, double);

/***************************************************************** 
 * constructs Givens plane rotation                              *
 * based on FORTRAN 77 routine from Linpack by J. Dongarra       *
 *****************************************************************/ 

void drotg(double *da, double *db, double *c, double *s)

{
   double r, roe, scale, z, temp1, temp2;

   roe = *db;
   temp1 = fabs(*da);
   temp2 = fabs(*db);
   if (temp1 > temp2) roe = *da;
   scale = temp1 + temp2;

   if (scale != ZERO) {
      temp1 = *da / scale;
      temp2 = *db / scale;
      r = scale * sqrt(temp1 * temp1 + temp2 * temp2);
      r *= fsign(ONE, roe);
      *c = *da / r;
      *s = *db / r;
   }
   else {
      *c = ONE;
      *s = ZERO;
      r = ZERO;
   }
   z = *s;

   temp1 = fabs(*c);
   if (temp1 > ZERO && temp1 <= *s) z = ONE / *c;

   *da = r;
   *db = z;
   return;
}
/************************************************************** 
 * Function forms the dot product of two vectors.      	      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 

double ddot(long n,double *dx,long incx,double *dy,long incy)

{
   long i;
   double dot_product;

   if (n <= 0 || incx == 0 || incy == 0) return(0.0);
   dot_product = 0.0;
   if (incx == 1 && incy == 1) 
      for (i=0; i < n; i++) dot_product += (*dx++) * (*dy++);
   else {
      if (incx < 0) dx += (-n+1) * incx;
      if (incy < 0) dy += (-n+1) * incy;
      for (i=0; i < n; i++) {
         dot_product += (*dx) * (*dy);
         dx += incx;
         dy += incy;
      }
   }
   return(dot_product);
}
/************************************************************** 
 * Function scales a vector by a constant.     		      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 

void dscal(long n,double da,double *dx,long incx)

{
   long i;

   if (n <= 0 || incx == 0) return;
   if (incx < 0) dx += (-n+1) * incx;
   for (i=0; i < n; i++) {
      *dx *= da;
      dx += incx;
   }
   return;
}
/************************************************************** 
 * Constant times a vector plus a vector     		      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 

void daxpy (long n,double da,double *dx,long incx,double *dy,long incy)

{
   long i;

   if (n <= 0 || incx == 0 || incy == 0 || da == 0.0) return;
   if (incx == 1 && incy == 1) 
      for (i=0; i < n; i++) {
	 *dy += da * (*dx++);
	 dy++;
      }
   else {
      if (incx < 0) dx += (-n+1) * incx;
      if (incy < 0) dy += (-n+1) * incy;
      for (i=0; i < n; i++) {
         *dy += da * (*dx);
         dx += incx;
         dy += incy;
      }
   }
   return;
}
#include <math.h>
/***********************************************************************
 *                                                                     *
 *				random()                               *
 *                        (double precision)                           *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   This is a translation of a Fortran-77 uniform random number
   generator.  The code is based  on  theory and suggestions  given in
   D. E. Knuth (1969),  vol  2.  The argument to the function should 
   be initialized to an arbitrary integer prior to the first call to 
   random.  The calling program should  not  alter  the value of the
   argument between subsequent calls to random.  Random returns values
   within the the interval (0,1).


   Arguments 
   ---------

   (input)
   iy	   an integer seed whose value must not be altered by the caller
	   between subsequent calls

   (output)
   random  a double precision random number between (0,1)

 ***********************************************************************/
double random(long *iy)

{
   static long m2 = 0;
   static long ia, ic, mic;
   static double halfm, s;

   /* If first entry, compute (max int) / 2 */
   if (!m2) {
      m2 = 1 << (8 * (int)sizeof(int) - 2); 
      halfm = m2;

      /* compute multiplier and increment for linear congruential 
       * method */
      ia = 8 * (long)(halfm * atan(1.0) / 8.0) + 5;
      ic = 2 * (long)(halfm * (0.5 - sqrt(3.0)/6.0)) + 1;
      mic = (m2-ic) + m2;

      /* s is the scale factor for converting to floating point */
      s = 0.5 / halfm;
   }

   /* compute next random number */
   *iy = *iy * ia;

   /* for computers which do not allow integer overflow on addition */
   if (*iy > mic) *iy = (*iy - m2) - m2;

   *iy = *iy + ic;

   /* for computers whose word length for addition is greater than
    * for multiplication */
   if (*iy / 2 > m2) *iy = (*iy - m2) - m2;
  
   /* for computers whose integer overflow affects the sign bit */
   if (*iy < 0) *iy = (*iy + m2) + m2;

   return((double)(*iy) * s);
}
/************************************************************** 
 * Function copies a vector x to a vector y	     	      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 

void dcopy(long n,double *dx,long incx,double *dy,long incy)

{
   long i;

   if (n <= 0 || incx == 0 || incy == 0) return;
   if (incx == 1 && incy == 1) 
      for (i=0; i < n; i++) *dy++ = *dx++;

   else {
      if (incx < 0) dx += (-n+1) * incx;
      if (incy < 0) dy += (-n+1) * incy;
      for (i=0; i < n; i++) {
         *dy = *dx;
         dx += incx;
         dy += incy;
      }
   }
   return;
}
