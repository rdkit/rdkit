/*************************************************************************
                           (c) Copyright 1993
                        University of Tennessee
                          All Rights Reserved                          
 *************************************************************************/
#include <stdio.h>
#include <errno.h>
#include <math.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include "bls2.h"

#define  ONE   1.0

#define  UNIX_CREAT

#ifdef UNIX_CREAT
#define PERMS 0664
#endif


double ddot(long, double *,long, double *, long);
void   opb(long, long, double *, double *);
void   opa(long, long, double *, double *);
void   dscal(long, double, double *,long);
void   daxpy(long, double, double *,long, double *, long);
float  timer();
long   imin(long, long);
long   validate(FILE *, long, long, long, long, long, long, double);
long   blklan2(FILE *, long, long, long, long, double *, double *,
	       long, long, double, double *, long, long *, long *,
	       long *, long*);

/***********************************************************************
 *                                                                     *
 *                        main()                                       *
 *     Sparse SVD Via Hybrid Block Lanczos Procedure for Equivalent    *
 *               A'A Eigensystems  (double precision)                  *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   This is a sample driver that invokes BLKLAN2 to compute singular 
   triplets of a large sparse matrix A.  In this test program, bls2,
   the Harwell-Boeing sparse matrix format is employed for accessing
   elements of the m by n sparse matrix A and its transpose.  Other 
   sparse matrix formats can be used, of course.  Approximate singular
   values of A and their residuals can be found in formatted output
   file blo2 and corresponding right and left singular vectors can be
   found in unformatted output file blv2.

   User supplied routines:  opa, opm, opb, timer
   
   opa(m, n, x, y) takes an n-vector x and returns the m-vector y = A * x.

   opb(m, n, x, y) takes an n-vector x and returns the m-vector y = A'A * x.  

   opm(m, n, nc, x, y) takes a block of nc vectors length n stored in two-
   dimensional array x and returns the block of n-vectors y = A'A * x. 

   User should edit timer() with an appropriate call to an longrinsic
   timing routine that returns elapsed user cpu time.


   External parameters
   -------------------

   Defined and documented in bls2.h


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
		If vtf is set to "TRUE", the unformatted output file blv2
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

   BLAS         ddot, daxpy, dscal
   USER         opa, opb, timer
   MISC         validate
   BLS2         blklan2


   BLS2 development
   ----------------

   BLS2 is a C translation of the Fortran-77 BLS2 from the SVDPACK
   library written by Michael W. Berry, University of Tennessee,
   Dept. of Computer Science, 107 Ayres Hall, Knoxville, TN, 37996-1301.
   This particular implementation is discussed in "Multiprocessor Sparse
   SVD Algorithms and Applications", Ph.D. thesis by M. Berry, University
   of Illinois at Urbana-Champaign, October 1990.

   Date written:  05 May 1992

   Theresa H. Do
   University of Tennessee
   Dept. of Computer Science
   107 Ayres Hall
   Knoxville, TN, 37996-1301
   longernet: tdo@cs.utk.edu

 ***********************************************************************/
  main() 
{
   double t0, exetime;
   long i, j, k, nrow, ncol, nnzero, nc, nb, nums, maxit, nn;
   long ncold, nbold, nk, size1, size2, indexu, indexv;
   long memory, vectors;
   double tol, *sing, *v, *u, *res, *tmpv, *tptr1;
   double tmp1, tmp2, xnorm;
   char title[73], name[41], vtf[6];
   char *in1, *in2, *out1, *out2;
   FILE *fp_in1, *fp_in2;
   FILE *fp_out1 = NULL; 
   long  fp_out2;
 
   in1 = "blp2";
   in2 = "matrix";
   out1 = "blo2";
   out2 = "blv2";

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
    title[73] = '\0';
    nn=imin(nrow,ncol);
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
       if ((fp_out2 = creat(out2, PERMS)) == -1) {
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
     * allocate memory						       *
     * pointr - column start array of harwell-boeing sparse matrix     *
     * 		format	       			          (ncol + 1)   *
     * rowind - row indices array of harwell-boeing format  (nnzero)   *
     * value  - nonzero values array of harwell-boeing sparse matrix   *
     *		format	       				    (nnzero)   *
     * sing   -                                               (nums)   *
     * res    -                                               (nums)   *
     * v      -                                        (ncol * nums)   *
     * ztemp  -                     (nrow * nc; assume nrow >= ncol)   *
     * (ztempp)                                                 (nc)   *
     *******************************************************************/
     size1 = sizeof(double) * (nnzero + nums * (ncol+2) + nrow * nc);
     size2 = sizeof(long) * (ncol + nnzero + 1);
     if (!(rowind = (long *)   malloc(size2))   ||
         !(value  = (double *) malloc(size1))   ||
	 !(ztempp = (double **) malloc(sizeof(double *) * nc)))    {
	 perror("MALLOC FAILED in MAIN()");
	 exit(errno);
     }
     tptr1 = value;

     /* calculate memory allocated for problem */
     memory = size1 + size2 + sizeof(double *) * nc;

     pointr = rowind + nnzero;
     tptr1 += nnzero;
     v      = tptr1;
     tptr1 += ncol * nums;
     sing   = tptr1;
     tptr1 += nums;
     res    = tptr1;
     tptr1 += nums;
     ztemp  = tptr1;
     j = 0;
     k = nc * nrow;
     for (i = 0; i < k; i += nrow) ztempp[j++] = &ztemp[i];
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
    if (blklan2(fp_out1, nnzero, nrow, ncol, nums, v, sing, ncold, 
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
    /*******************************************************************
     * allocate memory						       *
     * tmpv   -                                               (ncol)   *
     * u      -                                          (nrow * nk)   *
     *******************************************************************/
       size1 = sizeof(double) * (ncol + nk * nrow);
       memory += size1;
       if (!(tmpv = (double *) malloc(size1))) {
	   perror("SECOND MALLOC FAILED in MAIN()");
	   exit(errno);
       }
       u = tmpv + ncol;
       size1 = sizeof(double) * nrow;
       size2 = sizeof(double) * ncol;
       indexu = 0;
       indexv = 0;

       t0 = timer();
       for (i = 0; i < nk; i++) {

	  /* multiply by A'A */
          opb(nrow, ncol, &v[indexv], tmpv);
	  tmp1 = ddot(ncol, &v[indexv], 1, tmpv, 1);
	  daxpy(ncol, -tmp1, &v[indexv], 1, tmpv, 1);
	  tmp1 = sqrt(tmp1);
	  xnorm = sqrt(ddot(ncol, tmpv, 1, tmpv, 1));

	  /* multiply by A to get scaled left singular vectors */
	  opa(nrow, ncol, &v[indexv], &u[indexu]);
	  tmp2 = ONE / tmp1;
          dscal(nrow, tmp2, &u[indexu], 1);
	  res[i] = xnorm * tmp2;
	  sing[i] = tmp1;
          write (fp_out2, (char *)&u[indexu], size1);
          write (fp_out2, (char *)&v[indexv], size2);
	  indexu += nrow;
	  indexv += ncol;
       }
       exetime += timer() - t0;
    }
    /* write results to output file */
    fprintf(fp_out1,"\n");
    fprintf(fp_out1, " ... \n");
    fprintf(fp_out1, " ... HYBRID BLOCK LANCZOS [A^TA]\n");
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
           fprintf(fp_out1," ...... %3ld   %22.14E  (%11.2E)\n",i+1,sing[i],res[
i]);


    free(value);
    free(rowind);
    fclose(fp_in1);
    fclose(fp_in2);
    fclose(fp_out1);
    if (vectors) {
       free(tmpv);
       close(fp_out2);
    }
    exit(0);
}
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <fcntl.h>

extern double *tres, *y, *vv, *v0, *uvtmp, *pp, *qq;
extern double *alpha, *beta, *p, *q, *t, *z, *ztemp;
extern double **yp, **vvp, **uvtmpp, **ppp, **qqp, **ztempp;
extern long iconv, nn, iter;

#define   SEED      91827211
#define   MINBLKS   2
#define   TRANSP    1
#define   NTRANSP   0
#define   TRUE      1
#define   FALSE     0
#define   ZERO      0.0
#define   ONE       1.0

long   tql2(long, double *, double *, double **);
void   dgemv(long, long, long, double, double **, double *, double, double *);
void   dgemm2(long, long, long, long, long, double, double **, double **, 
              double, double **);
long   imin(long, long);
double random(long *);
long   validate(FILE *, long, long, long, long, long, long, double);
long   polong2(long *, long, long, double *, double *, double *, double);
void   block2(double **, double **, double **, long, long, long, long, long *);
void   orthg(long, long, long, double **, double **, double *);

/***********************************************************************
 *                                                                     *
 *                        blklan2()                                    *
 *                  Block Lanczos SVD Alogrithm                        *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   blklan2.c is designed to compute singular values and singular vectors
   of a large sparse matrix A.  The singular values (and singular vectors)
   of A are determined by the eigenvalues (and eigenvectors) of the matrix
   B, where B = A'A.

   The eigenvalues of B are the squares of the singular values of A, the
   eigenvectors of B correspond to the right singular vectors of A only.
   The left singular vectors of A are determined by 

		u = 1/sigma A*v,

   where {u, sigma, v} is a singular triplet of A.  This version of blklan2 
   is designed to approximate the ik-largest singular triplets of A. 
   This hybrid block Lanczos procedure consists of five phases:

   Phase 1:  Block Lanczos outer iteration to yield a symmetric block
             tridiagonal matrix S whose eigenvalues approximate those
	     of the matrix B = A'A.  Total (or complete)
	     re-orthogonalization is used here.

   Phase 2:  Lanczos method for tridiagonalizing the S matrix from
	     Phase 1 to yield the tridiagonal matrix T whose eigenvalues
	     also approximate those of B.  Total (or complete)
	     re-orthogonalization is used for this Lanczos recursion.  
	     A polong Lanczos method (single vector) is used if a blocksize
	     of 1 is encountered via normal deflation.

   Phase 3:  Apply an appropriate QL iteration to diagonalize T and
	     hence produce approximate eigenvalues (array alpha)
	     of matrix B = A'A and hence the squares of the singular
	     values of the original sparse matrix A.

   Phase 4:  Convergence test using a user-supplied residual tolerance.

   Phase 5:  Iteration restart with orthogonalization with respect
	     to any (all) converged eigenvectors of B = A'A.



   Arguments
   ---------

   (input)
   fp       polonger to output file
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
   eig     linear array containing the iko approximate singular values
   v       two-dimensional array containing the iko approximate right 
	      singular vectors corresponding to the approximate singular 
	      values in array eig
   res     linear array containing the iko residuals of the approximate 
	   singular triplets
   memory  memory storage needed in bytes


   External parameters
   -------------------

   Defined and documented in bls2.h


   Local parameters
   ----------------

    k        current # of desired triplets (exit when reduced to 0)
    k0       count of triplets found in current iteration
    nb       current block size
    nc       size of current subspace
    ns       number of blocks in current iteration
   

   Functions called
   --------------

   BLAS         dgemv, dgemm2, tql2, orthg
   MISC         validate, random, imin
   BLS2         polong2, block2


   NOTE:  Unlike Fortran, C stores arrays in row-major order.  For the
	  sake of efficiency, most matrices and vectors in the algorithm
	  are transposed, i.e., each vector and each column of a
	  matrix are stored row-wise.

 ***********************************************************************/

long blklan2(FILE *fp, long nnzero, long m, long n, long ik, double *v,
	     double *eig, long ic, long ib, double tol, double *res, 
	     long maxit, long *iko, long *ico, long *ibo, long *memory)

{
   long irand;
   long nbold, cont, ns, nb, nc, k, k0, length, memsize;
   long jj, ii, kk, ll, final, flag, i, j, *index;
   double **sp, **rp, **bigsp, **workptr2, **tempptr2;
   double *s, *r, *bigs, *workptr1, *tempptr1;
   double *ptr1, *ptr2;

   nb = ib;
   nc = ic;
   k  = ik;
   k0 = 0;
   ns = nc / nb;

   /* reset converged vector counter */
   iconv = 0;
   irand = (long) SEED;

   /* determine number of blocks for first iteration */
   if (ns < MINBLKS) {
      nb = nc / MINBLKS;
      ns = nc / nb;
   }
   /* return upon error on input */
   if (validate(fp, nnzero, m, n, *ico, *ibo, ik, tol))
      return(-1);

   /* allocate memory for linear arrays of doubles s, r, bigs, y, vv
    * v0, tres, uvtmp, z, p, q, t, alpha, beta, pp, qq, v */

   memsize = sizeof(double) * 
	     (nb * (nc + nc - nb + n + 1)     +
              nc * (nb + 2 * n + 3 + nc + nc) +
	      ik * (n + n)                    + 
	      4 * n);                   

   *memory += memsize;
 /***********************************************************************
  *         Allocate work area and initialize polongers                  *
  *   polonger              size                                         *
  * s     (sp)          nc      by nb                                   *
  * r     (rp)          (nc-nb) by nb                                   *
  * bigs  (bigsp)       nc      by (nb+1)                               *
  * y     (yp)          nb      by n                                    *
  * vv    (vvp)         nc      by n                                    *
  * v0                  ik      by n                                    *
  * tres                nb                                              *
  * uvtmp (uvtmpp)      (nc+ik) by n                                    *
  * z                   n                                               *
  * p                   n                                               *
  * q                   n                                               *
  * t                   n                                               *
  * alpha               nc                                              *
  * beta                nc                                              *
  * pp    (ppp)         nc      by nc                                   *
  * qq    (qqp)         nc      by nc                                   * 
  * v                   ik      by n  (allocated in bls2.c)             *
  * ztemp (ztempp)      nc      by n  (allocated in bls2.c)             *
  * index               nb                                              *
  ***********************************************************************/

   if (!(workptr1 = (double *)malloc(memsize)) ||
       !(index   = (long   *)malloc(sizeof(long) * nb))){
      perror("FIRST MALLOC FAILED in BLKLAN2()");
      exit(errno);
   }

   /* memory for linear array index */

   *memory += sizeof(long) * nb;

   /* allocate memory for arrays of polongers sp, rp, bigsp, yp, vvp,
    * uvtmpp, ppp, qqp.  This will allow the memory areas s, r, bigs, y,
    * vv, uvtmp, pp and qq to be addressed as 2-dimensional arrays */

   memsize = sizeof(double *) * (7 * nc + nb + ik);
   *memory += memsize;
   if (!(workptr2 = (double **)malloc(memsize))){
      perror("SECOND MALLOC FAILED in BLKLAN2()");
      exit(errno);
   }
   tempptr1  = workptr1;
   tempptr2  = workptr2;

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

   length    = n * nb;
   y         = tempptr1;
   tempptr1 += length;
   yp        = tempptr2;
   tempptr2 += nb;
   j = 0;
   for (i = 0; i < length; i += n) yp[j++] = &y[i];

   length    = n * nc;
   vv        = tempptr1;
   tempptr1 += length;
   vvp       = tempptr2;
   tempptr2 += nc;
   j = 0;
   for (i = 0; i < length; i += n) vvp[j++] = &vv[i];

   v0        = tempptr1;
   tempptr1 += n * ik;

   tres      = tempptr1;
   tempptr1 += nb;

   length    = n * (nc + ik);
   uvtmp     = tempptr1;
   tempptr1 += length;
   uvtmpp    = tempptr2;
   tempptr2 += nc + ik;
   j = 0;
   for (i = 0; i < length; i += n) uvtmpp[j++] = &uvtmp[i];

   z         = tempptr1;
   tempptr1 += n;

   p         = tempptr1;
   tempptr1 += n;

   q         = tempptr1;
   tempptr1 += n;

   t         = tempptr1;
   tempptr1 += n;

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
   length = n * nb;
   for (i = 0; i < length; i++) vv[i] = random(&irand);

   orthg(nb, 0, n, yp, vvp, ztemp); 
   iter = 0;
   while (iter < maxit) {
      nn = nb * ns;
      cont = TRUE;
      iter += 1;

      /*------------------------------------------------------------------*
       *           PHASE 1 and PHASE 2 (block algorithm)                  *
       *------------------------------------------------------------------*/

      if (nb > 1) 
	 block2(sp, rp, bigsp, m, n, nb, ns, &irand);
      else {
	 if (nbold != 1) k0 = 0;
	 /*------------------------------------------------------------------*
	  *                   PHASE 2A (polong algorithm)                     *
	  *------------------------------------------------------------------*/

	 cont = polong2(&k, m, n, v, res, eig, tol);
      }
      if (!cont) break;

      /*------------------------------------------------------------------*
       *                        PHASE 3                                   *
       *------------------------------------------------------------------*/

      /* solve symmetric tridiagonal EVP */
      for (i = 0; i < nn; i++) 
         for (j = 0; j < nn; j++) qqp[i][j] = ZERO;
      for (i = 0; i < nn; i++) qqp[i][i] = ONE;

      /* resort alpha's and rows of pp (need descending order) */
      if (tql2(nn, alpha, beta, qqp)) break;

      i = nn - 1;
      for (jj = 0; jj < nn; jj++) {
	 z[jj] = alpha[i];
	 for (ii = 0; ii < nn; ii++) uvtmpp[jj][ii] = qqp[i][ii];
	 i--;
      }
      for (jj = 0; jj < nn; jj++) {
	 alpha[jj] = z[jj];
	 for (ii = 0; ii < nn; ii++) qqp[jj][ii] = uvtmpp[jj][ii];
      }

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
	      eig[iconv + k0] = alpha[i];
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
         dgemm2(NTRANSP, NTRANSP, nn, nn, nn, ONE, qqp, ppp, ZERO, ztempp);
         dgemm2(NTRANSP, NTRANSP, nb + k0, n, nn, ONE, ztempp, vvp, ZERO, uvtmpp);
	 ptr1 = uvtmp;
	 length = (nb + k0) * n;
	 for (i = 0; i < length; i++) v[i] = *ptr1++;
      }
      else dgemv(TRANSP, nn, n, ONE, vvp, qq, ZERO, v);
      if (k0) {
         final = iconv + k0;
         ii = 0;
         for (jj = iconv; jj < final; jj++) {
	    ptr1 =  &v[n * index[ii++]];
	    ptr2 = &v0[n * jj];
            for (i = 0; i < n; i++) *ptr2++ = *ptr1++;
         }
         iconv = final;
         if (k <= 0) {
	    iter -= 1;
	    cont = FALSE;
         }
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
   ptr1 =  v;
   ptr2 = v0;
   for (j = 0; j < iconv; j++) 
      for (i = 0; i < n; i++) *ptr1++ = *ptr2++;
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

extern double *tres, *ztemp, *vv, *v0, *y, *uvtmp, *pp, *qq;
extern double **vvp, **uvtmpp, **yp, **ppp, **qqp;
extern double *p, *q, *t, *z, *alpha, *beta;
extern long iconv, nn, iter;

void   dsbmv(long, long, double, double **, double *, double, double *);
void   dgemm2(long, long, long, long, long, double, double **, double **,
            double, double **);
void   orthg(long, long, long, double **, double **, double *);
void   formbigs(long, long, double **, double **, double **);
double ddot(long, double *,long, double *, long);
void   daxpy(long, double, double *,long, double *, long);
double random(long *);
double enorm(long, double *);
void   opm(long, long, long, double **, double **);

/***********************************************************************
 *                                                                     *
 *                        block2()                                     *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   This function implements the first two phases of the hybrid block
   Lanczos procedure.  In the first phase, which is also known as the
   block Lanczos outer iteration, a symmetric block tridiagonal matrix S
   is formed.  The eigenvalues of the matrix S approximate those 
   of matrix B, where B = A'A and A is the original sparse matrix.  
   Total (or complete) re-orthogonalization is used.

   In the second phase, single-vector Lanczos tridiagonalization is used
   to reduce (preserving eigenvalues of S) the block matrix S to a 
   symmetric tridiagonal matrix T. 


   Arguments
   ---------

   (input)
   w, wp     work space
   sp        diagonal blocks (symmetric submatrices) of the
		symmetric block tridiagonal matrix S 
   rp        super-diagonal blocks (upper-triangular submatrices)
		of the symmetric block tridiagonal matrix S 
   bigsp     symmetric block tridiagonal matrix S
   m         row dimension of sparse matrix A
   n         column dimension of sparse matrix A
   nb        current block size
   ns        number of blocks in current iteration
   irand     seed for random number generator
   

   (output - globally defined)
   alpha     diagonal elements of symmetric tridiagonal matrix T (reduced 
                from matrix S)
   beta      off-diagonal elements of symmetric tridiagonal matrix T 
		(reduced from matrix S)
   tres      residuals of approximate eigenvalues determined from
	     a previous set of block Lanczos outer iterations.
   ppp       array of eigenvectors of S


   External parameters
   -------------------

   Defined and documented in bls2.h


   Functions called
   --------------

   BLAS         ddot, daxpy, enorm, dgemm2, orthg, dsbmv
   USER         opm
   MISC         random
   BLS2         formbigs

 ***********************************************************************/

void block2(double **sp, double **rp, double **bigsp, long m, long n, 
            long nb, long ns, long *irand)

{
   long jinc, nk, nj, i, j, k, blkptr;
   double *ptr, dum, pnm;

   for (i = 0; i < nn; i++)
      for (j = 0; j < nb; j++) sp[i][j] = ZERO;
   /* ns (number of blocks) is assumed to be at least 2 */
   nk = nn - nb;
   for (i = 0; i < nk; i++)
      for (j = 0; j < nb; j++) rp[i][j] = ZERO;

   opm(m, n, nb, vvp, yp);

   dgemm2(NTRANSP, TRANSP, nb, nb, n, ONE, vvp, yp, ZERO, sp);

   blkptr = 0;
   for (j = 1; j < ns; j++) {

      dgemm2(TRANSP, NTRANSP, nb, n, nb, -ONE, &sp[blkptr],
            &vvp[blkptr], ONE, yp);
      if (j > 1)
         dgemm2(NTRANSP, NTRANSP, nb, n, nb, -ONE, &rp[blkptr - nb],
               &vvp[blkptr - nb], ONE, yp);

      if (j == 1 && iter > 1)
         for (i = 0; i < nb; i++) tres[i] = enorm(n, yp[i]);

      if (iconv) {
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
      orthg(nb, nj, n, yp, uvtmpp, ztemp);

      for (k = 0; k < nb; k++) 
         for (i = k; i < nb; i++) rp[blkptr + k][i] = yp[i][k];

      jinc = blkptr + nb;
      ptr = vvp[jinc];
      for (k = nj; k < nj + nb; k++) 
         for (i = 0; i < n; i++) *ptr++ = uvtmpp[k][i]; 
  
      opm(m, n, nb, &vvp[jinc], yp);
      dgemm2(NTRANSP, TRANSP, nb, nb, n, ONE, &vvp[jinc], yp, ZERO, &sp[jinc]);
      blkptr += nb;
   }
   formbigs(nn, nb, rp, sp, bigsp);

   for (i = 0; i < nn; i++) p[i] = random(irand);
   pnm = enorm(nn,p);

   ptr = pp;
   for (i = 0; i < nn; i++) {
      p[i] /= pnm;
      *ptr++ = p[i];
      z[i] = ZERO;
   }
   for (j = 0; j < nn; j++) {
      dsbmv(nn, nb, ONE, bigsp, p, ONE, z);
      alpha[j] = ddot(nn, p, 1, z, 1);
      if (j == nn - 1) break;

      /* compute Z[j] := Z[j] - alpha[j] * P */
      daxpy(nn, -alpha[j], p, 1, z, 1);

      /* orthogonalize Z w.r.t. previous PP's */
      for (k = 0; k <= j; k++)
         for (i = 0; i < nn; i++) uvtmpp[k][i] = ppp[k][i];
      if (j) {
         ptr = uvtmpp[j + 1];
         for (i = 0; i < nn; i++) *ptr++ = z[i];
	 orthg(1, j+1, nn, yp, uvtmpp, ztemp);
	 ptr = uvtmpp[j + 1];
         for (i = 0; i < nn; i++) z[i] = *ptr++;
      }
      else {
         dum = -ddot(nn, uvtmp, 1, z, 1);
	 daxpy(nn, dum, uvtmp, 1, z, 1);
      }
      beta[j] = enorm(nn,z);
      if (beta[j] != ZERO) {

	 /* compute P[j+1] := Z[j] / beta[j] */
         ptr = ppp[j + 1];
         for (i = 0; i < nn; i++) {
	    t[i] = p[i];
	    p[i] = z[i] /beta[j];
	    *ptr++ = p[i];
	    z[i] = -beta[j] * t[i];
	 }
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

extern double *v0, *ztemp, *uvtmp, *vv, *alpha, *beta;
extern double *p, *q, *t, *z;
extern double **uvtmpp, **vvp, **yp;
extern long iconv, nn, iter;

double ddot(long, double *,long, double *, long);
double enorm(long, double *);
void   daxpy(long, double, double *,long, double *, long);
void   opb(long, long, double *, double *);
void   orthg(long, long, long, double **, double **, double *);

/***********************************************************************
 *                                                                     *
 *                        polong2()                                     *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   This function is a single-vector Lanczos tridiagonalization procedure
   (degenerate case of block size = 1 for block2.c) which is used when
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
   alpha   diagonal elements of symmetric tridiagonal matrix from the 
	   inner recursion
   beta    off-diagonal elements of symmetric tridiagonal matrix from the 
	   inner recursion

   External parameters
   -------------------

   Defined and documented in bls2.h


   Functions called
   --------------

   BLAS         ddot, enorm, daxpy, orthg
   USER         opb

 ***********************************************************************/

long polong2(long *k, long m, long n, 
            double *v, double *res, double *eig, double tol)

{
   double *ptr, dum, znm;
   long i, j, jj, convpj;

   for (j = 0; j < n; j++) p[j] = vv[j];
   for (j = 0; j < n; j++) z[j] = ZERO;

   for (j = 0; j < nn; j++) {
      opb(m, n, p, q);
      daxpy(n, ONE, q, 1, z, 1);
      alpha[j] = ddot(n, p, 1, z, 1);
      if (j == nn - 1) return(CONTINUE);

      /* compute Z[j] := Z[j] - alpha[j] * P */
      if (alpha[j] != ZERO) {
	 daxpy(n, -alpha[j], p, 1, z, 1);
	 if (!j) {
	    if ((znm = enorm(n, z)) <= tol) {
	       ptr = &v[iconv * n];
               for (i = 0; i < n; i++) *ptr++ = p[i];
	       eig[iconv] = alpha[0];
	       res[iconv] = znm;
               *k -= 1;
	       iter -= 1;
	       return(DONE);
            }
	 }

	 /* orthogonalize Z w.r.t. converged right S-vectors and previous 
	    VV's */
         convpj = iconv + j;
	 ptr = v0;
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
	    orthg(1, convpj + 1, n, yp, uvtmpp, ztemp);
	    ptr = uvtmpp[convpj + 1];
	    for (i = 0; i < n; i++) z[i] = *ptr++;
	 }
	 else {
	    dum = -ddot(n, uvtmp, 1, z, 1);
	    daxpy(n, dum, uvtmp, 1, z, 1);
	 }

	 /* compute beta[j] */
	 beta[j] = enorm(n,z);
	 if (beta[j] != ZERO) {

            /* compute P[j+1] := Z[j] / beta[j] */
	    ptr = vvp[j + 1];
	    for (i = 0; i < n; i++) {
	       t[i] = p[i];
	       p[i] = z[i] / beta[j];
	       *ptr++ = p[i];
	       z[i] = -beta[j] * t[i];
	    }
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
  fp_out1   polonger to output file
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
extern long mxvcount, mtxvcount;
extern long *pointr, *rowind;
extern double *value,*ztemp;

#define         ZERO         0.0

/**************************************************************
 *                                                            *
 * multiplication of matrix B by vector x, where B = A'A,     *
 * and A is nrow by ncol (nrow >> ncol) and is stored using   *
 * the Harwell-Boeing compressed column sparse matrix format. *
 * Hence, B is of order n = ncol.  y stores product vector.   *
 *                                                            *
 **************************************************************/

void opb(long nrow, long ncol, double *x, double *y)
{
   long i, j, end;
   
   mxvcount += 1;
   mtxvcount += 1;
   for (i = 0; i < ncol; i++) y[i] = ZERO;
   for (i = 0; i < nrow; i++) ztemp[i] = ZERO;

   for (i = 0; i < ncol; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++) 
	 ztemp[rowind[j]] += value[j] * (*x); 
      x++;
   }
   for (i = 0; i < ncol; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++) 
	 *y += value[j] * ztemp[rowind[j]];
      y++;
   }
   return;
}
extern long mxvcount, mtxvcount;
extern long *pointr, *rowind;
extern double *value, **ztempp;

#define       ZERO    0.0

/**************************************************************
 *                                                            *
 * multiplication of A'A by the transpose of an nc by n       *
 * matrix X.  A is nrow by ncol and is stored using the       *
 * Harwell-Boeing compressed column sparse matrix format.     *
 * The transpose of the n by nc product matrix Y is returned  *
 * in the two-dimensional array y.                            *
 *                                                            *
 *              Y = A'A * X'  and  y := Y'                    *
 *                                                            *
 **************************************************************/

void opm(long nrow, long ncol, long nc, double **x, double **y)
{
   long i, j, k, end;
   
   mxvcount  += nc;
   mtxvcount += nc;

   for (i = 0; i < nc; i++) 
      for (j = 0; j < nrow; j++) ztempp[i][j] = ZERO;
   for (i = 0; i < nc; i++) 
      for (j = 0; j < ncol; j++) y[i][j] = ZERO;

   /* multiply by sparse matrix A */
   for (i = 0; i < ncol; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++)
         for (k = 0; k < nc; k++)
	    ztempp[k][rowind[j]] += value[j] * x[k][i]; 
   }

   /* multiply by transpose of A (A') */
   for (k = 0; k < nc; k++)
      for (i = 0; i < ncol; i++) {
         end = pointr[i+1];
         for (j = pointr[i]; j < end; j++) 
	    y[k][i] += value[j] * ztempp[k][rowind[j]];
      }
   return;
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
   The resulting matrix Z is then factored longo the product of a
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
      fprintf(stderr, "%s %1d %s\n",
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
 *                         dgemm2()                                    *
 *                                                                     *
 * A C-translation of the level 3 BLAS routine DGEMM by Dongarra,      *
 * Duff, du Croz, and Hammarling (see LAPACK Users' Guide).            *
 * In this version, the arrays which store the matrices used in this   *
 * matrix-matrix multiplication are accessed as two-dimensional arrays.*
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   dgemm2() performs one of the matrix-matrix operations

	      C := alpha * op(A) * op(B) + beta * C,

   where op(X) = X or op(X) = X', alpha and beta are scalars, and A, B
   and C are matrices, with op(A) an m by k matrix, op(B) a k by n
   matrix and C an m by n matrix.


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
            leading m by k part of a must contain the matrix A. Otherwise,
	    the leading k by m part of a must contain  the matrix A.

   b        matrix B as a 2-dimensional array.  When transb = NTRANSP, the
            leading k by n part of a must contain the matrix B. Otherwise,
	    the leading n by k part of a must contain  the matrix B.

   beta     a scalar multiplier.  When beta is supplied as zero then C
	    need not be set on input.

   c        matrix C as a 2-dimensional array.  On entry, the leading
	    m by n part of c must contain the matrix C, except when
	    beta = 0.  In that case, c need not be set on entry. 
	    On exit, c is overwritten by the m by n matrix
	    (alpha * op(A) * op(B) + beta * C).

 ***********************************************************************/
void dgemm2(long transa, long transb, long m, long n, long k, 
            double alpha, double **a, double **b, double beta, double **c)
{
   long info;
   long i, j, l, nrowa, ncola, nrowb, ncolb;
   double temp, *atemp;

   info = 0;
   if      ( transa != TRANSP && transa != NTRANSP ) info = 1;
   else if ( transb != TRANSP && transb != NTRANSP ) info = 2;
   else if ( m < 0 ) 				     info = 3;
   else if ( n < 0 )				     info = 4;
   else if ( k < 0 )        			     info = 5;

   if (info) {
      fprintf(stderr, "%s %1d %s\n",
      "*** ON ENTRY TO DGEMM2, PARAMETER NUMBER",info,"HAD AN ILLEGAL VALUE");
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
   if (!m || !n || ((alpha == ZERO || !k) && beta == ONE))
      return;

   if (alpha == ZERO) {
      if (beta == ZERO) 
         for (i = 0; i < m; i++)
            for (j = 0; j < n; j++) c[i][j] = ZERO;

      else if (beta != ONE)
         for (i = 0; i < m; i++)
            for (j = 0; j < n; j++) c[i][j] *= beta;
      return;
   }

   if (beta == ZERO)
      for (i = 0; i < m; i++)
         for (j = 0; j < n; j++) c[i][j] = ZERO;

   else if (beta != ONE)
      for (i = 0; i < m; i++)
         for (j = 0; j < n; j++) c[i][j] *= beta;

   if (!transb) { 

      switch(transa) {

	 /* form C := alpha * A * B + beta * C */
	 case NTRANSP:  for(l = 0; l < nrowa; l++) {
                           atemp = *a++;
      		  	   for(j = 0; j < ncola; j++) {
	 	     	      temp = *atemp * alpha;
         	     	      for(i = 0; i < ncolb; i++) 
	    			 c[l][i] += temp * b[j][i];
	 	     	      atemp++;
      		  	   }
   	    	        }
			break;

	 /* form C := alpha * A' * B + beta * C */
	 case TRANSP:   for(l = 0; l < nrowa; l++) {
                           atemp = *a++;
      		  	   for(j = 0; j < ncola; j++) {
	 	     	      temp = *atemp * alpha;
         	     	      for(i = 0; i < ncolb; i++) 
	    			 c[j][i] += temp * b[l][i];
	 	     	      atemp++;
      		  	   }
   	    		}
			break;
      }
   }
   else { 
      switch(transa) {

	 /* form C := alpha * A * B' + beta * C */
	 case NTRANSP: for(l = 0; l < nrowa; l++) {
      		  	   for(j = 0; j < nrowb; j++) {
	 	     	      atemp = *a;
         	     	      for(i = 0; i < ncolb; i++) 
	    			 c[l][j] += (*atemp++) * alpha * b[j][i];
      		  	   }
			   a++;
   	    		}
			break;

	 /* form C := alpha * A' * B' + beta * C */
	 case TRANSP:   for(i = 0; i < ncola; i++) {
			   for (l = 0; l < nrowb; l++) {
      	       		      temp = ZERO;
	 		      for(j = 0; j < nrowa; j++) 
				 temp += a[j][i] * b[l][j];
                              c[i][l] += alpha * temp;
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
   computation of the unscaled sum of squares for the longermediate components. 
   The definitions of small, longermediate and large components depend on two
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
   double norm2, agiant, doublen, s1, s2, s3, xabs, x1max, x3max;
   long i;

   s1     = ZERO;
   s2     = ZERO;
   s3     = ZERO;
   x1max  = ZERO;
   x3max  = ZERO;
   doublen = (double)n;
   agiant = RGIANT / doublen;

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
#define		ZERO	 0.0
#define		ONE	 1.0

long imax(long, long);

/***********************************************************************
 *                                                                     *
 *                          dsbmv()                                    *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   The function performs the matrix-vector operation

	  y := alpha * A * y + beta * y,

   where alpha and beta are scalars, x and y are n-element vectors and
   A is an n by n symmetric band matrix, with k super-diagonals.

   Parameters
   ----------

   n         number of rows of matrix A; n must be at least 0.  Unchanged
	     upon exit.

   k         number of super-diagonals of matrix A

   a         2-dimensional array whose leading n by (k + 1) part must 
	     contain the upper triangular band part of the symmetric matrix,
	     supplied row by row, with the leading diagonal of the matrix 
	     in column (k + 1) of the array, the first super-diagonal 
	     starting at position 2 in column k, and so on.
	     The top left k by k triangle of the array A is not referenced.

   x         linear array of dimension of at least n.  Before entry,
	     x must contain the n elements of vector x.  Unchanged on exit.

   y         linear array of dimension of at least n.  Before entry,
	     y must contain the n elements of vector y.  On exit, y is
	     overwritten by the updated vector y.


   Functions called
   --------------

   MISC      imax  

 ***********************************************************************/

void dsbmv(long n, long k, double alpha, double **a, 
           double *x, double beta, double *y)

{
   long info, j, i, l;
   double *ptrtemp, temp1, temp2;

   info = 0;
   if ( n < 0 )                                      info = 1;
   else if ( k < 0 )                                 info = 2;

   if (info) {
      fprintf(stderr, "%s %1d %s\n",
      "*** ON ENTRY TO DSBMV, PARAMETER NUMBER",info,"HAD AN ILLEGAL VALUE");
      exit(info);
   }

   if (!n || (alpha == ZERO && beta == ONE))
      return;

   ptrtemp = y; 

   /* form y := beta * y */
   if (beta == ZERO) 
      for (i = 0; i < n; i++) *ptrtemp++ = ZERO;
   else if (beta != ONE) 
      for (i = 0; i < n; i++) *ptrtemp++ *= beta;

   if (alpha == ZERO) return;

   for (j = 0; j < n; j++) {
      temp1 = alpha * x[j];
      temp2 = ZERO;
      l = k - j;
      for (i = imax(0, j - k); i < j; i++) {
         y[i] = y[i] + temp1 * a[j][l+i];
         temp2 = temp2 + a[j][l+i] * x[i];
      }
      y[j] = y[j] + temp1 * a[j][k] + alpha * temp2;
   }
}
#include <math.h>
#define		ONE	1.0
#define		ZERO	0.0
double fsign(double, double);
double pythag(double, double);

/***********************************************************************
 *                                                                     *
 *				tql2()   			       *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   tql2() is a translation of a Fortran version of the Algol
   procedure TQL2, Num. Math. 11, 293-306(1968) by Dowdler, Martin, 
   Reinsch and Wilkinson.
   Handbook for Auto. Comp., vol.II-Linear Algebra, 227-240(1971).  

   This function finds the eigenvalues and eigenvectors of a symmetric
   tridiagonal matrix by the QL method.


   Arguments
   ---------

   (input)                                                             
   n      order of the symmetric tridiagonal matrix           
   d      contains the diagonal elements of the input matrix        
   e      contains the subdiagonal elements of the input matrix in its
            first n-1 positions.
   z      contains the identity matrix				    
                                                                   
   (output)                                                       
   d      contains the eigenvalues in ascending order.  if an error
            exit is made, the eigenvalues are correct but unordered for
            for indices 0,1,...,ierr.				   
   e      has been destroyed.					  
   z      contains orthonormal eigenvectors of the symmetric   
            tridiagonal (or full) matrix.  if an error exit is made,
            z contains the eigenvectors associated with the stored 
          eigenvalues.					
   ierr   set to zero for normal return, j if the j-th eigenvalue has
            not been determined after 30 iterations.		    


   Functions used
   --------------
   UTILITY	fsign
   MISC		pythag

 ***********************************************************************/

long tql2(long n, double *d, double *e, double **z)

{
   long j, last, l, l1, l2, m, i, k, iteration;
   double tst1, tst2, g, r, s, s2, c, c2, c3, p, f, h, el1, dl1;
   if (n == 1) return(0);
   f = ZERO;
   last = n - 1;
   tst1 = ZERO;
   e[last] = ZERO;

   for (l = 0; l < n; l++) {
      iteration = 0;
      h = fabs(d[l]) + fabs(e[l]);
      if (tst1 < h) tst1 = h;

      /* look for small sub-diagonal element */
      for (m = l; m < n; m++) {
	 tst2 = tst1 + fabs(e[m]);
	 if (tst2 == tst1) break;
      }
      if (m != l) {
	 while (iteration < 30) {
	    iteration += 1;

            /*  form shift */
	    l1 = l + 1;
	    l2 = l1 + 1;
	    g = d[l];
            p = (d[l1] - g) / (2.0 * e[l]);
	    r = pythag(p, ONE);
	    d[l] = e[l] / (p + fsign(r, p));
	    d[l1] = e[l] * (p + fsign(r, p));
	    dl1 = d[l1];
	    h = g - d[l];
	    if (l2 < n) 
	       for (i = l2; i < n; i++) d[i] -= h;
            f += h;

	    /* QL transformation */
	    p = d[m];
	    c = ONE;
	    c2 = c;
	    el1 = e[l1];
	    s = ZERO;
	    i = m - 1;
	    while (i >= l) {
	       c3 = c2;
	       c2 = c;
	       s2 = s;
	       g = c * e[i];
	       h = c * p;
	       r = pythag(p, e[i]);
	       e[i + 1] = s * r;
	       s = e[i] / r;
	       c = p / r;
	       p = c * d[i] - s * g;
	       d[i + 1]= h + s * (c * g + s * d[i]);

	       /*  form vector */
	       for (k = 0; k < n; k ++) {
	          h = z[i + 1][k];
	          z[i + 1][k] = s * z[i][k] + c * h;
	          z[i][k] = c * z[i][k] - s * h;
	       }
	       i--;
	    }
	    p = -s * s2 * c3 *el1 * e[l] / dl1;
	    e[l] = s * p;
	    d[l] = c * p;
	    tst2 = tst1 + fabs(e[l]);
	    if (tst2 <= tst1) break;
	    if (iteration == 30) 
	       return(l);
         }
      }
      d[l] += f;
   }

   /* order the eigenvalues */
   for (l = 1; l < n; l++) {
      i = l - 1;
      k = i;
      p = d[i];
      for (j = l; j < n; j++) {
	 if (d[j] < p) {
	    k = j;
	    p = d[j];
	 }
      }
      /* ...and corresponding eigenvectors */
      if (k != i) {
	 d[k] = d[i];
	 d[i] = p;
	  for (j = 0; j < n; j ++) {
	     p = z[i][j];
	     z[i][j] = z[k][j];
	     z[k][j] = p;
	  }
      }   
   }
   return(0);
}		/*...... end main ............................*/
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
#include <math.h>
double dmax(double, double);
double dmin(double, double);

/************************************************************** 
 *							      *
 * Function finds sqrt(a^2 + b^2) without overflow or         *
 * destructive underflow.				      *
 *							      *
 **************************************************************/ 
/************************************************************** 

   Funtions used
   -------------

   UTILITY	dmax, dmin

 **************************************************************/ 

double pythag(double a, double b)

{
   double p, r, s, t, u, temp;

   p = dmax(fabs(a), fabs(b));
   if (p != 0.0) {
      temp = dmin(fabs(a), fabs(b)) / p;
      r = temp * temp; 
      t = 4.0 + r;
      while (t != 4.0) {
	 s = r / t;
	 u = 1.0 + 2.0 * s;
	 p *= u;
	 temp = s / u;
	 r *= temp * temp;
	 t = 4.0 + r;
      }
   }
   return(p);
}
