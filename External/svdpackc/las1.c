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
#include "las1.h"

/* #define  UNIX_CREAT */

#ifdef UNIX_CREAT
#define PERMS 0664
#endif

long   landr(long, long, long, long, double, double, long, double,
	     double *, double *, double *);
void   dscal(long, double, double *,long);
void   daxpy(long, double, double *,long, double *, long);
double ddot(long, double *,long, double *, long);
void   opb(long, double *, double *);
void   write_data(long, long, double, double, long, double,
		  char *,char *, long, long, long);
long   check_parameters(long, long, long, double, double, long,
		        long);				
float  timer(void);

/***********************************************************************
 *                                                                     *
 *                        main()                                       *
 *     Sparse SVD Via Eigensystem of Equivalent 2-Cyclic Matrix        *
 *                  (double precision)                                 *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   This sample program uses landr to compute singular triplets of A via
   the equivalent symmetric eigenvalue problem                         

   B x = lambda x, where x' = (u',v'), lambda = +/- sigma, 
   where sigma is a singular value of A,
                                                                     
       [o  A]	
   B = [    ] , and A is m (nrow) by n (ncol) (nrow >> ncol),
       [A' o]
                                                                 
   so that {u, abs(lambda), v} is a singular triplet of A.        
   (A' = transpose of A)                                      
                                                            
   User supplied routines:  opb, store, timer              

   opb(nrow+ncol,x,y) takes an (m+n)-vector x and returns B*x in y.                                                     
   Based on operation flag isw, store(n,isw,j,s) stores/retrieves 
   to/from storage a vector of length n in s.                   
                                                               
   User should edit timer() with an appropriate call to an intrinsic
   timing routine that returns elapsed user cpu time.


   External parameters 
   -------------------

   Defined and documented in las1.h


   Local parameters 
   ----------------

  (input)
   endl     left end of interval containing unwanted eigenvalues of B
   endr     right end of interval containing unwanted eigenvalues of B
   kappa    relative accuracy of ritz values acceptable as eigenvalues
              of B
   r        work array
   n	    dimension of the eigenproblem for matrix B (nrow + ncol)
   maxprs   upper limit of desired number of singular triplets of A 
   lanmax   upper limit of desired number of Lanczos steps
   nnzero   number of nonzeros in A
   vectors  1 indicates both singular values and singular vectors are 
	      wanted and they can be found in output file lav1; 
	    0 indicates only singular values are wanted 
   		
  (output)
   ritz	    array of ritz values
   bnd      array of error bounds
   d        array of singular values
   memory   total memory allocated in bytes to solve B-eigenproblem


   Functions used
   --------------

   BLAS		daxpy, dscal, ddot
   USER		opb, timer
   MISC		write_data, check_parameters
   LAS1		landr


   Precision
   ---------

   All floating-point calculations use double precision;
   variables are declared as long and double.


   LAS1 development
   ----------------

   LAS1 is a C translation of the Fortran-77 LAS1 from the SVDPACK
   library written by Michael W. Berry, University of Tennessee,
   Dept. of Computer Science, 107 Ayres Hall, Knoxville, TN, 37996-1301

   31 Jan 1992:  Date written 

   Theresa H. Do
   University of Tennessee
   Dept. of Computer Science
   107 Ayres Hall
   Knoxville, TN, 37996-1301
   internet: tdo@cs.utk.edu

 ***********************************************************************/

main() 

{
   float t0, exetime;
   double endl, endr, kappa, tmp0, tmp1, tmp2, tmp3, xnorm;
   double *r, *ritz, *bnd, *d, *tptr1;
   long k, i, id, n, lanmax, maxprs, nnzero, size1, size2, memory, vectors;
   long *tptr3, hcount;
   char title[73], name[41], v[6];
   char *in1, *in2, *out1, *out2;
   FILE *fp_in1, *fp_in2;
 
   in1 = "lap1";
   in2 = "matrix";
   out1 = "lao1";
   out2 = "lav1";

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
    fscanf (fp_in1,"%s %ld %ld %lf %lf %s %lf", name, &lanmax,
            &maxprs, &endl, &endr, v, &kappa);

    if (!(strcmp(v, "TRUE"))) {
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

    n = nrow + ncol;

    /* write header of output file */
    write_data(lanmax, maxprs, endl, endr, vectors, kappa,
	       title, name, nrow, ncol, n);

    /* even though the validity of the parameters will be checked in the
     * SVD code, some parameter checking should also be done before
     * allocating memory to ensure that they are nonnegative */
    if (check_parameters(maxprs, lanmax, n, endl, endr, vectors, nnzero)) {
       fclose(fp_in1);
       fclose(fp_in2);
       fclose(fp_out1);
       if (vectors) close(fp_out2);
       exit(-1);
    }
    /*******************************************************************
     * allocate memory						       *
     * pointr - column start array of harwell-boeing sparse matrix     *
     *          format                                       (ncol+1)  *
     * rowind - row indices array of harwell-boeing format   (nnzero)  *
     * value  - nonzero values array of harwell-boeing sparse matrix   *
     *          format                                       (nnzero)  *
     * r      - work array                                        (n)  *
     * ritz   - array of ritz values                              (n)  *
     * bnd    - array of error bounds                             (n)  *
     * d      - array of approximate singular values of matrix A  (n)  *
     * a      - storage area for Lanczos vectors     (n * (lanmax+2))  *
     *******************************************************************/
	 size1 = sizeof(double) * (n * (lanmax + 6) + nnzero);
	 size2 = sizeof(long) * (ncol + nnzero + 1);
	 if (!(pointr = (long *)   malloc(size2))   ||
	     !(value  = (double *) malloc(size1))){
	     perror("MALLOC FAILED in MAIN()");
	     exit(errno);
         }
	 tptr1 = value;
	 
    /* calculate memory allocated for problem (including work space in landr */
    memory = size1 + size2 + sizeof(double) * (5 * n + 4 * lanmax + 1);	
     
     rowind = pointr + (ncol + 1);
     tptr1 += nnzero;
     r      = tptr1;
     tptr1 += n;
     ritz   = tptr1;
     tptr1 += n;
     bnd    = tptr1;
     tptr1 += n;
     d      = tptr1;
     tptr1 += n;
     a      = tptr1;
	 
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

    /* to get a random starting vector, the first n cells must be
     * initialize to zero */
    for (i = 0; i < n; i++) r[i] = 0.;

    exetime = timer();
    
    /* make a lanczos run; exit upon error */
    if(landr(n, lanmax, maxprs, nnzero, endl, endr, vectors, kappa,
       ritz, bnd, r)){
       free(value);
       free(pointr);
       fclose(fp_in1);
       fclose(fp_in2);
       fclose(fp_out1);
       if (vectors){
	  close(fp_out2);
          free(xv1);
          free(xv2);
       }
       exit(-1);
    }
      
    exetime = timer() - exetime;

    /* calculating singular vectors requires extra work space
     * (xv1, xv2, s) in landr() */
    if (vectors) memory += sizeof(double) * 
			   (n*(j+1) + n + (j+1)*(j+1));

    /* print error code if not zero */
    if (ierr)fprintf(fp_out1," ... RETURN FLAG = %9ld ...\n",ierr);

    /* ...Print ritz values and error bounds */
    fprintf(fp_out1,"\n");
    fprintf(fp_out1," ...... ALLOCATED MEMORY (BYTES)= %10.2E\n",(float)memory);
    fprintf(fp_out1," ...... LANSO EXECUTION TIME=%10.2E\n",exetime);
    fprintf(fp_out1," ...... \n");
    fprintf(fp_out1," ...... NUMBER OF LANCZOS STEPS = %3ld       NEIG = %3ld\n",j+1,neig);
    fprintf(fp_out1," ...... \n");
    fprintf(fp_out1," ......         COMPUTED RITZ VALUES  (ERROR BNDS)\n");
    fprintf(fp_out1," ...... \n");
    for (i = 0; i <= j; i++)
       fprintf(fp_out1," ...... %3ld   %22.14E  (%11.2E)\n",i+1,ritz[i],bnd[i]);

    /* compute residual error when singular values and vectors are	
     * computed */
    if (vectors) {
       t0 = timer();
       id = 0;

       for (i = 0; i < nsig; i++) {

	  /* multiply by 2-cyclic indefinite matrix */
	  opb(n, &xv1[id], xv2);
	  tmp0 = ddot(n, &xv1[id], 1, xv2, 1);
	  tmp1 = ddot(nrow, &xv1[id], 1, &xv1[id], 1);
	  tmp2 = ddot(ncol, &xv1[id + nrow], 1, &xv1[id + nrow], 1);
	  tmp3 = sqrt(tmp1 + tmp2);
	  daxpy(n, -tmp0, &xv1[id], 1, xv2, 1);
	  xnorm = ddot(n, xv2, 1, xv2, 1);
	  xnorm = sqrt(xnorm);
	  bnd[i] = xnorm / tmp3;
	  d[i] = fabs(tmp0);
	  id += n;
       }  
       exetime += (timer() - t0);
       hcount=mxvcount/2;
       fprintf(fp_out1," ...... \n");
       fprintf(fp_out1," ...... NO. MULTIPLICATIONS BY A  =%10ld\n", hcount);
       fprintf(fp_out1," ...... NO. MULT. BY TRANSPOSE(A) =%10ld\n", hcount);
       fprintf(fp_out1,"\n");
       fprintf(fp_out1," ...... LASVD EXECUTION TIME=%10.2E\n",exetime);
       fprintf(fp_out1," ...... \n");
       fprintf(fp_out1," ......         NSIG = %4ld\n", nsig);
       fprintf(fp_out1," ...... \n");
       fprintf(fp_out1," ......         COMPUTED S-VALUES     (RES. NORMS)\n");
       fprintf(fp_out1," ...... \n");
       for (i = 0; i < nsig; i++)
          fprintf(fp_out1," ...... %3ld   %22.14E  (%11.2E)\n",
		  i + 1, d[i], bnd[i]);
    }
    else {
       for (i = j; i >= 0; i--)
	  if (bnd[i] > kappa * fabs(ritz[i])) break;
       nsig = j - i;
       hcount=mxvcount/2;
       fprintf(fp_out1," ...... \n");
       fprintf(fp_out1," ...... NO. MULTIPLICATIONS BY A  =%10ld\n", hcount);
       fprintf(fp_out1," ...... NO. MULT. BY TRANSPOSE(A) =%10ld\n", hcount);
       fprintf(fp_out1,"\n");
       fprintf(fp_out1," ...... LASVD EXECUTION TIME=%10.2E\n", exetime);
       fprintf(fp_out1," ...... \n");
       fprintf(fp_out1," ......         NSIG = %4ld\n", nsig);
       fprintf(fp_out1," ...... \n");
       fprintf(fp_out1," ......         COMPUTED S-VALUES   (ERROR BNDS)\n");
       fprintf(fp_out1," ...... \n");

       k = j + 1 - nsig;
       for (i = 1 ; i <= nsig; i++)  {
          fprintf(fp_out1," ...... %3ld   %22.14E  (%11.2E)\n", 
                  i, ritz[k], bnd[k]);
          k++;
          }
   }

    free(value);
    free(pointr);
    fclose(fp_in1);
    fclose(fp_in2);
    fclose(fp_out1);
    if (vectors) {
       free(xv1);
       free(xv2);
       close(fp_out2);
    }
    exit(0);
}

extern long ncol,nrow;
extern char *error[];
extern FILE *fp_out1;
/***********************************************************************
 *								       *
 *		      check_parameters()			       *
 *								       *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------
   Function validates input parameters and returns error code (long)  

   Parameters 
   ----------
  (input)
   maxprs   upper limit of desired number of eigenpairs of B           
   lanmax   upper limit of desired number of lanczos steps             
   n        dimension of the eigenproblem for matrix B               
   endl     left end of interval containing unwanted eigenvalues of B
   endr     right end of interval containing unwanted eigenvalues of B
   vectors  1 indicates both eigenvalues and eigenvectors are wanted 
            and they can be found in lav1; 0 indicates eigenvalues only
   nnzero   number of nonzero elements in input matrix (matrix A)      
                                                                      
 ***********************************************************************/

long check_parameters(long maxprs, long lanmax, long n,
		      double endl, double endr, long vectors,
		      long nnzero) 
{
   long error_index, ncells;
   error_index = 0;

   /* assuming that nrow >= ncol... */
   if (ncol >= NMAX || nnzero > NZMAX) error_index = 1;
   else if (endl >= endr)  error_index = 2;
   else if (maxprs > lanmax)  error_index = 3;
   else if (n <= 0)  error_index = 4;
   else if (lanmax <= 0 || lanmax > n)  error_index = 5;
   else if (maxprs <= 0 || maxprs > lanmax)  error_index = 6;
   else {
       if (vectors) {
	  ncells = 6 * n + 4 * lanmax + 1 + lanmax * lanmax;
	  if (ncells > LMTNW) error_index = 7;
       }
       else {
	  ncells = 6 * n + 4 * lanmax + 1;
	  if (ncells > LMTNW) error_index = 8;
       }
   }
   if (error_index) fprintf(fp_out1, "%s\n", error[error_index]);
   return(error_index);
}

extern FILE *fp_out1;
/***********************************************************************
 *								       *
 *			  write_data()				       *
 *   Function writes out header of output file containing ritz values  *
 *								       *
 ***********************************************************************/

void write_data(long lanmax, long maxprs, double endl, double endr,
	        long vectors, double kappa, char *title,
		char *name, long nrow, long ncol, long n)

{
   fprintf(fp_out1, " ... \n");
   fprintf(fp_out1, " ... SOLVE THE CYCLIC   EIGENPROBLEM\n");
   fprintf(fp_out1, " ... NO. OF EQUATIONS          =%5ld\n", n);
   fprintf(fp_out1, " ... MAX. NO. OF LANCZOS STEPS =%5ld\n", lanmax);
   fprintf(fp_out1, " ... MAX. NO. OF EIGENPAIRS    =%5ld\n", maxprs);
   fprintf(fp_out1, " ... LEFT  END OF THE INTERVAL =%10.2E\n", endl);
   fprintf(fp_out1, " ... RIGHT END OF THE INTERVAL =%10.2E\n", endr);


 if (vectors) 
   fprintf(fp_out1, " ... WANT S-VECTORS?   [T/F]   =   T\n");
 else
   fprintf(fp_out1, " ... WANT S-VECTORS?   [T/F]   =   F\n");

 fprintf(fp_out1, " ... KAPPA                     =%10.2E\n", kappa);
 fprintf(fp_out1, " %s\n", title);
 fprintf(fp_out1, "           %s\n", name);
 fprintf(fp_out1, " ... NO. OF TERMS     (ROWS)   = %8ld\n", nrow);
 fprintf(fp_out1, " ... NO. OF DOCUMENTS (COLS)   = %8ld\n", ncol);
 fprintf(fp_out1, " ... ORDER OF MATRIX A         = %8ld\n", n);
 fprintf(fp_out1, " ... \n");
 return;
}
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <fcntl.h>
extern double eps, eps1, reps, eps34, *xv1, *xv2;
extern long j;
void   machar(long *, long *, long *, long *, long *);
long   check_parameters(long, long, long, double, double, long,
			long);
double dmax(double, double);
void   lanso(long, long, long, double, double, double *, double *,
	     double *[]);
void   ritvec(long, double, double *, double *, double *, double *,
              double *, double *);

/***********************************************************************
 *                                                                     *
 *				landr()				       *
 *        Lanczos algorithm with selective orthogonalization           *
 *                    Using Simon's Recurrence                         *
 *                       (double precision)                            *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   landr() is the LAS1 driver routine that, upon entry, 
     (1)  checks for the validity of input parameters of the 
	  B-eigenproblem 
     (2)  determines several machine constants
     (3)  makes a Lanczos run
     (4)  calculates B-eigenvectors (singular vectors of A) if requested
	    by user


   arguments
   ---------

   (input)
   n        dimension of the eigenproblem for matrix B                 
   lanmax   upper limit of desired number of Lanczos steps            
   maxprs   upper limit of desired number of eigenpairs              
   nnzero   number of nonzeros in matrix A
   endl     left end of interval containing unwanted eigenvalues of B
   endr     right end of interval containing unwanted eigenvalues of B
   vectors  1 indicates both eigenvalues and eigenvectors are wanted
              and they can be found in output file lav1; 
	    0 indicates only eigenvalues are wanted
   kappa    relative accuracy of ritz values acceptable as eigenvalues
	      of B (singular values of A)
   r        work array

   (output)
   j        number of Lanczos steps actually taken                     
   neig     number of ritz values stabilized                           
   ritz     array to hold the ritz values                              
   bnd      array to hold the error bounds


   External parameters
   -------------------

   Defined and documented in las1.h


   local parameters
   -------------------

   ibeta    radix for the floating-point representation
   it       number of base ibeta digits in the floating-point significand
   irnd     floating-point addition rounded or chopped
   machep   machine relative precision or round-off error
   negeps   largest negative integer
   wptr	    array of pointers each pointing to a work space


   Functions used
   --------------

   MISC         dmax, machar, check_parameters
   LAS1         ritvec, lanso

 ***********************************************************************/

long landr(long n, long lanmax, long maxprs, long nnzero, double endl,
	   double endr, long vectors, double kappa, double *ritz,
	   double *bnd, double *r)

{
   long i, size, ibeta, it, irnd, machep, negep;
   double *wptr[10], *tptr, *tptr2;

   /* data validation */
      if (check_parameters(maxprs, lanmax, n, endl, endr, vectors, nnzero))
	  return(-1);

   /* compute machine precision */
   machar(&ibeta, &it, &irnd, &machep, &negep);

   eps1 = eps * sqrt( (double) n );
   reps = sqrt(eps);
   eps34 = reps * sqrt(reps);


   size = 5*n + (lanmax*4 + 1);
   tptr = NULL;

   /* allocate work area and initialize pointers         *
    * pointer            name               size   	 *
    * wptr[0]             r                  n     	 *
    * wptr[1]             q		     n     	 *
    * wptr[2]             q_previous         n     	 *
    * wptr[3]             p		     n     	 *
    * wptr[4]             p_previous         n     	 *
    * wptr[5]             wrk                n     	 *
    * wptr[6]             alf              lanmax  	 *
    * wptr[7]             eta              lanmax  	 *
    * wptr[8]             oldeta           lanmax  	 *
    * wptr[9]             bet              lanmax+1	 */

   if (!(tptr = (double *) malloc(size * sizeof(double)))){
      perror("FIRST MALLOC FAILED in LANDR()");
      exit(errno);
   }

   tptr2 = tptr;
   wptr[0] = r;
   for (i = 1; i <= 5; i++) {
      wptr[i] = tptr;
      tptr += n;
   }
   for (i = 6; i <= 9; i++) {
      wptr[i] = tptr;
      tptr += lanmax;
   }

   lanso(n, lanmax, maxprs, endl, endr, ritz, bnd, wptr);

   /* compute eigenvectors */
   if (vectors) {
      if (!(xv1 = (double *) malloc(n * (j+1) * sizeof(double))) ||
	  !(xv2 = (double *) malloc(n * sizeof(double)))) {
	  perror("SECOND MALLOC FAILED in LANDR()");
	  exit(errno);
	 }
      kappa = dmax(fabs(kappa), eps34);
      ritvec(n, kappa, ritz, bnd, wptr[6], wptr[9], wptr[4], wptr[5]);
   }

   free(tptr2);
   return(0);
}

#define STORQ 1
#define RETRQ 2
#define STORP 3
#define RETRP 4

extern long ierr, j, nsig, fp_out2;
extern double *xv1;
void   dscal(long, double, double *,long);
void   dcopy(long, double *, long, double *, long);
void   daxpy(long, double, double *,long, double *, long);
void   store(long, long, long, double *);
void   imtql2(long, long, double *, double *, double *);
/***********************************************************************
 *                                                                     *
 *                        ritvec()                                     *
 * 	    Function computes the singular vectors of matrix A	       *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   This function is invoked by landr() only if eigenvectors of the cyclic
   eigenproblem are desired.  When called, ritvec() computes the 
   singular vectors of A and writes the result to an unformatted file.


   Parameters
   ----------

   (input)
   n	      dimension of the eigenproblem for matrix B
   kappa      relative accuracy of ritz values acceptable as 
		eigenvalues of B (singular values of A)
   ritz       array of ritz values
   bnd        array of error bounds
   alf        array of diagonal elements of the tridiagonal matrix T
   bet        array of off-diagonal elements of T
   w1, w2     work space
   fp_out2    pointer to unformatted output file
   j	      number of Lanczos iterations performed

   (output)
   xv1        array of eigenvectors of B (singular vectors of A)
   ierr	      error code
	      0 for normal return from imtql2()
	      k if convergence did not occur for k-th eigenvalue in 
		imtql2()
   nsig       number of accepted ritz values based on kappa

   (local)
   s	      work array which is initialized to the identity matrix of
	      order (j + 1) upon calling imtql2().  After the call, 
	      s contains the orthonormal eigenvectors of the symmetric 
	      tridiagonal matrix T

   Functions used
   --------------

   BLAS		dscal, dcopy, daxpy
   USER		store
   		imtql2

 ***********************************************************************/

void ritvec(long n, double kappa, double *ritz, double *bnd,
	    double *alf, double *bet, double *w1, double *w2)

{
   long js, jsq, i, k, size, id, id2, tmp;
   double *s;

   js = j + 1;
   jsq = js * js;
   size = sizeof(double) * n;

   if(!(s = (double *) malloc (jsq * sizeof(double)))) {
      perror("MALLOC FAILED in RITVEC()");
      exit(errno);
   }

   /* initialize s as an identity matrix */
   for (i = 0; i < jsq; i++) s[i] = 0.0;
   for (i = 0; i < jsq; i+= (js + 1)) s[i] = 1.0;
   dcopy(js, alf, 1, w1, -1);
   dcopy(j, &bet[1], 1, &w2[1], -1);

   /* compute eigenvalues of T */
   imtql2(js, js, w1, w2, s);
   if (ierr) return;

   /* on return w1 contains eigenvalues in ascending order and s 
    * contains the corresponding eigenvectors */

   write(fp_out2, (char *)&n, sizeof(n));
   write(fp_out2, (char *)&js, sizeof(js));
   write(fp_out2, (char *)&kappa, sizeof(kappa));
   nsig = 0;
   id = 0;
   id2 = jsq-js;
   for (k = 0; k < js; k++) {
      tmp = id2;
      if (bnd[k] <= kappa * fabs(ritz[k]) ) {
         for (i = 0; i < n; i++) w1[i] = 0.0;
         for (i = 0; i < js; i++) {
	    store(n, RETRQ, i, w2);
	    daxpy(n, s[tmp], w2, 1, w1, 1);
	    tmp -= js;
         }
         write(fp_out2, (char *)w1, size);

         /* store the w1 vector sequentially in array xv1 */
         for (i = 0; i < n; i++) xv1[id++] = w1[i];
         nsig += 1;
      }
      id2++;
   }
   free(s);
   return;
}
extern long ncol, nrow, mxvcount;
extern long *pointr, *rowind;
extern double *value;

/************************************************************** 
 * multiplication of 2-cyclic matrix B by a vector x, where   *
 *							      *
 * B = [0  A]						      *
 *     [A' 0] , where A is nrow by ncol (nrow >> ncol). Hence,*
 * B is of order n = nrow+ncol (y stores product vector).     *
 **************************************************************/ 
void opb(long n,double *x, double *y)

{
   long i, j, end;
   double *tmp;
   
   mxvcount += 2;
   for (i = 0; i < n; i++) y[i] = 0.0;

   tmp = &x[nrow]; 
   for (i = 0; i < ncol; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++) 
	 y[rowind[j]] += value[j] * (*tmp); 
      tmp++; 
   }
   for (i = nrow; i < n; i++) {
      end = pointr[i-nrow+1];
      for (j = pointr[i-nrow]; j < end; j++) 
	 y[i] += value[j] * x[rowind[j]];
   }
   return;
}
#include <stdio.h>
#include <math.h>
#define STORQ 1
#define RETRQ 2
#define STORP 3
#define RETRP 4
#define TRUE  1
#define FALSE 0

extern double rnm, anorm, tol, eps, eps1, reps, eps34;
extern long ierr, j, neig;
void   stpone(long, double *[]);
void   error_bound(long *, double, double, double *, double *);
void   lanczos_step(long, long, long, double *[], double *,
                    double *, double *, double *, long *, long *);
long   imin(long, long);
long   imax(long, long);
void   dsort2(long, long, double *, double *);
void   imtqlb(long n, double d[], double e[], double bnd[]);

/***********************************************************************
 *                                                                     *
 *                          lanso()                                    *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   Function determines when the restart of the Lanczos algorithm should 
   occur and when it should terminate.

   Arguments 
   ---------

   (input)
   n         dimension of the eigenproblem for matrix B
   lanmax    upper limit of desired number of lanczos steps           
   maxprs    upper limit of desired number of eigenpairs             
   endl      left end of interval containing unwanted eigenvalues
   endr      right end of interval containing unwanted eigenvalues
   ritz      array to hold the ritz values                       
   bnd       array to hold the error bounds                          
   wptr      array of pointers that point to work space:            
  	       wptr[0]-wptr[5]  six vectors of length n		
  	       wptr[6] array to hold diagonal of the tridiagonal matrix T
  	       wptr[9] array to hold off-diagonal of T	
  	       wptr[7] orthogonality estimate of Lanczos vectors at 
		 step j
 	       wptr[8] orthogonality estimate of Lanczos vectors at 
		 step j-1

   (output)
   j         number of Lanczos steps actually taken
   neig      number of ritz values stabilized
   ritz      array to hold the ritz values
   bnd       array to hold the error bounds
   ierr      (globally declared) error flag
	     ierr = 8192 if stpone() fails to find a starting vector
	     ierr = k if convergence did not occur for k-th eigenvalue
		    in imtqlb()
	     ierr = 0 otherwise


   Functions used
   --------------

   LAS		stpone, error_bound, lanczos_step
   MISC		dsort2
   UTILITY	imin, imax

 ***********************************************************************/

void lanso(long n, long lanmax, long maxprs, double endl,
	   double endr, double *ritz, double *bnd, double *wptr[])

{
   double *r, *alf, *eta, *oldeta, *bet, *wrk;
   long ll, first, last, ENOUGH, id1, id2, id3, i, l;

   r = wptr[0];
   alf = wptr[6];
   eta = wptr[7];
   oldeta = wptr[8];
   bet = wptr[9];
   wrk = wptr[5];
   j = 0;

   /* take the first step */
   stpone(n, wptr);
   if (!rnm || ierr) return;
   eta[0] = eps1;
   oldeta[0] = eps1;
   ll = 0;
   first = 1;
   last = imin(maxprs + imax(8,maxprs), lanmax);
   ENOUGH = FALSE;
   id1 = 0;
   while (id1 < maxprs && !ENOUGH) {
      if (rnm <= tol) rnm = 0.0;

      /* the actual lanczos loop */
      lanczos_step(n, first, last, wptr, alf, eta, oldeta, bet, &ll,
		   &ENOUGH);
      if (ENOUGH) j = j - 1;
      else j = last - 1;
      first = j + 1;
      bet[j+1] = rnm;

      /* analyze T */
      l = 0;
      for (id2 = 0; id2 < j; id2++) {
	 if (l > j) break;
         for (i = l; i <= j; i++) if (!bet[i+1]) break;
	 if (i > j) i = j;

	 /* now i is at the end of an unreduced submatrix */
	 dcopy(i-l+1, &alf[l],   1, &ritz[l],  -1);
	 dcopy(i-l,   &bet[l+1], 1, &wrk[l+1], -1);

	 imtqlb(i-l+1, &ritz[l], &wrk[l], &bnd[l]);

	 if (ierr) {
	    printf("IMTQLB FAILED TO CONVERGE (IERR = %d)\n",
		    ierr);
	    printf("L = %d    I = %d\n", l, i);
            for (id3 = l; id3 <= i; id3++) 
	       printf("%d  %lg  %lg  %lg\n",
	       id3, ritz[id3], wrk[id3], bnd[id3]);
	 }
         for (id3 = l; id3 <= i; id3++) 
	    bnd[id3] = rnm * fabs(bnd[id3]);
	 l = i + 1;
      }

      /* sort eigenvalues into increasing order */
      dsort2((j+1) / 2, j + 1, ritz, bnd);

      /* massage error bounds for very close ritz values */
      error_bound(&ENOUGH, endl, endr, ritz, bnd);

      /* should we stop? */
      if (neig < maxprs) {
	 if (!neig) last = first + 9;
	 else last = first + imax(3, 1 + ((j-5) * (maxprs-neig)) / neig);
	 last = imin(last, lanmax);
      }
      else ENOUGH = TRUE;
      ENOUGH = ENOUGH || first >= lanmax;
      id1++;
   }
   store(n, STORQ, j, wptr[1]);
   return;
}
#include <math.h>
#define STORQ 1
#define RETRQ 2
#define STORP 3
#define RETRP 4
#define TRUE  1
#define FALSE 0
#define MAXLL 2

extern double rnm, anorm, tol, eps, eps1, reps, eps34;
extern long ierr, j;
double ddot(long, double *,long, double *, long);
void   dscal(long, double, double *,long);
void   daxpy(long, double, double *,long, double *, long);
void   datx(long, double, double *,long, double *, long);
void   dcopy(long, double *, long, double *, long);
void   purge(long, long, double *, double *, double *, double *,
	     double *, double *, double *);
void   ortbnd(double *, double *, double *, double *);
double startv(long, double *[]);
void   store(long, long, long, double *);
long   imin(long, long);
long   imax(long, long);

/***********************************************************************
 *                                                                     *
 *			lanczos_step()                                 *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   Function embodies a single Lanczos step

   Arguments 
   ---------

   (input)
   n        dimension of the eigenproblem for matrix B
   first    start of index through loop				      
   last     end of index through loop				     
   wptr	    array of pointers pointing to work space		    
   alf	    array to hold diagonal of the tridiagonal matrix T
   eta      orthogonality estimate of Lanczos vectors at step j   
   oldeta   orthogonality estimate of Lanczos vectors at step j-1
   bet      array to hold off-diagonal of T                     
   ll       number of intitial Lanczos vectors in local orthog. 
              (has value of 0, 1 or 2)			
   enough   stop flag			

   Functions used
   --------------

   BLAS		ddot, dscal, daxpy, datx, dcopy
   USER		store
   LAS		purge, ortbnd, startv
   UTILITY	imin, imax

 ***********************************************************************/

void lanczos_step(long n, long first, long last, double *wptr[],
		  double *alf, double *eta, double *oldeta,
		  double *bet, long *ll, long *enough)

{
   double t, *mid;
   long i;

   for (j=first; j<last; j++) {
      mid     = wptr[2];
      wptr[2] = wptr[1];
      wptr[1] = mid;
      mid     = wptr[3];
      wptr[3] = wptr[4];
      wptr[4] = mid;

      store(n, STORQ, j-1, wptr[2]);
      if (j-1 < MAXLL) store(n, STORP, j-1, wptr[4]);
      bet[j] = rnm;

      /* restart if invariant subspace is found */
      if (!bet[j]) {
	 rnm = startv(n, wptr);
	 if (ierr) return;
	 if (!rnm) *enough = TRUE;
      }
      if (*enough) break;

      /* take a lanczos step */
      t = 1.0 / rnm;
      datx(n, t, wptr[0], 1, wptr[1], 1);
      dscal(n, t, wptr[3], 1);
      opb(n, wptr[3], wptr[0]);
      daxpy(n, -rnm, wptr[2], 1, wptr[0], 1);
      alf[j] = ddot(n, wptr[0], 1, wptr[3], 1);
      daxpy(n, -alf[j], wptr[1], 1, wptr[0], 1);

      /* orthogonalize against initial lanczos vectors */
      if (j <= MAXLL && (fabs(alf[j-1]) > 4.0 * fabs(alf[j])))
	 *ll = j;  
      for (i=0; i < imin(*ll, j-1); i++) {
	 store(n, RETRP, i, wptr[5]);
	 t = ddot(n, wptr[5], 1, wptr[0], 1);
	 store(n, RETRQ, i, wptr[5]);
         daxpy(n, -t, wptr[5], 1, wptr[0], 1);
	 eta[i] = eps1;
	 oldeta[i] = eps1;
      }

      /* extended local reorthogonalization */
      t = ddot(n, wptr[0], 1, wptr[4], 1);
      daxpy(n, -t, wptr[2], 1, wptr[0], 1);
      if (bet[j] > 0.0) bet[j] = bet[j] + t;
      t = ddot(n, wptr[0], 1, wptr[3], 1);
      daxpy(n, -t, wptr[1], 1, wptr[0], 1);
      alf[j] = alf[j] + t;
      dcopy(n, wptr[0], 1, wptr[4], 1);
      rnm = sqrt(ddot(n, wptr[0], 1, wptr[4], 1));
      anorm = bet[j] + fabs(alf[j]) + rnm;
      tol = reps * anorm;

      /* update the orthogonality bounds */
      ortbnd(alf, eta, oldeta, bet);

      /* restore the orthogonality state when needed */
      purge(n,*ll,wptr[0],wptr[1],wptr[4],wptr[3],wptr[5],eta,oldeta);
      if (rnm <= tol) rnm = 0.0;
   }
   return;
}
extern double rnm, eps, eps1, reps, eps34;
extern long j;
void   dswap(long, double *, long, double *, long);

/***********************************************************************
 *                                                                     *
 *                          ortbnd()                                   *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   Funtion updates the eta recurrence

   Arguments 
   ---------

   (input)
   alf      array to hold diagonal of the tridiagonal matrix T         
   eta      orthogonality estimate of Lanczos vectors at step j        
   oldeta   orthogonality estimate of Lanczos vectors at step j-1     
   bet      array to hold off-diagonal of T                          
   n        dimension of the eigenproblem for matrix B		    
   j        dimension of T					  
   rnm	    norm of the next residual vector			 
   eps1	    roundoff estimate for dot product of two unit vectors

   (output)
   eta      orthogonality estimate of Lanczos vectors at step j+1     
   oldeta   orthogonality estimate of Lanczos vectors at step j        


   Functions used
   --------------

   BLAS		dswap

 ***********************************************************************/

void ortbnd(double *alf, double *eta, double *oldeta, double *bet)

{
   long i;
   if (j < 1) return;
   if (rnm) {
      if (j > 1) {
	 oldeta[0] = (bet[1] * eta[1] + (alf[0]-alf[j]) * eta[0] -
		      bet[j] * oldeta[0]) / rnm + eps1;
      }
      for (i=1; i<=j-2; i++) 
	 oldeta[i] = (bet[i+1] * eta[i+1] + (alf[i]-alf[j]) * eta[i] +
		      bet[i] * eta[i-1] - bet[j] * oldeta[i])/rnm + eps1;
   }
   oldeta[j-1] = eps1;
   dswap(j, oldeta, 1, eta, 1);  
   eta[j] = eps1;
   return;
}
#include <math.h>
#define STORQ 1
#define RETRQ 2
#define STORP 3
#define RETRP 4
#define TRUE  1
#define FALSE 0

extern double tol, rnm, eps, eps1, reps, eps34;
extern long j;
void   store(long, long, long, double *);
void   daxpy(long, double, double *, long, double *, long);
void   dcopy(long, double *, long, double *, long);
long   idamax(long, double *, long);
double ddot(long, double *, long, double *, long);

/***********************************************************************
 *                                                                     *
 *				purge()                                *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   Function examines the state of orthogonality between the new Lanczos
   vector and the previous ones to decide whether re-orthogonalization 
   should be performed


   Arguments 
   ---------

   (input)
   n        dimension of the eigenproblem for matrix B		       
   ll       number of intitial Lanczos vectors in local orthog.       
   r        residual vector to become next Lanczos vector            
   q        current Lanczos vector			           
   ra       previous Lanczos vector
   qa       previous Lanczos vector
   wrk      temporary vector to hold the previous Lanczos vector
   eta      state of orthogonality between r and prev. Lanczos vectors 
   oldeta   state of orthogonality between q and prev. Lanczos vectors
   j        current Lanczos step				     

   (output)
   r	    residual vector orthogonalized against previous Lanczos 
	      vectors
   q        current Lanczos vector orthogonalized against previous ones


   Functions used
   --------------

   BLAS		daxpy,  dcopy,  idamax,  ddot
   USER		store

 ***********************************************************************/

void purge(long n, long ll, double *r, double *q, double *ra,  
	   double *qa, double *wrk, double *eta, double *oldeta)

{
   double t, tq, tr, reps1;
   long k, iteration, flag, i;

   if (j < ll+2) return; 

   k = idamax(j - (ll+1), &eta[ll], 1) + ll;
   if (fabs(eta[k]) > reps) {
      reps1 = eps1 / reps;
      iteration = 0;
      flag = TRUE;
      while (iteration < 2 && flag) {
	 if (rnm > tol) {

	    /* bring in a lanczos vector t and orthogonalize both 
	     * r and q against it */
	    tq = 0.0;
	    tr = 0.0;
            for (i = ll; i < j; i++) {
	       store(n,  RETRQ,  i,  wrk);
	       t   = -ddot(n, qa, 1, wrk, 1);
	       tq += fabs(t);
	       daxpy(n,  t,  wrk,  1,  q,  1);
	       t   = -ddot(n, ra, 1, wrk, 1);
	       tr += fabs(t);
	       daxpy(n, t, wrk, 1, r, 1);
	    }
	    dcopy(n, q, 1, qa, 1);
	    t   = -ddot(n, r, 1, qa, 1);
	    tr += fabs(t);
	    daxpy(n, t, q, 1, r, 1);
	    dcopy(n, r, 1, ra, 1);
	    rnm = sqrt(ddot(n, ra, 1, r, 1));
	    if (tq <= reps1 && tr <= reps1 * rnm) flag = FALSE;
	 }
	 iteration++;
      }
      for (i = ll; i <= j; i++) { 
	 eta[i] = eps1;
	 oldeta[i] = eps1;
      }
   }
   return;
}
#include <math.h>

extern double rnm, anorm, tol, eps, reps;
extern long j, ierr;
void   daxpy(long, double, double *,long, double *, long);
void   datx(long, double, double *,long, double *, long);
void   dcopy(long, double *, long, double *, long);
double ddot(long, double *,long, double *, long);
void   dscal(long, double, double *,long);
double startv(long, double *[]);
void   opb(long, double *, double *);
void   store(long, long, long, double *);

/***********************************************************************
 *                                                                     *
 *                         stpone()                                    *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   Function performs the first step of the Lanczos algorithm.  It also
   does a step of extended local re-orthogonalization.

   Arguments 
   ---------

   (input)
   n      dimension of the eigenproblem for matrix B

   (output)
   ierr   error flag
   wptr   array of pointers that point to work space that contains
	    wptr[0]             r[j]
	    wptr[1]             q[j]
	    wptr[2]             q[j-1]
	    wptr[3]             p
	    wptr[4]             p[j-1]
	    wptr[6]             diagonal elements of matrix T 


   Functions used
   --------------

   BLAS		daxpy, datx, dcopy, ddot, dscal
   USER		store, opb
   LAS		startv

 ***********************************************************************/

void stpone(long n, double *wrkptr[])

{
   double t, *alf;
   alf = wrkptr[6];

   /* get initial vector; default is random */
   rnm = startv(n, wrkptr);
   if (rnm == 0.0 || ierr != 0) return;

   /* normalize starting vector */
   t = 1.0 / rnm;
   datx(n, t, wrkptr[0], 1, wrkptr[1], 1);
   dscal(n, t, wrkptr[3], 1);

   /* take the first step */
   opb(n, wrkptr[3], wrkptr[0]);
   alf[0] = ddot(n, wrkptr[0], 1, wrkptr[3], 1);
   daxpy(n, -alf[0], wrkptr[1], 1, wrkptr[0], 1);
   t = ddot(n, wrkptr[0], 1, wrkptr[3], 1);
   daxpy(n, -t, wrkptr[1], 1, wrkptr[0], 1);
   alf[0] += t;
   dcopy(n, wrkptr[0], 1, wrkptr[4], 1);
   rnm = sqrt(ddot(n, wrkptr[0], 1, wrkptr[4], 1));
   anorm = rnm + fabs(alf[0]);
   tol = reps * anorm;
   return;
}
#include <stdio.h>
#include <math.h>
#define STORQ 1
#define RETRQ 2
#define STORP 3
#define RETRP 4
extern long j, ierr;
extern double eps;
double ddot(long, double *,long, double *, long);
void   daxpy(long, double, double *,long, double *, long);
void   dcopy(long, double *, long, double *, long);
static double _random(long *);
void   store(long, long, long, double *);
void   opb(long, double *, double *);

/***********************************************************************
 *                                                                     *
 *                         startv()                                    *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   Function delivers a starting vector in r and returns |r|; it returns 
   zero if the range is spanned, and ierr is non-zero if no starting 
   vector within range of operator can be found.

   Parameters 
   ---------

   (input)
   n      dimension of the eigenproblem matrix B
   wptr   array of pointers that point to work space
   j      starting index for a Lanczos run
   eps    machine epsilon (relative precision)

   (output)
   wptr   array of pointers that point to work space that contains
	  r[j], q[j], q[j-1], p[j], p[j-1]
   ierr   error flag (nonzero if no starting vector can be found)

   Functions used
   --------------

   BLAS		ddot, dcopy, daxpy
   USER		opb, store
   MISC		random

 ***********************************************************************/

double startv(long n, double *wptr[])

{
   double rnm2, *r, t;
   long irand;
   long id, i;

   /* get initial vector; default is random */
   rnm2 = ddot(n, wptr[0], 1, wptr[0], 1);
   irand = 918273 + j;
   r = wptr[0];
   for (id = 0; id < 3; id++) {
      if (id > 0 || j > 0 || rnm2 == 0) 
	 for (i = 0; i < n; i++) r[i] = _random(&irand);
      dcopy(n, wptr[0], 1, wptr[3], 1);

      /* apply operator to put r in range (essential if m singular) */
      opb(n, wptr[3], wptr[0]);
      dcopy(n, wptr[0], 1, wptr[3], 1);
      rnm2 = ddot(n, wptr[0], 1, wptr[3], 1);
      if (rnm2 > 0.0) break;
   }

   /* fatal error */
   if (rnm2 <= 0.0) {
      ierr = 8192;
      return(-1);
   }
   if (j > 0) {
      for (i = 0; i < j; i++) {
         store(n, RETRQ, i, wptr[5]);
	 t = -ddot(n, wptr[3], 1, wptr[5], 1);
	 daxpy(n, t, wptr[5], 1, wptr[0], 1);
      }

      /* make sure q[j] is orthogonal to q[j-1] */
      t = ddot(n, wptr[4], 1, wptr[0], 1);
      daxpy(n, -t, wptr[2], 1, wptr[0], 1);
      dcopy(n, wptr[0], 1, wptr[3], 1);
      t = ddot(n, wptr[3], 1, wptr[0], 1);
      if (t <= eps * rnm2) t = 0.0;
      rnm2 = t;
   }
   return(sqrt(rnm2));
}

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
static double _random(long *iy)

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
#include <math.h>
extern double tol, eps34, eps;
extern long j, neig;
long   idamax(long, double *, long);
double dmin(double, double);

/***********************************************************************
 *                                                                     *
 *			error_bound()                                  *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   Function massages error bounds for very close ritz values by placing 
   a gap between them.  The error bounds are then refined to reflect 
   this.


   Arguments 
   ---------

   (input)
   endl     left end of interval containing unwanted eigenvalues
   endr     right end of interval containing unwanted eigenvalues
   ritz     array to store the ritz values
   bnd      array to store the error bounds
   enough   stop flag


   Functions used
   --------------

   BLAS		idamax
   UTILITY	dmin

 ***********************************************************************/

void error_bound(long *enough, double endl, double endr, 
		 double *ritz, double *bnd)
{
   long mid, i;
   double gapl, gap;

   /* massage error bounds for very close ritz values */
   mid = idamax(j + 1, bnd, 1);

   for (i=((j+1) + (j-1)) / 2; i >= mid + 1; i -= 1)
      if (fabs(ritz[i-1] - ritz[i]) < eps34 * fabs(ritz[i])) 
         if (bnd[i] > tol && bnd[i-1] > tol) {
	    bnd[i-1] = sqrt(bnd[i] * bnd[i] + bnd[i-1] * bnd[i-1]);
	    bnd[i] = 0.0;
         }
	 

   for (i=((j+1) - (j-1)) / 2; i <= mid - 1; i +=1 ) 
      if (fabs(ritz[i+1] - ritz[i]) < eps34 * fabs(ritz[i])) 
	 if (bnd[i] > tol && bnd[i+1] > tol) {
	    bnd[i+1] = sqrt(bnd[i] * bnd[i] + bnd[i+1] * bnd[i+1]);
	    bnd[i] = 0.0;
         }

   /* refine the error bounds */
   neig = 0;
   gapl = ritz[j] - ritz[0];
   for (i = 0; i <= j; i++) {
      gap = gapl;
      if (i < j) gapl = ritz[i+1] - ritz[i];
      gap = dmin(gap, gapl);
      if (gap > bnd[i]) bnd[i] = bnd[i] * (bnd[i] / gap);
      if (bnd[i] <= 16.0 * eps * fabs(ritz[i])) {
	 neig += 1;
	 if (!*enough) *enough = endl < ritz[i] && ritz[i] < endr;
      }
   }   
   return;
}
#include <math.h>
#define		TRUE	1
#define		FALSE	0
extern long ierr;
double pythag(double, double);
double fsign(double, double);

/***********************************************************************
 *                                                                     *
 *				imtqlb()			       *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   imtqlb() is a translation of a Fortran version of the Algol
   procedure IMTQL1, Num. Math. 12, 377-383(1968) by Martin and 
   Wilkinson, as modified in Num. Math. 15, 450(1970) by Dubrulle.  
   Handbook for Auto. Comp., vol.II-Linear Algebra, 241-248(1971).  
   See also B. T. Smith et al, Eispack Guide, Lecture Notes in 
   Computer Science, Springer-Verlag, (1976).

   The function finds the eigenvalues of a symmetric tridiagonal
   matrix by the implicit QL method.


   Arguments 
   ---------

   (input)
   n      order of the symmetric tridiagonal matrix                   
   d      contains the diagonal elements of the input matrix           
   e      contains the subdiagonal elements of the input matrix in its
          last n-1 positions.  e[0] is arbitrary	             

   (output)
   d      contains the eigenvalues in ascending order.  if an error
            exit is made, the eigenvalues are correct and ordered for
            indices 0,1,...ierr, but may not be the smallest eigenvalues.
   e      has been destroyed.					    
   ierr   set to zero for normal return, j if the j-th eigenvalue has
            not been determined after 30 iterations.		    

   Functions used
   --------------

   UTILITY	fsign
   MISC		pythag

 ***********************************************************************/

void imtqlb(long n, double d[], double e[], double bnd[])

{
   long last, l, m, i, iteration;

   /* various flags */
   long exchange, convergence, underflow;	

   double b, test, g, r, s, c, p, f;

   if (n == 1) return;
   ierr = 0;
   bnd[0] = 1.0;
   last = n - 1;
   for (i = 1; i < n; i++) {
      bnd[i] = 0.0;
      e[i-1] = e[i];
   }
   e[last] = 0.0;
   for (l = 0; l < n; l++) {
      iteration = 0;
      while (iteration <= 30) {
	 for (m = l; m < n; m++) {
	    convergence = FALSE;
	    if (m == last) break;
	    else {
	       test = fabs(d[m]) + fabs(d[m+1]);
	       if (test + fabs(e[m]) == test) convergence = TRUE;
	    }
	    if (convergence) break;
	 }
	    p = d[l]; 
	    f = bnd[l]; 
	 if (m != l) {
	    if (iteration == 30) {
	       ierr = l;
	       return;
	    }
	    iteration += 1;
	    /*........ form shift ........*/
	    g = (d[l+1] - p) / (2.0 * e[l]);
	    r = pythag(g, 1.0);
	    g = d[m] - p + e[l] / (g + fsign(r, g));
	    s = 1.0;
	    c = 1.0;
	    p = 0.0;
	    underflow = FALSE;
	    i = m - 1;
	    while (underflow == FALSE && i >= l) {
	       f = s * e[i];
	       b = c * e[i];
	       r = pythag(f, g);
	       e[i+1] = r;
	       if (r == 0.0) underflow = TRUE;
	       else {
		  s = f / r;
		  c = g / r;
		  g = d[i+1] - p;
		  r = (d[i] - g) * s + 2.0 * c * b;
		  p = s * r;
		  d[i+1] = g + p;
		  g = c * r - b;
		  f = bnd[i+1];
		  bnd[i+1] = s * bnd[i] + c * f;
		  bnd[i] = c * bnd[i] - s * f;
		  i--;
	       }
	    }       /* end while (underflow != FALSE && i >= l) */
	    /*........ recover from underflow .........*/
	    if (underflow) {
	       d[i+1] -= p;
	       e[m] = 0.0;
	    }
	    else {
	       d[l] -= p;
	       e[l] = g;
	       e[m] = 0.0;
	    }
	 } 		       		   /* end if (m != l) */
	 else {

            /* order the eigenvalues */
	    exchange = TRUE;
	    if (l != 0) {
	       i = l;
	       while (i >= 1 && exchange == TRUE) {
	          if (p < d[i-1]) {
		     d[i] = d[i-1];
		     bnd[i] = bnd[i-1];
	             i--;
	          }
	          else exchange = FALSE;
	       }
	    }
	    if (exchange) i = 0;
	    d[i] = p;
	    bnd[i] = f; 
	    iteration = 31;
	 }
      }			       /* end while (iteration <= 30) */
   }				   /* end for (l=0; l<n; l++) */
   return;
}						  /* end main */
#include <math.h>
#define		TRUE	1
#define		FALSE	0
extern long ierr;
double fsign(double, double);
double pythag(double, double);

/***********************************************************************
 *                                                                     *
 *				imtql2()			       *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   imtql2() is a translation of a Fortran version of the Algol
   procedure IMTQL2, Num. Math. 12, 377-383(1968) by Martin and 
   Wilkinson, as modified in Num. Math. 15, 450(1970) by Dubrulle.  
   Handbook for Auto. Comp., vol.II-Linear Algebra, 241-248(1971).  
   See also B. T. Smith et al, Eispack Guide, Lecture Notes in 
   Computer Science, Springer-Verlag, (1976).

   This function finds the eigenvalues and eigenvectors of a symmetric
   tridiagonal matrix by the implicit QL method.


   Arguments
   ---------

   (input)                                                             
   nm     row dimension of the symmetric tridiagonal matrix           
   n      order of the matrix                                        
   d      contains the diagonal elements of the input matrix        
   e      contains the subdiagonal elements of the input matrix in its
            last n-1 positions.  e[0] is arbitrary	             
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

void imtql2(long nm, long n, double d[], double e[], double z[])

{
   long index, nnm, j, last, l, m, i, k, iteration, convergence, underflow;
   double b, test, g, r, s, c, p, f;
   if (n == 1) return;
   ierr = 0;
   last = n - 1;
   for (i = 1; i < n; i++) e[i-1] = e[i];
   e[last] = 0.0;
   nnm = n * nm;
   for (l = 0; l < n; l++) {
      iteration = 0;

      /* look for small sub-diagonal element */
      while (iteration <= 30) {
	 for (m = l; m < n; m++) {
	    convergence = FALSE;
	    if (m == last) break;
	    else {
	       test = fabs(d[m]) + fabs(d[m+1]);
	       if (test + fabs(e[m]) == test) convergence = TRUE;
	    }
	    if (convergence) break;
	 }
	 if (m != l) {

	    /* set error -- no convergence to an eigenvalue after
	     * 30 iterations. */     
	    if (iteration == 30) {
	       ierr = l;
	       return;
	    }
	    p = d[l]; 
	    iteration += 1;

	    /* form shift */
	    g = (d[l+1] - p) / (2.0 * e[l]);
	    r = pythag(g, 1.0);
	    g = d[m] - p + e[l] / (g + fsign(r, g));
	    s = 1.0;
	    c = 1.0;
	    p = 0.0;
	    underflow = FALSE;
	    i = m - 1;
	    while (underflow == FALSE && i >= l) {
	       f = s * e[i];
	       b = c * e[i];
	       r = pythag(f, g);
	       e[i+1] = r;
	       if (r == 0.0) underflow = TRUE;
	       else {
		  s = f / r;
		  c = g / r;
		  g = d[i+1] - p;
		  r = (d[i] - g) * s + 2.0 * c * b;
		  p = s * r;
		  d[i+1] = g + p;
		  g = c * r - b;

		  /* form vector */
		  for (k = 0; k < nnm; k += n) {
		     index = k + i;
		     f = z[index+1];
		     z[index+1] = s * z[index] + c * f;
		     z[index] = c * z[index] - s * f;
		  } 
		  i--;
	       }
	    }   /* end while (underflow != FALSE && i >= l) */
	    /*........ recover from underflow .........*/
	    if (underflow) {
	       d[i+1] -= p;
	       e[m] = 0.0;
	    }
	    else {
	       d[l] -= p;
	       e[l] = g;
	       e[m] = 0.0;
	    }
	 }
	 else break;
      }		/*...... end while (iteration <= 30) .........*/
   }		/*...... end for (l=0; l<n; l++) .............*/

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
	  for (j = 0; j < nnm; j += n) {
	     p = z[j+i];
	     z[j+i] = z[j+k];
	     z[j+k] = p;
	  }
      }   
   }
   return;
}		/*...... end main ............................*/
#include <math.h>
extern double eps;

/***********************************************************************
 *                                                                     *
 *				machar()			       *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   This function is a partial translation of a Fortran-77 subroutine 
   written by W. J. Cody of Argonne National Laboratory.
   It dynamically determines the listed machine parameters of the
   floating-point arithmetic.  According to the documentation of
   the Fortran code, "the determination of the first three uses an
   extension of an algorithm due to M. Malcolm, ACM 15 (1972), 
   pp. 949-951, incorporating some, but not all, of the improvements
   suggested by M. Gentleman and S. Marovich, CACM 17 (1974), 
   pp. 276-277."  The complete Fortran version of this translation is
   documented in W. J. Cody, "Machar: a Subroutine to Dynamically 
   Determine Determine Machine Parameters," TOMS 14, December, 1988.


   Parameters reported 
   -------------------

   ibeta     the radix for the floating-point representation       
   it        the number of base ibeta digits in the floating-point
               significand					 
   irnd      0 if floating-point addition chops		      
             1 if floating-point addition rounds, but not in the 
                 ieee style					
             2 if floating-point addition rounds in the ieee style
             3 if floating-point addition chops, and there is    
                 partial underflow				
             4 if floating-point addition rounds, but not in the
                 ieee style, and there is partial underflow    
             5 if floating-point addition rounds in the ieee style,
                 and there is partial underflow                   
   machep    the largest negative integer such that              
                 1.0+float(ibeta)**machep .ne. 1.0, except that 
                 machep is bounded below by  -(it+3)          
   negeps    the largest negative integer such that          
                 1.0-float(ibeta)**negeps .ne. 1.0, except that 
                 negeps is bounded below by  -(it+3)	       

 ***********************************************************************/

void machar(long *ibeta, long *it, long *irnd, long *machep, long *negep)

{

   double beta, betain, betah, a, b, ZERO, ONE, TWO, temp, tempa, temp1;
   long i, itemp;

   ONE = (double) 1;
   TWO = ONE + ONE;
   ZERO = ONE - ONE;

   a = ONE;
   temp1 = ONE;
   while (temp1 - ONE == ZERO) {
      a = a + a;
      temp = a + ONE;
      temp1 = temp - a;
   }

   b = ONE;
   itemp = 0;
   while (itemp == 0) {
      b = b + b;
      temp = a + b;
      itemp = (long)(temp - a);
   }
   *ibeta = itemp;
   beta = (double) *ibeta;

   *it = 0;
   b = ONE;
   temp1 = ONE;
   while (temp1 - ONE == ZERO) {
      *it = *it + 1;
      b = b * beta;
      temp = b + ONE;
      temp1 = temp - b;
   }
   *irnd = 0; 
   betah = beta / TWO; 
   temp = a + betah;
   if (temp - a != ZERO) *irnd = 1;
   tempa = a + beta;
   temp = tempa + betah;
   if ((*irnd == 0) && (temp - tempa != ZERO)) *irnd = 2;


   *negep = *it + 3;
   betain = ONE / beta;
   a = ONE;
   for (i = 0; i < *negep; i++) a = a * betain;
   b = a;
   temp = ONE - a;
   while (temp-ONE == ZERO) {
      a = a * beta;
      *negep = *negep - 1;
      temp = ONE - a;
   }
   *negep = -(*negep);

   *machep = -(*it) - 3;
   a = b;
   temp = ONE + a;
   while (temp - ONE == ZERO) {
      a = a * beta;
      *machep = *machep + 1;
      temp = ONE + a;
   }
   eps = a;
   return;
}
#include <stdio.h>
#define MAXLL 2
#define STORQ 1
#define RETRQ 2
#define STORP 3
#define RETRP 4
extern double *a;
void   dcopy(long, double *, long, double *, long);

/***********************************************************************
 *                                                                     *
 *                     store()                                         *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   store() is a user-supplied function which, based on the input
   operation flag, stores to or retrieves from memory a vector.


   Arguments 
   ---------

   (input)
   n       length of vector to be stored or retrieved
   isw     operation flag:
	     isw = 1 request to store j-th Lanczos vector q(j)
	     isw = 2 request to retrieve j-th Lanczos vector q(j)
	     isw = 3 request to store q(j) for j = 0 or 1
	     isw = 4 request to retrieve q(j) for j = 0 or 1
   s	   contains the vector to be stored for a "store" request 

   (output)
   s	   contains the vector retrieved for a "retrieve" request 

   Functions used
   --------------

   BLAS		dcopy

 ***********************************************************************/

void store(long n, long isw, long j, double *s)

{
   switch(isw) {
   case STORQ:	dcopy(n, s, 1, &a[(j+MAXLL) * n], 1);
		break;
   case RETRQ:	dcopy(n, &a[(j+MAXLL) * n], 1, s, 1);
		break;
   case STORP:	if (j >= MAXLL) {
		   fprintf(stderr,"store: (STORP) j >= MAXLL \n");
		   break;
		}
		dcopy(n, s, 1, &a[j*n], 1);
		break;
   case RETRP:	if (j >= MAXLL) {
		   fprintf(stderr,"store: (RETRP) j >= MAXLL \n");
		   break;
		}
   		dcopy(n, &a[j*n], 1, s, 1);
		break;
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
 * function scales a vector by a constant.	     	      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 

void datx(long n,double da,double *dx,long incx,double *dy,long incy)

{
   long i;

   if (n <= 0 || incx == 0 || incy == 0 || da == 0.0) return;
   if (incx == 1 && incy == 1) 
      for (i=0; i < n; i++) *dy++ = da * (*dx++);

   else {
      if (incx < 0) dx += (-n+1) * incx;
      if (incy < 0) dy += (-n+1) * incy;
      for (i=0; i < n; i++) {
         *dy = da * (*dx);
         dx += incx;
         dy += incy;
      }
   }
   return;
}
/********************************************************************* 
 * Function sorts array1 and array2 into increasing order for array1 *
 *********************************************************************/

void dsort2(long igap,long n,double *array1,double *array2)

{
    double temp;
    long i, j, index;

    if (!igap) return;
    else {
	for (i = igap; i < n; i++) {
	    j = i - igap;
	    index = i;
	    while (j >= 0 && array1[j] > array1[index]) {
		temp = array1[j];
		array1[j] = array1[index];
		array1[index] = temp;
		temp = array2[j];
		array2[j] = array2[index];
		array2[index] = temp;
	        j -= igap;
		index = j + igap;
	    }
	} 
    }
    dsort2(igap/2,n,array1,array2);
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

/***************************************************************** 
 * Function finds the index of element having max. absolute value*
 * based on FORTRAN 77 routine from Linpack by J. Dongarra       *
 *****************************************************************/ 

long idamax(long n,double *dx,long incx)

{
   long ix,i,imax;
   double dtemp, dmax;

   if (n < 1) return(-1);
   if (n == 1) return(0);
   if (incx == 0) return(-1);

   if (incx < 0) ix = (-n+1) * incx;
   else ix = 0;
   imax = ix;
   dx += ix;
   dmax = fabs(*dx);
   for (i=1; i < n; i++) {
      ix += incx;
      dx += incx;
      dtemp = fabs(*dx);
      if (dtemp > dmax) {
	 dmax = dtemp;
	 imax = ix;
      }
   }
   return(imax);
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
