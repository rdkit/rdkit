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
#include "sisg.h"
#include "sisc.h"

#define  UNIX_CREAT

#ifdef UNIX_CREAT
#define PERMS 0664
#endif

void   dscal(long, double, double *, long);
void   daxpy(long, double, double *,long, double *, long);
double ddot(long, double *,long, double *, long);
void   opa(long, double *, double *);
void   opb(long, double *, double * );
void   intros(long , long , long , double * , long);
void   write_data(long , long , long , long , long ,  long , long , double , char *, char *);
float  timer(void);
void   ritzit( long , long , long , double , void (*) (long , double *, double *), void (*)(long , long , long , double *, long ),long , double **, double *, double *, double *, double *, long *);
double fsign(double, double);
double dmax(double, double);
double dmin(double, double);
long   imin(long , long );
long   imax(long, long);
#define ABS(x) ((x) > (ZERO) ? (x) : (-(x)))



/***************************************************************
 *                                                             *
 *                    main()                                   *
 *    Sparse SVD via Eigensystem of them matrix A'A            *
 *              (double precision)                             *
 *                                                             *
 ***************************************************************/
/***************************************************************

 Description
 -----------

 this sample program uses ritzit to compute singular triplets of
 the matrix A via the equivalent symmetric eigenvalue problem

 B=A'A, A is m (nrow) by n (ncol)    (nrow >> ncol), and X'=[u', v'],

 so that {u, sqrt(lambda, v} is a singular triplet of A.
 (A' = transpose of A)

 user supplied routines: opb, opa

 opa( n, x, y) takes an n-vector x and should return A*x
                      in y
 opb( n, x, y) takes an n-vector x and should return B*x
                      in y.

 User should edit timer() with an appropriate call to an
 intrinsic timing routine that returns elapsed user cpu time.


 External parameters (constants)
 -------------------------------
 Defined and documented in sisg.h (sisc.h)
 All constants are in upper-case.

 Local parameters
 ----------------
 (input)

 em2		no. of  eigenpairs of B approximated
 numextra	number of extra vectors to carry
 km		max. number of iterations
 p		initial subspace dimension (= em2 + numextra)
 em2		no. of desired eigenpairs
 eps		tolerance for eigenpair residuals
 x		two dimensional array of iteration vectors
 cx,u   	work arrays

 (output)

 d		array of eigenvalues of B (squares of the singular
                values of A)
 f		array of residual norms
 u		left-singular vectors
 x		right-singular vectors
 
 Functions used
 --------------
 BLAS         daxpy, ddot
 USER         opb, timer, intros
 MISC         write_data 
 SIS2         ritzit

 Precision
 ---------
 All floating-point calculations use double precision;
 Variables are declared as long and double

 SIS2 development
 ----------------
 SIS2 is a C translation of the Fortran-77 SIS2 from the SVDPACK
 library written by Michael W. Berry, University of Tennessee,
 Dept of Computer Science, 107 Ayres hall, Knoxville, 
 TN 37996-1301

 Date Written: 21 Jun 1992

 Sowmini Varadhan
 University of Tennessee
 Dept. of Computer Science
 107 Ayres Hall
 Knoxville, TN 37996-1301
 
 internet: varadhan@cs.utk.edu

***************************************************************/


int main(int argc,char *argv[]) 
{
   float t0, exetime;
   long k, j, i, ii, n,  nnzero, size1, size2,  vectors;
   double *x[NSIG], *cx, *f, *d, *u  , *tptr ;
   long em, numextra, km, p, em2, imem;
   double eps, dsum, tmp1, tmp2, tmp3, tmp4, xnorm;
   char title[73], name[41], v[6];
   char in1[80], in2[80], out1[80], out2[80], out3[80];
   FILE *fp_in1, *fp_in2 ;
   long fp_out3;
 
   if(argc<3){
     fprintf(stderr,"Usage: sis2-rd parmfile matfile\n");
     exit(-1);
   }
   strncpy(in1,argv[1],80); /* input parameters    */;
   strncpy(in2,argv[2],80); /* matrix              */;
   strncpy(out1,"sio2",80); /* results             */;
   strncpy(out2,"sio5",80); /* info                */;
   strncpy(out3,"siv2",80); /* singular vectors    */;

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
    if (!(fp_out2 = fopen(out2, "w"))) { 
       printf("cannot open output file %s \n", out2);
       exit(-1);
    }
	
 
   /*write header of output file*/   

   fprintf(fp_out2,"\n\n ----------------------------------------\n");
   fprintf(fp_out2," INTERMEDIATE OUTPUT PARMS:\n\n");
   fprintf(fp_out2," M:=CURRENT DEGREE OF CHEBYSHEV POLYNOMIAL\n");
   fprintf(fp_out2," S:=NEXT ITERATION STEP\n");
   fprintf(fp_out2," G:=NUMBER OF ACCEPTED EIGENVALUES OF B\n");
   fprintf(fp_out2," H:=NUMBER OF ACCEPTED EIGEN VECTORS OF B\n");
   fprintf(fp_out2," F:=VECTOR OF ERROR FOR EIGENPAIRS OF B\n");
   fprintf(fp_out2," ---------------------------------------------\n");
   fprintf(fp_out2," M         S       G       H         F     \n");
   fprintf(fp_out2," ---------------------------------------------");

   /* read data */
    fscanf (fp_in2,"%72c%*s%*s%*s%ld%ld%ld%*d",
	    title, &nrow, &ncol, &nnzero);
    title[73] = '\0';

    /* skip data format line */
    fscanf(fp_in2,"%*s %*s %*s %*s");

    n = ncol; 
 
    if (nnzero>NZMAX)
    {  fprintf(stderr,"sorry, your matrix is too big (nnzero=%d)\n",nnzero);
       exit(-1);
    }

    fscanf(fp_in1,"%s%ld%ld%ld%lf%s",name, &em, &numextra, &km, &eps, v);
 
    if (!(strcmp(v, "TRUE"))) {
        vectors = 1;

#if  !defined UNIX_CREAT
       if ((fp_out3 = open(out3, O_CREAT | O_RDWR)) == -1) {
          printf("cannot open output file %s \n", out3);
          exit(-1);
       }
#else
       if ((fp_out3 = creat(out3, PERMS )) == -1) {
          printf("cannot open output file %s \n", out3);
          exit(-1);
       }
#endif


     }
    else vectors = 0;

    p= em + numextra;
    em2=em;
    /*******************************************************************
     * allocate memory						       *
     * pointr - column start array of harwell-boeing sparse matrix     *
     *          format                                       (ncol+1)  *
     * rowind - row indices array of harwell-boeing format   (nnzero)  *
     * value  - nonzero values array of harwell-boeing sparse matrix   *
     *          format                                       (nnzero)  *
     * x      - 2 dimensional array of  iteration vectors    (NSIG*n)  *
     * d      - array of eigenvalues of B                    (p)       *
                         (squares of the singular values of A)         *
     * f      - temporary storage array                      (p)       *
     * cx     - temporary storage array (control quantities) (n)       *
     * u      - work array                                   (nrow)    *
     *******************************************************************/
	 size1 = sizeof(double) * (nnzero + p* 2 + n  + nrow + NSIG * n);
	 size2 = sizeof(long) * (ncol + nnzero + 1);
	 if (!(pointr = (long *)   malloc(size2))   ||
	     !(value  = (double *) malloc(size1))){
	     perror("MALLOC FAILED in MAIN()");
	     exit(errno);
         }

         imem=size1 - (nnzero * sizeof(double));

         rowind=pointr + (ncol+1); 

         tptr=value + nnzero;
         for (ii=0;ii<NSIG; ii++){
          x[ii]=tptr;
	  tptr+=n;
         }
         d=tptr;
         tptr+=p;
         f=tptr;
         tptr+=p;
         cx=tptr;
         tptr+=n;
         u  =tptr;
	 
    /* read data */
    for (i = 0; i <= ncol; i++) fscanf(fp_in2, "%ld", &pointr[i]);
    
    for (i = 0; i < ncol; i++) pointr[i] -= 1;
 
    /* define last element of pointr in case it is not */ 
    pointr[i] = nnzero;
    for (i = 0; i < nnzero; i++) fscanf(fp_in2, "%ld", &rowind[i]);
    for (i = 0; i < nnzero; i++) rowind[i] -= 1;
    for (i = 0; i < nnzero; i++) fscanf(fp_in2, "%lf", &value[i]);


    exetime = timer();
    
    /*  call ritzit */
    ritzit(n,p,km,eps,opb,intros, em,x, d, f, cx, u, &imem);

    exetime = timer() - exetime;
    write_data(n,km,em2,em,p,imem,vectors,eps,title,name);
    t0=timer();

    for (j=0;j<em2;j++) {
      opb(n, &x[j][0], cx);
      tmp1=ddot(n,&(x[j][0]), 1, cx, 1);
      daxpy(n, -tmp1, &(x[j][0]), 1, cx, 1);
      tmp1=sqrt(tmp1);
      xnorm=sqrt(ddot(n, cx, 1, cx, 1));
      /* multiply by matrix A to get (scaled) left S-vector */
      opa(n, &x[j][0], u);
      tmp2=ONE/tmp1;
      dscal(nrow, tmp2, u, 1);
      xnorm=xnorm*tmp2;
      f[j]=xnorm;
      d[j]=tmp1;
      if (vectors) {
           write(fp_out3, (char *)u, nrow * sizeof(double));
           write(fp_out3, (char *)&x[j][0],n * sizeof(double));
     }
    }

   exetime=exetime + timer() -t0;

   fprintf(fp_out1,"\n ...... SISVD EXECUTION TIME= %10.2e\n",exetime);
   fprintf(fp_out1," ......\n ......");
   fprintf(fp_out1," MXV   = %12.ld\n ......    COMPUTED SINGULAR VALUES  (RESIDUAL NORMS)\n ......", mxvcount);

   for (ii=0;ii<em2;ii++)
      fprintf(fp_out1,"\n ...... %3.ld   %22.14e  (%11.2e)",
                                   ii+1,     d[ii],    f[ii]);

     
   fprintf(fp_out1,"\n");
   free(value);
   free(pointr);
   fclose(fp_in1);
   fclose(fp_in2);
   if (vectors) close(fp_out3);
   fclose(fp_out2);
   fclose(fp_out1);

}

/*********************************************************************/
void write_data(long n, long km, long em2, long em, long p,  long imem, long vectors, double eps, char *title, char *name)
{
 char svectors;
if (vectors) svectors='T'; else svectors='F';
  fprintf( fp_out1,"\n ...\n ... SOLVE THE [A^TA] EIGENPROBLEM\n");
  fprintf(fp_out1," ... NO. OF EQUATIONS           = %10.ld\n",n);
  fprintf( fp_out1," ... MAX. NO. OF ITERATIONS     = %10.ld\n ",km);
  fprintf(fp_out1,"... NO. OF DESIRED EIGENPAIRS  = %10.ld\n ",em2);
  fprintf(fp_out1,"... NO. OF APPROX. EIGENPAIRS  = %10.ld\n ",em);
  fprintf(fp_out1,"... INITIAL SUBSPACE DIM.      = %10.ld\n ",p);
  fprintf(fp_out1,"... FINAL   SUBSPACE DIM.      = %10.ld\n ",p-em);
  fprintf(fp_out1,"... MEM REQUIRED (BYTES)       = %10.ld\n",imem);
  fprintf(fp_out1," ... WANT S-VECTORS? [T/F]      = %10.c\n",svectors);
  fprintf(fp_out1," ... TOLERANCE                  = %10.2e\n ",eps);
  fprintf(fp_out1,"... NO. OF ITERATIONS TAKEN    = %10.ld\n ", ksteps);
  fprintf(fp_out1,"... MAX. DEG. CHEBYSHEV POLY.  = %10.ld\n ...\n ", maxm);
  fprintf(fp_out1,"%s\n %s\n ... NO. OF TERMS     (ROWS)    = %10.d\n",title, name, nrow);
  fprintf(fp_out1," ... NO. OF DOCUMENTS (COLS)    = %10.ld\n ",ncol);
  fprintf(fp_out1,"... ORDER OF MATRIX B          = %10.ld\n ...\n", n)
;

fflush(fp_out1);
}
/***********************************************************************/

void intros(long ks, long kg, long kh, double *f, long m)

{
  long l, i, j;

  l=kh+1;

  ksteps=ks;
  maxm=imax(maxm,m);

  fprintf(fp_out2,"\n%8.ld%8.ld%8.ld%8.ld%15.6e\n",m,ks,kg+1,kh+1,f[0]);
 
 if (l>=1)
  for (i=1;i<=l;i++){
   for (j=0;j<32;j++) fprintf(fp_out2," ");
   fprintf(fp_out2,"%15.6e\n",f[i]);
  }
}
  
#include <math.h>
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>

extern long mxvcount;
double ddot(long , double *, long, double *, long);
void dscal(long, double, double *, long);
void tred2(long, long, double **, double *, double *, double **);
long tql2(long, long, double *, double *, double **);
void daxpy(long, double, double *, long, double *, long);
double dmax(double, double);
long imax(long, long);

#include "sisc.h"

/****************************************************************
 *                                                              *
 *                        ritzit()                              *
 *                                                              *
 ****************************************************************/
/****************************************************************

 Description:
 ------------

 This subroutine is a translation of the Fortran-77 procedure
 RITZIT from the SVDPACK library written by Michael W. Berry,
 University of Tennessee, Dept. of Computer Science, 107 Ayres
 Hall, Knoxville TN 37919-1301

 This subroutine determines the absolutely largest eigenvalues
 and corresponding eigenvectors of a real symmetric matrix by
 simultaneous iteration.

 External parameters (constants)
 -------------------------------

 Defined and documented in sisg.h (sis1c.h)

 local parameters: 
 -----------------
 
 (input)

 n      the order of the matrix  whose eigenpairs are sought
        (matrix B for the SVD problem)
 kp     the number of simultaneous iteration vectors
 km     the maximum number of iteration steps to be performed. If
        starting values for the iteration vectors are available,
        km should be prefixed with a minus sign.
 eps    the tolerance for accepting eigenvectors
 opb    the name of the subroutine that defines the matrix B. opb
        is called with the parameters (n, u, w) and must 
        compute w=Bu without altering u for the SVD problem
 inf    the name of the subroutine that may be used to obtain 
        information or exert control during execution. inf is
        called with parameters (ks, kg, kh, f, m), where ks
        is the number of the next iteration step, kg is the
        number of already accepted eigenvectors, kh is the number
        already accepted eigenvalues, and f is the array of error
        quantities for the vectors of x. An element of f has the
        value 4.0 until the corresponding eigenvalue of the matrix
        B  has been accepted.
 kem    the number of eigenvalues and corresponding eigenvectors
         of matrix B desired. kem must be less than kp.
 x      contains, if km is negative, the starting values for the
        iteration vectors.

 (output)
 
 km     is unchanged
 kem    is reset to the number of eigenvalues and eigenvectors
        of matrix B actually accepted within the limit of km steps
 imem   the number of bytes needed for this invocation
 x      contains in its first kem columns orthonormal
        eigenvectors of the matrix B,  corresponding to the 
        eigenvalues in array d.  The remaining columns contain 
        approximations to further eigenvectors.
 d      contains in its first kem positions the absolutely
        largest eigenvalues of the matrix B. The remaining positions
        contain approximations to smaller eigenvalues.

 u, w, b, f, cx, x, s and  rq are temporary storage arrays.

 Functions used
 --------------

 BLAS:  ddot, dscal, daxpy
 EISP:  tred2, tql2, pythag
 MISC:  opb, inf, dmax, imax
 
*****************************************************************/
void ritzit( long n, long kp, long km, double eps, void (*opb) (long , double *, double * ), void (*inf)(long , long , long , double *, long ),long kem, double **x, double *d, double *f, double *cx, double *u, long *imem)

{
  long i, j, l, i1, l1, size, kg, kh, kz, kz1, kz2, ks, m, m1;
  long  jj, ii, ig, ip, ik, jp, k, flag1 ;  
  double  *w[NMAX], *rq, *b[NSIG], *s, *tptr, TmpRes;
  double xks, ee2, e, e2, e1, xkm, xk, xm1, t;

  /*******************************************************************
   * allocate memory for temporary storage arrays used in ritzit only*
   *                                                                 *
   *            rq - (kp)                                            *
   *             w - (NSIG * n)                                      *
   *             b - (NSIG * kp)                                     *
   *             s - (kp)                                            *
   *******************************************************************/

   size =    kp * (NSIG +2) + NSIG * n;
   if (!(rq=(double *) malloc(size * sizeof(double)))){
      perror("MALLOC FAILED IN RITZIT()");
      fprintf(stderr,"\t\t%d\n",size*sizeof(double));
      exit(errno);
   }
   tptr= rq + kp;
   for (i=0; i< NSIG; i++){
    w[i]= tptr;
    tptr+=n;
   }

   /*w=tptr; tptr+=(n * kp);*/
   for (i=0; i<NSIG; i++) {
    b[i]=tptr;
    tptr+=kp;
   };
   s=tptr;

   /*finished allocating memory*/


  *imem += size;
  *imem = sizeof(double) * (*imem);

  ee2 = ONE + 1.0e-1 * eps;
  e   = ZERO;
  kg  = -1;
  kh  = -1;
  kz  = 1367;
  kz1 = 0;
  kz2 = 0;
  ks  = 0;
  m   = 1;

  for (l=0; l<kp; l++){
    f[l] = 4.0e+0;
   cx[l] = ZERO;
   rq[l] = ZERO;
 } 
 
  if (km >= 0)    /*generate random initial iteration vectors*/
      for (j=0; j<kp; j++)
       for (l=0; l<n; l++){
         kz = (3125 * kz) % 65536;  
        x[j][l] = (double) (kz - 32768);
      }
  km =abs(km);
  l1 = 1;
  i1 = 1;
  jp=kp;
flag1=0;
#define NRow kp
#define NCol n

  /* extend orthonormalization to all kp rows of x 

     Variables used here:
 
     NCol                     : number of columns of x
     NRow                     : number of rows of x
     ii, ik, jp, k, l1        : indexing integers
     t, TmpRes                : double precision storage.


     at the end of this section of code, 

     x   : transpose(Q)
     b   : R

     where X = QR is the orthogonal factorization of the matrix X
     stored in array x

     invokes: ddot, daxpy, dscal (BLAS)
              sqrt               (math.h)
  */

    ik = l1 -1; 
  for (i=0; i< jp; i++){
    for (k=l1-1; k<jp; k++) b[i][k]= ZERO;
   
    if (i >= l1-1){
      ik=i+1;
      b[i][i]=sqrt(ddot(NCol, &x[i][0], 1 ,  &x[i][0], 1));
      t=ZERO;
      if (b[i][i] != ZERO) t=ONE / b[i][i];
      dscal(NCol, t, &(x[i][0]),  1);
     }

   for (ii=ik; ii< NRow; ii++){
     TmpRes=ZERO;
     for (jj=0; jj<NCol; jj++)
       TmpRes+= x[ii][jj] * x[i][jj];
     b[i][ii]=TmpRes;
  }
  
    for (k=ik; k<jp;k++) 
        daxpy(NCol, -b[i][k], &x[i][0], 1, &x[k][0], 1);

  } /* end for i */

  ig = 0;
  ip = kp - 1;
  
 /*statement 70  from original RITZIT begins here*/
while (1){

  flag1=1;  /* so that jacobi step is done at least once */

  while (flag1){
        /* jacobi step modified */
        for(k=ig; k< kp; k++){
          opb(n, &x[k][0], &w[0][0]);
          for (j=0; j<n; j++) x[k][j]=w[0][j];
         }
        l1=ig + 1;
        jp=kp;
        flag1=0; /*flag1 is set only if re-orthog needs to be done*/

     /* extend orthonormalization to all kp rows of x 

     Variables used here:
 
     NCol                     : number of columns of x
     NRow                     : number of rows of x
     ii, ik, jp, k, l1        : indexing integers
     t, TmpRes                : double precision storage.


     at the end of this section of code, 

     x   : transpose(Q)
     b   : R

     where X = QR is the orthogonal factorization of the matrix X
     stored in array x

     invokes: ddot, daxpy, dscal (BLAS)
              sqrt               (math.h)
  */

    ik = l1 -1; 
  for (i=0; i< jp; i++){
    for (k=l1-1; k<jp; k++) b[i][k]= ZERO;
   
    if (i >= l1-1){
      ik=i+1;
      b[i][i]=sqrt(ddot(NCol, &x[i][0], 1 ,  &x[i][0], 1));
      t=ZERO;
      if (b[i][i] != ZERO) t=ONE / b[i][i];
      dscal(NCol, t, &(x[i][0]),  1);
     }

   for (ii=ik; ii< NRow; ii++){
     TmpRes=ZERO;
     for (jj=0; jj<NCol; jj++)
       TmpRes+= x[ii][jj] * x[i][jj];
     b[i][ii]=TmpRes;
  }
  
    for (k=ik; k<jp;k++) 
        daxpy(NCol, -b[i][k], &x[i][0], 1, &x[k][0], 1);

  } /* end for i */

       if ( ks<=0 ) {
         /* measures against unhappy choice of initial vectors*/
         for (k=0;k<kp; k++)
           if (b[k][k] ==ZERO)
              for (j=0;j<n;j++){
                 kz= (3125 * kz) % 65536;
                 x[k][j] = (double)(kz-32768);
                 flag1=1;
               }
       }
     if (flag1){
       l1=1;
       ks=1; /*we dont want to re-initialize x[][] again*/

   /* extend orthonormalization to all kp rows of x 

     Variables used here:
 
     NCol                     : number of columns of x
     NRow                     : number of rows of x
     ii, ik, jp, k, l1        : indexing integers
     t, TmpRes                : double precision storage.


     at the end of this section of code, 

     x   : transpose(Q)
     b   : R

     where X = QR is the orthogonal factorization of the matrix X
     stored in array x

     invokes: ddot, daxpy, dscal (BLAS)
              sqrt               (math.h)
  */

    ik = l1 -1; 
  for (i=0; i< jp; i++){
    for (k=l1-1; k<jp; k++) b[i][k]= ZERO;
   
    if (i >= l1-1){
      ik=i+1;
      b[i][i]=sqrt(ddot(NCol, &x[i][0], 1 ,  &x[i][0], 1));
      t=ZERO;
      if (b[i][i] != ZERO) t=ONE / b[i][i];
      dscal(NCol, t, &(x[i][0]),  1);
     }

   for (ii=ik; ii< NRow; ii++){
     TmpRes=ZERO;
     for (jj=0; jj<NCol; jj++)
       TmpRes+= x[ii][jj] * x[i][jj];
     b[i][ii]=TmpRes;
  }
  
    for (k=ik; k<jp;k++) 
        daxpy(NCol, -b[i][k], &x[i][0], 1, &x[k][0], 1);

  } /* end for i */
      } 
    } /*end while flag1 */

  for (k= ig; k<kp; k++)
    for (l=k; l<kp; l++){
      t = ZERO;

      for (i=l; i<kp; i++)
        t+= b[k][i] * b[l][i];

      /* negate matrix to reverse eigenvalue ordering */
      b[k][l] = -t;
    }
  j=kp - kg - 1;

  tred2(ig, kp, b, d, u, b);
  ii=tql2(ig, kp, d, u, b);


  for (k=ig; k< kp; k++)
   d[k]=sqrt(dmax(-d[k], ZERO));
  

 for (j=ig; j<kp; j++)
  for (k=0; k<n; k++) 
   w[j][k]=ZERO;

   for (j=ig; j<kp; j++)
     for (l=ig; l<kp; l++){
       TmpRes=b[j][l];
       for (k=0; k<n; k++)
         w[j][k] += TmpRes  * x[l][k];}   /* w is going to be transposed as 
                                         compared to the fortran version */
     

  for (j=ig; j<kp; j++)
   for (k=0; k<n; k++)
      x[j][k]=w[j][k]; 

  xks=(double)(++ks);
  if (d[kp-1] > e) e=d[kp-1 ];

  /* randomization */
  if (kz1<3) {
    for (j=0; j<n; j++) {
       kz= (3125 * kz) % 65536;
       x[kp-1][j] = (double) (kz -32768);
    }
    l1=kp;
    jp=kp;

  /* extend orthonormalization to all kp rows of x 

     Variables used here:
 
     NCol                     : number of columns of x
     NRow                     : number of rows of x
     ii, ik, jp, k, l1        : indexing integers
     t, TmpRes                : double precision storage.


     at the end of this section of code, 

     x   : transpose(Q)
     b   : R

     where X = QR is the orthogonal factorization of the matrix X
     stored in array x

     invokes: ddot, daxpy, dscal (BLAS)
              sqrt               (math.h)
  */

    ik = l1 -1; 
  for (i=0; i< jp; i++){
    for (k=l1-1; k<jp; k++) b[i][k]= ZERO;
   
    if (i >= l1-1){
      ik=i+1;
      b[i][i]=sqrt(ddot(NCol, &x[i][0], 1 ,  &x[i][0], 1));
      t=ZERO;
      if (b[i][i] != ZERO) t=ONE / b[i][i];
      dscal(NCol, t, &(x[i][0]),  1);
     }

   for (ii=ik; ii< NRow; ii++){
     TmpRes=ZERO;
     for (jj=0; jj<NCol; jj++)
       TmpRes+= x[ii][jj] * x[i][jj];
     b[i][ii]=TmpRes;
  }
  
    for (k=ik; k<jp;k++) 
        daxpy(NCol, -b[i][k], &x[i][0], 1, &x[k][0], 1);

  } /* end for i */

  }
  
  /*compute control quantities cx */
  for (k=ig; k<ip;k++){
    t=(d[k] - e)*(d[k] + e);
    if (t <= ZERO) cx[k]=ZERO; 
    else 
           if (e ==ZERO) cx[k] = 1.0e+3 + log(d[k]);
           else  cx[k]= log( (d[k]+sqrt(t)) / e);
  }
    
    /*acceptance test for eigenvalues including adjustment of 
      kem and kh such that 
           d[kem] > e, 
           d[kh]  > e, and,
           d[kem] does not oscillate strongly */

  for (k=ig; k<kem; k++)
    if ((d[k] <= e) ||
	((kz1>1) && (d[k] <= 9.99e-1 * rq[k]))){
      kem=k-1;
      break;
    }

  if (!kem) break; 

  k=kh+1;

  while ((d[k] != ZERO) && ( d[k] <= ee2 * rq[k])){
    kh=k++;
  }

  do {
    if (d[k] <= e) kh= k- 1;
    --k;
  }
  while (k > kem-1);

  /*acceptance test for eigenvectors */
  l= kg;
  e2= ZERO;
 
  for (k=ig; k<ip; k++){
    if (k == l+1){
      /*check for nested eigenvalues */
      l=k;
      l1=k;
      if (k != ip){
	ik= k + 1;
	s[0] = 5.0e-1 / xks;
	t= ONE / (double)(ks * m);
	for (j=ik; j<ip; j++){
	  if ((cx[j] * (cx[j] + s[0]) + t) <= (cx[j-1] * cx[j-1]))
	    break;
	  l=j;
	}
      }
    }
    if (l > kh) {l = l1; break;}
  
    opb(n, &x[k][0], &w[0][0]);
    s[0]= ZERO;
      
    for (j=0; j<=l; j++)
      if (fabs(d[j] - d[k]) < 1.0e-2 * d[k]){
	t=ZERO;
 
	for (i=0; i<n; i++) 
	  t+= w[0][i] * x[j][i];

	for (i=0; i<n; i++) {
	  w[0][i] =w[0][i] - t * x[j][i];
	}

	s[0]+= t*t;
      }
  
    t=ZERO;

    for (i=0; i<n; i++) 
      t+= w[0][i] * w[0][i];

    if (s[0] != ZERO) t=sqrt(t/(s[0] + t)); 
    else t=ONE;
   
    if (t > e2) e2=t;

    if (k==l){
      /* test for acceptance of group of eigenvectors*/
      if ((l>=kem-1) &&
	  (d[kem] * f[kem] < eps * (d[kem] -e)))
	kg=kem;

      if (e2 < f[l])
	for (j=l1; j<=l; j++) f[j]=e2;

      if ((l <= kem) &&
	  (d[l] * f[l] < eps * (d[l]-e)))
	kg=l;

      ig=kg+1;
    }
  }

  /* statements 570 to 660  from original RITZIT*/

   if (e<= 4.0e-2 * d[0]) {
     m=1;
     k=1;
   }
   else {
     e2=2.0e0/e;
     e1=5.1e-1 * e2;
     k= 2 * imax((long)(4.0e0/cx[0]), 1);
     if (m>k) m=k;
   }
   /*reduce kem if convergence would be too slow */
   xkm=(double)km;
   if ((f[kem-1] != ZERO) &&
       (xks < 9.0e-1 *xkm)){
         xk=(double)k;
         s[0]=xk * cx[kem-1];
         if (s[0] < 5.0e-2) t=5.0e-1 * s[0] *cx[kem-1];
         else t=cx[kem-1] + log(5.0e-1 * (ONE + exp(-2.0e0 * s[0])))/xk;
         s[0]=log(d[kem-1] * f[kem-1]/(eps * (d[kem-1] -e)))/t;
         if ((xkm - xks) * xkm < s[0] * xks) kem--;
   }
   inf(ks,kg,kh,f,m);

  if ((kg >= kem-1) ||
      (ks >= km))      break;
 
  for (k=ig; k<kp; k++)  rq[k] = d[k];

  do {
      /*statements 680-700 */
      if (ks + m > km) {
           kz2=-1;
           if (m >1) m = 2* ((km -ks +1)/2);
      }
      else
          m1=m;

      /*shortcut last intermediate block if all error quantities f are
        sufficiently small */
      if (l >= kem){
        s[0]= d[kem-1] * f[kem-1]/(eps *(d[kem-1] -e));
        t= s[0] * s[0] - ONE;
        if (t <= ZERO) break;
        s[0] = log(s[0] + sqrt(t))/(cx[kem-1] -cx[kh+1]);
        m1=2 * (long)(5.0e-1 * s[0] + 1.01e0);
        if (m1<=m) kz2=-1;
        else m1=m;
     }
     xm1=(double) m1;
     
                                   /*chebyshev iteration */
     if (m==1)
       for (k=ig; k<kp; k++){
       opb(n, &x[k][0], &w[0][0]);
       for (i=0; i<n; i++) 
           x[k][i] = w[0][i];
       }
     else                                /*degree != ONE */
       for (k=ig; k<kp; k++){
          opb(n, &x[k][0], &w[0][0]);

          for (i=0; i<n; i++) u[i]=e1 * w[0][i];

          opb(n, u, &w[0][0]);

          for (i=0; i<n; i++) x[k][i]= e2 *w[0][i] - x[k][i];

          if (m1>=4)
               for (j=3; j<m1; j+=2){
                 opb(n, &x[k][0], &w[0][0]);
 
                 for (i=0; i<n; i++)
                    u[i]=e2 * w[0][i] - u[i];

                 opb(n, u, &w[0][0]);

                 for (i=0; i<n; i++)
                     x[k][i] = e2 * w[0][i] - x[k][i];
                }
       }

       
       l1=ig+1;

    /* extend orthonormalization to all kp rows of x

     Variables used here:
 
     NCol                     : number of columns of x
     NRow                     : number of rows of x
     ii, ik, jp, k, l1        : indexing integers
     t, TmpRes                : double precision storage.


     at the end of this section of code, 

     x   : transpose(Q)
     b   : R

     where X = QR is the orthogonal factorization of the matrix X
     stored in array x

     invokes: ddot, daxpy, dscal (BLAS)
              sqrt               (math.h)
  */

    ik = l1 -1; 
  for (i=0; i< jp; i++){
    for (k=l1-1; k<jp; k++) b[i][k]= ZERO;
   
    if (i >= l1-1){
      ik=i+1;
      b[i][i]=sqrt(ddot(NCol, &x[i][0], 1 ,  &x[i][0], 1));
      t=ZERO;
      if (b[i][i] != ZERO) t=ONE / b[i][i];
      dscal(NCol, t, &(x[i][0]),  1);
     }

   for (ii=ik; ii< NRow; ii++){
     TmpRes=ZERO;
     for (jj=0; jj<NCol; jj++)
       TmpRes+= x[ii][jj] * x[i][jj];
     b[i][ii]=TmpRes;
  }
  
    for (k=ik; k<jp;k++) 
        daxpy(NCol, -b[i][k], &x[i][0], 1, &x[k][0], 1);

  } /* end for i */

      /*discounting the error quantities F */
      if (kg!=kh){
         if (m ==1) 
              for (k=ig; k<=kh; k++) f[k]= f[k] * (d[kh+1]/d[k]);
         else {
               t=exp(-xm1 * cx[kh+1]);
               for (k=ig; k<=kh; k++){
                   s[0]=exp(-xm1 * (cx[k]-cx[kh+1]));
                   f[k] = s[0] * f[k] * (ONE + t*t)/(ONE + (s[0] *t)*(s[0]*t));
               }
         } 
      }
      ks=ks+m1;
      kz2=kz2 - m1;
     
   } /*possible repetition of intermediate steps*/
   while (kz2>=0);

   kz1++;
   kz2 = 2 * kz1;
   m = 2 * m; 

}  /* end while kz2<=0 ... go to 70*/

/*statement 900 of original RITZIT program begins here*/

kem = kg;
l1=1; 
jp=kp-1;

  /* extend orthonormalization to all kp rows of x 

     Variables used here:
 
     NCol                     : number of columns of x
     NRow                     : number of rows of x
     ii, ik, jp, k, l1        : indexing integers
     t, TmpRes                : double precision storage.


     at the end of this section of code, 

     x   : transpose(Q)
     b   : R

     where X = QR is the orthogonal factorization of the matrix X
     stored in array x

     invokes: ddot, daxpy, dscal (BLAS)
              sqrt               (math.h)
  */

    ik = l1 -1; 
  for (i=0; i< jp; i++){
    for (k=l1-1; k<jp; k++) b[i][k]= ZERO;
   
    if (i >= l1-1){
      ik=i+1;
      b[i][i]=sqrt(ddot(NCol, &x[i][0], 1 ,  &x[i][0], 1));
      t=ZERO;
      if (b[i][i] != ZERO) t=ONE / b[i][i];
      dscal(NCol, t, &(x[i][0]),  1);
     }

   for (ii=ik; ii< NRow; ii++){
     TmpRes=ZERO;
     for (jj=0; jj<NCol; jj++)
       TmpRes+= x[ii][jj] * x[i][jj];
     b[i][ii]=TmpRes;
  }
  
    for (k=ik; k<jp;k++) 
        daxpy(NCol, -b[i][k], &x[i][0], 1, &x[k][0], 1);

  } /* end for i */

/*statements 920 up to last card of the original RITZIT program */
 
for (k=0; k<ip; k++)
  for (i=k; i<ip; i++) b[k][i]=ZERO;

for (k=0; k<ip; k++) {
   opb(n, &x[k][0], &w[0][0]);

   for (i=0; i<=k; i++)
        for (l=0; l<n; l++)
          b[i][k]=b[i][k] - x[i][l] * w[0][l];
   }

tred2(0, ip, b, d, u, b);
tql2(0, ip, d, u, b);

/*reordering of eigenvalues and eigenvectors according to the magnitudes 
  of the former */

for (i=0; i< ip; i++)
  if (i!=ip-1){ 
       k=i;
       t=d[i];
       ii=i+1;
       for (j=ii; j<ip; j++)
         if (fabs(d[j]) > fabs(t)){
             k=j; 
             t=d[j];
         }
       if (k!=i) {
          d[k] = d[i];
          d[i] = t;

          for (j=0; j<ip; j++){
            s[j]= b[i][j];
            b[i][j] = b[k][j];
            b[k][j] = s[j];
          }
      }
   d[i]=-d[i];
   }

  for (i=0; i<ip; i++)
   for(j=0; j<n; j++)
    w[i][j]=ZERO;

  for (i=0; i<ip; i++)
   for (k=0; k<ip; k++){
    TmpRes=b[i][k];
    for (j=0; j<n; j++)
      w[i][j]+=TmpRes * x[k][j];
   }

     for (i=0; i<ip; i++)
        for (j=0; j<n; j++)
           x[i][j] = w[i][j];

     d[kp-1]=e;

  /* free memory at the end of ritzit*/
  free(rq);
  return;

}  /*end ritzit*/

#include "sisc.h"

extern long ncol, nrow, mxvcount;
extern long *pointr, *rowind;
extern double *value ;

/**************************************************************
 * multiplication of matrix A by a vector x, where            *
 *                                                            *
 * A is nrow by ncol (nrow >> ncol)                           *
 * y stores product vector                                    *
 **************************************************************/
void opa(long n,double *x, double *y)

{
   long i, j, end;

   mxvcount += 1;
   for (i = 0; i < nrow; i++) y[i] = ZERO;

/* multiply by sparse C */
   for (i = 0; i < ncol; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++)
	y[rowind[j]] += value[j] * x[i];
   }

   return;
}
#include "sisc.h"

extern long ncol, nrow, mxvcount;
extern long *pointr, *rowind;
extern double *value ;

/************************************************************** 
 * multiplication of matrix B by a vector x, where            *
 *							      *
 * B =  A'A, where A is nrow by ncol (nrow >> ncol)           *
 * Hence, B is of order n:=ncol                               *
 * y stores product vector                                    *
 **************************************************************/ 
void opb(long n,double *x, double *y)

{
   long i, j, end;
   double *ztemp;
 

  ztemp=(double *) malloc(nrow * sizeof(double));
   mxvcount += 2;

   for (i = 0; i < ncol; i++) y[i] = ZERO;

   for (i=0; i<nrow; i++) ztemp[i] = ZERO;

/* multiply by sparse C */
   for (i = 0; i < ncol; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++)
	ztemp[rowind[j]] += value[j] * x[i];
   }

/*multiply by sparse C' */
   for (i=0;i<ncol; i++){
      end=pointr[i+1];
      for (j=pointr[i]; j< end; j++)
         y[i] += value[j] * ztemp[rowind[j]];
   }

   free(ztemp);
   return;
}
#include <math.h>
double fsign(double, double);
/***********************************************************************
 *                                                                     *
 *                              tred2()                                *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

  Description
  -----------
  
  tred2() is a translation of the algol procedure TRED2, Num. Math. 11, 
  181-195 (1968) by Martin, Reinsch, and Wikinson.  Handbook for Auto.
  Comp., Vol. II- Linear Algebra, 212-226 (1971)

  This subroutine reduces a real symmetric matrix to a symmetric
  tridiagonal matrix using and accumulating orthogonal similarity
  transformations.

  Arguments
  ---------

  (input)
  offset index of the leading element of the matrix to be
         tridiagonalized. The matrix tridiagonalized should be 
         stored in a[offset:n-1, offset:n-1]

  n	 order of the matrix

  a	 contains the real symmetric input matrix. Only the upper
	 triangle of the matrix need be supplied

  (output)
  d	 contains the diagonal elements of the tridiagonal matrix.
  
  e	 contains the subdiagonal elements of the tridiagonal matrix
	 in its first n-1 positions.

  z	 contains the orthogonal transformation matrix produced in the
	 reduction.

  a and z may coincide. If distinct, a is unaltered.

  Functions used:
  UTILITY: fsign

***********************************************************************/

void tred2(long offset, long n, double **a, double *d, double *e, double **z)
{
 long jj,ii,i,j,k,l, jp1;
 double h, scale, f, g,  hh, tmp;

 

 for (i=offset;i<n;i++) 
  { 
   for (j=i;j<n;j++)
     {
      z[j][i]=a[i][j];   /*fix this later?.. the rest of the routine 
                           assumes that z has the lower triangular part
                           of the symmetric matrix */
     }
   d[i]=a[i][n-1];
  }


  if (n==1) 
   {
    for (i=offset;i<n;i++)
     {
       d[i]=z[n-1][i];
       z[n-1][i]=ZERO;
     }
    z[n-1][n-1]=ONE;
    e[1]=ZERO;
    return;
   }

  /*for i = n step -1 until 2 do*/

  for (ii=3;ii<n+2-offset;ii++)
   {
     i=n+2-ii;
     l=i-1;
     h=ZERO; 
     scale=ZERO;

    /*scale row (algol tol then not needed)*/
     if (l>=1)
       for (k=offset;k<=l;k++)
        {
         scale+= fabs(d[k]);
        }
	
    if ((scale==ZERO)||(l<1))
     {
      e[i]=d[l];
      for (j=offset;j<=l;j++)
          {
            d[j]=z[l][j];
            z[i][j]=ZERO;
            z[j][i]=ZERO;
          }
     }
   else                   /*scale <> ZERO */
     {
       for (k=offset;k<=l;k++)
        {
         d[k]=d[k]/scale;
         h+=d[k]*d[k];
        }


       f=d[l];
       g=-fsign(sqrt(h), f);
       e[i]=scale * g;
       h-=f*g;
       d[l]=f-g;
   
       /* form A*u */
  
       for (j=offset; j<=l; j++)
          e[j]=ZERO;
          
          for (j=offset;j<=l;j++)
            {
             f=d[j];
             z[j][i]=f;
             g= e[j] + z[j][j] * f;
             
             jp1= j + 1;
   
             if (l >= jp1) 
                 {
                  for (k=jp1; k<=l; k++)
                   {
                     g+= z[k][j] * d[k];
                     e[k] += z[k][j] * f;
                   }
                 };
             e[j]=g;
           }

       /* form P */
 
       f= ZERO;
 
       for (j=offset; j<=l; j++)
        {
          e[j]=e[j]/h;
          f+= e[j] * d[j];
        }

       hh= f/ (h+h);
  
       /* form Q */
  
      for (j=offset; j<=l; j++)
       e[j] -= hh * d[j];

      /* form reduced A */

      for (j=offset; j<=l; j++)
       {
         f= d[j];
         g = e[j];

         for (k=j; k<=l; k++)
          z[k][j]= z[k][j] - f * e[k] - g * d[k];

         d[j]=z[l][j];
         z[i][j]=ZERO;
       }
    }  /* end scale <> zero */

    d[i]=h;
   }   /* end for ii */
   /*accumulation of transformation matrices */

   for (i=offset + 1;i<n;i++)
    {
     l=i-1;
     z[n-1][l] = z[l][l];
     z[l][l] = ONE;
     h=d[i];

     if (h != ZERO) 
       {
        for (k=offset; k<=l; k++)
          d[k]= z[k][i]/h;

        for (j=offset; j<=l; j++)
         {
           g= ZERO;
           
           for (k=offset;k<=l; k++)
            g+= z[k][i]*z[k][j];

	   for (k=offset;k<=l;k++)
            z[k][j] -= g * d[k];
         }
       }
       for (k=offset;k<=l;k++) z[k][i]=ZERO;
     }
  
     for (i=offset;i<n;i++)
       {
        d[i]=z[n-1][i];
        z[n-1][i]=ZERO;
       }
     z[n-1][n-1]=ONE;
     e[0]=ZERO;

/*preparation for tql2.c.. reorder e[]*/
for (i=1+offset;i<n;i++) e[i-1]=e[i]; 

/*preparation for tql2.c.. z has to be transposed for 
  tql2 to give correct eigenvectors */
for (ii=offset; ii<n; ii++)
 for (jj=ii; jj<n; jj++)
 {
   tmp=z[ii][jj];
  z[ii][jj]=z[jj][ii];
  z[jj][ii]=tmp;
 }

     return;
}
           
        
#include <math.h>
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
   offset the index of the leading element  of the input(full) matrix
          to be factored.
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
	  (return value)


   Functions used
   --------------
   UTILITY	fsign
   MISC		pythag

 ***********************************************************************/

long tql2(long offset, long n, double *d, double *e, double **z)

{
   long j, last, l, l1, l2, m, i, k, iteration;
   double tst1, tst2, g, r, s, s2, c, c2, c3, p, f, h, el1, dl1;
   if (n == 1) return(0);
   f = ZERO;
   last = n - 1;
   tst1 = ZERO;
   e[last] = ZERO;

   for (l = offset; l < n; l++) {
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
	       for (k = offset; k < n; k ++) {
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
   for (l = 1+offset; l < n; l++) {
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
	  for (j = offset; j < n; j ++) {
	     p = z[i][j];
	     z[i][j] = z[k][j];
	     z[k][j] = p;
	  }
      }   
   }
   return(0);
}		/*...... end main ............................*/
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
/************************************************************** 
 * Function forms the dot product of two vectors.      	      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 

extern long flag;

double ddot( long n,double *dx,long incx,double *dy,long incy)

{
   long i;
   double dot_product;

   if (n <= 0 || incx == 0 || incy == 0) return(0.0);
   dot_product = 0.0;
   if (incx == 1 && incy == 1) 
      for (i=0; i < n; i++) 
      dot_product += (*dx++) * (*dy++);
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
