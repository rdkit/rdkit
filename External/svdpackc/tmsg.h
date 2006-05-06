/**************************************************************
 * Sparse SVD Via Trace Minimization Procedure for Equivalent *
 * 2-Cyclic and A'A Eigensystems.                             *
 *                                                            *
 * Global variables and common areas used by tms1 and its     *
 * procedures.                                                *
 **************************************************************/

long    ierr,           /* error flag                         */
        ncol,           /* number of columns of A             */
        nrow;           /* number of rows of A                */

/**************************************************************
 * pointers to areas holding input matrix which is stored in  *
 * harwell-boeing format.                                     *
 **************************************************************/

long    *pointr = NULL, /* pointer to column start array      */
        *rowind = NULL; /* pointer to row indices array       */
double  *value  = NULL; /* pointer to nonzero values array    */

double  alpha,          /* 1-norm of the matrix A             */
        tol,
        eps;            /* positive machine epsilon           */

long    *iw,**iwork,*lwork;

double  *w1,*w2,*w3,*w4,*w5,*yy;
double  **work1, **work2, **work3, **work4, **work5, **y;
double  *z, *v2, *r, *r2, *pp, *z2, *zz; /* temporary variables  */
        
FILE    *fp_out1 = NULL;/* output file pointers               */
long    fp_out2; /* output file pointers                      */

char    *error[10] = { /*error messages used by function      *
                        *check_parameters                     */
          NULL,
          "SORRY, YOUR MATRIX IS TOO BIG ",
          NULL,
          "*** P CANNOT EXCEED MAXI ***",
          "*** N = NROW + NCOL MUST BE GREATER THAN ZERO ***",
          "*** MAXI (NUMBER OF MIN. ITERATIONS) IS INVALID ***",
          "*** P (NUMBER OF EIGENPAIRS DESIRED) IS INVALID ***",
          NULL,
          NULL,
          NULL};
