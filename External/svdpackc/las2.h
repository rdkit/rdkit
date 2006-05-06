/**************************************************************
 * Sparse svd via eigensystem of A'A matrix   		      *
 * The equivalent symmetric eigenvalue problem:               *
 *                                                            *
 *  B x = lambda x, where x' = (u',v'), lambda = sigma**2,    *
 *                                                            *
 *  B = A'A, and A is m (nrow) by n (ncol) (nrow >> ncol),    *
 *							      *
 *  so that {u, sqrt(lambda), v} is a singular triplet of A.  *
 *  (A' = transpose of A)				      *
 *                                                            *
 * global variables and common areas used by las2 and its     *
 * procedures.                                                *
 **************************************************************/

#define LMTNW   600000  /* max. size of working area allowed  */
#define NMAX    3000    /* bound on ncol, order of A          */
#define NZMAX   100000  /* bound on number of nonzeros in a   */

long    ierr,           /* error flag                         */
	j,              /* number of lanczos steps taken      */
	neig,           /* number of ritz values stabilized   */
	nsig,           /* number of accepted ritz values     *
			 * based on kappa (relative accuracy) */
    	ncol,           /* number of columns of A             */
    	nrow,           /* number of rows of A                */
	mxvcount = 0;

/**************************************************************
 * pointers to areas holding input matrix which is stored in  *
 * harwell-boeing format.                                     *
 **************************************************************/
long    *pointr = NULL, /* pointer to column start array      */
   	*rowind = NULL; /* pointer to row indices array       */
double  *value = NULL;  /* pointer to nonzero values array    */

double  rnm,            /* norm of the next residual vector   */
	anorm,
	tol,
	eps,            /* positive machine epsilon           */
	eps1,           /* roundoff estimate for dot product  *
			 * of two unit vector                 */
 	reps,
	eps34;

double  *xv1 = NULL,    /* temp arrays needed for computing   */
	*xv2 = NULL,    /* singular vectors                   */
	*ztemp = NULL,

        *a = NULL;      /* pointer to area used by user-      *
			 * supplied procedure store and holds *
		  	 * lanczos vectors                    */

FILE	*fp_out1 = NULL;/* output file pointers               */
long	 fp_out2;

char	*error[10] = {  /* error messages used by function    *
			 * check_parameters                   */
	    NULL,
	  " SORRY, YOUR MATRIX IS TOO BIG ",
	  " ***** ENDL MUST BE LESS THAN ENDR *****",
	  " ***** MAXPRS CANNOT EXCEED LANMAX *****",
	  " ***** N = NROW + NCOL MUST BE GREATER THAN ZERO *****",
	  " ***** LANMAX (NUMBER OF LANCZOS STEPS) IS INVALID *****",
	  " ***** MAXPRS (NUMBER OF IEGENPAIRS DESIRED) IS INVALID *****",
	  " ***** 6*N+4*LANMAX+1 + LANMAX*LANMAX CANNOT EXCEED NW *****",
	  " ***** 6*N+4*LANMAX+1 CANNOT EXCEED NW *****",
	    NULL};
