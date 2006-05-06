/****************************************************************
 * Sparse SVD Via Hybrid Block Lanczos Procedure for Equivalent *
 * A'A Eigensystems.                                            *
 *                                                              *
 * Global variables and common areas used by bls2 and its       *
 * functions.                                                   *
 ****************************************************************/

long   *pointr ,       /* pointer to column start array         */
       *rowind ;       /* pointer to row indices array          */
double  *value  ;      /* pointer to nonzero values array       */

long   mxvcount,       /* matrix-vector multiplications counter */ 
       mtxvcount,      /* transposed matrix-vector mult. counter*/
       iconv,          /* converged vector counter              */
       nn,             /* current subspace size                 */
       iter;           /* iteration counter                     */


/* memory for the following arrays are allocated in bls2.c and  *
 * blklan2.c (see blklan2.c for dimensions of these arrays)     */

double *alpha,         /* diagonal elements of bidiagonal matrix*
		        * from single vector Lanczos (inner)    *
			* recursion                             */
       *beta,          /* super-diagonals of bidiagonal matrix  *
			* from single vector Lanczos (inner)    *
			* recursion                             */
       *p,             /* work array                            */
       *q,             /* work array                            */
       *t,             /* work array                            */
       *z;             /* work array                            */
double *tres,          /* temporary residual array              */
       *y,             /* work array                            */
       *vv,            /* work array                            */
       *v0,            /* array of converged right S-vectors    */
       *uvtmp,         /* temporary work space                  */
       *pp,            /* left S-vectors of block upper bi-     * 
			* diagonal matrix from outer recursion  */ 
       *qq,            /* right S-vectors of block upper bi-    *
			* diagonal matrix from outer recursion  */
       *ztemp;         /* temporary work space                  */

double **yp,           /* corresponding 2-dimensional           */ 
       **vvp,          /* ...array representation of the above  */
       **uvtmpp,       /* ...linear arrays                      */
       **ppp, 
       **qqp,
       **ztempp;
