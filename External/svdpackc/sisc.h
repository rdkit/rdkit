/****************************************************************
 * Sparse SVD Via Subspace Iteraton Procedure for Equivalent    *
 * 2-Cyclic and A'A Eigensystems.                               *
 *                                                              *
 * Constants used by sis1 and sis2.                             *
 ****************************************************************/

#define NMAX  25000     /*Cyclic eigensystem:
                        upper bound on {number of rows (nrow) + number 
                        of columns (ncol) } for the sparse matrix A 
                        stored in Harwell-Boeing format. The input 
			file "matrix" must have (nrow + ncol) <= NMAX */

#define NZMAX 400000   /*NZMAX is an upper bound on the number of 
                        non-zeroes for the sparse matrix A stored 
                        in "matrix", i.e., nnzero <= NZMAX */

#define NSIG  2000      /*NSIG is an upper bound on the number of 
                        desired singular triplets of the sparse 
                        matrix A, i.e., em2 <= NSIG*/
#define ZERO  0.0
#define ONE   1.0
