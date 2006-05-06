/****************************************************************
 * Sparse SVD Via Subspace Iteraton Procedure for Equivalent    *
 * 2-Cyclic and A'A Eigensystems.                               *
 *                                                              *
 * Global variables and common areas used by sis1, sis2,and     *
 * their functions.                                             *
 ****************************************************************/

/************************************************************** 
 * pointers to areas holding input matrix which is stored in  *
 * harwell-boeing format.			              *	
 **************************************************************/
long	*pointr = NULL, /* pointer to column start array      */
	*rowind = NULL; /* pointer to row indices array	      */
double  *value = NULL;  /* pointer to nonzero values array    */

FILE	*fp_out1 = NULL, *fp_out2;
long nrow, ncol, mxvcount, maxm = 0, ksteps = 0 ;
double alpha;
