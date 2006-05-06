#ifndef SVDLIB_H
#define SVDLIB_H

#ifndef FALSE
#  define FALSE 0
#endif
#ifndef TRUE
#  define TRUE  1
#endif

/******************************** Structures *********************************/
typedef struct smat *SMat;
typedef struct dmat *DMat;
typedef struct svdrec *SVDRec;

/* Harwell-Boeing sparse matrix. */
struct smat {
  long rows;
  long cols;
  long vals;     /* Total non-zero entries. */
  long *pointr;  /* For each col (plus 1), index of first non-zero entry. */
  long *rowind;  /* For each nz entry, the row index. */
  double *value; /* For each nz entry, the value. */
};

/* Row-major dense matrix.  Rows are consecutive vectors. */
struct dmat {
  long rows;
  long cols;
  double **value; /* Accessed by [row][col]. Free value[0] and value to free.*/
};

struct svdrec {
  int d;      /* Dimensionality (rank) */
  DMat Ut;    /* Transpose of left singular vectors. (d by m)
                 The vectors are the rows of Ut. */
  double *S;  /* Array of singular values. (length d) */
  DMat Vt;    /* Transpose of right singular vectors. (d by n)
                 The vectors are the rows of Vt. */
};


/******************************** Variables **********************************/

/* Version info */
extern char *SVDVersion;

/* How verbose is the package: 0, 1 (default), 2 */
extern long SVDVerbosity;

/* Counter(s) used to track how much work is done in computing the SVD. */
enum svdCounters {SVD_MXV, SVD_COUNTERS};
extern long SVDCount[SVD_COUNTERS];
extern void svdResetCounters(void);

enum svdFileFormats {SVD_F_STH, SVD_F_ST, SVD_F_SB, SVD_F_DT, SVD_F_DB};
/*
File formats:
SVD_F_STH: sparse text, SVDPACK-style
SVD_F_ST:  sparse text, SVDLIB-style
SVD_F_DT:  dense text
SVD_F_SB:  sparse binary
SVD_F_DB:  dense binary
*/

/* True if a file format is sparse: */
#define SVD_IS_SPARSE(format) ((format >= SVD_F_STH) && (format <= SVD_F_SB))


/******************************** Functions **********************************/

/* Creates an empty dense matrix. */
extern DMat svdNewDMat(int rows, int cols);
/* Frees a dense matrix. */
extern void svdFreeDMat(DMat D);

/* Creates an empty sparse matrix. */
SMat svdNewSMat(int rows, int cols, int vals);
/* Frees a sparse matrix. */
void svdFreeSMat(SMat S);

/* Creates an empty SVD record. */
SVDRec svdNewSVDRec(void);
/* Frees an svd rec and all its contents. */
void svdFreeSVDRec(SVDRec R);

/* Converts a sparse matrix to a dense one (without affecting former) */
DMat svdConvertStoD(SMat S);
/* Converts a dense matrix to a sparse one (without affecting former) */
SMat svdConvertDtoS(DMat D);


/* Writes an array to a file. */
extern void svdWriteDenseArray(double *a, int n, char *filename, char binary);

/* Loads a matrix file (in various formats) into a sparse matrix. */
extern SMat svdLoadSparseMatrix(char *filename, int format);
/* Loads a matrix file (in various formats) into a dense matrix. */
extern DMat svdLoadDenseMatrix(char *filename, int format);

/* Writes a dense matrix to a file in a given format. */
extern void svdWriteDenseMatrix(DMat A, char *filename, int format);
/* Writes a sparse matrix to a file in a given format. */
extern void svdWriteSparseMatrix(SMat A, char *filename, int format);


/* Performs the las2 SVD algorithm and returns the resulting Ut, S, and Vt. */
extern SVDRec svdLAS2(SMat A, long lanmax, long maxprs, double end[2], 
                      double kappa);


#endif /* SVDLIB_H */
