#include <stdio.h>
#include <stdlib.h>
#include "svdlib.h"
#include "svdutil.h"

char *SVDVersion = "1.01";
long SVDVerbosity = 1;
long SVDCount[SVD_COUNTERS];

void svdResetCounters(void) {
  int i;
  for (i = 0; i < SVD_COUNTERS; i++)
    SVDCount[i] = 0;
}

/********************************* Allocation ********************************/

/* Row major order.  Rows are vectors that are consecutive in memory.  Matrix
   is initialized to empty. */
DMat svdNewDMat(int rows, int cols) {
  int i;
  DMat D = malloc(sizeof(struct dmat));
  if (!D) {perror("svdNewDMat"); return NULL;}
  D->rows = rows;
  D->cols = cols;

  D->value = (double **) malloc(rows * sizeof(double *));
  if (!D->value) {SAFE_FREE(D); return NULL;}

  D->value[0] = (double *) calloc(rows * cols, sizeof(double));
  if (!D->value[0]) {SAFE_FREE(D->value); SAFE_FREE(D); return NULL;}

  for (i = 1; i < rows; i++) D->value[i] = D->value[i-1] + cols;
  return D;
}

void svdFreeDMat(DMat D) {
  if (!D) return;
  SAFE_FREE(D->value[0]);
  SAFE_FREE(D->value);
  free(D);
}


SMat svdNewSMat(int rows, int cols, int vals) {
  SMat S = calloc(1, sizeof(struct smat));
  if (!S) {perror("svdNewSMat"); return NULL;}
  S->rows = rows;
  S->cols = cols;
  S->vals = vals;
  S->pointr = longArray(cols + 1, TRUE, "svdNewSMat: pointr");
  if (!S->pointr) {svdFreeSMat(S); return NULL;}
  S->rowind = longArray(vals, FALSE, "svdNewSMat: rowind");
  if (!S->rowind) {svdFreeSMat(S); return NULL;}
  S->value  = doubleArray(vals, FALSE, "svdNewSMat: value");
  if (!S->value)  {svdFreeSMat(S); return NULL;}
  return S;
}

void svdFreeSMat(SMat S) {
  if (!S) return;
  SAFE_FREE(S->pointr);
  SAFE_FREE(S->rowind);
  SAFE_FREE(S->value);
  free(S);
}


/* Creates an empty SVD record */
SVDRec svdNewSVDRec(void) {
  SVDRec R = calloc(1, sizeof(struct svdrec));
  if (!R) {perror("svdNewSVDRec"); return NULL;}
  return R;
}

/* Frees an svd rec and all its contents. */
void svdFreeSVDRec(SVDRec R) {
  if (!R) return;
  if (R->Ut) svdFreeDMat(R->Ut);
  if (R->S) SAFE_FREE(R->S);
  if (R->Vt) svdFreeDMat(R->Vt);
  free(R);
}


/**************************** Conversion *************************************/

/* Converts a sparse matrix to a dense one (without affecting the former) */
DMat svdConvertStoD(SMat S) {
  int i, c;
  DMat D = svdNewDMat(S->rows, S->cols);
  if (!D) {
    error("svdConvertStoD: failed to allocate D");
    return NULL;
  }
  for (i = 0, c = 0; i < S->vals; i++) {
    while (S->pointr[c + 1] <= i) c++;
    D->value[S->rowind[i]][c] = S->value[i];
  }
  return D;
}

/* Converts a dense matrix to a sparse one (without affecting the dense one) */
SMat svdConvertDtoS(DMat D) {
  SMat S;
  int i, j, n;
  for (i = 0, n = 0; i < D->rows; i++)
    for (j = 0; j < D->cols; j++)
      if (D->value[i][j] != 0) n++;
  
  S = svdNewSMat(D->rows, D->cols, n);
  if (!S) {
    error("svdConvertDtoS: failed to allocate S");
    return NULL;
  }
  for (j = 0, n = 0; j < D->cols; j++) {
    S->pointr[j] = n;
    for (i = 0; i < D->rows; i++)
      if (D->value[i][j] != 0) {
        S->rowind[n] = i;
        S->value[n] = D->value[i][j];
        n++;
      }
  }
  S->pointr[S->cols] = S->vals;
  return S;
}


/**************************** Input/Output ***********************************/

void svdWriteDenseArray(double *a, int n, char *filename, char binary) {
  int i;
  FILE *file = writeFile(filename, FALSE);
  if (!file){ 
    error("svdWriteDenseArrayFile: failed to write %s", filename);
    return;
  }
  if (binary) {
    writeBinInt(file, n);
    for (i = 0; i < n; i++)
      writeBinFloat(file, (float) a[i]);
  } else {
    fprintf(file, "%d\n", n);
    for (i = 0; i < n; i++)
      fprintf(file, "%g\n", a[i]);
  }
  closeFile(file);
}

/* File format has a funny header, then first entry index per column, then the
   row for each entry, then the value for each entry.  Indices count from 1.
   Assumes A is initialized. */
static SMat svdLoadSparseTextHBFile(FILE *file) {
  long i, x, rows, cols, vals;
  SMat S;
  if (fscanf(file, "%*72c%*s %*s %*s %ld %ld %ld %*d %*s %*s %*s %*s", 
             &rows, &cols, &vals) != 3) {
    error("svdLoadSparseTextHBFile: bad file format");
    return NULL;
  }
  
  S = svdNewSMat(rows, cols, vals);
  if (!S) return NULL;

  /* Read column pointers. */
  for (i = 0; i <= S->cols; i++) {
    if (fscanf(file, " %ld", &x) != 1) {
      error("svdLoadSparseTextHBFile: error reading pointr %d", i);
      return NULL;
    }
    S->pointr[i] = x - 1;
  }
  S->pointr[S->cols] = S->vals;
  
  /* Read row indices. */
  for (i = 0; i < S->vals; i++) {
    if (fscanf(file, " %ld", &x) != 1) {
      error("svdLoadSparseTextHBFile: error reading rowind %d", i);
      return NULL;
    }
    S->rowind[i] = x - 1;
  }
  for (i = 0; i < S->vals; i++) 
    if (fscanf(file, " %lf", S->value + i) != 1) {
      error("svdLoadSparseTextHBFile: error reading value %d", i);
      return NULL;
    }
  return S;
}

static void svdWriteSparseTextHBFile(SMat S, FILE *file) {
  int i;
  fprintf(file, "%-72sX\nTransposed\nrra %ld %ld %ld 0\nX X X X\n",
          "SVDLIBC", S->rows, S->cols, S->vals);
  for (i = 0; i <= S->cols; i++)
    fprintf(file, "%ld%s", S->pointr[i] + 1, (((i+1) % 20) == 0) ? "\n" : " ");
  fprintf(file, "\n");
  for (i = 0; i < S->vals; i++)
    fprintf(file, "%ld%s", S->rowind[i] + 1, (((i+1) % 20) == 0) ? "\n" : " ");
  fprintf(file, "\n");
  for (i = 0; i < S->vals; i++)
    fprintf(file, "%g%s", S->value[i], (((i+1) % 20) == 0) ? "\n" : " ");
  fprintf(file, "\n");
}


static SMat svdLoadSparseTextFile(FILE *file) {
  long c, i, n, v, rows, cols, vals;
  SMat S;
  if (fscanf(file, " %ld %ld %ld", &rows, &cols, &vals) != 3) {
    error("svdLoadSparseTextFile: bad file format");
    return NULL;
  }

  S = svdNewSMat(rows, cols, vals);
  if (!S) return NULL;
  
  for (c = 0, v = 0; c < cols; c++) {
    if (fscanf(file, " %ld", &n) != 1) {
      error("svdLoadSparseTextFile: bad file format");
      return NULL;
    }
    S->pointr[c] = v;
    for (i = 0; i < n; i++, v++) {
      if (fscanf(file, " %ld %lf", S->rowind + v, S->value + v) != 2) {
        error("svdLoadSparseTextFile: bad file format");
        return NULL;
      }
    }
  }
  S->pointr[cols] = vals;
  return S;
}

static void svdWriteSparseTextFile(SMat S, FILE *file) {
  int c, v;
  fprintf(file, "%ld %ld %ld\n", S->rows, S->cols, S->vals);
  for (c = 0, v = 0; c < S->cols; c++) {
    fprintf(file, "%ld\n", S->pointr[c + 1] - S->pointr[c]);
    for (; v < S->pointr[c+1]; v++)
      fprintf(file, "%ld %g\n", S->rowind[v], S->value[v]);
  }
}


static SMat svdLoadSparseBinaryFile(FILE *file) {
  int rows, cols, vals, n, c, i, v, r, e = 0;
  float f;
  SMat S;
  e += readBinInt(file, &rows);
  e += readBinInt(file, &cols);
  e += readBinInt(file, &vals);
  if (e) {
    error("svdLoadSparseBinaryFile: bad file format");
    return NULL;
  }

  S = svdNewSMat(rows, cols, vals);
  if (!S) return NULL;
  
  for (c = 0, v = 0; c < cols; c++) {
    if (readBinInt(file, &n)) {
      error("svdLoadSparseBinaryFile: bad file format");
      return NULL;
    }
    S->pointr[c] = v;
    for (i = 0; i < n; i++, v++) {
      e += readBinInt(file, &r);
      e += readBinFloat(file, &f);
      if (e) {
        error("svdLoadSparseBinaryFile: bad file format");
        return NULL;
      }
      S->rowind[v] = r;
      S->value[v] = f;
    }
  }
  S->pointr[cols] = vals;
  return S;
}

static void svdWriteSparseBinaryFile(SMat S, FILE *file) {
  int c, v;
  writeBinInt(file, (int) S->rows);
  writeBinInt(file, (int) S->cols);
  writeBinInt(file, (int) S->vals);
  for (c = 0, v = 0; c < S->cols; c++) {
    writeBinInt(file, (int) (S->pointr[c + 1] - S->pointr[c]));
    for (; v < S->pointr[c+1]; v++) {
      writeBinInt(file, (int) S->rowind[v]);
      writeBinFloat(file, (float) S->value[v]);
    }
  }
}


static DMat svdLoadDenseTextFile(FILE *file) {
  long rows, cols, i, j;
  DMat D;
  if (fscanf(file, " %ld %ld", &rows, &cols) != 2) {
    error("svdLoadDenseTextFile: bad file format");
    return NULL;
  }

  D = svdNewDMat(rows, cols);
  if (!D) return NULL;

  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++) {
      if (fscanf(file, " %lf", &(D->value[i][j])) != 1) {
        error("svdLoadDenseTextFile: bad file format");
        return NULL;
      }
    }
  return D;
}

static void svdWriteDenseTextFile(DMat D, FILE *file) {
  int i, j;
  fprintf(file, "%ld %ld\n", D->rows, D->cols);
  for (i = 0; i < D->rows; i++)
    for (j = 0; j < D->cols; j++) 
      fprintf(file, "%g%c", D->value[i][j], (j == D->cols - 1) ? '\n' : ' ');
}


static DMat svdLoadDenseBinaryFile(FILE *file) {
  int rows, cols, i, j, e = 0;
  float f;
  DMat D;
  e += readBinInt(file, &rows);
  e += readBinInt(file, &cols);
  if (e) {
    error("svdLoadDenseBinaryFile: bad file format");
    return NULL;
  }

  D = svdNewDMat(rows, cols);
  if (!D) return NULL;

  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++) {
      if (readBinFloat(file, &f)) {
        error("svdLoadDenseBinaryFile: bad file format");
        return NULL;
      }
      D->value[i][j] = f;
    }
  return D;
}

static void svdWriteDenseBinaryFile(DMat D, FILE *file) {
  int i, j;
  writeBinInt(file, (int) D->rows);
  writeBinInt(file, (int) D->cols);
  for (i = 0; i < D->rows; i++)
    for (j = 0; j < D->cols; j++) 
      writeBinFloat(file, (float) D->value[i][j]);
}


SMat svdLoadSparseMatrix(char *filename, int format) {
  SMat S = NULL;
  DMat D = NULL;
  FILE *file = fatalReadFile(filename);
  switch (format) {
  case SVD_F_STH: 
    S = svdLoadSparseTextHBFile(file);
    break;
  case SVD_F_ST:
    S = svdLoadSparseTextFile(file);
    break;
  case SVD_F_SB:
    S = svdLoadSparseBinaryFile(file);
    break;
  case SVD_F_DT:
    D = svdLoadDenseTextFile(file);
    break;
  case SVD_F_DB:
    D = svdLoadDenseBinaryFile(file);
    break;
  default: error("svdLoadSparseMatrix: unknown format %d", format);
  }
  closeFile(file);
  if (D) {
    S = svdConvertDtoS(D);
    svdFreeDMat(D);
  }
  return S;
}

DMat svdLoadDenseMatrix(char *filename, int format) {
  SMat S = NULL;
  DMat D = NULL;
  FILE *file = fatalReadFile(filename);
  switch (format) {
  case SVD_F_STH: 
    S = svdLoadSparseTextHBFile(file);
    break;
  case SVD_F_ST:
    S = svdLoadSparseTextFile(file);
    break;
  case SVD_F_SB:
    S = svdLoadSparseBinaryFile(file);
    break;
  case SVD_F_DT:
    D = svdLoadDenseTextFile(file);
    break;
  case SVD_F_DB:
    D = svdLoadDenseBinaryFile(file);
    break;
  default: error("svdLoadSparseMatrix: unknown format %d", format);
  }
  closeFile(file);
  if (S) {
    D = svdConvertStoD(S);
    svdFreeSMat(S);
  }
  return D;
}

void svdWriteSparseMatrix(SMat S, char *filename, int format) {
  DMat D = NULL;
  FILE *file = writeFile(filename, FALSE);
  if (!file) {
    error("svdWriteSparseMatrix: failed to write file %s\n", filename);
    return;
  }
  switch (format) {
  case SVD_F_STH: 
    svdWriteSparseTextHBFile(S, file);
    break;
  case SVD_F_ST:
    svdWriteSparseTextFile(S, file);
    break;
  case SVD_F_SB:
    svdWriteSparseBinaryFile(S, file);
    break;
  case SVD_F_DT:
    D = svdConvertStoD(S);
    svdWriteDenseTextFile(D, file);
    break;
  case SVD_F_DB:
    D = svdConvertStoD(S);
    svdWriteDenseBinaryFile(D, file);
    break;
  default: error("svdLoadSparseMatrix: unknown format %d", format);
  }
  closeFile(file);
  if (D) svdFreeDMat(D);
}

void svdWriteDenseMatrix(DMat D, char *filename, int format) {
  SMat S = NULL;
  FILE *file = writeFile(filename, FALSE);
  if (!file) {
    error("svdWriteDenseMatrix: failed to write file %s\n", filename);
    return;
  }
  switch (format) {
  case SVD_F_STH: 
    S = svdConvertDtoS(D);
    svdWriteSparseTextHBFile(S, file);
    break;
  case SVD_F_ST:
    S = svdConvertDtoS(D);
    svdWriteSparseTextFile(S, file);
    break;
  case SVD_F_SB:
    S = svdConvertDtoS(D);
    svdWriteSparseBinaryFile(S, file);
    break;
  case SVD_F_DT:
    svdWriteDenseTextFile(D, file);
    break;
  case SVD_F_DB:
    svdWriteDenseBinaryFile(D, file);
    break;
  default: error("svdLoadSparseMatrix: unknown format %d", format);
  }
  closeFile(file);
  if (S) svdFreeSMat(S);
}
