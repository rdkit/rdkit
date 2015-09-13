
#include "JacobiEigenSolver.h"


#define TOLERANCE 1e-6

namespace RDNumeric {
  namespace EigenSolvers{

    unsigned int jacobiEigenSolver (const DoubleSymmMatrix &matrix,
                                    DoubleVector &eigenValues,
                                    DoubleSquareMatrix &eigenVectors,
                                    unsigned int maxIter) {
      unsigned int dim=matrix.numRows();
      PRECONDITION(dim==matrix.numCols(), "Matrix is not square Matrix.");
      PRECONDITION(dim==eigenValues.size(), "EigenValues Vector has wrong length.");
      PRECONDITION(dim==eigenVectors.numRows() && dim==eigenVectors.numCols(), "EigenVectors Matrix has wrong dimensions.");

      double offDiagNorm, diagNorm;
      double b, dma, q, t, c, s;
      double atemp, vtemp, dtemp;

      double mat[dim*dim];
      for (unsigned int j = 0; j < dim; j++) {
        for (unsigned int i = 0; i <= j; i++) {
          mat[j*dim+i] = matrix.getVal(i,j);
          mat[i*dim+j] = mat[j*dim+i];
        }
      }

      double *eVal = eigenValues.getData(),
          *eVec = eigenVectors.getData();

      // initialize the eigen vector to Identity
      for (unsigned int j = 0; j < dim; j++) {
        for (unsigned int i = 0; i < dim; i++) {
          eVec[j*dim+i] = 0.0;
        }
        eVec[j*dim+j] = 1.0;
        eVal[j] = mat[j*dim+j];
      }

      for (unsigned int l = 0; l < maxIter; l++) {
        diagNorm = 0.0;
        offDiagNorm = 0.0;
        for (unsigned int j = 0; j < dim; j++) {
          diagNorm += fabs(eVal[j]);
          for (unsigned int i = 0; i < j; i++) {
            offDiagNorm += fabs(mat[j*dim+i]);
          }
        }
        if((offDiagNorm/diagNorm) <= TOLERANCE) { //converged
          for (unsigned int j = 0; j < dim-1; j++) {
            unsigned int k = j;
            dtemp = eVal[k];
            for (unsigned int i = j+1; i < dim; i++) {
              if(eVal[i] < dtemp) {
                k = i;
                dtemp = eVal[k];
              }
            }

            if(k > j) {
              eVal[k] = eVal[j];
              eVal[j] = dtemp;
              for (unsigned int i = 0; i < dim; i++) {
                dtemp = eVec[k*dim+i];
                eVec[k*dim+i] = eVec[j*dim+i];
                eVec[j*dim+i] = dtemp;
              }
            }
          }
          return l+1;
        }

        for (unsigned int j = 1; j < dim; j++) {
          for (unsigned int i = 0; i < j; i++) {
            b = mat[j*dim+i];
            if(fabs(b) > 0.0) {
              dma = eVal[j] - eVal[i];
              if((fabs(dma) + fabs(b)) <=  fabs(dma)) {
                t = b / dma;
              }
              else {
                q = 0.5 * dma / b;
                t = 1.0/(fabs(q) + sqrt(1.0+q*q));
                if(q < 0.0) {
                  t = -t;
                }
              }
              c = 1.0/sqrt(t * t + 1.0);
              s = t * c;
              mat[j*dim+i] = 0.0;
              for (unsigned int k = 0; k < i; k++) {
                atemp = c * mat[i*dim+k] - s * mat[j*dim+k];
                mat[j*dim+k] = s * mat[i*dim+k] + c * mat[j*dim+k];
                mat[i*dim+k] = atemp;
              }
              for (unsigned int k = i+1; k < j; k++) {
                atemp = c * mat[k*dim+i] - s * mat[j*dim+k];
                mat[j*dim+k] = s * mat[k*dim+i] + c * mat[j*dim+k];
                mat[k*dim+i] = atemp;
              }
              for (unsigned int k = j+1; k < dim; k++) {
                atemp = c * mat[k*dim+i] - s * mat[k*dim+j];
                mat[k*dim+j] = s * mat[k*dim+i] + c * mat[k*dim+j];
                mat[k*dim+i] = atemp;
              }
              for (unsigned int k = 0; k < dim; k++) {
                vtemp = c * eVec[i*dim+k] - s * eVec[j*dim+k];
                eVec[j*dim+k] = s * eVec[i*dim+k] + c * eVec[j*dim+k];
                eVec[i*dim+k] = vtemp;
              }
              dtemp = c*c*eVal[i] + s*s*eVal[j] - 2.0*c*s*b;
              eVal[j] = s*s*eVal[i] + c*c*eVal[j] +  2.0*c*s*b;
              eVal[i] = dtemp;
            }  /* end if */
          } /* end for i */
        } /* end for j */
      } /* end for l */
      return maxIter;
    }
  }
}
