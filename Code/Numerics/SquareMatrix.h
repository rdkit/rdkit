//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_SQUARE_MATRIX_H__
#define __RD_SQUARE_MATRIX_H__

#include "Matrix.h"

namespace RDNumeric {
  template <typename TYPE> class SquareMatrix : public Matrix<TYPE> {
  public:
    //! brief Square matrix of size N
    SquareMatrix() {};

    explicit SquareMatrix(unsigned int N) : Matrix<TYPE>(N, N) {};
          
    SquareMatrix(unsigned int N, TYPE val) : Matrix<TYPE>(N, N, val) {};

    SquareMatrix(unsigned int N, typename Matrix<TYPE>::DATA_SPTR data) : Matrix<TYPE>(N, N, data) {};

    //inline unsigned int size() const {
    //  return d_nRows;
    //};

    virtual SquareMatrix<TYPE>& operator*=(TYPE scale) {
      Matrix<TYPE>::operator*=(scale);
      return *this;
    }

    //! In place matrix multiplication
    virtual SquareMatrix<TYPE> & operator*=(const SquareMatrix<TYPE> & B) {
      CHECK_INVARIANT(this->d_nCols == B.numRows(), "Size mismatch during multiplication");

      const TYPE *bData = B.getData();
      TYPE *newData = new TYPE[this->d_dataSize];
      unsigned int i, j, k;
      unsigned int idA, idAt, idC, idCt, idB;
      TYPE* data = this->d_data.get();
      for (i = 0; i < this->d_nRows; i++) {
        idA = i*this->d_nRows;
        idC = idA;
        for (j = 0; j < this->d_nCols; j++) {
          idCt = idC + j;
          newData[idCt] = (TYPE)(0.0);
          for (k = 0; k < this->d_nCols; k++) {
            idAt = idA + k;
            idB = k*this->d_nRows + j;
            newData[idCt] += (data[idAt]*bData[idB]);
          }
        }
      }
      boost::shared_array<TYPE>  tsptr(newData);
      this->d_data.swap(tsptr);
      return (*this);
    }

    //! In place matrix transpose
    virtual SquareMatrix<TYPE> &transposeInplace() {
      unsigned int i,j;
      unsigned int id1, id1t, id2;
      TYPE temp;
      TYPE *data = this->d_data.get();
      for (i = 1; i < this->d_nRows; i++) {
        id1 = i*this->d_nCols;
        for (j = 0; j < i; j++) {
          
          id1t = id1 + j;
          id2 = j*this->d_nCols + i;
          temp = data[id1t];
          data[id1t] = data[id2];
          data[id2] = temp;
        }
      }
      return (*this);
    }

  };
  typedef SquareMatrix<double> DoubleSquareMatrix;
}

#endif
    
