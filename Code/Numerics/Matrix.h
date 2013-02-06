//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_MATRIX_H__
#define __RD_MATRIX_H__

#include <RDGeneral/Invariant.h>
#include "Vector.h"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <boost/smart_ptr.hpp>

//#ifndef INVARIANT_SILENT_METHOD
//#define INVARIANT_SILENT_METHOD
//#endif

namespace RDNumeric {
  
  //! A matrix class for general, non-square matrices
  template <class TYPE> class Matrix {
  public:

    typedef boost::shared_array<TYPE> DATA_SPTR;

    //! Initialize with a size.
    Matrix(unsigned int nRows, unsigned int nCols) : 
      d_nRows(nRows), d_nCols(nCols), d_dataSize(nRows*nCols) {
      TYPE *data = new TYPE[d_dataSize];
      memset(static_cast<void *>(data),0,d_dataSize*sizeof(TYPE));
      d_data.reset(data);
    };

    //! Initialize with a size and default value.
    Matrix(unsigned int nRows, unsigned int nCols, TYPE val) :
      d_nRows(nRows), d_nCols(nCols), d_dataSize(nRows*nCols) {
      TYPE *data = new TYPE[d_dataSize];
      unsigned int i;
      for (i = 0; i < d_dataSize; i++) {
        data[i] = val;
      }
      d_data.reset(data);
    }

    //! Initialize from a pointer.
    /*!
      <b>NOTE:</b> this does not take ownership of the data,
      if you delete the data externally, this Matrix will be sad.
    */
    Matrix(unsigned int nRows, unsigned int nCols, DATA_SPTR data) :
      d_nRows(nRows), d_nCols(nCols), d_dataSize(nRows*nCols) {
      d_data = data;
    }

    //! copy constructor
    /*! We make a copy of the other vector's data.
     */
    Matrix(const Matrix<TYPE> &other) :
      d_nRows(other.numRows()), d_nCols(other.numCols()), d_dataSize(d_nRows*d_nCols) {
      TYPE *data = new TYPE[d_dataSize];
      const TYPE *otherData = other.getData();
      memcpy(static_cast<void *>(data), static_cast<const void *>(otherData),
	     d_dataSize*sizeof(TYPE));
      d_data.reset(data);
    }

    virtual ~Matrix() {
    }

    //! returns the number of rows
    inline unsigned int numRows() const {
      return d_nRows;
    }
    
    //! returns the number of columns
    inline unsigned int numCols() const {
      return d_nCols;
    }

    inline unsigned int getDataSize() const {
      return d_dataSize;
    }

    //! returns a particular element of the matrix
    inline virtual TYPE getVal(unsigned int i, unsigned int j) const {
      PRECONDITION(i<d_nRows,"bad index");
      PRECONDITION(j<d_nCols,"bad index");
      unsigned int id = i*d_nCols + j;
      return d_data[id];
    }

    //! sets a particular element of the matrix
    inline virtual void setVal(unsigned int i, unsigned int j, TYPE val) {
      PRECONDITION(i<d_nRows,"bad index");
      PRECONDITION(j<d_nCols,"bad index");
      unsigned int id = i*d_nCols + j;
      
      d_data[id] = val;
    }

    //! returns a copy of a row of the matrix
    inline virtual void getRow(unsigned int i, Vector<TYPE> &row) const { 
      PRECONDITION(i<d_nRows,"bad index");
      PRECONDITION(d_nCols == row.size(), "");
      unsigned int id = i*d_nCols;
      TYPE *rData  = row.getData(); 
      TYPE *data = d_data.get();
      memcpy(static_cast<void *>(rData), static_cast<void *>(&data[id]), d_nCols*sizeof(TYPE));
      
    }
     
    //! returns a copy of a column of the matrix
    inline virtual void getCol(unsigned int i, Vector<TYPE> &col) const { 
      PRECONDITION(i<d_nCols,"bad index");
      PRECONDITION(d_nRows == col.size(), "");
      unsigned int j,id;
      TYPE *rData  = col.getData(); 
      TYPE *data = d_data.get();
      for (j = 0; j < d_nRows; j++) {
        id = j*d_nCols + i;
        rData[j] = data[id];
      }
    }

    //! returns a pointer to our data array
    inline TYPE *getData() {
      return d_data.get();
    }
    
    //! returns a const pointer to our data array
    inline const TYPE *getData() const {
      return d_data.get();
    }


    //! Copy operator.
    /*! We make a copy of the other Matrix's data.
     */
    
    Matrix<TYPE>& assign(const Matrix<TYPE> &other) {
      PRECONDITION(d_nRows == other.numRows(), "Num rows mismatch in matrix copying");
      PRECONDITION(d_nCols == other.numCols(), "Num cols mismatch in matrix copying");
      const TYPE *otherData = other.getData();
      TYPE *data = d_data.get();
      memcpy(static_cast<void *>(data), static_cast<const void *>(otherData), d_dataSize*sizeof(TYPE));
      return *this;
    }

    //! Matrix addition.
    /*! Perform a element by element addition of other Matrix to this Matrix
     */
    virtual Matrix<TYPE>& operator+=(const Matrix<TYPE> &other) {
      PRECONDITION(d_nRows == other.numRows(), "Num rows mismatch in matrix addition");
      PRECONDITION(d_nCols == other.numCols(), "Num cols mismatch in matrix addition");
      const TYPE *oData = other.getData();
      unsigned int i;
      TYPE *data = d_data.get();
      for (i = 0; i < d_dataSize; i++) {
        data[i] += oData[i];
      }
      return *this;
    }

    //! Matrix subtraction.
    /*! Perform a element by element subtraction of other Matrix from this Matrix
     */
    virtual Matrix<TYPE>& operator-=(const Matrix<TYPE> &other) {
      PRECONDITION(d_nRows == other.numRows(), "Num rows mismatch in matrix addition");
      PRECONDITION(d_nCols == other.numCols(), "Num cols mismatch in matrix addition");
      const TYPE *oData = other.getData();
      unsigned int i;
      TYPE *data = d_data.get();
      for (i = 0; i < d_dataSize; i++) {
        data[i] -= oData[i];
      }
      return *this;
    }

    //! Multiplication by a scalar
    virtual Matrix<TYPE>& operator*=(TYPE scale) {
      unsigned int i;
      TYPE *data = d_data.get();
      for (i = 0; i < d_dataSize; i++) {
        data[i] *= scale;
      }
      return *this;
    }

    //! division by a scalar
    virtual Matrix<TYPE>& operator/=(TYPE scale) {
      unsigned int i;
      TYPE *data = d_data.get();
      for (i = 0; i < d_dataSize; i++) {
        data[i] /= scale;
      }
      return *this;
    }

    //! copies the transpose of this Matrix into another, returns the result
    /*!
      \param transpose the Matrix to store the results

      \return the transpose of this matrix.
         This is just a reference to the argument.

     */
    virtual Matrix<TYPE>& transpose(Matrix<TYPE> &transpose) const {
      unsigned int tRows = transpose.numRows();
      unsigned int tCols = transpose.numCols();
      PRECONDITION(d_nCols == tRows, "Size mismatch during transposing");
      PRECONDITION(d_nRows == tCols, "Size mismatch during transposing");
      unsigned int i, j;
      unsigned int idA, idAt, idT;
      TYPE *tData = transpose.getData(); 
      TYPE *data = d_data.get();
      for (i = 0; i < d_nRows; i++) {
        idA = i*d_nCols;
        for (j = 0; j < d_nCols; j++) {
          idAt = idA + j;
          idT = j*tCols + i;
          tData[idT] = data[idAt];
        }
      }
      return transpose;
    }

  protected:
    Matrix() : d_nRows(0), d_nCols(0), d_dataSize(0), d_data(){} ;
    unsigned int d_nRows;
    unsigned int d_nCols;
    unsigned int d_dataSize;
    DATA_SPTR d_data;
   
  private:
    Matrix<TYPE>& operator=(const Matrix<TYPE> &other);
  };

  //! Matrix multiplication
  /*!
    Multiply a Matrix A with a second Matrix B 
    so the result is C = A*B
    
    \param A  the the first Matrix used in the multiplication
    \param B  the Matrix by which to multiply
    \param C  Matrix to use for the results
    
    \return the results of multiplying A by B.
    This is just a reference to C.
  */
  template <class TYPE>
    Matrix<TYPE>& multiply(const Matrix<TYPE> &A, const Matrix<TYPE> &B, 
                           Matrix<TYPE> &C)  {
    unsigned int aRows = A.numRows();
    unsigned int aCols = A.numCols();
    unsigned int cRows = C.numRows(); 
    unsigned int cCols = C.numCols();
    unsigned int bRows = B.numRows();
    unsigned int bCols = B.numCols();
    CHECK_INVARIANT(aCols == bRows, "Size mismatch during multiplication");
    CHECK_INVARIANT(aRows == cRows, "Size mismatch during multiplication");
    CHECK_INVARIANT(bCols == cCols, "Size mismatch during multiplication");
    
    // we have the sizes check do do the multiplication
    TYPE *cData = C.getData();
    const TYPE *bData = B.getData();
    const TYPE *aData = A.getData();
    unsigned int i, j, k;
    unsigned int idA, idAt, idB, idC, idCt;
    for (i = 0; i < aRows; i++) {
      idC = i*cCols;
      idA = i*aCols;
      for (j = 0; j < cCols; j++) {
        idCt = idC + j;
        cData[idCt] = (TYPE)0.0;
        for (k = 0; k < aCols; k++) {
          idAt = idA + k;
          idB = k*bCols + j;
          cData[idCt] += (aData[idAt]*bData[idB]);
        }
      }
    }
    return C;
  };

  //! Matrix-Vector multiplication
  /*!
    Multiply a Matrix A with a Vector x
    so the result is y = A*x
    
    \param A  the matrix used in the multiplication
    \param x  the Vector by which to multiply
    \param y  Vector to use for the results
    
    \return the results of multiplying x by this
    This is just a reference to y.
  */
  template <class TYPE>
    Vector<TYPE>& multiply(const Matrix<TYPE> &A, const Vector<TYPE> &x, 
                           Vector<TYPE> &y) {
    unsigned int aRows = A.numRows();
    unsigned int aCols = A.numCols();
    unsigned int xSiz = x.size();
    unsigned int ySiz = y.size();
    CHECK_INVARIANT(aCols == xSiz, "Size mismatch during multiplication");
    CHECK_INVARIANT(aRows == ySiz, "Size mismatch during multiplication");
    unsigned int i, j;
    unsigned int idA, idAt;
    const TYPE *xData = x.getData();
    const TYPE *aData = A.getData();
    TYPE *yData = y.getData();
    for (i = 0; i < aRows; i++) {
      idA = i*aCols;
      yData[i] = (TYPE)(0.0);
      for (j = 0; j < aCols; j++) {
        idAt = idA + j;
        yData[i] += (aData[idAt]*xData[j]);
      }
    }
    return y;
  };

  typedef Matrix<double> DoubleMatrix;
};

//! ostream operator for Matrix's
template <class TYPE > std::ostream & operator<<(std::ostream& target, 
                                                 const RDNumeric::Matrix<TYPE> &mat) {
  unsigned int nr = mat.numRows();
  unsigned int nc = mat.numCols();
  target << "Rows: " << mat.numRows() << " Columns: " << mat.numCols() << "\n";

  unsigned int i, j;
  for (i = 0; i < nr; i++) {
    for (j = 0; j < nc; j++) {
      target << std::setw(7) << std::setprecision(3) << mat.getVal(i, j);
    }
    target << "\n";
  }
  return target;
}

#endif
