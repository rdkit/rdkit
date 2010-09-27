//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_SYMM_MATRIX_H__
#define __RD_SYMM_MATRIX_H__

#include "Matrix.h"
#include "SquareMatrix.h"
#include <cstring>
#include <boost/smart_ptr.hpp>

//#ifndef INVARIANT_SILENT_METHOD
//#define INVARIANT_SILENT_METHOD
//#endif
namespace RDNumeric {
  //! A symmetric matrix class
  /*! 
    The data is stored as the lower triangle, so
     A[i,j] = data[i*(i+1) + j] when i >= j and
     A[i,j] = data[j*(j+1) + i] when i < j
  */
  template <class TYPE> class SymmMatrix {
  public:
    typedef boost::shared_array<TYPE> DATA_SPTR;

    explicit SymmMatrix(unsigned int N) : 
      d_size(N), d_dataSize(N*(N+1)/2)  {
      TYPE *data = new TYPE[d_dataSize];
      memset(static_cast<void *>(data),0,d_dataSize*sizeof(TYPE));
      d_data.reset(data);
    }

    SymmMatrix(unsigned int N, TYPE val) : 
      d_size(N), d_dataSize(N*(N+1)/2)  {
      TYPE *data = new TYPE[d_dataSize];
      unsigned int i;
      for (i = 0; i < d_dataSize; i++) {
        data[i] = val;
      }
      d_data.reset(data);
    }
    
    SymmMatrix(unsigned int N, DATA_SPTR data) :
      d_size(N), d_dataSize(N*(N+1)/2)  {
      d_data = data;
    }
    
    SymmMatrix(const SymmMatrix<TYPE> &other) :
      d_size(other.numRows()), d_dataSize(other.getDataSize())  {
      TYPE *data = new TYPE[d_dataSize];
      const TYPE *otherData = other.getData();

      memcpy(static_cast<void *>(data), static_cast<const void *>(otherData),
	     d_dataSize*sizeof(TYPE));
      d_data.reset(data);
    }

    ~SymmMatrix() {}
    
    //! returns the number of rows
    inline unsigned int numRows() const {
      return d_size;
    }

    //! returns the number of columns
    inline unsigned int numCols() const {
      return d_size;
    }

    inline unsigned int getDataSize() const {
      return d_dataSize;
    }

    void setToIdentity() {
      TYPE *data = d_data.get();
      memset(static_cast<void *>(data), 0, d_dataSize*sizeof(TYPE));
      for (unsigned int i = 0; i < d_size; i++) {
        data[i*(i+3)/2] = (TYPE)1.0;
      }
    }

    TYPE getVal(unsigned int i, unsigned int j) const {
      RANGE_CHECK(0, i, d_size-1);
      RANGE_CHECK(0, j, d_size-1);
      unsigned int id;
      if (i >= j) {
        id = i*(i+1)/2 + j;
      } else {
        id = j*(j+1)/2 + i;
      }
      return d_data[id];
    }

    void setVal(unsigned int i, unsigned int j, TYPE val) {
      RANGE_CHECK(0, i, d_size-1);
      RANGE_CHECK(0, j, d_size-1);
      unsigned int id;
      if (i >= j) {
        id = i*(i+1)/2 + j;
      } else {
        id = j*(j+1)/2 + i;
      }
      d_data[id] = val;
    }

    void getRow(unsigned int i, Vector<TYPE> &row) { 
      CHECK_INVARIANT(d_size == row.size(), "");
      TYPE *rData  = row.getData(); 
      TYPE *data = d_data.get();
      for (unsigned int j = 0; j < d_size; j++) {
	unsigned int id;
        if (j <= i) {
          id = i*(i+1)/2 + j;
        } else {
          id = j*(j+1)/2 + i;
        }
        rData[j] = data[id];
      }
    }
     
    void getCol(unsigned int i, Vector<TYPE> &col) { 
      CHECK_INVARIANT(d_size == col.size(), "");
      TYPE *rData  = col.getData();
      TYPE *data = d_data.get();
      for (unsigned int j = 0; j < d_size; j++) {
	unsigned int id;
        if (i <= j) {
          id = j*(j+1)/2 + i;
        } else {
          id = i*(i+1)/2 + j;
        }
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

    SymmMatrix<TYPE>& operator*=(TYPE scale) {
      TYPE *data = d_data.get();
      for (unsigned int i = 0; i < d_dataSize; i++) {
        data[i] *= scale;
      }
      return *this;
    }

    SymmMatrix<TYPE>& operator/=(TYPE scale) {
      TYPE *data = d_data.get();
      for (unsigned int i = 0; i < d_dataSize; i++) {
        data[i] /= scale;
      }
      return *this;
    }

    SymmMatrix<TYPE>& operator+=(const SymmMatrix<TYPE> &other) {
      CHECK_INVARIANT(d_size == other.numRows(), "Sizes don't match in the addition");
      const TYPE *oData = other.getData();
      TYPE *data = d_data.get();
      for (unsigned int i = 0; i < d_dataSize; i++) {
        data[i] += oData[i];
      }
      return *this;
    }

    SymmMatrix<TYPE>& operator-=(const SymmMatrix<TYPE> &other) {
      CHECK_INVARIANT(d_size == other.numRows(), "Sizes don't match in the addition");
      const TYPE *oData = other.getData();
      TYPE *data = d_data.get();
      for (unsigned int i = 0; i < d_dataSize; i++) {
        data[i] -= oData[i];
      }
      return *this;
    }

    //! in-place matrix multiplication
    SymmMatrix<TYPE>& operator*=(const SymmMatrix<TYPE> &B) {
      CHECK_INVARIANT(d_size == B.numRows(), "Size mismatch during multiplication");
      TYPE *cData = new TYPE[d_dataSize];
      const TYPE *bData = B.getData();
      TYPE *data = d_data.get();
      for (unsigned int i = 0; i < d_size; i++) {
        unsigned int idC = i*(i+1)/2;
        for (unsigned int j = 0; j < i+1; j++) {
          unsigned int idCt = idC + j;
          cData[idCt] = (TYPE)0.0;
          for (unsigned int k = 0; k < d_size; k++) {
	    unsigned int idA,idB;
            if (k <= i) {
              idA = i*(i+1)/2 + k;
            } else {
              idA = k*(k+1)/2 + i;
            } 
            if (k <= j) {
              idB = j*(j+1)/2 + k;
            } else {
              idB = k*(k+1)/2 + j;
            }
            cData[idCt] += (data[idA]*bData[idB]);
          }
        }
      }
      
      for (unsigned int i = 0; i < d_dataSize; i++) {
        data[i] = cData[i];
      }
      delete [] cData;
      return (*this);
    }

    /* Transpose will basically return a copy of itself
     */
    SymmMatrix<TYPE>& transpose(SymmMatrix<TYPE> &transpose) const { 
      CHECK_INVARIANT(d_size == transpose.numRows(), "Size mismatch during transposing");
      TYPE *tData = transpose.getData(); 
      TYPE *data = d_data.get();
      for (unsigned int i = 0; i < d_dataSize; i++) {
        tData[i] = data[i];
      }
      return transpose;
    }

    SymmMatrix<TYPE> &transposeInplace() {
      // nothing to be done we are symmetric
      return (*this);
    }

  protected: 
    
    SymmMatrix() : d_size(0), d_dataSize(0), d_data(0){};
    unsigned int d_size;
    unsigned int d_dataSize;
    DATA_SPTR d_data;

  private:
    SymmMatrix<TYPE>& operator=(const SymmMatrix<TYPE> &other);
  };
  
  //! SymmMatrix-SymmMatrix multiplication 
  /*!
    Multiply SymmMatrix A with a second SymmMatrix B 
    and write the result to C = A*B

    \param A  the first SymmMatrix 
    \param B  the second SymmMatrix to multiply 
    \param C  SymmMatrix to use for the results
    
    \return the results of multiplying A by B.
    This is just a reference to C.
    
    This method is reimplemented here for efficiency reasons
    (we basically don't want to use getter and setter functions)
    
  */
  template <class TYPE>
    SymmMatrix<TYPE>& multiply(const SymmMatrix<TYPE> &A,
			       const SymmMatrix<TYPE> &B, 
                               SymmMatrix<TYPE> &C) {
    unsigned int aSize = A.numRows();
    CHECK_INVARIANT(B.numRows() == aSize, "Size mismatch in matric multiplication");
    CHECK_INVARIANT(C.numRows() == aSize, "Size mismatch in matric multiplication");
    TYPE *cData = C.getData();
    const TYPE *aData = A.getData();
    const TYPE *bData = B.getData();
    for (unsigned int i = 0; i < aSize; i++) {
      unsigned int idC = i*(i+1)/2;
      for (unsigned int j = 0; j < i+1; j++) {
        unsigned int idCt = idC + j;
        cData[idCt] = (TYPE)0.0;
        for (unsigned int k = 0; k < aSize; k++) {
	  unsigned int idA,idB;
          if (k <= i) {
            idA = i*(i+1)/2 + k;
          } else {
            idA = k*(k+1)/2 + i;
          } 
          if (k <= j) {
            idB = j*(j+1)/2 + k;
          } else {
            idB = k*(k+1)/2 + j;
          }
          cData[idCt] += (aData[idA]*bData[idB]);
        }
      }
    }
    return C;
  }

  //! SymmMatrix-Vector multiplication
  /*!
    Multiply a SymmMatrix A with a Vector x
    so the result is y = A*x
    
    \param A  the SymmMatrix for multiplication 
    \param x  the Vector by which to multiply
    \param y  Vector to use for the results
    
    \return the results of multiplying x by A
    This is just a reference to y.
    
    This method is reimplemented here for efficiency reasons
    (we basically don't want to use getter and setter functions)
    
  */
  template <class TYPE>
    Vector<TYPE>& multiply(const SymmMatrix<TYPE> &A, const Vector<TYPE> &x, 
                                   Vector<TYPE> &y) {
    unsigned int aSize = A.numRows();
    CHECK_INVARIANT(aSize == x.size(), "Size mismatch during multiplication");
    CHECK_INVARIANT(aSize == y.size(), "Size mismatch during multiplication");
    const TYPE *xData = x.getData();
    const TYPE *aData = A.getData();
    TYPE *yData = y.getData();
    for (unsigned int i = 0; i < aSize; i++) {
      yData[i] = (TYPE)(0.0);
      unsigned int idA = i*(i+1)/2;
      for (unsigned int j = 0; j < i+1; j++) {
        //idA = i*(i+1)/2 + j;
        yData[i] += (aData[idA]*xData[j]);
	idA++;
      }
      idA--;
      for (unsigned int j = i+1; j < aSize; j++) {
        //idA = j*(j+1)/2 + i;
	idA += j;
        yData[i] += (aData[idA]*xData[j]);
      }
    }
    return y;
  }

  typedef SymmMatrix<double> DoubleSymmMatrix;
  typedef SymmMatrix<int> IntSymmMatrix;
  typedef SymmMatrix<unsigned int> UintSymmMatrix;
}

//! ostream operator for Matrix's
template <class TYPE > std::ostream & operator<<(std::ostream& target, 
                                                 const RDNumeric::SymmMatrix<TYPE> &mat) {
  unsigned int nr = mat.numRows();
  unsigned int nc = mat.numCols();
  target << "Rows: " << mat.numRows() << " Columns: " << mat.numCols() << "\n";

  for (unsigned int i = 0; i < nr; i++) {
    for (unsigned int j = 0; j < nc; j++) {
      target << std::setw(7) << std::setprecision(3) << mat.getVal(i, j);
    }
    target << "\n";
  }
  return target;
}

#endif
    
