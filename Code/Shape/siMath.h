/*******************************************************************************
siMath.h - Shape-it
 
Copyright 2012 by Silicos-it, a division of Imacosi BVBA

This file is part of Shape-it.

	Shape-it is free software: you can redistribute it and/or modify
	it under the terms of the GNU Lesser General Public License as published 
	by the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	Shape-it is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public License
	along with Shape-it.  If not, see <http://www.gnu.org/licenses/>.

Shape-it is linked against OpenBabel version 2.

	OpenBabel is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation version 2 of the License.

***********************************************************************/



#ifndef __SILICOS_PHARAO_SIMATH_H__
#define __SILICOS_PHARAO_SIMATH_H__



// General
#include <vector>
#include <math.h>
#include <stdlib.h>

// OpenBabel

// Shape-it
#ifndef INF
#define INF HUGE_VAL
#endif

#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef TAU
#define TAU 1e-12
#endif

#ifndef min
template < class T > inline T min(T x, T y)
{
    return (x < y) ? x : y;
}
#endif

#ifndef max
template < class T > inline T max(T x, T y)
{
    return (x > y) ? x : y;
}
#endif

#ifndef sign
template < class T > inline T sign(const T & a, const T & b)
{
    return (b >= 0.0) ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}
#endif




namespace SiMath {




    inline double triangle(const double &a, const double &b) {
	double A(fabs(a)), B(fabs(b));
	if (A > B) {
	    return A * sqrt(1.0 + (B / A) * (B / A));
	} else if (B == 0) {
	    return 0;
	}
	return B * sqrt(1.0 + (A / B) * (A / B));
    }



    class Vector {
      private:
	unsigned int _n;	///< Number of data points in vector
	 std::vector < double >_pVector;	///< std::vector to hold all values

      public:
	 Vector():_n(0), _pVector(0) {
	};			///< Empty vector
	Vector(const unsigned int n):_n(n), _pVector(n) {
	};			///< vector of length n, no initial value
	Vector(const unsigned int n, const double &v):_n(n), _pVector(n, v) {
	};			///< vector of length n, initial constant value v
	Vector(const unsigned int n, const double *);	///< Copy of data stored in an array double[n]
	Vector(const std::vector < double >&);	///< Copy of data stored in std::vector<double>

	Vector(const Vector &);

	~Vector();

	void clear();
	void reset(unsigned int);

	void resize(unsigned int);

	double getValueAt(const unsigned int);
	double getValueAt(const unsigned int) const;

	double max() const;
	double max(unsigned int &) const;
	double min() const;
	double min(unsigned int &) const;
	double sum() const;
	double mean() const;
	double stDev() const;
	double stDev(double m) const;
	unsigned int size() const {
	    return _n;
	};

	Vector & operator=(const Vector &);	///< copy assignment, resets the size of the Vector if needed
	Vector & operator=(const double &);	///< set all elements in Vector to constant value
	Vector & operator+=(const double &);	///< add constant value to all elements in Vector
	Vector & operator+=(const Vector &);	///< add full Vector element-wise
	Vector & operator-=(const double &);	///< subtract constant value to all elements in Vector
	Vector & operator-=(const Vector &);	///< subtract full Vector element-wise
	Vector & operator*=(const double &);	///< multiply all elements with a constant value
	Vector & operator*=(const Vector &);	///< multiply full Vector element-wise 
	Vector & operator/=(const double &);	///< divide all elements with a constant value
	Vector & operator/=(const Vector &);	///< divide full Vector element-wise 
	Vector & operator-();	///< change sign of all elements in Vector
	Vector operator+(const Vector &) const;	///< operator to write C = A + B
	Vector operator-(const Vector &) const;	///< operator to write C = A - B
	Vector operator*(const Vector &) const;	///< operator to write C = A * B
	Vector operator/(const Vector &) const;	///< operator to write C = A / B

	bool operator==(const Vector &) const;	///< check if two vectors are the same, which is only true if all elements are the same
	bool operator!=(const Vector &) const;	///< check if two vectors are different

	inline double &operator[] (const unsigned int i) {
	    return _pVector[i];
	};			///< set i-th element from vector
	inline double operator[] (const unsigned int i) const {
	    return _pVector[i];
	};			///< get i-th element from vector (const implementation)

	void swap(const unsigned int, const unsigned int);

	double dotProd(const Vector &);

	const double *getArrayPointer() const {
	    return &(_pVector[0]);
	};			///< direct access to the data
	double *getArrayPointer() {
	    return &(_pVector[0]);
	};			///< direct access to the data
    };



    class Matrix {
      private:
	unsigned int _nRows;
	unsigned int _nCols;
	double **_pMatrix;

      public:
	 Matrix():_nRows(0), _nCols(0), _pMatrix(NULL) {
	};
	Matrix(const unsigned int, const unsigned int);
	Matrix(const unsigned int, const unsigned int, const double &);
	Matrix(const Matrix &);

	Matrix(const unsigned int, const unsigned int, const Vector &);

	~Matrix();

	void reset(const unsigned int, const unsigned int);
	void clear();

	inline unsigned int nbrRows() const {
	    return _nRows;
	};
	inline unsigned int nbrColumns() const {
	    return _nCols;
	};

	double getValueAt(const unsigned int, const unsigned int);
	const double getValueAt(const unsigned int,
				const unsigned int) const;
	Vector getRow(const unsigned int) const;
	Vector getColumn(const unsigned int) const;

	void setValueAt(const unsigned int, const unsigned int, double);
	void setRow(const unsigned int, Vector &);
	void setColumn(const unsigned int, Vector &);

	void swapRows(unsigned int, unsigned int);
	void swapColumns(unsigned int, unsigned int);
	Matrix transpose(void);


	Matrix & operator=(const Matrix &);	///< copy assignment, resets the size of the matrix if needed
	Matrix & operator=(const double &);	///< set all elements in matrix to constant value
	Matrix & operator+=(const double &);	///< add constant value to all elements in matrix
	Matrix & operator+=(const Matrix &);	///< add full matrix element-wise
	Matrix & operator-=(const double &);	///< subtract constant value from all elements in matrix
	Matrix & operator-=(const Matrix &);	///< subtract full matrix element-wise
	Matrix & operator*=(const double &);	///< multiply all elements with a constant value
	Matrix & operator*=(const Matrix &);	///< multiply full matrix element-wise 
	Matrix & operator/=(const double &);	///< divide all elements with a constant value
	Matrix & operator/=(const Matrix &);	///< divide full matrix element-wise 
	Matrix & operator-();	///< change sign of all elements in matrix

	Matrix operator+(const Matrix &) const;	///< add two matrices element by element and store the result in a new matrix
	Matrix operator-(const Matrix &) const;	///< substract two matrices element by element and store the result in a new matrix
	Matrix operator*(const Matrix &) const;	///< multiply two matrices element by element and store the result in a new matrix
	Matrix operator/(const Matrix &) const;	///< divide two matrices element by element and store the result in a new matrix

	inline double *operator[] (const unsigned int i) {
	    return _pMatrix[i];
	};
	inline const double *operator[] (const unsigned int i) const {
	    return _pMatrix[i];
	};
    };



    Vector rowProduct(const Matrix & A, const Vector & U);
    Vector colProduct(const Vector & U, const Matrix & A);



    class SVD {
      public:

	SVD(const Matrix &, bool bU = true, bool bV = true);

	Vector getSingularValues() {
	    return _S;
	};
	Matrix getSingularMatrix();

	Matrix getU() {
	    return _U;
	};
	Matrix getV() {
	    return _V;
	};

	double norm2() {
	    return _S[0];
	};

	double cond() {
	    return _S[0] / _S[_S.size() - 1];
	};

	int rank();

      private:

	int _m;			///< number of rows
	int _n;			///< number of columns
	Matrix _U;		///< Left singular vectors
	Matrix _V;		///< Right singular vectors
	Vector _S;		///< Singular values

	bool _computeV;		///< Check if V should be computed
	bool _computeU;		///< Check if U should be computed

    };



    double randD(double, double);



};				// end of namespace SiMath



#endif
