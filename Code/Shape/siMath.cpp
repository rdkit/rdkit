/*******************************************************************************
siMath.cpp - Shape-it
 
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



#include <Shape/siMath.h>



using namespace SiMath;



Vector::Vector(const unsigned int n, const double *v):_n(n), _pVector(n)
{
    for (unsigned int i = 0; i < _n; ++i)
	_pVector[i] = v[i];
}



Vector::Vector(const std::vector < double >&v):_n(v.size()), _pVector(_n)
{
    for (unsigned int i = 0; i < _n; ++i)
	_pVector[i] = v[i];
}



Vector::Vector(const Vector & v):_n(v._n), _pVector(_n)
{
    for (unsigned int i = 0; i < _n; ++i)
	_pVector[i] = v._pVector[i];
}



Vector::~Vector()
{
    _pVector.clear();
}



void
 Vector::clear()
{
    _pVector.clear();
    _n = 0;
}



void Vector::reset(unsigned int n)
{
    if (_n != n)		// only reset the vector itself if the new size is larger
	_pVector.resize(n);
    _n = n;
    for (unsigned int i = 0; i < _n; ++i)
	_pVector[i] = 0;
}



void Vector::resize(unsigned int n)
{
    if (_n != n)
	_pVector.resize(n);
    _n = n;
}




double Vector::getValueAt(const unsigned int i)
{
    return _pVector[i];
}



double Vector::getValueAt(const unsigned int i) const
{
    return _pVector[i];
}



double Vector::max() const
{
    double d = _pVector[0];
    for (unsigned int i = 1; i < _n; ++i) {
	if (_pVector[i] > d) {
	    d = _pVector[i];
	}
    }
    return d;
}



double Vector::max(unsigned int &index) const
{
    double d = _pVector[0];
    for (unsigned int i = 1; i < _n; ++i) {
	if (_pVector[i] > d) {
	    d = _pVector[i];
	    index = i;
	}
    }
    return d;
}



double Vector::min() const
{
    double d = _pVector[0];
    for (unsigned int i = 1; i < _n; ++i) {
	if (_pVector[i] < d) {
	    d = _pVector[i];
	}
    }
    return d;
}



double Vector::min(unsigned int &index) const
{
    double d = _pVector[0];
    for (unsigned int i = 1; i < _n; ++i) {
	if (_pVector[i] > d) {
	    d = _pVector[i];
	    index = i;
	}
    }
    return d;
}



double Vector::sum() const
{
    double m(0.0);
    for (unsigned int i = 0; i < _n; ++i)
	m += _pVector[i];
    return m;
}



double Vector::mean() const
{
    double m(0.0);
    for (unsigned int i = 0; i < _n; ++i)
	m += _pVector[i];
    return m / _n;
}



double Vector::stDev() const
{
    double m(0.0);
    for (unsigned int i = 0; i < _n; ++i)
	m += _pVector[i];
    double s(0.0);
    for (unsigned int i = 0; i < _n; ++i)
	s += (m - _pVector[i]) * (m - _pVector[i]);
    return sqrt(s / (_n - 1));
}



double Vector::stDev(double m) const
{
    double s(0.0);
    for (unsigned int i = 0; i < _n; ++i)
	s += (m - _pVector[i]) * (m - _pVector[i]);
    return sqrt(s / (_n - 1));
}



Vector & Vector::operator=(const Vector & src)
{
    if (_n != src._n) {
	_n = src._n;
	_pVector.resize(_n);
    }
    for (unsigned int i = 0; i < _n; ++i)
	_pVector[i] = src._pVector[i];
    return *this;
}



Vector & Vector::operator=(const double &v)
{
    for (unsigned int i = 0; i < _n; ++i)
	_pVector[i] = v;
    return *this;
}



Vector & Vector::operator+=(const double &v)
{
    for (unsigned int i = 0; i < _n; ++i)
	_pVector[i] += v;
    return *this;
}



Vector & Vector::operator+=(const Vector & V)
{
    for (unsigned int i = 0; i < _n; ++i)
	_pVector[i] += V._pVector[i];
    return *this;
}



Vector & Vector::operator-=(const double &v)
{
    for (unsigned int i = 0; i < _n; ++i)
	_pVector[i] -= v;
    return *this;
}



Vector & Vector::operator-=(const Vector & V)
{
    for (unsigned int i = 0; i < _n; ++i)
	_pVector[i] -= V._pVector[i];
    return *this;
}



Vector & Vector::operator*=(const double &v)
{
    for (unsigned int i = 0; i < _n; ++i)
	_pVector[i] *= v;
    return *this;
}



Vector & Vector::operator*=(const Vector & V)
{
    for (unsigned int i = 0; i < _n; ++i)
	_pVector[i] *= V._pVector[i];
    return *this;
}



Vector & Vector::operator/=(const double &v)
{
    for (unsigned int i = 0; i < _n; ++i)
	_pVector[i] /= v;
    return *this;
}



Vector & Vector::operator/=(const Vector & V)
{
    for (unsigned int i = 0; i < _n; ++i)
	_pVector[i] /= V._pVector[i];
    return *this;
}



Vector & Vector::operator-()
{
    for (unsigned int i = 0; i < _n; ++i)
	_pVector[i] = -_pVector[i];
    return *this;
}



Vector Vector::operator+(const Vector & V) const
{
    Vector r(_n);
    for (unsigned int i = 0; i < _n; ++i)
	r[i] = _pVector[i] + V._pVector[i];
    return r;
}



Vector Vector::operator-(const Vector & V) const
{
    Vector r(_n);
    for (unsigned int i = 0; i < _n; ++i)
	r[i] = _pVector[i] - V._pVector[i];
    return r;
}



Vector Vector::operator*(const Vector & V) const
{
    Vector r(_n);
    for (unsigned int i = 0; i < _n; ++i)
	r[i] = _pVector[i] * V._pVector[i];
    return r;
}



Vector Vector::operator/(const Vector & V) const
{
    Vector r(_n);
    for (unsigned int i = 0; i < _n; ++i)
	r[i] = _pVector[i] / V._pVector[i];
    return r;
}



bool Vector::operator==(const Vector & V) const
{
    for (unsigned int i = 0; i < _n; ++i) {
	if (_pVector[i] != V._pVector[i])
	    return false;
    }
    return true;
}



bool Vector::operator!=(const Vector & V) const
{
    for (unsigned int i = 0; i < _n; ++i) {
	if (_pVector[i] != V._pVector[i])
	    return true;
    }
    return false;
}



double Vector::dotProd(const Vector & v)
{
    double d(0.0);
    for (unsigned int i = 0; i < _n; ++i) {
	d += _pVector[i] * v[i];
    }
    return d;
}



void Vector::swap(const unsigned int i, const unsigned int j)
{
    double dummy = _pVector[i];
    _pVector[i] = _pVector[j];
    _pVector[j] = dummy;
    return;
}



Matrix::Matrix(const unsigned int n, const unsigned int m):
_nRows(n), _nCols(m), _pMatrix(0)
{
    if (n && m) {
	double *dummy = new double[n * m];	// data
	_pMatrix = new double *[n];	// row pointers
	for (unsigned int i = 0; i < n; ++i) {
	    _pMatrix[i] = dummy;
	    dummy += m;
	}
    }
}



Matrix::Matrix(const unsigned int n, const unsigned int m,
	       const double &v):_nRows(n), _nCols(m), _pMatrix(0)
{
    if (n && m) {
	double *dummy = new double[n * m];
	_pMatrix = new double *[n];
	for (unsigned int i = 0; i < n; ++i) {
	    _pMatrix[i] = dummy;
	    dummy += m;
	}
	for (unsigned int i = 0; i < n; ++i)
	    for (unsigned int j = 0; j < m; ++j)
		_pMatrix[i][j] = v;
    }
}



Matrix::Matrix(const unsigned int n, const unsigned int m,
	       const Vector & vec):_nRows(n), _nCols(m), _pMatrix(0)
{
    double *dummy(new double[n * m]);
    _pMatrix = new double *[n];
    for (unsigned int i = 0; i < n; ++i) {
	_pMatrix[i] = dummy;
	dummy += m;
    }
    for (unsigned int i = 0; i < n; ++i) {
	for (unsigned int j = 0; j < m; ++j) {
	    _pMatrix[i][j] = vec[i * m + j];
	}
    }
}



Matrix::Matrix(const Matrix & src):_nRows(src._nRows),
_nCols(src._nCols), _pMatrix(0)
{
    if (_nRows && _nCols) {
	double *dummy(new double[_nRows * _nCols]);
	_pMatrix = new double *[_nRows];
	for (unsigned int i = 0; i < _nRows; ++i) {
	    _pMatrix[i] = dummy;
	    dummy += _nCols;
	}
	for (unsigned int i = 0; i < _nRows; ++i)
	    for (unsigned int j = 0; j < _nCols; ++j)
		_pMatrix[i][j] = src[i][j];
    }
}



Matrix::~Matrix()
{
    if (_pMatrix != NULL) {
	if (_pMatrix[0] != NULL)
	    delete[](_pMatrix[0]);
	delete[](_pMatrix);
    }
    _pMatrix = NULL;
}



double
 Matrix::getValueAt(const unsigned int i, const unsigned int j)
{
    return _pMatrix[i][j];
}



const double Matrix::getValueAt(const unsigned int i, const unsigned int j) const
{
    return _pMatrix[i][j];
}



Vector Matrix::getRow(const unsigned int i) const
{
    Vector v(_nCols);
    for (unsigned int j = 0; j < _nCols; ++j)
	v[j] = _pMatrix[i][j];
    return v;
}



Vector Matrix::getColumn(const unsigned int i) const
{
    Vector v(_nRows);
    for (unsigned int j = 0; j < _nRows; ++j)
	v[j] = _pMatrix[j][i];

    return v;
}



inline void
    Matrix::setValueAt(const unsigned int i, const unsigned int j,
		       double v)
{
    _pMatrix[i][j] = v;
}



void
 Matrix::setRow(const unsigned int i, Vector & src)
{
    for (unsigned int j = 0; j < _nCols; ++j)
	_pMatrix[i][j] = src[j];
}



void Matrix::setColumn(const unsigned int i, Vector & src)
{
    for (unsigned int j = 0; j < _nRows; ++j)
	_pMatrix[j][i] = src[j];
}



Matrix & Matrix::operator=(const Matrix & M)
{
    // check dimensions
    if (_nRows != M.nbrRows() || _nCols != M.nbrColumns()) {
	if (_nRows && _pMatrix != 0) {
	    // delete old matrix
	    if (_nCols && _pMatrix[0] != NULL)
		delete[]_pMatrix[0];
	    delete[]_pMatrix;
	}
	_pMatrix = NULL;
	// create a new matrix
	_nRows = M.nbrRows();
	_nCols = M.nbrColumns();
	_pMatrix = new double *[_nRows];
	_pMatrix[0] = new double[_nRows * _nCols];
	for (unsigned int i = 1; i < _nRows; ++i)
	    _pMatrix[i] = _pMatrix[i - 1] + _nCols;
    }
    // fill in all new values       
    for (unsigned int i = 0; i < _nRows; ++i)
	for (unsigned int j = 0; j < _nCols; ++j)
	    _pMatrix[i][j] = M[i][j];
    return *this;
}



Matrix & Matrix::operator=(const double &v)
{
    for (unsigned int i = 0; i < _nRows; ++i)
	for (unsigned int j = 0; j < _nCols; ++j)
	    _pMatrix[i][j] = v;
    return *this;
}



Matrix & Matrix::operator+=(const double &v)
{
    for (int i = 0; i < _nRows; i++)
	for (int j = 0; j < _nCols; j++)
	    _pMatrix[i][j] += v;
    return *this;
}



Matrix & Matrix::operator+=(const Matrix & M)
{
    for (unsigned int i = 0; i < _nRows; ++i)
	for (unsigned int j = 0; j < _nCols; ++j)
	    _pMatrix[i][j] += M[i][j];
    return *this;
}



Matrix & Matrix::operator-=(const double &v)
{
    for (int i = 0; i < _nRows; i++)
	for (int j = 0; j < _nCols; j++)
	    _pMatrix[i][j] -= v;
    return *this;
}



Matrix & Matrix::operator-=(const Matrix & M)
{
    for (unsigned int i = 0; i < _nRows; ++i)
	for (unsigned int j = 0; j < _nCols; ++j)
	    _pMatrix[i][j] -= M[i][j];
    return *this;
}



Matrix & Matrix::operator*=(const double &v)
{
    for (unsigned int i = 0; i < _nRows; ++i)
	for (unsigned int j = 0; j < _nCols; ++j)
	    _pMatrix[i][j] *= v;
    return *this;
}



Matrix & Matrix::operator*=(const Matrix & M)
{
    for (unsigned int i = 0; i < _nRows; ++i)
	for (unsigned int j = 0; j < _nCols; ++j)
	    _pMatrix[i][j] *= M[i][j];
    return *this;
}



Matrix & Matrix::operator/=(const double &v)
{
    for (unsigned int i = 0; i < _nRows; ++i)
	for (unsigned int j = 0; j < _nCols; ++j)
	    _pMatrix[i][j] /= v;
    return *this;
}



Matrix & Matrix::operator/=(const Matrix & M)
{
    for (unsigned int i = 0; i < _nRows; ++i)
	for (unsigned int j = 0; j < _nCols; ++j)
	    _pMatrix[i][j] /= M[i][j];
    return *this;
}



Matrix & Matrix::operator-()
{
    for (unsigned int i = 0; i < _nRows; ++i)
	for (unsigned int j = 0; j < _nCols; ++j)
	    _pMatrix[i][j] = -_pMatrix[i][j];

    return *this;
}



Matrix Matrix::operator+(const Matrix & M) const
{
    Matrix B(M);
    for (unsigned int i = 0; i < _nRows; ++i)
	for (unsigned int j = 0; j < _nCols; ++j)
	    B[i][j] = _pMatrix[i][j] + M[i][j];
    return B;
}



Matrix Matrix::operator-(const Matrix & M) const
{
    Matrix B(M);
    for (unsigned int i = 0; i < _nRows; ++i)
	for (unsigned int j = 0; j < _nCols; ++j)
	    B[i][j] = _pMatrix[i][j] - M[i][j];
    return B;
}



Matrix Matrix::operator*(const Matrix & M) const
{
    Matrix B(M);
    for (unsigned int i = 0; i < _nRows; ++i)
	for (unsigned int j = 0; j < _nCols; ++j)
	    B[i][j] = _pMatrix[i][j] * M[i][j];
    return B;
}



Matrix Matrix::operator/(const Matrix & M) const
{
    Matrix B(M);
    for (unsigned int i = 0; i < _nRows; ++i)
	for (unsigned int j = 0; j < _nCols; ++j)
	    B[i][j] = _pMatrix[i][j] / M[i][j];
    return B;
}



void Matrix::swapRows(unsigned int i, unsigned int j)
{
    double dummy;
    for (unsigned int k = 0; k < _nCols; ++k)	// loop over all columns
    {
	dummy = _pMatrix[i][k];	// store original element at [i,k]
	_pMatrix[i][k] = _pMatrix[j][k];	// replace [i,k] with [j,k]
	_pMatrix[j][k] = dummy;	// replace [j,k] with element originally at [i,k] 
    }
    return;
}



void Matrix::swapColumns(unsigned int i, unsigned int j)
{
    double dummy;
    for (unsigned int k = 0; k < _nRows; ++k)	// loop over all rows
    {
	dummy = _pMatrix[k][i];	// store original element at [k,i]
	_pMatrix[k][i] = _pMatrix[k][j];	// replace [k,i] with [k,j]
	_pMatrix[k][j] = dummy;	// replace [k,j] with element orignally at [k,i]
    }
    return;
}



void Matrix::reset(const unsigned int r, const unsigned int c)
{
    // check dimensions
    if (_nRows != r || _nCols != c) {
	if (_nRows != 0 && _nCols != 0 && _pMatrix != 0) {
	    // delete old matrix
	    if (_pMatrix[0] != NULL)
		delete[]_pMatrix[0];
	    delete[]_pMatrix;
	}
	// create a new matrix
	_nRows = r;
	_nCols = c;
	if (_nRows == 0 || _nCols == 0) {
	    _pMatrix = NULL;
	    return;
	}
	_pMatrix = new double *[_nRows];
	_pMatrix[0] = new double[_nRows * _nCols];
	for (unsigned int i = 1; i < _nRows; ++i)
	    _pMatrix[i] = _pMatrix[i - 1] + _nCols;
    }
    // fill in all new values       
    for (unsigned int i = 0; i < _nRows; ++i)
	for (unsigned int j = 0; j < _nCols; ++j)
	    _pMatrix[i][j] = 0;

}



void Matrix::clear()
{
    // delete old matrix
    if (_pMatrix != NULL) {
	if (_pMatrix[0] != NULL)
	    delete[]_pMatrix[0];
	delete[]_pMatrix;
    }
    _pMatrix = NULL;
    _nRows = 0;
    _nCols = 0;
}



Matrix Matrix::transpose(void)
{
    Matrix T(_nCols, _nRows);
    for (unsigned int i(0); i < _nRows; ++i) {
	for (unsigned int j(0); j < _nCols; ++j) {
	    T[j][i] = _pMatrix[i][j];
	}
    }
    return T;
}



SiMath::Vector
    SiMath::rowProduct(const SiMath::Matrix & A, const SiMath::Vector & U)
{
    Vector v(A.nbrRows(), 0.0);

    for (unsigned int i = 0; i < A.nbrRows(); ++i) {
	double s(0.0);
	for (unsigned int j = 0; j < A.nbrColumns(); ++j) {
	    s += A[i][j] * U[j];
	}
	v[i] = s;
    }
    return v;
}



SiMath::Vector
    SiMath::colProduct(const SiMath::Vector & U, const SiMath::Matrix & A)
{
    Vector v(A.nbrColumns(), 0.0);
    for (unsigned int i = 0; i < A.nbrColumns(); ++i) {
	double s(0.0);
	for (unsigned int j = 0; j < A.nbrRows(); ++j) {
	    s += U[j] * A[j][i];
	}
	v[i] = s;
    }
    return v;
}

double SiMath::randD(double a, double b)
{
    double d(a);
    d += (b - a) * ((double) rand() / RAND_MAX);
    return d;
}
