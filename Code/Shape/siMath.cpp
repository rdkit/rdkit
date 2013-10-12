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



Vector::Vector(const unsigned int n, const double * v) :
	_n(n), 
	_pVector(n)
{
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] = v[i];
}



Vector::Vector(const std::vector<double> & v) :
	_n(v.size()), 
	_pVector(_n)
{
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] = v[i];
}



Vector::Vector(const Vector & v) : 
	_n(v._n), 
	_pVector(_n)
{
	for ( unsigned int i=0; i<_n; ++i)
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



void
Vector::reset(unsigned int n)
{
	if ( _n !=  n ) // only reset the vector itself if the new size is larger
		_pVector.resize(n);
	_n = n;
	for ( unsigned int i=0; i<_n; ++i )
		_pVector[i] = 0;
}



void
Vector::resize(unsigned int n)
{
	if ( _n !=  n )
		_pVector.resize(n);
	_n = n;
}




double
Vector::getValueAt(const unsigned int i)
{
	return _pVector[i];
}

 	
   
double 
Vector::getValueAt(const unsigned int i) const
{
   return _pVector[i];
}


 	
double 
Vector::max() const
{
	double d = _pVector[0];
	for ( unsigned int i=1; i<_n; ++i)
	{
		if ( _pVector[i] > d )
		{
			d = _pVector[i];
		}
	}
	return d;
}



double
Vector::max(unsigned int & index) const
{
	double d = _pVector[0];
	for ( unsigned int i=1; i<_n; ++i)
	{
		if ( _pVector[i] > d )
		{
			d = _pVector[i];
			index = i;
		}
	}
	return d;
}



double 
Vector::min() const
{
	double d = _pVector[0];
	for ( unsigned int i=1; i<_n; ++i)
	{
		if ( _pVector[i] < d )
		{
			d = _pVector[i];
		}
	}
	return d;
}



double
Vector::min(unsigned int & index) const 
{
	double d = _pVector[0];
	for ( unsigned int i=1; i<_n; ++i)
	{
		if ( _pVector[i] > d )
		{
			d = _pVector[i];
			index = i;
		}
	}
	return d;
}



double 
Vector::sum() const
{
	double m(0.0);
	for ( unsigned int i=0; i<_n; ++i)
		m += _pVector[i];
	return m;
}



double 
Vector::mean() const
{
	double m(0.0);
	for ( unsigned int i=0; i<_n; ++i)
		m += _pVector[i];
	return m/_n;
}



double 
Vector::stDev() const
{
	double m(0.0);
	for ( unsigned int i=0; i<_n; ++i)
		m += _pVector[i];
	double s(0.0);
	for ( unsigned int i=0; i<_n; ++i)
		s += (m-_pVector[i])*(m-_pVector[i]);
	return sqrt(s/(_n-1));
}



double
Vector::stDev(double m) const
{
	double s(0.0);
	for ( unsigned int i=0; i<_n; ++i)
		s += (m-_pVector[i])*(m-_pVector[i]);
	return sqrt(s/(_n-1));
}



Vector & 
Vector::operator=(const Vector & src)
{
	if ( _n != src._n )
	{	
		_n = src._n;
		_pVector.resize(_n);
	}
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] = src._pVector[i];
	return *this;
}



Vector & 
Vector::operator= (const double & v)
{
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] = v;
	return *this;
}



Vector &
Vector::operator+= (const double & v)
{
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] += v;
	return *this;
}



Vector & 
Vector::operator+= (const Vector & V)
{
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] += V._pVector[i];
	return *this;
}



Vector &
Vector::operator-= (const double & v)
{
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] -= v;
	return *this;
}



Vector & 
Vector::operator-= (const Vector & V)
{
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] -= V._pVector[i];
	return *this;
}



Vector &
Vector::operator*= (const double & v)
{
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] *= v;
	return *this;
}



Vector & 
Vector::operator*= (const Vector & V)
{
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] *= V._pVector[i];
	return *this;
}



Vector &
Vector::operator/= (const double & v)
{
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] /= v;
	return *this;
}



Vector &
Vector::operator/= (const Vector & V)
{
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] /= V._pVector[i];
	return *this;
}

 
   
Vector &
Vector::operator- ()
{
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] = -_pVector[i];
	return *this;
}


  
Vector 
Vector::operator+ (const Vector & V) const
{
	Vector r(_n);
	for ( unsigned int i=0; i<_n; ++i)
		r[i] = _pVector[i] + V._pVector[i];
	return r;
}



Vector 
Vector::operator- (const Vector & V) const
{
	Vector r(_n);
	for ( unsigned int i=0; i<_n; ++i)
		r[i] = _pVector[i] - V._pVector[i];
	return r;
}



Vector 
Vector::operator* (const Vector & V) const
{
	Vector r(_n);
	for ( unsigned int i=0; i<_n; ++i)
		r[i] = _pVector[i] * V._pVector[i];
	return r;
}



Vector 
Vector::operator/ (const Vector & V) const
{
	Vector r(_n);
	for ( unsigned int i=0; i<_n; ++i)
		r[i] = _pVector[i] / V._pVector[i];
	return r;
}



bool  
Vector::operator== (const Vector & V) const
{
	for ( unsigned int i=0; i<_n; ++i)
	{	
		if ( _pVector[i] != V._pVector[i] ) return false;		
	}
	return true;	
}



bool  
Vector::operator!= (const Vector & V) const
{
	for ( unsigned int i=0; i<_n; ++i)
	{	
		if ( _pVector[i] != V._pVector[i] ) return true;
	}
	return false;	
}



double 
Vector::dotProd(const Vector & v){
	double d(0.0);
	for (unsigned int i=0; i<_n; ++i )
	{
		d += _pVector[i] * v[i];
	}
	return d;
}



void 
Vector::swap(const unsigned int i, const unsigned int j)
{
   double dummy = _pVector[i];
   _pVector[i] = _pVector[j];
   _pVector[j] = dummy;
   return;
}



Matrix::Matrix(const unsigned int n, const unsigned int m) :
	_nRows(n),
	_nCols(m),
	_pMatrix(0)
{
	if ( n && m )
	{
		double* dummy = new double[n*m];  // data
		_pMatrix = new double*[n];             // row pointers
		for ( unsigned int i=0; i<n; ++i)
		{
			_pMatrix[i] = dummy;
			dummy += m;
		}
	}
}



Matrix::Matrix(const unsigned int n, const unsigned int m, const double & v) :
	_nRows(n),
	_nCols(m),
	_pMatrix(0)
{
	if ( n && m )
	{
		double * dummy = new double[n*m];
		_pMatrix = new double*[n];
		for ( unsigned int i=0; i<n; ++i)
		{
			_pMatrix[i] = dummy;
			dummy += m;
		}
		for ( unsigned int i=0; i<n; ++i)
			for ( unsigned int j=0; j<m; ++j )
				_pMatrix[i][j] = v; 
	}
}



Matrix::Matrix(const unsigned int n, const unsigned int m, const Vector & vec) :
	_nRows(n),
	_nCols(m),
	_pMatrix(0)
{
	double * dummy(new double[n*m]);
	_pMatrix = new double*[n];
	for(unsigned int i=0 ; i<n ; ++i) {
		_pMatrix[i] = dummy;
		dummy += m;
	}
	for(unsigned int i=0 ; i<n ; ++i) {
		for(unsigned int j=0 ; j<m ; ++j) {
			_pMatrix[i][j] = vec[i*m+j];
		}
	}
}



Matrix::Matrix(const Matrix & src) : 
	_nRows(src._nRows),
	_nCols(src._nCols),
	_pMatrix(0)
{
	if ( _nRows && _nCols )
	{
		double * dummy(new double[_nRows * _nCols]);
		_pMatrix = new double*[_nRows];
		for ( unsigned int i=0; i<_nRows; ++i)
		{
			_pMatrix[i] = dummy;
			dummy += _nCols;
		}
		for ( unsigned int i=0; i<_nRows; ++i)
			for ( unsigned int j=0; j<_nCols; ++j )
				_pMatrix[i][j] = src[i][j]; 
	}
}



Matrix::~Matrix()
{
	if ( _pMatrix != NULL )
	{
		if ( _pMatrix[0] != NULL )
			delete[] (_pMatrix[0]);
		delete[] (_pMatrix);
	}
	_pMatrix = NULL;
}



double
Matrix::getValueAt(const unsigned int i, const unsigned int j)
{
	return _pMatrix[i][j];
}



const double 
Matrix::getValueAt(const unsigned int i, const unsigned int j) const
{
	return _pMatrix[i][j];
}



Vector
Matrix::getRow(const unsigned int i) const
{
	Vector v(_nCols); 
	for ( unsigned int j=0; j<_nCols; ++j )
		v[j] = _pMatrix[i][j];
	return v;
}



Vector
Matrix::getColumn(const unsigned int i) const
{
	Vector v(_nRows);
	for ( unsigned int j=0; j<_nRows; ++j )
		v[j] = _pMatrix[j][i];
	
	return v;
}



inline void
Matrix::setValueAt(const unsigned int i, const unsigned int j, double v)
{
	_pMatrix[i][j] = v;
}



void 
Matrix::setRow(const unsigned int i,Vector & src)
{
	for ( unsigned int j=0; j<_nCols; ++j )
		_pMatrix[i][j] = src[j];
}



void 
Matrix::setColumn(const unsigned int i, Vector & src)
{
	for ( unsigned int j=0; j<_nRows; ++j )
		_pMatrix[j][i] = src[j];
}



Matrix &
Matrix::operator= (const Matrix & M)
{
	// check dimensions
	if ( _nRows != M.nbrRows() ||	_nCols != M.nbrColumns() )
	{
		if ( _nRows && _pMatrix != 0 )
		{
			// delete old matrix
			if ( _nCols && _pMatrix[0] != NULL )
				delete[] _pMatrix[0];
			delete[] _pMatrix;
		}
		_pMatrix = NULL;
		// create a new matrix
		_nRows = M.nbrRows();
		_nCols = M.nbrColumns();
		_pMatrix = new double*[_nRows];
		_pMatrix[0] = new double[_nRows*_nCols];
		for (unsigned int i=1; i<_nRows; ++i)
			_pMatrix[i] = _pMatrix[i-1] + _nCols;
	}
	// fill in all new values	
	for (unsigned int i=0; i<_nRows; ++i)
		for (unsigned int j=0; j<_nCols; ++j)
			_pMatrix[i][j] = M[i][j];
	return *this;
}



Matrix & 
Matrix::operator= (const double & v) 
{
	for ( unsigned int i=0; i<_nRows; ++i)
		for ( unsigned int j=0; j<_nCols; ++j)
			_pMatrix[i][j] = v;
	return *this;
}



Matrix & 
Matrix::operator+= (const double & v)
{
	for ( int i=0; i<_nRows; i++)
				for ( int j=0; j<_nCols; j++)
					_pMatrix[i][j] += v;
	return *this;
}



Matrix & 
Matrix::operator+= (const Matrix & M)
{
	for ( unsigned int i=0; i<_nRows; ++i)
		for ( unsigned int j=0; j<_nCols; ++j)
			_pMatrix[i][j] += M[i][j];
	return *this;
}



Matrix & 
Matrix::operator-= (const double & v)
{
	for ( int i=0; i<_nRows; i++)
		for ( int j=0; j<_nCols; j++)
			_pMatrix[i][j] -= v;
	return *this;
}



Matrix & 
Matrix::operator-= (const Matrix & M)
{
	for ( unsigned int i=0; i<_nRows; ++i)
		for ( unsigned int j=0; j<_nCols; ++j)
			_pMatrix[i][j] -= M[i][j];
	return *this;
}



Matrix& 
Matrix::operator*= (const double & v)
{
	for ( unsigned int i=0; i<_nRows; ++i )
				for ( unsigned int j=0; j<_nCols; ++j )
					_pMatrix[i][j] *= v;
	return *this;
}



Matrix&
Matrix::operator*= (const Matrix & M)
{
	for ( unsigned int i=0; i<_nRows; ++i)
		for ( unsigned int j=0; j<_nCols; ++j)
			_pMatrix[i][j] *= M[i][j];
	return *this;
}



Matrix& 
Matrix::operator/= (const double & v)
{
	for ( unsigned int i=0; i<_nRows; ++i )
		for ( unsigned int j=0; j<_nCols; ++j )
			_pMatrix[i][j] /= v;
	return *this;
}



Matrix & 
Matrix::operator/= (const Matrix & M)
{
	for ( unsigned int i=0; i<_nRows; ++i)
		for ( unsigned int j=0; j<_nCols; ++j)
			_pMatrix[i][j] /= M[i][j];
	return *this;
}



Matrix & 
Matrix::operator- ()
{
	for ( unsigned int i=0; i<_nRows; ++i)
		for ( unsigned int j=0; j<_nCols; ++j)
			_pMatrix[i][j] = -_pMatrix[i][j];
	
	return *this;
}



Matrix
Matrix::operator+ (const Matrix & M) const
{
	Matrix B(M);
	for ( unsigned int i=0; i<_nRows; ++i)
		for ( unsigned int j=0; j<_nCols; ++j)
			B[i][j] = _pMatrix[i][j] + M[i][j];
	return B;
}



Matrix
Matrix::operator- (const Matrix & M) const
{
	Matrix B(M);
	for ( unsigned int i=0; i<_nRows; ++i)
		for ( unsigned int j=0; j<_nCols; ++j)
			B[i][j] = _pMatrix[i][j] - M[i][j];
	return B;
}



Matrix
Matrix::operator* (const Matrix & M) const
{
	Matrix B(M);
	for ( unsigned int i=0; i<_nRows; ++i)
		for ( unsigned int j=0; j<_nCols; ++j)
			B[i][j] = _pMatrix[i][j] * M[i][j];
	return B;
}



Matrix
Matrix::operator/ (const Matrix & M) const
{
	Matrix B(M);
	for ( unsigned int i=0; i<_nRows; ++i)
		for ( unsigned int j=0; j<_nCols; ++j)
			B[i][j] = _pMatrix[i][j] / M[i][j];
	return B;
}



void 
Matrix::swapRows(unsigned int i, unsigned int j)
{
	double dummy;
	for (unsigned int k=0; k<_nCols; ++k)         // loop over all columns
	{
		dummy = _pMatrix[i][k];            // store original element at [i,k]
		_pMatrix[i][k] = _pMatrix[j][k];   // replace [i,k] with [j,k]
		_pMatrix[j][k] = dummy;            // replace [j,k] with element originally at [i,k] 
	}
	return;
}



void 
Matrix::swapColumns(unsigned int i, unsigned int j)
{
	double dummy;
	for (unsigned int k=0; k<_nRows; ++k)         // loop over all rows
	{
		dummy = _pMatrix[k][i];            // store original element at [k,i]
		_pMatrix[k][i] = _pMatrix[k][j];   // replace [k,i] with [k,j]
		_pMatrix[k][j] = dummy;            // replace [k,j] with element orignally at [k,i]
	}
	return;
}



void
Matrix::reset (const unsigned int r, const unsigned int c)
{
	// check dimensions
	if ( _nRows != r ||	_nCols != c )
	{
		if ( _nRows != 0 && _nCols != 0 && _pMatrix != 0 )
		{
			// delete old matrix
			if ( _pMatrix[0] != NULL )
				delete[] _pMatrix[0];
			delete[] _pMatrix;
		}
		// create a new matrix
		_nRows = r;
		_nCols = c;
		if ( _nRows == 0 || _nCols == 0 )
		{
			_pMatrix = NULL;
			return;
		}
		_pMatrix = new double*[_nRows];
		_pMatrix[0] = new double[_nRows*_nCols];
		for ( unsigned int i=1; i<_nRows; ++i)
			_pMatrix[i] = _pMatrix[i-1] + _nCols;
	}
	// fill in all new values	
	for ( unsigned int i=0; i<_nRows; ++i)
		for ( unsigned int j=0; j<_nCols; ++j)
			_pMatrix[i][j] = 0;
	
}



void 
Matrix::clear()
{
	// delete old matrix
	if ( _pMatrix != NULL )
	{
		if ( _pMatrix[0] != NULL )
			delete[] _pMatrix[0];
		delete[] _pMatrix;
	}
	_pMatrix = NULL;
	_nRows = 0;
	_nCols = 0;
}



Matrix 
Matrix::transpose(void)
{	
	Matrix T(_nCols, _nRows);
	for (unsigned int i(0); i < _nRows; ++i)
	{
		for (unsigned int j(0); j < _nCols; ++j)
		{
         T[j][i] = _pMatrix[i][j]; 
		}
	}
	return T;
}



SiMath::Vector
SiMath::rowProduct(const SiMath::Matrix & A, const SiMath::Vector & U)
{
	Vector v(A.nbrRows(),0.0);
	
	for ( unsigned int i=0; i<A.nbrRows(); ++i)
	{
		double s(0.0);
		for ( unsigned int j=0; j<A.nbrColumns(); ++j)
		{
			s += A[i][j] * U[j];
		}
		v[i] = s;
	}
	return v;  
}



SiMath::Vector
SiMath::colProduct(const SiMath::Vector & U, const SiMath::Matrix & A)
{
	Vector v(A.nbrColumns(),0.0);
	for ( unsigned int i=0; i<A.nbrColumns(); ++i)
	{
		double s(0.0);
		for ( unsigned int j=0; j<A.nbrRows(); ++j)
		{
			s += U[j] * A[j][i];
		}
		v[i] = s;
	}
	return v;  
}



SVD::SVD(const Matrix& Aorig, bool bU, bool bV) :
	_m(Aorig.nbrRows()),
	_n(Aorig.nbrColumns()),
	_U(),
	_V(),
	_S(0),
	_computeV(bV),
	_computeU(bU)
{
	// dimensionality of the problem
	int nu = min(_m,_n);
	int nct = min(_m-1, _n);
	int nrt = max(0, std::min(_n-2,_m));

	// define the dimensions of the internal matrices and vetors
	_S.reset(min(_m+1, _n)); 
	
	if ( _computeU )
		_U.reset(_m,nu);

	if ( _computeV )
		_V.reset(_n,_n);

	// local working vectors
	Vector e(_n);
	Vector work(_m);
	
	// make a copy of A to do the computations on
	Matrix Acopy(Aorig);

	// loop indices
	int i=0, j=0, k=0;
	
	// Reduce A to bidiagonal form, storing the diagonal elements
	// in _S and the super-diagonal elements in e.
	
	for (k = 0; k < max(nct,nrt); k++)
	{
		if (k < nct) 
		{
			// Compute the transformation for the k-th column and place the k-th diagonal in _S[k].
			_S[k] = 0;
			for (i = k; i < _m; i++) {
				_S[k] = triangle(_S[k], Acopy[i][k]);
			}
			if (_S[k] != 0.0) {
				if (Acopy[k][k] < 0.0) {
					_S[k] = -_S[k];
				}
				for (i = k; i < _m; i++) {
					Acopy[i][k] /= _S[k];
				}
				Acopy[k][k] += 1.0;
			}
			_S[k] = -_S[k];
		}
		for (j = k+1; j < _n; j++) 
		{
			if ((k < nct) && (_S[k] != 0.0))
			{
				// Apply the transformation to Acopy
				double t = 0;
				for (i = k; i < _m; i++) {
					t += Acopy[i][k]*Acopy[i][j];
				}
				t = -t/Acopy[k][k];
				for (i = k; i < _m; i++) {
					Acopy[i][j] += t*Acopy[i][k];
				}
			}
			
			// Place the k-th row of A into e for the subsequent calculation of the row transformation.
			e[j] = Acopy[k][j];
		}
		
		// Place the transformation in _U for subsequent back multiplication.
		if ( _computeU & (k < nct) ) 
		{			
			for (i = k; i < _m; i++) 
			{
				_U[i][k] = Acopy[i][k];
			}
		}
		
		if ( k < nrt )
		{
			// Compute the k-th row transformation and place the k-th super-diagonal in e[k].
			// Compute 2-norm without under/overflow.
			e[k] = 0.0;
			for (i = k+1; i < _n; i++) {
				e[k] = triangle(e[k],e[i]);
			}
			if (e[k] != 0.0) 
			{
				if (e[k+1] < 0.0) { // switch sign
					e[k] = -e[k]; 
				}
				for (i = k+1; i < _n; i++) { // scale
					e[i] /= e[k];
				}
				e[k+1] += 1.0;
			}
			e[k] = -e[k]; 
			if ((k+1 < _m) & (e[k] != 0.0))
			{
				// Apply the transformation.
				
				for (i = k+1; i < _m; i++) {
					work[i] = 0.0;
				}
				for (j = k+1; j < _n; j++) {
					for (i = k+1; i < _m; i++) {
						work[i] += e[j]*Acopy[i][j];
					}
				}
				for (j = k+1; j < _n; j++) {
					double t = -e[j]/e[k+1];
					for (i = k+1; i < _m; i++) {
						Acopy[i][j] += t*work[i];
					}
				}
			}
			
			
			// Place the transformation in _V for subsequent back multiplication.			
			if ( _computeV ) 
			{
				for (i = k+1; i < _n; i++) 
				{
					_V[i][k] = e[i];
				}
			}
		}
	}
	
	// Set up the final bidiagonal matrix of order p.
	int p = min(_n,_m+1);
	if (nct < _n) {
		_S[nct] = Acopy[nct][nct];
	}
	if (_m < p) {
		_S[p-1] = 0.0;
	}
	if (nrt+1 < p) {
		e[nrt] = Acopy[nrt][p-1];
	}
	e[p-1] = 0.0;
	
	// If required, generate U.
	if ( _computeU ) 
	{
		for (j = nct; j < nu; j++) {
			for (i = 0; i < _m; i++) {
				_U[i][j] = 0.0;
			}
			_U[j][j] = 1.0;
		}
		for (k = nct-1; k >= 0; k--) 
		{
			if (_S[k] != 0.0) 
			{
				for (j = k+1; j < nu; j++) 
				{
					double t = 0;
					for (i = k; i < _m; i++) {
						t += _U[i][k]*_U[i][j];
					}
					t = -t/_U[k][k];
					for (i = k; i < _m; i++) {
						_U[i][j] += t*_U[i][k];
					}
				}
				for (i = k; i < _m; i++ ) {
					_U[i][k] = -_U[i][k];
				}
				_U[k][k] = 1.0 + _U[k][k];
				for (i = 0; i < k-1; i++) {
					_U[i][k] = 0.0;
				}
			} 
			else
			{
				for (i = 0; i < _m; i++) {
					_U[i][k] = 0.0;
				}
				_U[k][k] = 1.0;
			}
		}
	}
	
	// If required, generate _V.
	if ( _computeV ) 
	{
		for (k = _n-1; k >= 0; k--) 
		{
			if ((k < nrt) & (e[k] != 0.0)) 
			{
				for (j = k+1; j < nu; j++) {
					double t = 0;
					for (i = k+1; i < _n; i++) {
						t += _V[i][k]*_V[i][j];
					}
					t = -t/_V[k+1][k];
					for (i = k+1; i < _n; i++) {
						_V[i][j] += t*_V[i][k];
					}
				}
			}
			for (i = 0; i < _n; i++) {
				_V[i][k] = 0.0;
			}
			_V[k][k] = 1.0;
		}
	}
	
	// Main iteration loop for the singular values.
	int pp = p-1;
	int iter = 0;
	double eps = pow(2.0,-52.0);
	while (p > 0) {
		k=0;
		unsigned int mode=0;
		
		// Here is where a test for too many iterations would go.
		// This section of the program inspects for negligible elements in the s and e arrays.  
		// On completion the variables mode and k are set as follows.
		
		// mode = 1     if s(p) and e[k-1] are negligible and k<p
		// mode = 2     if s(k) is negligible and k<p
		// mode = 3     if e[k-1] is negligible, k<p, and s(k), ..., s(p) are not negligible (qr step).
		// mode = 4     if e(p-1) is negligible (convergence).
		for (k = p-2; k >= -1; k--)
		{
			if (k == -1) {
				break;
			}
			if (fabs(e[k]) <= eps*(fabs(_S[k]) + fabs(_S[k+1]))) {
				e[k] = 0.0;
				break;
			}
		}
		if ( k == p-2 ) 
		{
			mode = 4;
		} 
		else 
		{
			int ks(p-1); // start from ks == p-1
			for ( ; ks >= k; ks--) {
				if (ks == k) {
					break;
				}
				double t = ( (ks != p) ? fabs(e[ks]) : 0.0) + ( (ks != k+1) ? fabs(e[ks-1]) : 0.0);
				if (fabs(_S[ks]) <= eps*t)  
				{
					_S[ks] = 0.0;
					break;
				}
			}
			if (ks == k) {
				mode = 3;
			} else if (ks == p-1) {
				mode = 1;
			} else {
				mode = 2;
				k = ks;
			}
		}
		k++;
				
		// Perform the task indicated by the selected mode.		
		switch ( mode ) {
			
			case 1: 
			{ 			// Deflate negligible _S[p]
				double f = e[p-2];
				e[p-2] = 0.0;
				for (j = p-2; j >= k; j--) 
				{
					double t = SiMath::triangle(_S[j],f);
					double cs = _S[j]/t;
					double sn = f/t;
					_S[j] = t;
					if (j != k) {
						f = -sn*e[j-1];
						e[j-1] = cs*e[j-1];
					}
					
					// update V 
					if ( _computeV ) 
					{
						for (i = 0; i < _n; i++) 
						{
							t = cs*_V[i][j] + sn*_V[i][p-1];
							_V[i][p-1] = -sn*_V[i][j] + cs*_V[i][p-1];
							_V[i][j] = t;
						}
					}
				}
			}
			break; // end case 1
				
			case 2: 
			{ // Split at negligible _S[k]
				double f = e[k-1];
				e[k-1] = 0.0;
				for (j = k; j < p; j++) 
				{
					double t = triangle(_S[j],f);
					double cs = _S[j]/t;
					double sn = f/t;
					_S[j] = t;
					f = -sn*e[j];
					e[j] = cs*e[j];
					
					if ( _computeU ) 
					{
						for (i = 0; i < _m; i++) {
							t = cs*_U[i][j] + sn*_U[i][k-1];
							_U[i][k-1] = -sn*_U[i][j] + cs*_U[i][k-1];
							_U[i][j] = t;
						}
					}
				}
			}
			break; // end case 2 
				
			case 3: 
			{ // Perform one qr step.
				
				// Calculate the shift.
				double scale = max(max(max(max(fabs(_S[p-1]),fabs(_S[p-2])),fabs(e[p-2])),fabs(_S[k])),fabs(e[k]));
				double sp = _S[p-1]/scale;
				double spm1 = _S[p-2]/scale;
				double epm1 = e[p-2]/scale;
				double sk = _S[k]/scale;
				double ek = e[k]/scale;
				double b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
				double c = (sp*epm1)*(sp*epm1);
				double shift = 0.0;
				if ((b != 0.0) || (c != 0.0)) {
					shift = sqrt(b*b + c);
					if (b < 0.0) {
						shift = -shift;
					}
					shift = c/(b + shift);
				}
				double f = (sk + sp)*(sk - sp) + shift;
				double g = sk*ek;
				
				// Chase zeros.
				
				for (j = k; j < p-1; j++) 
				{
					double t = SiMath::triangle(f,g);
					double cs = f/t;
					double sn = g/t;
					if (j != k) {
						e[j-1] = t;
					}
					f = cs*_S[j] + sn*e[j];
					e[j] = cs*e[j] - sn*_S[j];
					g = sn*_S[j+1];
					_S[j+1] = cs*_S[j+1];
					
					if ( _computeV ) 
					{
						for (i = 0; i < _n; i++) 
						{
							t = cs*_V[i][j] + sn*_V[i][j+1];
							_V[i][j+1] = -sn*_V[i][j] + cs*_V[i][j+1];
							_V[i][j] = t;
						}
					}
					t = SiMath::triangle(f,g);
					cs = f/t;
					sn = g/t;
					_S[j] = t;
					f = cs*e[j] + sn*_S[j+1];
					_S[j+1] = -sn*e[j] + cs*_S[j+1];
					g = sn*e[j+1];
					e[j+1] = cs*e[j+1];
					
					if ( _computeU && (j < _m-1) ) 
					{
						for (i = 0; i < _m; i++) 
						{
							t = cs*_U[i][j] + sn*_U[i][j+1];
							_U[i][j+1] = -sn*_U[i][j] + cs*_U[i][j+1];
							_U[i][j] = t;
						}
					}
				}
				e[p-2] = f;
				iter++;
			}
				break; // end case 3
				
				// convergence step				
			case 4: 
			{ 
				
				// Make the singular values positive.
				if (_S[k] <= 0.0) 
				{
					_S[k] = (_S[k] < 0.0 ) ? -_S[k] : 0.0;

					if ( _computeV )
					{
						for (i = 0; i <= pp; i++) {
							_V[i][k] = -_V[i][k];
						}
					}
				}
				
				// Order the singular values.
				while (k < pp) 
				{
					if (_S[k] >= _S[k+1])
						break;

					// swap values and columns if necessary
					_S.swap(k,k+1);
					
					if ( _computeV && (k < _n-1)) 
						_V.swapColumns(k,k+1);
					
					if ( _computeU && (k < _m-1) )
						_U.swapColumns(k,k+1);

					k++;
				}
				iter = 0;
				p--;
			}
				break; // end case 4
		} 
	}
}


Matrix 
SVD::getSingularMatrix() 
{
	unsigned int n = _S.size();
	Matrix A(n,n,0.0);
	// set diagonal elements
	for (int i = 0; i < n; i++) 
	{
		A[i][i] = _S[i];
	}
	
	return A;
}



int 
SVD::rank() 
{
	double eps = pow(2.0,-52.0);
	double tol = max(_m,_n) * _S[0] * eps;
	int r = 0;
	for (int i = 0; i < _S.size(); i++) {
		if (_S[i] > tol) {
			r++;
		}
	}
	return r;
}



double 
SiMath::randD(double a, double b)
{
	double d(a);
	d += (b-a) * ((double)rand()/RAND_MAX);
	return d;
}


