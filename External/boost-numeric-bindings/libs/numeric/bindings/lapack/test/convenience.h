#ifndef CONVENIENCE_H
#define CONVENIENCE_H

#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <iomanip>
#include <sstream>

#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/matrix_traits.hpp>
#include <boost/numeric/bindings/traits/vector_traits.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

namespace ublas = boost::numeric::ublas;
namespace traits = boost::numeric::bindings::traits;

// single precision typedefs
typedef float freal_t;
typedef ublas::matrix<freal_t, ublas::column_major> fmat_t;
typedef ublas::vector<freal_t> fvec_t;

// double precision typedefs
typedef double dreal_t;
typedef ublas::matrix<dreal_t, ublas::column_major> dmat_t;
typedef ublas::vector<dreal_t> dvec_t;

// single precision complex typedefs
typedef std::complex<float> fcmplx_t;
typedef ublas::matrix<fcmplx_t, ublas::column_major> fcmat_t;
typedef ublas::vector<fcmplx_t> fcvec_t;

// double precision complex typedefs
typedef std::complex<double> dcmplx_t;
typedef ublas::matrix<dcmplx_t, ublas::column_major> dcmat_t;
typedef ublas::vector<dcmplx_t> dcvec_t;

// matrix/vector test sizes
const int row_size = 3;
const int col_size = 3;
const int row_range = 2;
const int col_range = 2;

//////////////////////////////////////////////////////////////////
//
//	Helper functions and structs to aid with testing 
//
/////////////////////////////////////////////////////////////////
template <typename StreamType, typename MatType>
void matrix_print(StreamType& oss, const std::string& name, const MatType& mat)
{
	const int m = traits::matrix_size1(mat);
	const int n = traits::matrix_size2(mat);

	oss << name << std::endl;
	for (int i=0; i < m; ++i)
	{
		for (int j=0; j < n; ++j)
		{
			oss << mat(i,j) << std::setw(10);
		}
		oss << std::setw(0) << std::endl;
	}
}

template <typename StreamType, typename VecType>
void vector_print(StreamType& oss, const std::string& name, const VecType& vec)
{
	const int m = traits::vector_size(vec);

	oss << name << std::endl;
	for (int i=0; i < m; ++i)
	{
		oss << vec(i) << std::endl;
	}
}

// structs to create matrices for testing
template <typename T>
struct MatrixGenerator {};

template <>
struct MatrixGenerator<fmat_t>
{
	typedef fmat_t Result;
	inline Result operator() (size_t m, size_t n)
	{
		Result mat(row_size, col_size);
		mat(0,0) = 8.;    
		mat(0,1) = 1.;
		mat(0,2) = 6.;
		mat(1,0) = 3.;
		mat(1,1) = 5.;
		mat(1,2) = 7.;
		mat(2,0) = 4.;
		mat(2,1) = 9.;
		mat(2,2) = 2.;

		return Result(ublas::project(mat, ublas::range(0, m), ublas::range(0, n)));
	}
};

template <>
struct MatrixGenerator<dmat_t>
{
	typedef dmat_t Result;
	inline Result operator() (size_t m, size_t n)
	{
		Result mat(row_size, col_size);
		mat(0,0) = 8.;    
		mat(0,1) = 1.;
		mat(0,2) = 6.;
		mat(1,0) = 3.;
		mat(1,1) = 5.;
		mat(1,2) = 7.;
		mat(2,0) = 4.;
		mat(2,1) = 9.;
		mat(2,2) = 2.;

		return Result(ublas::project(mat, ublas::range(0, m), ublas::range(0, n)));
	}
};

template <>
struct MatrixGenerator<fcmat_t>
{
	typedef fcmat_t Result;
	inline Result operator() (size_t m, size_t n)
	{
		typedef Result::value_type val_t;

		Result mat(row_size , col_size);
		mat(0,0) = val_t(35.,1.);
		mat(0,1) = val_t(6.,26.);
		mat(0,2) = val_t(19.,24.);
		mat(1,0) = val_t(3.,32.);
		mat(1,1) = val_t(7.,21.);
		mat(1,2) = val_t(23.,25.);
		mat(2,0) = val_t(31.,9.);
		mat(2,1) = val_t(2.,22.);
		mat(2,2) = val_t(27.,20.);

		return Result(ublas::project(mat, ublas::range(0, m), ublas::range(0, n)));
	}
};

template <>
struct MatrixGenerator<dcmat_t>
{
	typedef dcmat_t Result;
	inline Result operator() (size_t m, size_t n)
	{
		typedef Result::value_type val_t;

		Result mat(row_size , col_size);
		mat(0,0) = val_t(35.,1.);
		mat(0,1) = val_t(6.,26.);
		mat(0,2) = val_t(19.,24.);
		mat(1,0) = val_t(3.,32.);
		mat(1,1) = val_t(7.,21.);
		mat(1,2) = val_t(23.,25.);
		mat(2,0) = val_t(31.,9.);
		mat(2,1) = val_t(2.,22.);
		mat(2,2) = val_t(27.,20.);

		return Result(ublas::project(mat, ublas::range(0, m), ublas::range(0, n)));
	}
};


// structs to create vectors for testing
template <typename T>
struct VectorGenerator {};

template <>
struct VectorGenerator<fvec_t>
{
	typedef fvec_t Result;
	inline Result operator() (size_t m)
	{
		typedef Result::value_type val_t;
		Result v(m);

		for (size_t i=0; i < m; ++i)
		{
			val_t val = val_t(i+1);
			v(i) = val;
		}

		return v;
	}
};

template <>
struct VectorGenerator<dvec_t>
{
	typedef dvec_t Result;
	inline Result operator() (size_t m)
	{
		typedef Result::value_type val_t;
		Result v(m);

		for (size_t i=0; i < m; ++i)
		{
			val_t val = val_t(i+1);
			v(i) = val;
		}

		return v;
	}
};

template <>
struct VectorGenerator<fcvec_t>
{
	typedef fcvec_t Result;
	inline Result operator() (size_t m)
	{
		typedef Result::value_type val_t;
		Result v(m);

		for (size_t i=0; i < m; ++i)
		{
			val_t::value_type val = val_t::value_type(i);
			v(i) = val_t(val+1, val+1);
		}

		return v;
	}
};

template <>
struct VectorGenerator<dcvec_t>
{
	typedef dcvec_t Result;
	inline Result operator() (size_t m)
	{
		typedef Result::value_type val_t;
		Result v(m);

		for (size_t i=0; i < m; ++i)
		{
			val_t::value_type val = val_t::value_type(i);
			v(i) = val_t(val+1, val+1);
		}

		return v;
	}
};

#endif