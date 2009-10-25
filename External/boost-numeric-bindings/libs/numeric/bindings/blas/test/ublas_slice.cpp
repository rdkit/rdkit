//
//  Copyright Markus Rickert 2008
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include <algorithm>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/bindings/blas/blas.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_vector2.hpp>

int
main(int argc, char** argv)
{
	// a * b' = C ; a' * b = d
	{
		boost::numeric::ublas::vector<double> a(3);
		for (std::size_t i = 0; i < a.size(); ++i) a(i) = i;
		std::cout << "a=" << a << std::endl;
		
		boost::numeric::ublas::vector<double> b(3);
		for (std::size_t i = 0; i < b.size(); ++i) b(i) = i;
		std::cout << "b=" << b << std::endl;
		
		boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> c(3, 3);
		boost::numeric::bindings::blas::gemm(
			boost::numeric::bindings::traits::NO_TRANSPOSE,
			boost::numeric::bindings::traits::TRANSPOSE,
			1.0, a, b, 0.0, c
		);
		std::cout << "C=" << c << std::endl;
		
		boost::numeric::ublas::vector<double> d(1);
		boost::numeric::bindings::blas::gemm(
			boost::numeric::bindings::traits::TRANSPOSE,
			boost::numeric::bindings::traits::NO_TRANSPOSE,
			1.0, a, b, 0.0, d
		);
		std::cout << "d=" << d << std::endl;
	}
	
	std::cout << std::endl;
	
	// a * b' = C ; a' * b = d
	{
		boost::numeric::ublas::bounded_vector<double, 3> a;
		for (std::size_t i = 0; i < a.size(); ++i) a(i) = i;
		std::cout << "a=" << a << std::endl;
		
		boost::numeric::ublas::bounded_vector<double, 3> b;
		for (std::size_t i = 0; i < b.size(); ++i) b(i) = i;
		std::cout << "b=" << b << std::endl;
		
		boost::numeric::ublas::bounded_matrix<double, 3, 3, boost::numeric::ublas::column_major> c;
		boost::numeric::bindings::blas::gemm(
			boost::numeric::bindings::traits::NO_TRANSPOSE,
			boost::numeric::bindings::traits::TRANSPOSE,
			1.0, a, b, 0.0, c
		);
		std::cout << "C=" << c << std::endl;
		
		boost::numeric::ublas::bounded_vector<double, 1> d;
		boost::numeric::bindings::blas::gemm(
			boost::numeric::bindings::traits::TRANSPOSE,
			boost::numeric::bindings::traits::NO_TRANSPOSE,
			1.0, a, b, 0.0, d
		);
		std::cout << "d=" << d << std::endl;
	}
	
	std::cout << std::endl;
	
	// A * B = C
	{
		boost::numeric::ublas::bounded_matrix<double, 4, 3, boost::numeric::ublas::column_major> a;
		for (std::size_t i = 0; i < a.size1(); ++i) for (std::size_t j = 0; j < a.size2(); ++j) a(i, j) = i * a.size2() + j;
		std::cout << "A=" << a << std::endl;
		
		boost::numeric::ublas::bounded_matrix<double, 3, 4, boost::numeric::ublas::column_major> b;
		for (std::size_t i = 0; i < b.size1(); ++i) for (std::size_t j = 0; j < b.size2(); ++j) b(i, j) = i * b.size2() + j;
		std::cout << "B=" << b << std::endl;
		
		boost::numeric::ublas::bounded_matrix<double, 4, 4, boost::numeric::ublas::column_major> c;
		boost::numeric::bindings::blas::gemm(a, b, c);
		std::cout << "C=" << c << std::endl;
	}
	
	std::cout << std::endl;
	
	// A[0:3;0:2] * B[0:2;0:3] = C
	{
		boost::numeric::ublas::bounded_matrix<double, 4, 3, boost::numeric::ublas::column_major> a;
		for (std::size_t i = 0; i < a.size1(); ++i) for (std::size_t j = 0; j < a.size2(); ++j) a(i, j) = i * a.size2() + j;
		std::cout << "A=" << a << std::endl;
		
		boost::numeric::ublas::matrix_range<
			boost::numeric::ublas::bounded_matrix<double, 4, 3, boost::numeric::ublas::column_major>
		> a2 = boost::numeric::ublas::subrange(a, 0, 3, 0, 2);
		std::cout << "A2=" << a2 << std::endl;
		
		boost::numeric::ublas::bounded_matrix<double, 3, 4, boost::numeric::ublas::column_major> b;
		for (std::size_t i = 0; i < b.size1(); ++i) for (std::size_t j = 0; j < b.size2(); ++j) b(i, j) = i * b.size2() + j;
		std::cout << "B=" << b << std::endl;
		
		boost::numeric::ublas::matrix_range<
			boost::numeric::ublas::bounded_matrix<double, 3, 4, boost::numeric::ublas::column_major>
		> b2 = boost::numeric::ublas::subrange(b, 0, 2, 0, 3);
		std::cout << "B2=" << b2 << std::endl;
		
		boost::numeric::ublas::bounded_matrix<double, 4, 4, boost::numeric::ublas::column_major> c;
		std::fill(c.data().begin(), c.data().end(), 0.0);
		boost::numeric::ublas::matrix_range<
			boost::numeric::ublas::bounded_matrix<double, 4, 4, boost::numeric::ublas::column_major>
		> c2 = boost::numeric::ublas::subrange(c, 0, 3, 0, 3);
		boost::numeric::bindings::blas::gemm(a2, b2, c2);
		std::cout << "C2=" << c2 << std::endl;
		std::cout << "C=" << c << std::endl;
	}
	
	std::cout << std::endl;
	
	// a + b = b ; b - a = b
	{
		boost::numeric::ublas::bounded_vector<double, 3> a;
		for (std::size_t i = 0; i < a.size(); ++i) a(i) = i;
		std::cout << "a=" << a << std::endl;
		
		boost::numeric::ublas::bounded_vector<double, 3> b;
		for (std::size_t i = 0; i < b.size(); ++i) b(i) = i;
		std::cout << "b=" << b << std::endl;
		
		boost::numeric::bindings::blas::axpy(1.0, a, b);
		std::cout << "b=" << b << std::endl;
		
		boost::numeric::bindings::blas::axpy(-1.0, a, b);
		std::cout << "b=" << b << std::endl;
	}
	
	std::cout << std::endl;
	
	// b + c = c ; c - b = c
	{
		boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> a(5, 5);
		for (std::size_t i = 0; i < a.size1(); ++i) for (std::size_t j = 0; j < a.size2(); ++j) a(i, j) = i * a.size2() + j;
		std::cout << "A=" << a << std::endl;
		
		boost::numeric::ublas::matrix_vector_range<
			boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major>
		> b(a, boost::numeric::ublas::range(1, 4), boost::numeric::ublas::range(0, 3));
		std::cout << "b=" << b << std::endl;
		
		boost::numeric::ublas::matrix_vector_slice<
			boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major>
		> c(a, boost::numeric::ublas::slice(0, 1, 3), boost::numeric::ublas::slice(3, 0, 3));
		std::cout << "c=" << c << std::endl;
		
		boost::numeric::bindings::blas::axpy(1.0, b, c);
		std::cout << "c=" << c << std::endl;
		
		boost::numeric::bindings::blas::axpy(-1.0, b, c);
		std::cout << "c=" << c << std::endl;
	}
	
	std::cout << std::endl;
	
	// b + c = c ; c - b = c
	{
		boost::numeric::ublas::bounded_matrix<double, 5, 5, boost::numeric::ublas::column_major> a;
		for (std::size_t i = 0; i < a.size1(); ++i) for (std::size_t j = 0; j < a.size2(); ++j) a(i, j) = i * a.size2() + j;
		std::cout << "A=" << a << std::endl;
		
		boost::numeric::ublas::matrix_vector_range<
			boost::numeric::ublas::bounded_matrix<double, 5, 5, boost::numeric::ublas::column_major>
		> b(a, boost::numeric::ublas::range(1, 4), boost::numeric::ublas::range(0, 3));
		std::cout << "b=" << b << std::endl;
		
		boost::numeric::ublas::matrix_vector_slice<
			boost::numeric::ublas::bounded_matrix<double, 5, 5, boost::numeric::ublas::column_major>
		> c(a, boost::numeric::ublas::slice(0, 1, 3), boost::numeric::ublas::slice(3, 0, 3));
		std::cout << "c=" << c << std::endl;
		
		boost::numeric::bindings::blas::axpy(1.0, b, c);
		std::cout << "c=" << c << std::endl;
		
		boost::numeric::bindings::blas::axpy(-1.0, b, c);
		std::cout << "c=" << c << std::endl;
	}
	
	return 0;
}
