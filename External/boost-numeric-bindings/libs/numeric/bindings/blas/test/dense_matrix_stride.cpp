//
//  Copyright Markus Rickert 2008
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/bindings/traits/c_array2.hpp>
#include <boost/numeric/bindings/traits/dense_traits.hpp>
#include <boost/numeric/bindings/traits/std_valarray2.hpp>
#include <boost/numeric/bindings/traits/std_vector2.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
//#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_vector2.hpp>

int
main(int argc, char** argv)
{
	// matrix
	boost::numeric::ublas::matrix<double> a(3, 4);
	for (std::size_t i = 0; i < a.size1(); ++i) for (std::size_t j = 0; j < a.size2(); ++j) a(i, j) = i * a.size2() + j;
	std::cout << "A=";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride1(a) << ",";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride2(a) << std::endl;
	std::cout << a << std::endl;
	
	// matrix matrix_range
	boost::numeric::ublas::matrix_range<
		boost::numeric::ublas::matrix<double>
	> a2 = boost::numeric::ublas::subrange(a, 1, 3, 0, 2);
	std::cout << "A2=";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride1(a2) << ",";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride2(a2) << std::endl;
	std::cout << a2 << std::endl;
	
	// matrix column_major
	boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> b(3, 4);
	std::cout << "B=";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride1(b) << ",";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride2(b) << std::endl;
	
	// bounded_matrix
	boost::numeric::ublas::bounded_matrix<double, 3, 4> c;
	std::cout << "C=";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride1(c) << ",";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride2(c) << std::endl;
	
	// bounded_matrix matrix_range
	boost::numeric::ublas::matrix_range<
		boost::numeric::ublas::bounded_matrix<double, 3, 4>
	> c2 = boost::numeric::ublas::subrange(c, 0, 1, 0, 2);
	std::cout << "C2=";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride1(c2) << ",";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride2(c2) << std::endl;
	
	// bounded_matrix column_major
	boost::numeric::ublas::bounded_matrix<double, 3, 4, boost::numeric::ublas::column_major> d;
	std::cout << "D=";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride1(d) << ",";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride2(d) << std::endl;
	
	// vector
	boost::numeric::ublas::vector<double> e(3);
	std::cout << "e=";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride1(e) << ",";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride2(e) << std::endl;
	
	// vector vector_range
	boost::numeric::ublas::vector_range<
		boost::numeric::ublas::vector<double>
	> f = boost::numeric::ublas::subrange(e, 0, 2);
	std::cout << "f=";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride1(f) << ",";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride2(f) << std::endl;
	
	// bounded_vector
	boost::numeric::ublas::bounded_vector<double, 3> h;
	std::cout << "h=";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride1(h) << ",";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride2(h) << std::endl;
	
	// bounded_vector vector_range
	boost::numeric::ublas::vector_range<
		boost::numeric::ublas::bounded_vector<double, 3>
	> i = boost::numeric::ublas::subrange(h, 0, 2);
	std::cout << "i=";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride1(i) << ",";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride2(i) << std::endl;
	
	// std::vector
	std::vector<double> k(3);
	std::cout << "k=";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride1(k) << ",";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride2(k) << std::endl;
	
	// std::valarray
	std::valarray<double> l(3);
	std::cout << "l=";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride1(l) << ",";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride2(l) << std::endl;
	
	// array
	double m[3];
	std::cout << "m=";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride1(m) << ",";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride2(m) << std::endl;
	
	// matrix
	boost::numeric::ublas::matrix<double> n(3, 1);
	std::cout << "N=";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride1(n) << ",";
	std::cout << boost::numeric::bindings::traits::dense_matrix_stride2(n) << std::endl;
	std::cout << boost::numeric::bindings::traits::matrix_stride1(n) << ",";
	std::cout << boost::numeric::bindings::traits::matrix_stride2(n) << std::endl;
	
	return 0;
}
