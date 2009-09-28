
#include "convenience.h"
#include <boost/numeric/bindings/lapack/gelsd.hpp>

// set to 1 for more output
//#define VERBOSE_OUTPUT 0

// set to 1 to write test output to file, otherwise outputs to console
#define OUTPUT_TO_FILE 0

// determines which tests to run
#define TEST_SQUARE 1
#define TEST_UNDERDETERMINED 0
#define TEST_OVERDETERMINED 0
#define TEST_MULTIPLE_SOLUTION_VECTORS 0

// determines if optimal, minimal, or both workspaces are used for testing
#define USE_OPTIMAL_WORKSPACE 1
#define USE_MINIMAL_WORKSPACE 0 

namespace lapack = boost::numeric::bindings::lapack;

// test function declarations
template <typename StreamType, typename MatrType, typename VecType>
int test_square_gelsd(StreamType& oss);

template <typename StreamType, typename MatType, typename VecType>
int test_under_gelsd(StreamType& oss);

template <typename StreamType, typename MatType, typename VecType>
int test_over_gelsd(StreamType& oss);

template <typename StreamType, typename MatType, typename VecType>
int test_multiple_gelsd(StreamType& oss);

template <typename StreamType, typename MatType, typename VecType>
int test_transpose_gel(StreamType& oss, const char& trans);


int main()
{
	// stream for test output
	typedef std::ostringstream stream_t;
	stream_t oss;

#if TEST_SQUARE
	oss << "Start Square Matrix Least Squares Tests" << std::endl;
	oss << "Testing sgelsd" << std::endl;
	if (test_square_gelsd<stream_t, fmat_t, fvec_t>(oss) == 0)
	{
		oss << "sgelsd passed." << std::endl;
	}
	oss << "End sgelsd tests" << std::endl;
	oss << std::endl;
	oss << "Testing dgelsd" << std::endl;
	if (test_square_gelsd<stream_t, dmat_t, dvec_t>(oss) == 0)
	{
		oss << "dgelsd passed." << std::endl;
	}
	oss << "End dgelsd tests" << std::endl;
	oss << std::endl;
	oss << "Testing cgelsd" << std::endl;
	if (test_square_gelsd<stream_t, fcmat_t, fcvec_t>(oss) == 0)
	{
		oss << "cgelsd passed." << std::endl;
	}
	oss << "End cgelsd tests" << std::endl;
	oss << std::endl;
	oss << "Testing zgelsd" << std::endl;
	if (test_square_gelsd<stream_t, fcmat_t, fcvec_t>(oss) == 0)
	{
		oss << "zgelsd passed." << std::endl;
	}
	oss << "End zgelsd tests" << std::endl;
	oss << std::endl;
	oss << "End Square Matrix Least Squares Tests" << std::endl;
#endif

#if TEST_UNDERDETERMINED
	oss << std::endl;
	oss << "Start Under-determined Matrix Least Squares Test" << std::endl;
	oss << "Testing sgelsd" << std::endl;
	if (test_under_gelsd<stream_t, fmat_t, fvec_t>(oss) == 0)
	{
		oss << "sgelsd passed." << std::endl;
	}
	oss << "End sgelsd tests" << std::endl;
	oss << std::endl;
	oss << "Testing dgelsd" << std::endl;
	if (test_under_gelsd<stream_t, dmat_t, dvec_t>(oss) == 0)
	{
		oss << "dgelsd passed." << std::endl;
	}
	oss << "End dgelsd tests" << std::endl;
	oss << std::endl;
	oss << "Testing cgelsd" << std::endl;
	if (test_under_gelsd<stream_t, fcmat_t, fcvec_t>(oss) == 0)
	{
		oss << "cgelsd passed." << std::endl;
	}
	oss << "End cgelsd tests" << std::endl;
	oss << std::endl;
	oss << "Testing zgelsd" << std::endl;
	if (test_under_gelsd<stream_t, fcmat_t, fcvec_t>(oss) == 0)
	{
		oss << "zgelsd passed." << std::endl;
	}
	oss << "End zgelsd tests" << std::endl;
	oss << std::endl;
	oss << "End Underdetermined Matrix Least Squares Tests" << std::endl;
#endif

#if TEST_OVERDETERMINED
	oss << std::endl;
	oss << "Start Overdetermined Matrix Least Squares Test" << std::endl;
	oss << "Testing sgelsd" << std::endl;
	if (test_over_gelsd<stream_t, fmat_t, fvec_t>(oss) == 0)
	{
		oss << "sgelsd passed." << std::endl;
	}
	oss << "End sgelsd tests" << std::endl;
	oss << std::endl;
	oss << "Testing dgelsd" << std::endl;
	if (test_over_gelsd<stream_t, dmat_t, dvec_t>(oss) == 0)
	{
		oss << "dgelsd passed." << std::endl;
	}
	oss << "End dgelsd tests" << std::endl;
	oss << std::endl;
	oss << "Testing cgelsd" << std::endl;
	if (test_over_gelsd<stream_t, fcmat_t, fcvec_t>(oss) == 0)
	{
		oss << "cgelsd passed." << std::endl;
	}
	oss << "End cgelsd tests" << std::endl;
	oss << std::endl;
	oss << "Testing zgelsd" << std::endl;
	if (test_over_gelsd<stream_t, fcmat_t, fcvec_t>(oss) == 0)
	{
		oss << "zgelsd passed." << std::endl;
	}
	oss << "End zgelsd tests" << std::endl;
	oss << std::endl;
	oss << "End Overdetermined Matrix Least Squares Test" << std::endl;
#endif

#if TEST_MULTIPLE_SOLUTION_VECTORS
	oss << std::endl;
	oss << "Start Multiple Solution Vectors Least Squares Test" << std::endl;
	oss << "Testing sgelsd" << std::endl;
	if (test_multiple_gelsd<stream_t, fmat_t, fvec_t>(oss) == 0)
	{
		oss << "sgelsd passed." << std::endl;
	}
	oss << "End sgelsd tests" << std::endl;
	oss << std::endl;
	oss << "Testing dgelsd" << std::endl;
	if (test_multiple_gelsd<stream_t, dmat_t, dvec_t>(oss) == 0)
	{
		oss << "dgelsd passed." << std::endl;
	}
	oss << "End dgelsd tests" << std::endl;
	oss << std::endl;
	oss << "Testing cgelsd" << std::endl;
	if (test_multiple_gelsd<stream_t, fcmat_t, fcvec_t>(oss) == 0)
	{
		oss << "cgelsd passed." << std::endl;
	}
	oss << "End cgelsd tests" << std::endl;
	oss << std::endl;
	oss << "Testing zgelsd" << std::endl;
	if (test_multiple_gelsd<stream_t, fcmat_t, fcvec_t>(oss) == 0)
	{
		oss << "zgelsd passed." << std::endl;
	}
	oss << "End zgelsd tests" << std::endl;
	oss << std::endl;
	oss << "End Multiple Solution Vectors Least Squares Test" << std::endl;
#endif

#if OUTPUT_TO_FILE
	// Finished testing
	std::cout << std::endl;
	std::cout << "Tests Completed." << std::endl;
	std::cout << std::endl;

	std::string filename;
	std::cout << "Enter filename to write test results: ";
	std::getline(std::cin, filename);

	std::ofstream testFile(filename.c_str());

	if (testFile)
	{
		testFile << oss.str();
		testFile.close();
	}
#else
	std::cout << oss.str() << std::endl;

	// Finished testing
	std::cout << std::endl;
	std::cout << "Tests Completed." << std::endl;
	std::cout << std::endl;
#endif

	// wait for user to finish
	std::string done;
	std::cout << "Press Enter to exit";
	std::getline(std::cin, done);

}

// tests square system (m-by-n where m == n)
template <typename StreamType, typename MatType, typename VecType>
int test_square_gelsd(StreamType& oss)
{
	// return value
	int err = 0;

	// square matrix test
	MatType mat(MatrixGenerator<MatType>()(row_size, col_size));
	VecType vec(VectorGenerator<VecType>()(row_size));
	
	const int m = traits::matrix_size1(mat);
	const int n = traits::matrix_size2(mat);

#if USE_OPTIMAL_WORKSPACE
	MatType optimalmat(mat);
	VecType optimalvec(vec);
	err += lapack::gelsd(optimalmat, optimalvec, lapack::optimal_workspace());
	VecType optimalanswer(ublas::project(optimalvec, ublas::range(0, n)));
	VecType optimal_check = ublas::prod(mat, optimalanswer);
#endif
#if USE_MINIMAL_WORKSPACE
	MatType minimalmat(mat);
	VecType minimalvec(vec);
	err += lapack::gelsd(minimalmat, minimalvec, lapack::minimal_workspace());
	VecType minimalanswer(ublas::project(minimalvec, ublas::range(0, n)));
	VecType minimal_check = ublas::prod(mat, minimalanswer);
#endif

	matrix_print(oss, "A", mat);
	oss << std::endl;
	vector_print(oss, "B", vec);
	oss << std::endl;

#if USE_OPTIMAL_WORKSPACE
	vector_print(oss, "optimal workspace x", optimalanswer);
	oss << std::endl;
#endif
#if USE_MINIMAL_WORKSPACE
	vector_print(oss, "minimal workspace x", minimalanswer);
	oss << std::endl;
#endif
#if USE_OPTIMAL_WORKSPACE
	// check A*x=B
	vector_print(oss, "optimal A*x=B", optimal_check);
	oss << std::endl;
#endif
#if USE_MINIMAL_WORKSPACE
	vector_print(oss, "minimal A*x=B", minimal_check);
	oss << std::endl;
#endif

	return err;
}

// tests overdetermined system (m-by-n where m < n)
template <typename StreamType, typename MatType, typename VecType>
int test_under_gelsd(StreamType& oss)
{
	// return value
	int err = 0;

	// under-determined matrix test
	MatType mat(MatrixGenerator<MatType>()(row_range, col_size));
	VecType vec(VectorGenerator<VecType>()(row_size));

	const int m = traits::matrix_size1(mat);
	const int n = traits::matrix_size2(mat);

#if USE_OPTIMAL_WORKSPACE
	MatType optimalmat(mat);
	VecType optimalvec(vec);
	err += lapack::gelsd(optimalmat, optimalvec, lapack::optimal_workspace());
	VecType optimalanswer(ublas::project(optimalvec, ublas::range(0, n)));
	VecType optimal_check = ublas::prod(mat, optimalanswer);
#endif
#if USE_MINIMAL_WORKSPACE
	MatType minimalmat(mat);
	VecType minimalvec(vec);
	err += lapack::gelsd(minimalmat, minimalvec, lapack::minimal_workspace());
	VecType minimalanswer(ublas::project(minimalvec, ublas::range(0, n)));
	VecType minimal_check = ublas::prod(mat, minimalanswer);
#endif

	matrix_print(oss, "A", mat);
	oss << std::endl;
	vector_print(oss, "B", vec);
	oss << std::endl;

#if USE_OPTIMAL_WORKSPACE
	vector_print(oss, "optimal workspace x", optimalanswer);
	oss << std::endl;
#endif
#if USE_MINIMAL_WORKSPACE
	vector_print(oss, "minimal workspace x", minimalanswer);
	oss << std::endl;
#endif
#if USE_OPTIMAL_WORKSPACE
	// check A*x=B
	vector_print(oss, "optimal A*x=B", optimal_check);
	oss << std::endl;
#endif
#if USE_MINIMAL_WORKSPACE
	vector_print(oss, "minimal A*x=B", minimal_check);
	oss << std::endl;
#endif

	return err;
}

// tests overdetermined system (m-by-n where m > n)
template <typename StreamType, typename MatType, typename VecType>
int test_over_gelsd(StreamType& oss)
{
	// return value
	int err = 0;

	// overdetermined matrix test
	MatType mat(MatrixGenerator<MatType>()(row_size, col_range));
	VecType vec(VectorGenerator<VecType>()(row_size));

	const int m = traits::matrix_size1(mat);
	const int n = traits::matrix_size2(mat);

#if USE_OPTIMAL_WORKSPACE
	MatType optimalmat(mat);
	VecType optimalvec(vec);
	err += lapack::gelsd(optimalmat, optimalvec, lapack::optimal_workspace());
	VecType optimalanswer(ublas::project(optimalvec, ublas::range(0, n)));
	VecType optimal_check = ublas::prod(mat, optimalanswer);
#endif
#if USE_MINIMAL_WORKSPACE
	MatType minimalmat(mat);
	VecType minimalvec(vec);
	err += lapack::gelsd(minimalmat, minimalvec, lapack::minimal_workspace());
	VecType minimalanswer(ublas::project(minimalvec, ublas::range(0, n)));
	VecType minimal_check = ublas::prod(mat, minimalanswer);
#endif

	matrix_print(oss, "A", mat);
	oss << std::endl;
	vector_print(oss, "B", vec);
	oss << std::endl;

#if USE_OPTIMAL_WORKSPACE
	vector_print(oss, "optimal workspace x", optimalanswer);
	oss << std::endl;
#endif
#if USE_MINIMAL_WORKSPACE
	vector_print(oss, "minimal workspace x", minimalanswer);
	oss << std::endl;
#endif
#if USE_OPTIMAL_WORKSPACE
	// check A*x=B
	vector_print(oss, "optimal A*x=B", optimal_check);
	oss << std::endl;
#endif
#if USE_MINIMAL_WORKSPACE
	vector_print(oss, "minimal A*x=B", minimal_check);
	oss << std::endl;
#endif

	return err;
}

// tests multiple solution vectors stored column-wise in B for equation A*x=B
template <typename StreamType, typename MatType, typename VecType>
int test_multiple_gelsd(StreamType& oss)
{
	// return value
	int err = 0;

	// multiple solutions vectors test
	MatType mat(MatrixGenerator<MatType>()(row_size, col_size));
	MatType vec(mat.size1(), 2);
	ublas::column(vec, 0) = VectorGenerator<VecType>()(mat.size1());
	ublas::column(vec, 1) = VectorGenerator<VecType>()(mat.size1());

	const int m = traits::matrix_size1(mat);
	const int n = traits::matrix_size2(mat);
	const int nrhs = traits::matrix_size2(vec);

#if USE_OPTIMAL_WORKSPACE
	MatType optimalmat(mat);
	MatType optimalvec(vec);
	err += lapack::gelsd(optimalmat, optimalvec, lapack::optimal_workspace());
	MatType optimalanswer(ublas::project(optimalvec, ublas::range(0, n), ublas::range(0, nrhs)));
	MatType optimal_check = ublas::prod(mat, optimalanswer);
#endif
#if USE_MINIMAL_WORKSPACE
	MatType minimalmat(mat);
	MatType minimalvec(vec);
	err += lapack::gelsd(minimalmat, minimalvec, lapack::minimal_workspace());
	MatType minimalanswer(ublas::project(minimalvec, ublas::range(0, n), ublas::range(0, nrhs)));
	MatType minimal_check = ublas::prod(mat, minimalanswer);
#endif

	matrix_print(oss, "A", mat);
	oss << std::endl;
	matrix_print(oss, "B", vec);
	oss << std::endl;

#if USE_OPTIMAL_WORKSPACE
	matrix_print(oss, "optimal workspace x", optimalanswer);
	oss << std::endl;
#endif
#if USE_MINIMAL_WORKSPACE
	matrix_print(oss, "minimal workspace x", minimalanswer);
	oss << std::endl;
#endif
#if USE_OPTIMAL_WORKSPACE
	// check A*x=B
	matrix_print(oss, "optimal A*x=B", optimal_check);
	oss << std::endl;
#endif
#if USE_MINIMAL_WORKSPACE
	matrix_print(oss, "minimal A*x=B", minimal_check);
	oss << std::endl;
#endif

	return err;
}
