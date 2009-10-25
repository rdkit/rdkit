//
// Copyright Jesse Manning 2007
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_GELSD_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_GELSD_HPP

#include <algorithm>

#include <boost/numeric/bindings/traits/type.hpp>
#include <boost/numeric/bindings/traits/traits.hpp>
#include <boost/numeric/bindings/traits/type_traits.hpp>
#include <boost/numeric/bindings/lapack/lapack.h>
#include <boost/numeric/bindings/lapack/workspace.hpp>
#include <boost/numeric/bindings/traits/detail/array.hpp>
#include <boost/numeric/bindings/traits/detail/utils.hpp>
#include <boost/numeric/bindings/lapack/ilaenv.hpp>

// included to implicitly convert a vector to an nx1 matrix
// so that it is compatible with lapack binding
#include <boost/numeric/bindings/traits/ublas_vector2.hpp> 


#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK 
#  include <boost/static_assert.hpp>
#  include <boost/type_traits.hpp>
#endif

namespace boost { namespace numeric { namespace bindings {

	namespace lapack {

		namespace detail {

			inline void gelsd(const int m, const int n, const int nrhs, 
							  float *a, const int lda, float *b, const int ldb, 
							  float *s, const float rcond, int *rank, float *work,
							  const int lwork, int *iwork, int *info)
			{
				LAPACK_SGELSD(&m, &n, &nrhs, a, &lda, b, &ldb, s, 
							  &rcond, rank, work, &lwork, iwork, info);
			}

			inline void gelsd(const int m, const int n, const int nrhs, 
							  double *a, const int lda, double *b, const int ldb, 
							  double *s, const double rcond, int *rank, double *work,
							  const int lwork, int *iwork, int *info)
			{
				LAPACK_DGELSD(&m, &n, &nrhs, a, &lda, b, &ldb, s, 
							  &rcond, rank, work, &lwork, iwork, info);
			}

			inline void gelsd(const int m, const int n, const int nrhs, 
							  traits::complex_f *a, const int lda, traits::complex_f *b, 
							  const int ldb, float *s, const float rcond, int *rank, 
							  traits::complex_f *work, const int lwork, float *rwork, 
							  int *iwork, int *info)
			{
				LAPACK_CGELSD(&m, &n, &nrhs, traits::complex_ptr(a), 
							  &lda, traits::complex_ptr(b), &ldb, s, 
							  &rcond, rank, traits::complex_ptr(work), 
							  &lwork, rwork, iwork, info);
			}

			inline void gelsd(const int m, const int n, const int nrhs, 
							  traits::complex_d *a, const int lda, traits::complex_d *b, 
							  const int ldb, double *s, const double rcond, int *rank, 
							  traits::complex_d *work, const int lwork, double *rwork, 
							  int *iwork, int *info)
			{
				LAPACK_ZGELSD(&m, &n, &nrhs, traits::complex_ptr(a), 
							  &lda, traits::complex_ptr(b), &ldb, s, 
							  &rcond, rank, traits::complex_ptr(work), 
							  &lwork, rwork, iwork, info);
			}

			// gelsd for real type
			template <typename MatrA, typename MatrB, typename VecS, typename Work>
			inline int gelsd(MatrA& A, MatrB& B, VecS& s, Work& work)
			{
				typedef typename MatrA::value_type val_t;
				typedef typename traits::type_traits<val_t>::real_type real_t;

				const int m = traits::matrix_size1(A);
				const int n = traits::matrix_size2(A);
				const int nrhs = traits::matrix_size2(B);
				const int maxmn = std::max(m, n);
				const int minmn = std::min(m, n);

				// sanity checks
				assert(m >= 0 && n >= 0);
				assert(nrhs >= 0);
				assert(traits::leading_dimension(A) >= std::max(1, m));
				assert(traits::leading_dimension(B) >= std::max(1, maxmn));
				assert(traits::vector_size(work) >= 1);
				assert(traits::vector_size(s) >= std::max(1, minmn));

				int info;
				const real_t rcond = -1;	// use machine precision
				int rank;

				// query for maximum size of subproblems
				const int smlsiz = ilaenv(9, "GELSD", "");
				const int nlvl = static_cast<int>(((std::log(static_cast<float>(minmn))/std::log(2.f))/ (smlsiz+1)) + 1);

				traits::detail::array<int> iwork(3*minmn*nlvl + 11*minmn);

				detail::gelsd(traits::matrix_size1(A),
							  traits::matrix_size2(A),
							  traits::matrix_size2(B),
							  traits::matrix_storage(A),
							  traits::leading_dimension(A),
							  traits::matrix_storage(B),
							  traits::leading_dimension(B),
							  traits::vector_storage(s),
							  rcond,
							  &rank,
							  traits::vector_storage(work),
							  traits::vector_size(work),
							  traits::vector_storage(iwork),
							  &info);

				return info;
			}

			// gelsd for complex type
			template <typename MatrA, typename MatrB, typename VecS, 
						typename Work, typename RWork>
			inline int gelsd(MatrA& A, MatrB& B, VecS& s, Work& work, RWork& rwork)
			{
				typedef typename MatrA::value_type val_t;
				typedef typename traits::type_traits<val_t>::real_type real_t;

				const int m = traits::matrix_size1(A);
				const int n = traits::matrix_size2(A);
				const int nrhs = traits::matrix_size2(B);
				const int maxmn = std::max(m, n);
				const int minmn = std::min(m, n);

				// sanity checks
				assert(m >= 0 && n >= 0);
				assert(nrhs >= 0);
				assert(traits::leading_dimension(A) >= std::max(1, m));
				assert(traits::leading_dimension(B) >= std::max(1, maxmn));
				assert(traits::vector_size(work) >= 1);
				assert(traits::vector_size(s) >= std::max(1, minmn));

				int info;
				const real_t rcond = -1;	// use machine precision
				int rank;

				// query for maximum size of subproblems
				const int smlsiz = ilaenv(9, "GELSD", "");
				const int nlvl = static_cast<int>(((std::log(static_cast<float>(minmn))/std::log(2.f))/ (smlsiz+1)) + 1);

				traits::detail::array<int> iwork(3*minmn*nlvl + 11*minmn);

				detail::gelsd(traits::matrix_size1(A),
					traits::matrix_size2(A),
					traits::matrix_size2(B),
					traits::matrix_storage(A),
					traits::leading_dimension(A),
					traits::matrix_storage(B),
					traits::leading_dimension(B),
					traits::vector_storage(s),
					rcond,
					&rank,
					traits::vector_storage(work),
					traits::vector_size(work),
					traits::vector_storage(rwork),
					traits::vector_storage(iwork),
					&info);

				return info;
			}

			template <int N>
			struct Gelsd { };

			// specialization for gelsd real flavors (sgelsd, dgelsd)
			template <>
			struct Gelsd<1>
			{
				template <typename MatrA, typename MatrB, typename VecS>
				inline int operator() (MatrA& A, MatrB& B, VecS& s, minimal_workspace) const
				{
					typedef typename MatrA::value_type val_t;

					const int m = traits::matrix_size1(A);
					const int n = traits::matrix_size2(A);
					const int nrhs = traits::matrix_size2(B);

					const int minmn = std::min(m, n);			// minmn = m < n ? m : n
					const int maxmn = std::max(m, n);			// maxmn = m > n ? m : n
					const int maxmnr = std::max(maxmn, nrhs);	// maxmnr = maxmn > nrhs ? maxmn : nrhs

					// query for maximum size of subproblems
					const int smlsiz = ilaenv(9, "GELSD", "");
					const int nlvl = static_cast<int>(((std::log(static_cast<float>(minmn))/std::log(2.f))/ (smlsiz+1)) + 1);

					const int lwork = 12*minmn + 2*minmn*smlsiz + 8*minmn*nlvl + 
							     	  minmn*nrhs + (smlsiz+1)*(smlsiz+1);

					traits::detail::array<val_t> work(lwork);

					return gelsd(A, B, s, work);
				}

				template <typename MatrA, typename MatrB, typename VecS>
				inline int operator() (MatrA& A, MatrB& B, VecS& s, optimal_workspace) const
				{
					typedef typename MatrA::value_type val_t;
					typedef typename traits::type_traits<val_t>::real_type real_t;

					const int m = traits::matrix_size1(A);
					const int n = traits::matrix_size2(A);
					const int nrhs = traits::matrix_size2(B);

					const int minmn = std::min(m, n);			// minmn = m < n ? m : n
					const int maxmn = std::max(m, n);			// maxmn = m > n ? m : n
					const int maxmnr = std::max(maxmn, nrhs);	// maxmnr = maxmn > nrhs ? maxmn : nrhs

					val_t temp_work;
					int temp_iwork;

					const real_t rcond = -1;
					int rank;
					int info;

					// query for optimal workspace size
					detail::gelsd(traits::matrix_size1(A),
								  traits::matrix_size2(A),
								  traits::matrix_size2(B),
								  traits::matrix_storage(A),
								  traits::leading_dimension(A),
								  traits::matrix_storage(B),
								  traits::leading_dimension(B),
								  traits::vector_storage(s),
								  rcond,
								  &rank,
								  &temp_work,	//traits::vector_storage(work),
								  -1,			//traits::vector_size(work),
								  &temp_iwork,
								  &info);

					assert(info == 0);

					const int lwork = traits::detail::to_int(temp_work);

					traits::detail::array<val_t> work(lwork);

					return gelsd(A, B, s, work);
				}

				template <typename MatrA, typename MatrB, typename VecS, typename Work>
				inline int operator() (MatrA& A, MatrB& B, VecS& s, detail::workspace1<Work>& workspace) const
				{
					const int m = traits::matrix_size1(A);
					const int n = traits::matrix_size2(A);
					const int nrhs = traits::matrix_size2(B);

					const int minmn = std::min(m, n);			// minmn = m < n ? m : n
					const int maxmn = std::max(m, n);			// maxmn = m > n ? m : n
					const int maxmnr = std::max(maxmn, nrhs);	// maxmnr = maxmn > nrhs ? maxmn : nrhs

					return gelsd(A, B, s, workspace.w_);
				}
			};

			// specialization for gelsd (cgelsd, zgelsd)
			template <>
			struct Gelsd<2>
			{
				template <typename MatrA, typename MatrB, typename VecS>
				inline int operator() (MatrA& A, MatrB& B, VecS& s, minimal_workspace) const
				{
					typedef typename MatrA::value_type val_t;
					typedef typename traits::type_traits<val_t>::real_type real_t;

					const int m = traits::matrix_size1(A);
					const int n = traits::matrix_size2(A);
					const int nrhs = traits::matrix_size2(B);

					const int minmn = std::min(m, n);			// minmn = m < n ? m : n
					const int maxmn = std::max(m, n);			// maxmn = m > n ? m : n
					const int maxmnr = std::max(maxmn, nrhs);	// maxmnr = maxmn > nrhs ? maxmn : nrhs

					// query for maximum size of subproblems
					const int smlsiz = ilaenv(9, "GELSD", "");
					const int nlvl = static_cast<int>(((std::log(static_cast<float>(minmn))/std::log(2.f))/ (smlsiz+1)) + 1);

					traits::detail::array<val_t> work(2*minmn + minmn*nrhs);

					const int rwork_size = 10*minmn + 2*minmn*smlsiz + 8*minmn*nlvl + 
										   3*smlsiz*nrhs + (smlsiz+1)*(smlsiz+1);

					traits::detail::array<real_t> rwork(std::max(1, rwork_size));

					return gelsd(A, B, s, work, rwork);
				}

				template <typename MatrA, typename MatrB, typename VecS>
				inline int operator() (MatrA& A, MatrB& B, VecS& s, optimal_workspace) const
				{
					typedef typename MatrA::value_type val_t;
					typedef typename traits::type_traits<val_t>::real_type real_t;

					const int m = traits::matrix_size1(A);
					const int n = traits::matrix_size2(A);
					const int nrhs = traits::matrix_size2(B);

					const int minmn = std::min(m, n);			// minmn = m < n ? m : n
					const int maxmn = std::max(m, n);			// maxmn = m > n ? m : n
					const int maxmnr = std::max(maxmn, nrhs);	// maxmnr = maxmn > nrhs ? maxmn : nrhs

					val_t temp_work;
					real_t temp_rwork;
					int temp_iwork;

					const real_t rcond = -1;
					int rank;
					int info;

					// query for optimal workspace size
					detail::gelsd(traits::matrix_size1(A),
								  traits::matrix_size2(A),
								  traits::matrix_size2(B),
								  traits::matrix_storage(A),
								  traits::leading_dimension(A),
								  traits::matrix_storage(B),
								  traits::leading_dimension(B),
								  traits::vector_storage(s),
								  rcond,
								  &rank,
								  &temp_work,	//traits::vector_storage(work),
								  -1,			//traits::vector_size(work),
								  &temp_rwork,
								  &temp_iwork,
								  &info);

					assert(info == 0);

					const int lwork = traits::detail::to_int(temp_work);

					traits::detail::array<val_t> work(lwork);

					// query for maximum size of subproblems
					const int smlsiz = ilaenv(9, "GELSD", "");
					const int nlvl = static_cast<int>(((std::log(static_cast<float>(minmn))/std::log(2.f))/ (smlsiz+1)) + 1);

					const int rwork_size = 10*minmn + 2*minmn*smlsiz + 8*minmn*nlvl + 
											3*smlsiz*nrhs + (smlsiz+1)*(smlsiz+1);

					traits::detail::array<real_t> rwork(std::max(1, rwork_size));

					return gelsd(A, B, s, work, rwork);
				}

				template <typename MatrA, typename MatrB, typename VecS, typename Work, typename RWork>
				inline int operator() (MatrA& A, MatrB& B, VecS& s, detail::workspace2<Work, RWork>& workspace) const
				{
					const int m = traits::matrix_size1(A);
					const int n = traits::matrix_size2(A);
					const int nrhs = traits::matrix_size2(B);

					const int minmn = std::min(m, n);			// minmn = m < n ? m : n
					const int maxmn = std::max(m, n);			// maxmn = m > n ? m : n
					const int maxmnr = std::max(maxmn, nrhs);	// maxmnr = maxmn > nrhs ? maxmn : nrhs

					return gelsd(A, B, s, workspace.w_, workspace.wr_);
				}
			};

		} // detail

		// gelsd
		// Parameters:
		//	A:			matrix of coefficients
		//	B:			matrix of solutions (stored column-wise)
		//	s:			vector to store singular values on output, length >= max(1, min(m,n))
		//  workspace:	either optimal, minimal, or user supplied
		//
		template <typename MatrA, typename MatrB, typename VecS, typename Work>
		inline int gelsd(MatrA& A, MatrB& B, VecS& s, Work& workspace)
		{
			typedef typename MatrA::value_type val_t;

			return detail::Gelsd<n_workspace_args<val_t>::value>() (A, B, s, workspace);
		}

		// gelsd, no singular values are returned
		// Parameters:
		//	A:			matrix of coefficients
		//	B:			matrix of solutions (stored column-wise)
		//	workspace:	either optimal, minimal, or user supplied
		//
		template <typename MatrA, typename MatrB, typename Work>
		inline int gelsd(MatrA& A, MatrB& B, Work& workspace)
		{
			typedef typename MatrA::value_type val_t;
			typedef typename traits::type_traits<val_t>::real_type real_t;

			const int m = traits::matrix_size1(A);
			const int n = traits::matrix_size2(A);

			const int s_size = std::max(1, std::min(m,n));
			traits::detail::array<real_t> s(s_size);

			return detail::Gelsd<n_workspace_args<val_t>::value>() (A, B, s, workspace);
		}

	} // lapack

}}}

#endif
