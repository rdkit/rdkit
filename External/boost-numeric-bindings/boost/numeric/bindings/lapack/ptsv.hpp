/*
 * 
 * Copyright (c) Karl Meerbergen 2008
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 */

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_PTSV_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_PTSV_HPP

#include <boost/numeric/bindings/traits/type_traits.hpp>
#include <boost/numeric/bindings/traits/traits.hpp>
#include <boost/numeric/bindings/lapack/lapack.h>

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK 
#  include <boost/static_assert.hpp>
#endif 

#include <cassert>

namespace boost { namespace numeric { namespace bindings { 

  namespace lapack {

    /////////////////////////////////////////////////////////////////////
    //
    // system of linear equations A * X = B
    // with A tridiagonal symmetric or Hermitian positive definite matrix
    //
    /////////////////////////////////////////////////////////////////////

    /*
     * ptsv() computes the solution to a system of linear equations
     * A*X = B, where A is an N-by-N Hermitian positive definite tridiagonal
     * matrix, and X and B are N-by-NRHS matrices.
     *
     * A is factored as A = L*D*L**T, and the factored form of A is then
     * used to solve the system of equations.
     */

    namespace detail {

      inline 
      void ptsv ( int const n, int const nrhs,
                 float* d, float* e, float* b, int const ldb, int* info) 
      {
        LAPACK_SPTSV (&n, &nrhs, d, e, b, &ldb, info);
      }

      inline 
      void ptsv ( int const n, int const nrhs,
                 double* d, double* e, double* b, int const ldb, int* info) 
      {
        LAPACK_DPTSV (&n, &nrhs, d, e, b, &ldb, info);
      }

      inline 
      void ptsv ( int const n, int const nrhs,
                 float* d, traits::complex_f* e, traits::complex_f* b, int const ldb, 
                 int* info) 
      {
        LAPACK_CPTSV (&n, &nrhs, d, traits::complex_ptr(e), traits::complex_ptr(b), &ldb, info);
      }

      inline 
      void ptsv ( int const n, int const nrhs,
                 double* d, traits::complex_d* e, traits::complex_d* b, int const ldb, 
                 int* info) 
      {
        LAPACK_ZPTSV (&n, &nrhs, d, traits::complex_ptr(e), traits::complex_ptr(b), &ldb, info);
      }

    }

    template <typename D, typename E, typename B>
    inline int ptsv( D& d, E& e, B& b ) {
      int const n = traits::vector_size(d) ;
      assert( n==traits::vector_size(e)+1 ) ;
      assert( n==traits::matrix_num_rows(b) ) ;

      int info ;
      detail::ptsv( n, traits::matrix_num_columns (b)
                  , traits::vector_storage(d)
                  , traits::vector_storage(e)
                  , traits::matrix_storage(b)
                  , traits::leading_dimension(b)
                  , &info
                  ) ;
      return info ;
    } // ptsv()


    /*
     * pttrf() computes the L * D * L^H factorization of a Hermitian
     * positive definite tridiagonal matrix A.  The factorization may also
     * be regarded as having the form A = U^H * D *U.
     */

    namespace detail {

      inline 
      void pttrf ( int const n, float* d, float* e, int* info) {
        LAPACK_SPTTRF ( &n, d, e, info) ;
      }

      inline 
      void pttrf ( int const n, double* d, double* e, int* info) {
        LAPACK_DPTTRF ( &n, d, e, info);
      }

      inline 
      void pttrf ( int const n, float* d, traits::complex_f* e, int* info) 
      {
        LAPACK_CPTTRF ( &n, d, traits::complex_ptr(e), info);
      }

      inline 
      void pttrf ( int const n, double* d, traits::complex_d* e, int* info) 
      {
        LAPACK_ZPTTRF ( &n, d, traits::complex_ptr(e), info);
      }

    }

    template <typename D, typename E>
    inline
    int pttrf (D& d, E& e) {
      int const n = traits::vector_size (d);
      assert (n == traits::vector_size (e) + 1);
      int info; 
      detail::pttrf ( n, traits::vector_storage(d), traits::vector_storage(e), &info);
      return info; 
    }


    /*
     * pttrs() solves a tridiagonal system of the form
     *   A * X = B
     * using the factorization A = U^H * D * U or A = L * D * L^H computed by pttrf().
     * D is a diagonal matrix specified in the vector D, U (or L) is a unit
     * bidiagonal matrix whose superdiagonal (subdiagonal) is specified in
     * the vector E, and X and B are N by NRHS matrices.
     */

    namespace detail {

      inline 
      void pttrs (char const uplo, int const n, int const nrhs,
                  float const* d, float const* e, float* b, int const ldb, int* info) 
      {
        LAPACK_SPTTRS (&n, &nrhs, d, e, b, &ldb, info);
      }

      inline 
      void pttrs (char const uplo, int const n, int const nrhs,
                  double const* d, double const* e, double* b, int const ldb, int* info) 
      {
        LAPACK_DPTTRS (&n, &nrhs, d, e, b, &ldb, info);
      }

      inline 
      void pttrs (char const uplo, int const n, int const nrhs,
                  float const* d, 
                  traits::complex_f const* e, 
                  traits::complex_f* b, int const ldb, int* info) 
      {
        LAPACK_CPTTRS (&uplo, &n, &nrhs, d, 
                       traits::complex_ptr (e), 
                       traits::complex_ptr (b), &ldb, info);
      }

      inline 
      void pttrs (char const uplo, int const n, int const nrhs,
                  double const* d, 
                  traits::complex_d const* e, 
                  traits::complex_d* b, int const ldb, int* info) 
      {
        LAPACK_ZPTTRS (&uplo, &n, &nrhs, d, 
                       traits::complex_ptr (e), 
                       traits::complex_ptr (b), &ldb, info);
      }

    }

    template <typename D, typename E, typename MatrB>
    inline
    int pttrs (char uplo, D const& d, E const& e, MatrB& b) {
      int const n = traits::vector_size (d);
      assert (n == traits::vector_size (e) + 1);
      assert (n == traits::matrix_num_rows (b));
      
      int info; 
      detail::pttrs (uplo, n, traits::matrix_num_columns (b),
                     traits::vector_storage (d), 
                     traits::vector_storage (e), 
                     traits::matrix_storage (b), 
                     traits::leading_dimension (b), 
                     &info);
      return info; 
    } // pttrs()

  }

}}}

#endif 
