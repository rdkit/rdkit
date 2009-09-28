//
// Copyright Vardan Akopian 2007
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_GBSV_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_GBSV_HPP

#include <boost/numeric/bindings/traits/type_traits.hpp>
#include <boost/numeric/bindings/traits/traits.hpp>
#include <boost/numeric/bindings/lapack/lapack.h>
#include <boost/numeric/bindings/traits/ublas_banded.hpp>
#include <boost/numeric/bindings/traits/ublas_vector2.hpp>

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK 
#  include <boost/static_assert.hpp>
#  include <boost/type_traits/is_same.hpp>
#endif 


namespace boost { namespace numeric { namespace bindings { 

  namespace lapack {


    namespace detail {
      inline 
      void gbtrf (int const n, int const m, int const kl, int const ku,
                  double* ab, int const ldab, int* ipiv, int* info) 
      {
        LAPACK_DGBTRF (&n, &m, &kl, &ku, ab, &ldab, ipiv, info);
      }
    }

    template <typename MatrA, typename IVec>
    inline
    int gbtrf (MatrA& a, IVec& ipiv) {

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK 
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<MatrA>::matrix_structure, 
        traits::banded_t
      >::value)); 
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<MatrA>::ordering_type, 
        traits::row_major_t
      >::value)); 
#endif 

      int const n = traits::matrix_size1 (a);
      int const m = traits::matrix_size2 (a); 
      assert (traits::vector_size (ipiv) == (m < n ? m : n));

      // if the matrix has kl lower and ku upper diagonals, then we should have
      // allocated kl lower and kl+ku upper diagonals
      int const kl = traits::matrix_lower_bandwidth (a);
      int const ku = traits::matrix_upper_bandwidth (a) - kl;
      int const ld = traits::leading_dimension (a);

      assert(ku >= 0);

      int info; 
      detail::gbtrf (n, m, kl, ku,
                     traits::matrix_storage (a), 
		     ld,
                     traits::vector_storage (ipiv),  
                     &info);
      return info; 
    }


    namespace detail {
      inline 
      void gbtrs (char const trans, int const n, int const kl, int const ku, int const m,
                  double const* ab, int const ldab, int const* ipiv,
		  double* b, int const ldb, int* info) 
      {
        LAPACK_DGBTRS (&trans, &n, &kl, &ku, &m, ab, &ldab, ipiv, b, &ldb, info);
      }
    }


    template <typename MatrA, typename MatrB, typename IVec>
    inline
    int gbtrs (char const trans, MatrA const& a, IVec const& ipiv, MatrB& b) 
    {
      assert (trans == 'N' || trans == 'T' || trans == 'C'); 

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK 
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<MatrA>::matrix_structure, 
        traits::banded_t
      >::value)); 
#endif 

      int const n = traits::matrix_size1 (a);
      assert (n == traits::matrix_size2 (a)); 
      assert (n == traits::matrix_size1 (b)); 
      assert (n == traits::vector_size (ipiv)); 

      // if the matrix has kl lower and ku upper diagonals, then we should have
      // allocated kl lower and kl+ku upper diagonals
      int const kl = traits::matrix_lower_bandwidth (a);
      int const ku = traits::matrix_upper_bandwidth (a) - kl;
      int const ld = traits::leading_dimension (a);

      assert(ku >= 0);

      int info; 
      detail::gbtrs (trans, n, kl, ku, traits::matrix_size2 (b), 
#ifndef BOOST_NO_FUNCTION_TEMPLATE_ORDERING
                     traits::matrix_storage (a), 
#else
                     traits::matrix_storage_const (a), 
#endif 
                     ld,
#ifndef BOOST_NO_FUNCTION_TEMPLATE_ORDERING
                     traits::vector_storage (ipiv),  
#else
                     traits::vector_storage_const (ipiv),  
#endif
                     traits::matrix_storage (b),
                     traits::leading_dimension (b),
                     &info);
      return info; 
    }


  }

}}}

#endif 
