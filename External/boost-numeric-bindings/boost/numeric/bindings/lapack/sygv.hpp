//
// Copyright Fabien Dekeyser, Quoc-Cuong Pham 2008
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_SYGV_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_SYGV_HPP

#include <boost/numeric/bindings/lapack/hegv.hpp>

namespace boost { namespace numeric { namespace bindings { 
  namespace lapack {

    template <typename A, typename B, typename W, typename Work>
    int sygv (int itype, char jobz, char uplo, A& a, B& b, W& w, Work work = optimal_workspace()) {

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK 
      typedef typename A::value_type                               value_type ;
      typedef typename traits::type_traits< value_type >::real_type real_type ;
      BOOST_STATIC_ASSERT((boost::is_same<value_type, real_type>::value));
#endif

      return hegv (itype, jobz, uplo, a, b, w, work);
    }
  }

}}}

#endif
