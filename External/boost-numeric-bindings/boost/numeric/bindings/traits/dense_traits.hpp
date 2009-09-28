//
//  Copyright Markus Rickert 2008
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_TRAITS_DENSE_TRAITS_H
#define BOOST_NUMERIC_BINDINGS_TRAITS_DENSE_TRAITS_H

#include <boost/numeric/bindings/traits/config.hpp> 

#ifndef BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 

#include <boost/numeric/bindings/traits/traits.hpp>
#include <boost/numeric/bindings/traits/detail/dense_ordering.hpp>

namespace boost { namespace numeric { namespace bindings { namespace traits {

  template <typename M>
  inline
  std::ptrdiff_t 
  dense_matrix_stride1 (M& m) { 
    return detail::dense_ordering< typename matrix_traits<M>::ordering_type >::stride1 (m); 
  }
  template <typename M>
  inline
  std::ptrdiff_t 
  dense_matrix_stride2 (M& m) { 
    return detail::dense_ordering< typename matrix_traits<M>::ordering_type >::stride2 (m); 
  }

}}}}  

#else // BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 

#error with your compiler dense matrices cannot be used in bindings

#endif // BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 

#endif // BOOST_NUMERIC_BINDINGS_TRAITS_DENSE_TRAITS_H
