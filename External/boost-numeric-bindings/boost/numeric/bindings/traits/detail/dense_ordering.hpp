//
//  Copyright Markus Rickert 2008
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_TRAITS_DETAIL_DENSE_ORDERING_H
#define BOOST_NUMERIC_BINDINGS_TRAITS_DETAIL_DENSE_ORDERING_H

#include <boost/numeric/ublas/fwd.hpp> 

namespace boost { namespace numeric { namespace bindings { namespace traits {

  namespace detail {
    
    template <typename StOrdTag>
    struct dense_ordering {};
    
    template<> 
    struct dense_ordering<row_major_t> {
      typedef row_major_t type; 
      
      template <typename M>
      static std::ptrdiff_t stride1( M const& m ) {
        return leading_dimension (m) ;
      }
      
      template <typename M>
      static std::ptrdiff_t stride2( M const& m ) {
        return 1 ;
      }
    };
    
    template<> 
    struct dense_ordering<column_major_t> {
      typedef column_major_t type; 
      
      template <typename M>
      static std::ptrdiff_t stride1( M const& m ) {
        return 1 ;
      }
      
      template <typename M>
      static std::ptrdiff_t stride2( M const& m ) {
        return leading_dimension (m) ;
      }
    };
    
  }

}}}}

#endif // BOOST_NUMERIC_BINDINGS_TRAITS_DETAIL_DENSE_ORDERING_H
