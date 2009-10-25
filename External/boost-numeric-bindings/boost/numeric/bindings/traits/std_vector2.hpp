/*
 * 
 * Copyright (c) 2002, 2003 Kresimir Fresl, Toon Knapen and Karl Meerbergen
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * KF acknowledges the support of the Faculty of Civil Engineering, 
 * University of Zagreb, Croatia.
 *
 */

#ifndef BOOST_NUMERIC_BINDINGS_TRAITS_STD_VECTOR2_H
#define BOOST_NUMERIC_BINDINGS_TRAITS_STD_VECTOR2_H

#include <boost/numeric/bindings/traits/matrix_traits.hpp>

#ifndef BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 

#include <vector>


namespace boost { namespace numeric { namespace bindings { namespace traits {

  // std::vector<> treated as matrix (nx1)
  template <typename T, typename Alloc, typename V>
  struct matrix_detail_traits<std::vector<T, Alloc>, V> 
  {
#ifndef BOOST_NUMERIC_BINDINGS_NO_SANITY_CHECK
    BOOST_STATIC_ASSERT( 
      (boost::is_same< 
         std::vector<T, Alloc>, 
         typename boost::remove_const<V>::type 
       >::value) );
#endif

    typedef std::vector<T, Alloc> identifier_type;
    typedef V matrix_type; 
    typedef general_t matrix_structure; 
    typedef column_major_t ordering_type; 

    typedef T value_type; 
    typedef typename default_vector_traits< V, T >::pointer pointer; 

    static pointer storage (matrix_type& v) {
      return vector_traits<matrix_type>::storage (v); 
    }
    static std::ptrdiff_t num_rows (matrix_type& v) { return v.size(); } 
    static std::ptrdiff_t num_columns (matrix_type&) { return 1; }
    static std::ptrdiff_t storage_size (matrix_type& v) { return v.size(); }
//    static std::ptrdiff_t stride1 (matrix_type& v) { return vector_traits<V>::stride (v); } 
//    static std::ptrdiff_t stride2 (matrix_type&) { return 1; }
    static std::ptrdiff_t leading_dimension (matrix_type& v) { return v.size(); }
  }; 

}}}}  

#endif // BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 

#endif // BOOST_NUMERIC_BINDINGS_TRAITS_STD_VECTOR2_H
