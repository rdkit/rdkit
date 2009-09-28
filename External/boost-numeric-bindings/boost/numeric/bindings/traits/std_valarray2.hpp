/*
 * 
 * Copyright (c) 2002, 2003 Kresimir Fresl, Toon Knapen and Karl Meerbergen
 * Copyright (c) 2008 Markus Rickert
 *
 * Permission to copy, modify, use and distribute this software 
 * for any non-commercial or commercial purpose is granted provided 
 * that this license appear on all copies of the software source code.
 *
 * Authors assume no responsibility whatsoever for its use and makes 
 * no guarantees about its quality, correctness or reliability.
 *
 * KF acknowledges the support of the Faculty of Civil Engineering, 
 * University of Zagreb, Croatia.
 *
 */

#ifndef BOOST_NUMERIC_BINDINGS_TRAITS_STD_VALARRAY2_H
#define BOOST_NUMERIC_BINDINGS_TRAITS_STD_VALARRAY2_H

#include <boost/numeric/bindings/traits/config.hpp> 

#ifndef BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 

#include <boost/numeric/bindings/traits/matrix_traits.hpp>
#include <valarray>


namespace boost { namespace numeric { namespace bindings { namespace traits {

  // std::valarray<> treated as matrix (nx1)
  template <typename T, typename V>
  struct matrix_detail_traits<std::valarray<T>, V> 
  {
#ifndef BOOST_NUMERIC_BINDINGS_NO_SANITY_CHECK
    BOOST_STATIC_ASSERT( 
      (boost::is_same< 
         std::valarray<T>, 
         typename boost::remove_const<V>::type 
       >::value) );
#endif

    typedef std::valarray<T> identifier_type;
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

#else // BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 

#error with your compiler std::valarray<> cannot be used in bindings

#endif // BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 

#endif // BOOST_NUMERIC_BINDINGS_TRAITS_STD_VALARRAY2_H
