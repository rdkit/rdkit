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

#ifndef BOOST_NUMERIC_BINDINGS_TRAITS_UBLAS_VECTOR_AS_MATRIX_H
#define BOOST_NUMERIC_BINDINGS_TRAITS_UBLAS_VECTOR_AS_MATRIX_H

#include <boost/numeric/bindings/traits/config.hpp> 

#ifndef BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS

#ifndef BOOST_UBLAS_HAVE_BINDINGS
#  include <boost/numeric/bindings/traits/ublas_vector.hpp>
#endif 
#include <boost/numeric/bindings/traits/matrix_traits.hpp>

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK
#  include <boost/static_assert.hpp>
#  include <boost/type_traits/same_traits.hpp>
#  include <boost/mpl/if.hpp> 
#endif


namespace boost { namespace numeric { namespace bindings { namespace traits {

  // ublas::vector<> treated as matrix (nx1)
  template <typename T, typename ArrT, typename V>
  struct matrix_detail_traits<boost::numeric::ublas::vector<T, ArrT>, V> 
  {
#ifndef BOOST_NUMERIC_BINDINGS_NO_SANITY_CHECK
    BOOST_STATIC_ASSERT( 
      (boost::is_same< 
         boost::numeric::ublas::vector<T, ArrT>, 
         typename boost::remove_const<V>::type 
       >::value) );
#endif

    typedef boost::numeric::ublas::vector<T, ArrT> identifier_type;
    typedef V matrix_type; 
    typedef general_t matrix_structure; 
    typedef column_major_t ordering_type; 

    typedef T value_type; 
    typedef typename detail::generate_const<V,T>::type* pointer; 

    static pointer storage (matrix_type& v) {
      typedef typename detail::generate_const<V,ArrT>::type array_type;
      return vector_traits<array_type>::storage (v.data()); 
    }
    static std::ptrdiff_t num_rows (matrix_type& v) { return v.size(); } 
    static std::ptrdiff_t num_columns (matrix_type&) { return 1; }
    static std::ptrdiff_t storage_size (matrix_type& v) { return v.size(); }
//    static std::ptrdiff_t stride1 (matrix_type& v) { return vector_traits<V>::stride (v); } 
//    static std::ptrdiff_t stride2 (matrix_type&) { return 1; }
    static std::ptrdiff_t leading_dimension (matrix_type& v) { return v.size(); }
  }; 


  // ublas::vector_range<> treated as matrix (nx1)
  template <typename T, typename V>
  struct matrix_detail_traits<boost::numeric::ublas::vector_range<T>, V> 
  {
#ifndef BOOST_NUMERIC_BINDINGS_NO_SANITY_CHECK
    BOOST_STATIC_ASSERT( 
      (boost::is_same< 
         boost::numeric::ublas::vector_range<T>, 
         typename boost::remove_const<V>::type 
       >::value) );
#endif

    typedef boost::numeric::ublas::vector_range<T> identifier_type;
    typedef V matrix_type; 
    typedef general_t matrix_structure; 
    typedef column_major_t ordering_type; 

    typedef typename T::value_type value_type; 
    typedef typename detail::generate_const<V,value_type>::type* pointer; 

    static pointer storage (matrix_type& v) {
      return vector_traits<V>::storage (v); 
    }
    static std::ptrdiff_t num_rows (matrix_type& v) { return v.size(); } 
    static std::ptrdiff_t num_columns (matrix_type&) { return 1; }
    static std::ptrdiff_t storage_size (matrix_type& v) { return v.size(); }
//    static std::ptrdiff_t stride1 (matrix_type& v) { return vector_traits<V>::stride (v); } 
//    static std::ptrdiff_t stride2 (matrix_type&) { return 1; }
    static std::ptrdiff_t leading_dimension (matrix_type& v) { return v.size(); }
  }; 


#ifndef BOOST_NUMERIC_BINDINGS_FORTRAN 

  // (undocumented) ublas::c_vector<>
  template <typename T, std::size_t N, typename V>
  struct matrix_detail_traits<boost::numeric::ublas::c_vector<T,N>, V> 
  {
#ifndef BOOST_NUMERIC_BINDINGS_NO_SANITY_CHECK
    BOOST_STATIC_ASSERT( 
      (boost::is_same< 
         boost::numeric::ublas::c_vector<T,N>, 
         typename boost::remove_const<V>::type 
       >::value) );
#endif

    typedef boost::numeric::ublas::c_vector<T,N> identifier_type;
    typedef V matrix_type; 
    typedef general_t matrix_structure; 
    typedef row_major_t ordering_type; // consistent with c_matrix<> 

    typedef T value_type; 
    typedef typename detail::generate_const<V,T>::type* pointer; 

    static pointer storage (matrix_type& v) { return v.data(); }
    static std::ptrdiff_t num_rows (matrix_type&) { return 1; } 
    static std::ptrdiff_t num_columns (matrix_type& v) { return v.size(); }
//    static std::ptrdiff_t stride1 (matrix_type& v) { return vector_traits<V>::stride (v); } 
//    static std::ptrdiff_t stride2 (matrix_type&) { return 1; }
    static std::ptrdiff_t storage_size (matrix_type&) { return N; }
    static std::ptrdiff_t leading_dimension (matrix_type&) { return N; }
  }; 

#endif // BOOST_NUMERIC_BINDINGS_FORTRAN 


  // ublas::bounded_vector<> treated as matrix (nx1)
  template <typename T, std::size_t N, typename V>
  struct matrix_detail_traits<boost::numeric::ublas::bounded_vector<T, N>, V> 
  {
#ifndef BOOST_NUMERIC_BINDINGS_NO_SANITY_CHECK
    BOOST_STATIC_ASSERT( 
      (boost::is_same< 
         boost::numeric::ublas::bounded_vector<T, N>, 
         typename boost::remove_const<V>::type 
       >::value) );
#endif

    typedef boost::numeric::ublas::bounded_vector<T, N> identifier_type;
    typedef V matrix_type; 
    typedef general_t matrix_structure; 
    typedef column_major_t ordering_type; 

    typedef T value_type; 
    typedef typename detail::generate_const<V,T>::type* pointer; 

    static pointer storage (matrix_type& v) {
      typedef typename detail::generate_const<V,typename identifier_type::array_type>::type array_type;
      return vector_traits<array_type>::storage (v.data()); 
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

#error with your compiler ublas::vector<> cannot be used as matrix

#endif // BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS

#endif // BOOST_NUMERIC_BINDINGS_TRAITS_UBLAS_VECTOR_AS_MATRIX_H
