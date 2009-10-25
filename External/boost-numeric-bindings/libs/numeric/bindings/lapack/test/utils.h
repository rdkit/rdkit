
#ifndef F_UTILS_H
#define F_UTILS_H

#include <stddef.h>
#include <iostream>
#include <boost/numeric/bindings/traits/traits.hpp>


///////////////////////////////
// vectors

// element access:
template <typename V>
struct vct_access_traits {
  typedef typename 
  boost::numeric::bindings::traits::vector_traits<V>::value_type val_t;
  typedef val_t& ref_t; 
  static ref_t elem (V& v, size_t i) { return v[i]; }
};

template <typename V>
struct vct_access_traits<V const> {
  typedef typename 
  boost::numeric::bindings::traits::vector_traits<V>::value_type val_t;
  typedef val_t ref_t; 
  static ref_t elem (V const& v, size_t i) { return v[i]; }
};

template <typename V>
inline
typename vct_access_traits<V>::ref_t elem_v (V& v, size_t i) {
  return vct_access_traits<V>::elem (v, i); 
}

// initialization:
struct ident {
  size_t operator() (size_t i) const { return i; } 
}; 
struct iplus1 {
  size_t operator() (size_t i) const { return i + 1; } 
}; 
template <typename T>
struct const_val {
  T val;
  const_val (T v = T()) : val (v) {}
  T operator() (size_t) const { return val; } 
  T operator() (size_t, size_t) const { return val; } 
};  
struct kpp {
  size_t val; 
  kpp (size_t v = 0) : val (v) {}
  size_t operator() (size_t) { 
    size_t tmp = val;
    ++val; 
    return tmp; 
  } 
  size_t operator() (size_t, size_t) { 
    size_t tmp = val;
    ++val; 
    return tmp; 
  } 
}; 
template <typename T>
struct times_plus { 
  T s, a, b; 
  times_plus (T ss, T aa = T(), T bb = T()) : s (ss), a (aa), b (bb) {} 
  T operator() (size_t i) const { 
    return s * (T) i + a; 
  } 
  T operator() (size_t i, size_t j) const { 
    return s * ((T) i + a) + (T) j + b; 
  } 
}; 

template <typename F, typename V>
void init_v (V& v, F f = F()) {
  size_t sz 
    = boost::numeric::bindings::traits::vector_traits<V>::size (v);
  for (std::size_t i = 0; i < sz; ++i) 
    elem_v (v, i) = f (i); 
}

// printing: 
template <typename V>
void print_v (V const& v, char const* ch = 0) {
  if (ch)
    std::cout << ch << ": "; 
  size_t sz 
    = boost::numeric::bindings::traits::vector_traits<V const>::size (v);
  for (std::size_t i = 0; i < sz; ++i)
    std::cout << elem_v (v, i) << " ";
  std::cout << std::endl; 
}


/////////////////////////////////////
// matrices 

// element access: 
template <typename M>
struct matr_access_traits {
  typedef typename 
  M::reference ref_t;
  //boost::numeric::bindings::traits::matrix_traits<M>::value_type val_t;
  //typedef val_t& ref_t; 
  static ref_t elem (M& m, size_t i, size_t j) { return m (i, j); }
};

template <typename M>
struct matr_access_traits<M const> {
  typedef typename 
  boost::numeric::bindings::traits::matrix_traits<M>::value_type val_t;
  typedef val_t ref_t; 
  static ref_t elem (M const& m, size_t i, size_t j) { return m (i, j); }
};

template <typename M>
inline
typename matr_access_traits<M>::ref_t elem_m (M& m, size_t i, size_t j) {
  return matr_access_traits<M>::elem (m, i, j); 
}

// initialization: 
struct rws {
  size_t operator() (size_t i, size_t) const { return i; } 
}; 
struct rws1 {
  size_t operator() (size_t i, size_t) const { return i + 1; } 
}; 
struct cls {
  size_t operator() (size_t, size_t j) const { return j; } 
}; 
struct cls1 {
  size_t operator() (size_t, size_t j) const { return j + 1; } 
}; 

template <typename F, typename M>
void init_m (M& m, F f = F()) {
  size_t sz1 
    = boost::numeric::bindings::traits::matrix_traits<M>::num_rows (m);
  size_t sz2
    = boost::numeric::bindings::traits::matrix_traits<M>::num_columns (m);
  for (std::size_t i = 0; i < sz1; ++i) 
    for (std::size_t j = 0; j < sz2; ++j) 
      elem_m (m, i, j) = f (i, j); 
}

template <typename M>
void init_symm (M& m, char uplo = 'f') {
  size_t n 
    = boost::numeric::bindings::traits::matrix_traits<M>::num_rows (m);
  for (size_t i = 0; i < n; ++i) {
    elem_m (m, i, i) = n;
    for (size_t j = i + 1; j < n; ++j) {
      if (uplo == 'u' || uplo == 'U')
        elem_m (m, i, j) = n - (j - i);
      else if (uplo == 'l' || uplo == 'L')
        elem_m (m, j, i) = n - (j - i);
      else 
        elem_m (m, i, j) = elem_m (m, j, i) = n - (j - i);
    }
  }
}

// printing: 
template <typename M>
void print_m (M const& m, char const* ch = 0) {
  if (ch)
    std::cout << ch << ":\n"; 
  size_t sz1 
    = boost::numeric::bindings::traits::matrix_traits<M const>::num_rows (m);
  size_t sz2
    = boost::numeric::bindings::traits::matrix_traits<M const>::num_columns (m);
  for (std::size_t i = 0 ; i < sz1 ; ++i) {
    for (std::size_t j = 0 ; j < sz2 ; ++j) 
      std::cout << elem_m (m, i, j) << " ";
    std::cout << std::endl; 
  }
  //  std::cout << std::endl; 
}

template <typename M>
void print_m_data (M const& m, char const* ch = 0) {
  if (ch)
    std::cout << ch << " data:\n"; 
  size_t sz = 
    boost::numeric::bindings::traits::matrix_traits<M const>
    ::storage_size (m); 
  typename 
    boost::numeric::bindings::traits::matrix_traits<M const>::pointer st = 
    boost::numeric::bindings::traits::matrix_traits<M const>::storage (m); 
  for (std::size_t i = 0 ; i < sz ; ++i, ++st) 
      std::cout << *st << " ";
  std::cout << std::endl; 
}

#endif
