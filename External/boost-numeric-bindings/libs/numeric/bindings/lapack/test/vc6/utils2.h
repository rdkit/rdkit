
#ifndef F_UTILS_H
#define F_UTILS_H

#include <stddef.h>
#include <iostream>
#include <boost/numeric/bindings/traits/traits.hpp>


///////////////////////////////
// vectors

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
    = boost::numeric::bindings::traits::vector_size (v);
  for (std::size_t i = 0; i < sz; ++i) 
    v [i] = f (i); 
}

// printing: 
template <typename V>
void print_v (V const& v, char const* ch = 0) {
  if (ch)
    std::cout << ch << ": "; 
  size_t sz 
    = boost::numeric::bindings::traits::vector_size (v);
  for (std::size_t i = 0; i < sz; ++i)
    std::cout << v [i] << " ";
  std::cout << std::endl; 
}


/////////////////////////////////////
// matrices 

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
    = boost::numeric::bindings::traits::matrix_size1 (m);
  size_t sz2
    = boost::numeric::bindings::traits::matrix_size2 (m);
  for (std::size_t i = 0; i < sz1; ++i) 
    for (std::size_t j = 0; j < sz2; ++j) 
      m (i, j) = f (i, j); 
}

template <typename M>
void init_symm (M& m, char uplo = 'f') {
  size_t n 
    = boost::numeric::bindings::traits::matrix_size1 (m);
  for (size_t i = 0; i < n; ++i) {
    m (i, i) = n;
    for (size_t j = i + 1; j < n; ++j) {
      if (uplo == 'u' || uplo == 'U')
        m (i, j) = n - (j - i);
      else if (uplo == 'l' || uplo == 'L')
        m (j, i) = n - (j - i);
      else 
        m (i, j) = m (j, i) = n - (j - i);
    }
  }
}

// printing: 
template <typename M>
void print_m (M const& m, char const* ch = 0) {
  if (ch)
    std::cout << ch << ":\n"; 
  size_t sz1 
    = boost::numeric::bindings::traits::matrix_size1 (m);
  size_t sz2
    = boost::numeric::bindings::traits::matrix_size2 (m);
  for (std::size_t i = 0 ; i < sz1 ; ++i) {
    for (std::size_t j = 0 ; j < sz2 ; ++j) 
      std::cout << m (i, j) << " ";
    std::cout << std::endl; 
  }
  //  std::cout << std::endl; 
}

#endif
