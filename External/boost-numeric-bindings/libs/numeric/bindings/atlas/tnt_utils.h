
#ifndef F_TNT_UTILS_H
#define F_TNT_UTILS_H

#include "utils.h" 

template <typename T>
struct vct_access_traits<TNT::Fortran_Array1D<T> > {
  typedef T& ref_t; 
  typedef TNT::Fortran_Array1D<T> v_t;
  static ref_t elem (v_t& v, size_t i) { return v (i + 1); }
};
template <typename T>
struct vct_access_traits<TNT::Fortran_Array1D<T> const> {
  typedef T ref_t; 
  typedef TNT::Fortran_Array1D<T> v_t;
  static ref_t elem (v_t const& v, size_t i) { return v (i + 1); 
  }
};

template <typename T>
struct matr_access_traits<TNT::Array2D<T> > {
  typedef T& ref_t; 
  typedef TNT::Array2D<T> m_t;
  static ref_t elem (m_t& m, size_t i, size_t j) { 
#ifndef BOOST_MSVC
    return m[i][j];
#else
    return static_cast<T**> (m)[i][j];
#endif 
  }
};
template <typename T>
struct matr_access_traits<TNT::Array2D<T> const> {
  typedef T ref_t; 
  typedef TNT::Array2D<T> m_t;
  static ref_t elem (m_t const& m, size_t i, size_t j) { return m[i][j]; }
};

template <typename T>
struct matr_access_traits<TNT::Fortran_Array2D<T> > {
  typedef T& ref_t; 
  typedef TNT::Fortran_Array2D<T> m_t;
  static ref_t elem (m_t& m, size_t i, size_t j) { return m (i + 1, j + 1); }
};
template <typename T>
struct matr_access_traits<TNT::Fortran_Array2D<T> const> {
  typedef T ref_t; 
  typedef TNT::Fortran_Array2D<T> m_t;
  static ref_t elem (m_t const& m, size_t i, size_t j) { 
    return m (i + 1, j + 1); 
  }
};

#endif 
