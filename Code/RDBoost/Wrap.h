//
// Copyright (c) 2003-2006 greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _RD_WRAP_H_
#define _RD_WRAP_H_

//
// Generic Wrapper utility functionality
//
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "list_indexing_suite.hpp"
#include <vector>
#include <list>
#include <iostream>
#include "Exceptions.h"

namespace python = boost::python;

void throw_index_error(int key);  //!< construct and throw an \c IndexError
void throw_value_error(const std::string err); //!< construct and throw a \c ValueError
void translate_index_error(IndexErrorException const&e);
void translate_value_error(ValueErrorException const&e);

//! \brief Registers a templated converter for returning \c vectors of a 
//!        particular type.
//! This should be used instead of calling \c vector_to_python<T>()
//!    directly because this will catch the appropriate exception if
//!    the specified converter has already been registered.
template <typename T>
void RegisterVectorConverter(bool noproxy=false){
  std::string name = "_vect";
  name += typeid(T).name();

  if(noproxy){
    python::class_< std::vector<T> >(name.c_str())
      .def(python::vector_indexing_suite< std::vector<T>, 1>());
  } else {
    python::class_< std::vector<T> >(name.c_str())
      .def(python::vector_indexing_suite< std::vector<T> >());
  }
}


//! \brief Registers a templated converter for returning \c lists of a 
//!        particular type.
//! This should be used instead of calling \c list_to_python<T>()
//!    directly because this will catch the appropriate exception if
//!    the specified converter has already been registered.
template <typename T>
void RegisterListConverter(bool noproxy=false){
  std::string name = "_list";
  name += typeid(T).name();

  if(noproxy){
    python::class_< std::list<T> >(name.c_str())
      .def(python::list_indexing_suite< std::list<T>, 1 >());
  } else {
    python::class_< std::list<T> >(name.c_str())
      .def(python::list_indexing_suite< std::list<T> >());

  }
}
#endif
