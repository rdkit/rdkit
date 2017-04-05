//
// Copyright (c) 2003-2008 greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_WRAP_H_
#define _RD_WRAP_H_

#include <RDGeneral/Invariant.h>

#include <RDGeneral/BoostStartInclude.h>
//
// Generic Wrapper utility functionality
//
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <memory>
// ha ha... we don't actually want to use this, but we need to keep support old
// C++ versions:
#if __cplusplus < 201103L
#define rdk_auto_ptr std::auto_ptr
#else
template <typename T>
using rdk_auto_ptr = std::unique_ptr<T>;
#endif

// code for windows DLL handling taken from
// http://www.boost.org/more/separate_compilation.html
#include <boost/config.hpp>

#ifdef BOOST_HAS_DECLSPEC  // defined in config system
// we need to import/export our code only if the user has specifically
// asked for it by defining either BOOST_ALL_DYN_LINK if they want all boost
// libraries to be dynamically linked, or RDKIT_WRAP_DYN_LINK
// if they want just this one to be dynamically liked:
#if defined(BOOST_ALL_DYN_LINK) || defined(RDKIT_WRAP_DYN_LINK)
// export if this is our own source, otherwise import:
#ifdef RDKIT_WRAP_SOURCE
#define RDKIT_WRAP_DECL __declspec(dllexport)
#else
#define RDKIT_WRAP_DECL __declspec(dllimport)
#endif  // RDKIT_WRAP_SOURCE
#endif  // DYN_LINK
#endif  // BOOST_HAS_DECLSPEC
//
// if RDKIT_WRAP_DECL isn't defined yet define it now:
#ifndef RDKIT_WRAP_DECL
#define RDKIT_WRAP_DECL
#endif

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/cstdint.hpp>
#include "list_indexing_suite.hpp"
#include <RDGeneral/BoostEndInclude.h>

#include <vector>
#include <list>
#include <iostream>
#include <RDGeneral/Exceptions.h>

namespace python = boost::python;

RDKIT_WRAP_DECL void throw_index_error(
    int key);  //!< construct and throw an \c IndexError
RDKIT_WRAP_DECL void throw_value_error(
    const std::string err);  //!< construct and throw a \c ValueError
RDKIT_WRAP_DECL void translate_index_error(IndexErrorException const &e);
RDKIT_WRAP_DECL void translate_value_error(ValueErrorException const &e);

#ifdef INVARIANT_EXCEPTION_METHOD
RDKIT_WRAP_DECL void throw_runtime_error(
    const std::string err);  //!< construct and throw a \c ValueError
RDKIT_WRAP_DECL void translate_invariant_error(Invar::Invariant const &e);
#endif
//! \brief Registers a templated converter for returning \c vectors of a
//!        particular type.
//! This should be used instead of calling \c vector_to_python<T>()
//!    directly because this will catch the appropriate exception if
//!    the specified converter has already been registered.
template <typename T>
void RegisterVectorConverter(bool noproxy = false) {
  std::string name = "_vect";
  name += typeid(T).name();

  if (noproxy) {
    python::class_<std::vector<T> >(name.c_str())
        .def(python::vector_indexing_suite<std::vector<T>, 1>());
  } else {
    python::class_<std::vector<T> >(name.c_str())
        .def(python::vector_indexing_suite<std::vector<T> >());
  }
}

//! \brief Registers a templated converter for returning \c lists of a
//!        particular type.
//! This should be used instead of calling \c list_to_python<T>()
//!    directly because this will catch the appropriate exception if
//!    the specified converter has already been registered.
template <typename T>
void RegisterListConverter(bool noproxy = false) {
  std::string name = "_list";
  name += typeid(T).name();

  if (noproxy) {
    python::class_<std::list<T> >(name.c_str())
        .def(python::list_indexing_suite<std::list<T>, 1>());
  } else {
    python::class_<std::list<T> >(name.c_str())
        .def(python::list_indexing_suite<std::list<T> >());
  }
}

template <typename T>
rdk_auto_ptr<std::vector<T> > pythonObjectToVect(const python::object &obj,
                                                 T maxV) {
  rdk_auto_ptr<std::vector<T> > res;
  if (obj) {
    res.reset(new std::vector<T>);
    python::stl_input_iterator<T> beg(obj), end;
    while (beg != end) {
      T v = *beg;
      if (v >= maxV) {
        throw_value_error("list element larger than allowed value");
      }
      res->push_back(v);
      ++beg;
    }
  }
  return res;
}
template <typename T>
rdk_auto_ptr<std::vector<T> > pythonObjectToVect(const python::object &obj) {
  rdk_auto_ptr<std::vector<T> > res;
  if (obj) {
    res.reset(new std::vector<T>);
    unsigned int nFrom = python::extract<unsigned int>(obj.attr("__len__")());
    for (unsigned int i = 0; i < nFrom; ++i) {
      T v = python::extract<T>(obj[i]);
      res->push_back(v);
    }
  }
  return res;
}
template <typename T>
void pythonObjectToVect(const python::object &obj, std::vector<T> &res) {
  if (obj) {
    res.clear();
    python::stl_input_iterator<T> beg(obj), end;
    while (beg != end) {
      T v = *beg;
      res.push_back(v);
      ++beg;
    }
  }
}

// Quiet warnings on GCC
#if defined(__GNUC__) || defined(__GNUG__)
#define RDUNUSED __attribute__((__unused__))
#else
#define RDUNUSED
#endif

#ifdef RDK_THREADSAFE_SSS
// Release the Global Interpreter lock at certain places
//  on construction - release the lock
//  on destruction - grab the lock
//  no entry into the python interpreter can be performed
//   between releasing and grabbing the lock
class RDUNUSED NOGIL {
 public:
  inline NOGIL() { m_thread_state = PyEval_SaveThread(); }

  inline ~NOGIL() {
    PyEval_RestoreThread(m_thread_state);
    m_thread_state = NULL;
  }

 private:
  PyThreadState *m_thread_state;
};
#else
// Never release the lock when not compiling thread-safe
struct RDUNUSED NOGIL {};
#endif

// -------------------
// This block was adapted from this mailing list post by Matthew Scouten:
// https://mail.python.org/pipermail/cplusplus-sig/2009-May/014505.html
// Matthew credits Hans Meine in his post.
template <class T>
inline PyObject *managingPyObject(T *p) {
  return typename python::manage_new_object::apply<T *>::type()(p);
}

template <class Copyable>
python::object generic__copy__(python::object copyable) {
  Copyable *newCopyable(
      new Copyable(python::extract<const Copyable &>(copyable)));
  python::object result(
      python::detail::new_reference(managingPyObject(newCopyable)));

  python::extract<python::dict>(result.attr("__dict__"))().update(
      copyable.attr("__dict__"));

  return result;
}

template <class Copyable>
python::object generic__deepcopy__(python::object copyable, python::dict memo) {
  python::object copyMod = python::import("copy");
  python::object deepcopy = copyMod.attr("deepcopy");

  Copyable *newCopyable(
      new Copyable(python::extract<const Copyable &>(copyable)));
  python::object result(
      python::detail::new_reference(managingPyObject(newCopyable)));

  // HACK: copyableId shall be the same as the result of id(copyable) in Python
  // -
  // please tell me that there is a better way! (and which ;-p)
  size_t copyableId = (size_t)(copyable.ptr());
  memo[copyableId] = result;

  python::extract<python::dict>(result.attr("__dict__"))().update(deepcopy(
      python::extract<python::dict>(copyable.attr("__dict__"))(), memo));
  return result;
}
// -------------------

#endif
