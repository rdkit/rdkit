//
// Copyright (c) 2003-2008 greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef _RD_WRAP_H_
#define _RD_WRAP_H_

#include <RDGeneral/Invariant.h>

#include <RDGeneral/BoostStartInclude.h>
//
// Generic Wrapper utility functionality
//
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/dynamic_bitset.hpp>
#include <memory>

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <cstdint>
#include "list_indexing_suite.hpp"
#include <RDGeneral/BoostEndInclude.h>

#include <vector>
#include <list>
#include <iostream>
#include <RDGeneral/Exceptions.h>

namespace python = boost::python;

RDKIT_RDBOOST_EXPORT void throw_index_error(
    int key);  //!< construct and throw an \c IndexError
RDKIT_RDBOOST_EXPORT void throw_value_error(
    const std::string err);  //!< construct and throw a \c ValueError
RDKIT_RDBOOST_EXPORT void throw_key_error(
    const std::string key);  //!< construct and throw a \c KeyError
RDKIT_RDBOOST_EXPORT void translate_index_error(IndexErrorException const &e);
RDKIT_RDBOOST_EXPORT void translate_value_error(ValueErrorException const &e);
RDKIT_RDBOOST_EXPORT void translate_key_error(KeyErrorException const &e);

#ifdef INVARIANT_EXCEPTION_METHOD
RDKIT_RDBOOST_EXPORT void throw_runtime_error(
    const std::string err);  //!< construct and throw a \c ValueError
RDKIT_RDBOOST_EXPORT void translate_invariant_error(Invar::Invariant const &e);
#endif
//! \brief Registers a templated converter for returning \c vectors of a
//!        particular type.
//! This should be used instead of calling \c vector_to_python<T>()
//!    directly because this will catch the appropriate exception if
//!    the specified converter has already been registered.
template <typename T>
void RegisterVectorConverter(const char *name, bool noproxy = false) {
  if (noproxy) {
    python::class_<std::vector<T>>(name).def(
        python::vector_indexing_suite<std::vector<T>, 1>());
  } else {
    python::class_<std::vector<T>>(name).def(
        python::vector_indexing_suite<std::vector<T>>());
  }
}

template <typename T>
void RegisterVectorConverter(bool noproxy = false) {
  std::string name = "_vect";
  name += typeid(T).name();

  RegisterVectorConverter<T>(name.c_str(), noproxy);
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
    python::class_<std::list<T>>(name.c_str())
        .def(python::list_indexing_suite<std::list<T>, 1>());
  } else {
    python::class_<std::list<T>>(name.c_str())
        .def(python::list_indexing_suite<std::list<T>>());
  }
}

template <typename T>
std::unique_ptr<std::vector<T>> pythonObjectToVect(const python::object &obj,
                                                   T maxV) {
  std::unique_ptr<std::vector<T>> res;
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
std::unique_ptr<std::vector<T>> pythonObjectToVect(const python::object &obj) {
  std::unique_ptr<std::vector<T>> res;
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

RDKIT_RDBOOST_EXPORT boost::dynamic_bitset<> pythonObjectToDynBitset(const python::object &obj,
                                                   boost::dynamic_bitset<>::size_type maxV);


// Quiet warnings on GCC
#if defined(__GNUC__) || defined(__GNUG__)
#define RDUNUSED __attribute__((__unused__))
#else
#define RDUNUSED
#endif

class PyGILStateHolder {
public:
  PyGILStateHolder() :
    d_gstate(PyGILState_Ensure()) {}
  ~PyGILStateHolder() {
    PyGILState_Release(d_gstate);
  }
private:
  PyGILState_STATE d_gstate;
};

#ifdef RDK_THREADSAFE_SSS
// Release the Global Interpreter lock at certain places
//  on construction - release the lock
//  on destruction - grab the lock
//  no entry into the python interpreter can be performed
//   between releasing and grabbing the lock
class RDKIT_RDBOOST_EXPORT RDUNUSED NOGIL {
 public:
  inline NOGIL() { m_thread_state = PyEval_SaveThread(); }

  inline ~NOGIL() {
    PyEval_RestoreThread(m_thread_state);
    m_thread_state = nullptr;
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

/// Awesome StackOverflow response:
/// http://stackoverflow.com/questions/15842126/feeding-a-python-list-into-a-function-taking-in-a-vector-with-boost-python
/// I know a lot more about how boost works.
/// @brief Type that allows for registration of conversions from
///        python iterable types.
struct iterable_converter {
  /// @note Registers converter from a python iterable type to the
  ///       provided type.
  template <typename Container>
  iterable_converter &from_python() {
    boost::python::converter::registry::push_back(
        &iterable_converter::convertible,
        &iterable_converter::construct<Container>,
        boost::python::type_id<Container>());

    // Support chaining.
    return *this;
  }

  /// @brief Check if PyObject is iterable.
  static void *convertible(PyObject *object) {
    return PyObject_GetIter(object) ? object : nullptr;
  }

  /// @brief Convert iterable PyObject to C++ container type.
  ///
  /// Container Concept requirements:
  ///
  ///   * Container::value_type is CopyConstructable.
  ///   * Container can be constructed and populated with two iterators.
  ///     I.e. Container(begin, end)
  template <typename Container>
  static void construct(
      PyObject *object,
      boost::python::converter::rvalue_from_python_stage1_data *data) {
    namespace python = boost::python;
    // Object is a borrowed reference, so create a handle indicting it is
    // borrowed for proper reference counting.
    python::handle<> handle(python::borrowed(object));

    // Obtain a handle to the memory block that the converter has allocated
    // for the C++ type.
    typedef python::converter::rvalue_from_python_storage<Container>
        storage_type;
    void *storage = reinterpret_cast<storage_type *>(data)->storage.bytes;

    typedef python::stl_input_iterator<typename Container::value_type> iterator;

    // Allocate the C++ type into the converter's memory block, and assign
    // its handle to the converter's convertible variable.  The C++
    // container is populated by passing the begin and end iterators of
    // the python object to the container's constructor.
    new (storage) Container(iterator(python::object(handle)),  // begin
                            iterator());                       // end
    data->convertible = storage;
  }
};

#endif
