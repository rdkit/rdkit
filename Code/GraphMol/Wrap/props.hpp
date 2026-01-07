//  Copyright (c) 2015-2018, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#ifndef RD_WRAPPED_PROPS_H
#define RD_WRAPPED_PROPS_H

#include <RDBoost/python.h>
#include <RDBoost/pyint_api.h>
#include <RDBoost/Wrap.h>
#include <RDGeneral/Dict.h>
#include <algorithm>
#include <boost/algorithm/string.hpp>

namespace RDKit {

template <class T>
inline const char *GetTypeName() {
  return "unregistered C++ type";
}

template <>
inline const char *GetTypeName<double>() {
  return "a double value";
}
template <>
inline const char *GetTypeName<int>() {
  return "an integer value";
}
template <>
inline const char *GetTypeName<unsigned int>() {
  return "an unsigned integer value";
}
template <>
inline const char *GetTypeName<bool>() {
  return "a True or False value";
}

template <class T, class U>
bool AddToDict(const U &ob, boost::python::dict &dict, const std::string &key) {
  T res;
  try {
    if (ob.getPropIfPresent(key, res)) {
      dict[key] = res;
    }
  } catch (std::bad_any_cast &) {
    return false;
  }
  return true;
}

const std::string getPropsAsDictDocString =
    "Returns a dictionary populated with properties.\n"
    "When possible, string values will be converted to integers or doubles (trimming if necessary)\n"

    " n.b. Some properties are not able to be converted to python "
    "types.\n\n"
    "  ARGUMENTS:\n"
    "    - includePrivate: (optional) toggles inclusion of private "
    "properties in the result set.\n"
    "                      Defaults to False.\n"
    "    - includeComputed: (optional) toggles inclusion of computed "
    "properties in the result set.\n"
    "                      Defaults to False.\n\n"
    "  RETURNS: a dictionary\n";

template <class T>
boost::python::dict GetPropsAsDict(const T &obj, bool includePrivate,
                                   bool includeComputed,
                                   bool autoConvertStrings = true) {
  boost::python::dict dict;
  STR_VECT keys = obj.getPropList(includePrivate, includeComputed);

  for (const auto &key : keys) {
    if (key == "__computedProps") {
      continue;
    }

    bool found = false;
    // Try each type with exception handling
    try {
      int v;
      if (obj.getPropIfPresent(key, v)) {
        dict[key] = v;
        found = true;
      }
    } catch (std::bad_any_cast &) {}
    if (found) continue;

    try {
      double v;
      if (obj.getPropIfPresent(key, v)) {
        dict[key] = v;
        found = true;
      }
    } catch (std::bad_any_cast &) {}
    if (found) continue;

    try {
      bool v;
      if (obj.getPropIfPresent(key, v)) {
        dict[key] = v;
        found = true;
      }
    } catch (std::bad_any_cast &) {}
    if (found) continue;

    try {
      unsigned int v;
      if (obj.getPropIfPresent(key, v)) {
        dict[key] = v;
        found = true;
      }
    } catch (std::bad_any_cast &) {}
    if (found) continue;

    try {
      float v;
      if (obj.getPropIfPresent(key, v)) {
        dict[key] = v;
        found = true;
      }
    } catch (std::bad_any_cast &) {}
    if (found) continue;

    // Try vectors before string to avoid implicit vector->string conversion
    try {
      std::vector<double> v;
      if (obj.getPropIfPresent(key, v)) {
        dict[key] = v;
        found = true;
      }
    } catch (std::bad_any_cast &) {}
    if (found) continue;

    try {
      std::vector<float> v;
      if (obj.getPropIfPresent(key, v)) {
        dict[key] = v;
        found = true;
      }
    } catch (std::bad_any_cast &) {}
    if (found) continue;

    try {
      std::vector<int> v;
      if (obj.getPropIfPresent(key, v)) {
        dict[key] = v;
        found = true;
      }
    } catch (std::bad_any_cast &) {}
    if (found) continue;

    try {
      std::vector<unsigned int> v;
      if (obj.getPropIfPresent(key, v)) {
        dict[key] = v;
        found = true;
      }
    } catch (std::bad_any_cast &) {}
    if (found) continue;

    try {
      std::vector<std::string> v;
      if (obj.getPropIfPresent(key, v)) {
        dict[key] = v;
        found = true;
      }
    } catch (std::bad_any_cast &) {}
    if (found) continue;

    // Try string last to avoid catching vector->string conversions
    try {
      std::string v;
      if (obj.getPropIfPresent(key, v)) {
        if (autoConvertStrings) {
          auto trimmed = v;
          boost::trim(trimmed);
          int iconv;
          if (boost::conversion::try_lexical_convert(trimmed, iconv)) {
            dict[key] = iconv;
            found = true;
          } else {
            double dconv;
            if (boost::conversion::try_lexical_convert(trimmed, dconv)) {
              dict[key] = dconv;
              found = true;
            }
          }
        }
        if (!found) {
          dict[key] = v;
          found = true;
        }
      }
    } catch (std::bad_any_cast &) {}
  }

  return dict;
}

static PyObject *rawPy(python::object &&pyobj) {
  Py_INCREF(pyobj.ptr());
  return pyobj.ptr();
}

template <class T>
PyObject *rawPy(T &&thing) {
  return rawPy(python::object(thing));
}

template <class RDOb, class T>
PyObject* GetProp(const RDOb *ob, const std::string &key) {
  T res;
  try {
    if (!ob->getPropIfPresent(key, res)) {
      PyErr_SetString(PyExc_KeyError, key.c_str());
      return nullptr;
    }
  } catch (const std::exception &e) {
    auto msg = std::string("key `") + key +
                              "` exists but does not result in " +
                              GetTypeName<T>() + " reason: " + e.what();
    PyErr_SetString(PyExc_ValueError, msg.c_str());
    return nullptr;
  }
  return rawPy(std::move(res));
}

template <class RDOb>
python::object autoConvertString(const RDOb *ob, const std::string &key) {
  int ivalue;
  double dvalue;
  std::string svalue;

  if (ob->getPropIfPresent(key, ivalue))
    return python::object(ivalue);
  else if (ob->getPropIfPresent(key, dvalue))
    return python::object(dvalue);
  else if (ob->getPropIfPresent(key, svalue))
    return python::object(svalue);

  return python::object();
}

template <class RDOb>
PyObject *GetPyProp(const RDOb *obj, const std::string &key, bool autoConvert) {
  // Try different types with exception handling
  try {
    int v;
    if (obj->getPropIfPresent(key, v)) {
      return rawPy(v);
    }
  } catch (std::bad_any_cast &) {}

  try {
    double v;
    if (obj->getPropIfPresent(key, v)) {
      return rawPy(v);
    }
  } catch (std::bad_any_cast &) {}

  try {
    bool v;
    if (obj->getPropIfPresent(key, v)) {
      return rawPy(v);
    }
  } catch (std::bad_any_cast &) {}

  try {
    unsigned int v;
    if (obj->getPropIfPresent(key, v)) {
      return rawPy(v);
    }
  } catch (std::bad_any_cast &) {}

  try {
    float v;
    if (obj->getPropIfPresent(key, v)) {
      return rawPy(v);
    }
  } catch (std::bad_any_cast &) {}

  // Try vectors before string to avoid implicit vector->string conversion
  try {
    std::vector<double> v;
    if (obj->getPropIfPresent(key, v)) {
      return rawPy(v);
    }
  } catch (std::bad_any_cast &) {}

  try {
    std::vector<float> v;
    if (obj->getPropIfPresent(key, v)) {
      return rawPy(v);
    }
  } catch (std::bad_any_cast &) {}

  try {
    std::vector<int> v;
    if (obj->getPropIfPresent(key, v)) {
      return rawPy(v);
    }
  } catch (std::bad_any_cast &) {}

  try {
    std::vector<unsigned int> v;
    if (obj->getPropIfPresent(key, v)) {
      return rawPy(v);
    }
  } catch (std::bad_any_cast &) {}

  try {
    std::vector<std::string> v;
    if (obj->getPropIfPresent(key, v)) {
      return rawPy(v);
    }
  } catch (std::bad_any_cast &) {}

  // Try string last to avoid catching vector->string conversions
  try {
    std::string v;
    if (obj->getPropIfPresent(key, v)) {
      if (autoConvert) {
        auto trimmed = v;
        boost::trim(trimmed);
        int iconv;
        if (boost::conversion::try_lexical_convert(trimmed, iconv)) {
          return rawPy(iconv);
        }
        double dconv;
        if (boost::conversion::try_lexical_convert(trimmed, dconv)) {
          return rawPy(dconv);
        }
      }
      return rawPy(v);
    }
  } catch (std::bad_any_cast &) {}

  // Property not found
  PyErr_SetString(PyExc_KeyError, key.c_str());
  return nullptr;
}

// Return policy for functions that directly return a PyObject* and
// are fully responsible for setting the Python error state.
struct return_pyobject_passthrough {
  template <class T>
  struct apply {
    struct type {
      static bool convertible() { return true; }

      PyObject *operator()(PyObject *inner) const { return inner; }
#ifndef BOOST_PYTHON_NO_PY_SIGNATURES
      PyTypeObject const *get_pytype() const {
        return boost::python::converter::expected_pytype_for_arg<
            T>::get_pytype();
      }
#endif
    };
  };
};

template <class RDOb>
int MolHasProp(const RDOb &mol, const std::string &key) {
  int res = mol.hasProp(key);
  return res;
}

template <class RDOb, class T>
void MolSetProp(const RDOb &mol, const std::string &key, const T &val,
                bool computed = false) {
  mol.setProp(key, val, computed);
}

template <class RDOb>
void MolClearProp(const RDOb &mol, const std::string &key) {
  if (!mol.hasProp(key)) {
    return;
  }
  mol.clearProp(key);
}

template <class RDOb>
void MolClearComputedProps(const RDOb &mol) {
  mol.clearComputedProps();
}
}  // namespace RDKit

#endif
