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

namespace RDKit {

template <class T>
inline const char *GetTypeName() {
  // PRECONDITION(0, "Unregistered c++ type");
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
  } catch (boost::bad_any_cast &) {
    return false;
  }
  return true;
}

template <class T>
boost::python::dict GetPropsAsDict(const T &obj, bool includePrivate,
                                   bool includeComputed) {
  boost::python::dict dict;
  // precedence double, int, unsigned, std::vector<double>,
  // std::vector<int>, std::vector<unsigned>, string
  STR_VECT keys = obj.getPropList(includePrivate, includeComputed);
  for (size_t i = 0; i < keys.size(); ++i) {
    if (AddToDict<int>(obj, dict, keys[i])) {
      continue;
    }
    if (AddToDict<unsigned int>(obj, dict, keys[i])) {
      continue;
    }
    if (AddToDict<bool>(obj, dict, keys[i])) {
      continue;
    }
    if (AddToDict<double>(obj, dict, keys[i])) {
      continue;
    }
    if (AddToDict<std::vector<int>>(obj, dict, keys[i])) {
      continue;
    }
    if (AddToDict<std::vector<unsigned int>>(obj, dict, keys[i])) {
      continue;
    }
    if (AddToDict<std::vector<double>>(obj, dict, keys[i])) {
      continue;
    }
    if (AddToDict<std::vector<std::string>>(obj, dict, keys[i])) {
      continue;
    }
    if (AddToDict<std::string>(obj, dict, keys[i])) {
      continue;
    }
  }
  return dict;
}

template <class RDOb, class T>
T GetProp(RDOb *ob, const char *key) {
  T res;
  try {
    if (!ob->getPropIfPresent(key, res)) {
      PyErr_SetString(PyExc_KeyError, key);
      throw python::error_already_set();
    }
    return res;
  } catch (const std::exception &e) {
    throw ValueErrorException(std::string("key `") + key +
                              "` exists but does not result in " +
                              GetTypeName<T>() + " reason: " + e.what());
  }

  return res;
}

template <class RDOb>
int MolHasProp(const RDOb &mol, const char *key) {
  int res = mol.hasProp(key);
  // std::cout << "key: "  << key << ": " << res << std::endl;
  return res;
}

template <class RDOb, class T>
void MolSetProp(const RDOb &mol, const char *key, const T &val,
                bool computed = false) {
  mol.setProp(key, val, computed);
}

template <class RDOb>
void MolClearProp(const RDOb &mol, const char *key) {
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
