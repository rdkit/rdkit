//  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
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
#include <RDGeneral/export.h>
#ifndef RDKIT_RDANY_H
#define RDKIT_RDANY_H
#include <RDGeneral/BoostStartInclude.h>
#include <boost/any.hpp>
#include <boost/utility.hpp>
#include <boost/lexical_cast.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include "LocaleSwitcher.h"
#include "RDValue.h"
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
namespace RDKit {

// RDValue does not dynamically create POD types (kind of like
//  cdiggins::any)  However, it doesn't use RTTI type info
//  directly, it uses a companion short valued type
//  to determine what to do.
// For unregistered types, it falls back to boost::any.
//  The Size of an RDAny is (sizeof(double) + sizeof(short) == 10 bytes)
//
//   For the sake of compatibility, errors throw boost::bad_any_cast
//
// Examples:
//
//   RDAny v(2.);
//   v = 1;
//   std::vector<double> d;
//   v == d;
//   v.asDoubleVect().push_back(4.)
//   rdany_cast<std::vector<double>(v).push_back(4.)
//
//   Falls back to boost::any for non registered types
//   v = boost::shared_ptr<ROMol>(new ROMol(m));
//

// Safe container for RDValue -- cleans up memory and copy constructs
struct RDAny {
  RDValue m_value;

  RDAny() : m_value() {}
  template <class T>
  RDAny(const T &d) : m_value(d) {}
  /*
  explicit RDAny(bool v) : m_value(v) {}
  template <class T>
  explicit RDAny(std::vector<T> *v) : m_value(v) {}
  template <class T>
  explicit RDAny(const boost::shared_ptr<T> &v) : m_value(v) {}
  */
  RDAny(const RDAny &rhs) { copy_rdvalue(m_value, rhs.m_value); }

  ~RDAny() { RDValue::cleanup_rdvalue(m_value); }

  // For easy of use:
  //   RDAny v;
  //   v = 2.0;
  //   v = std::string("foo...");

  RDAny &operator=(const RDAny &rhs) {
    copy_rdvalue(m_value, rhs.m_value);
    return *this;
  }

  RDAny &operator=(float d) {
    RDValue::cleanup_rdvalue(m_value);
    m_value = RDValue(d);
    return *this;
  }

  RDAny &operator=(int d) {
    RDValue::cleanup_rdvalue(m_value);
    m_value = RDValue(d);
    return *this;
  }

  RDAny &operator=(unsigned int d) {
    RDValue::cleanup_rdvalue(m_value);
    m_value = RDValue(d);
    return *this;
  }

  RDAny &operator=(bool d) {
    RDValue::cleanup_rdvalue(m_value);
    m_value = RDValue(d);
    return *this;
  }

  RDAny &operator=(const std::string &d) {
    RDValue::cleanup_rdvalue(m_value);
    m_value = RDValue(d);
    return *this;
  }

  RDAny &operator=(const std::vector<double> &d) {
    RDValue::cleanup_rdvalue(m_value);
    m_value = RDValue(d);
    return *this;
  }

  RDAny &operator=(const std::vector<float> &d) {
    RDValue::cleanup_rdvalue(m_value);
    m_value = RDValue(d);
    return *this;
  }

  RDAny &operator=(const std::vector<int> &d) {
    RDValue::cleanup_rdvalue(m_value);
    m_value = RDValue(d);
    return *this;
  }

  RDAny &operator=(const std::vector<unsigned int> &d) {
    RDValue::cleanup_rdvalue(m_value);
    m_value = d;
    return *this;
  }

  RDAny &operator=(const std::vector<std::string> &d) {
    RDValue::cleanup_rdvalue(m_value);
    m_value = RDValue(d);
    return *this;
  }

  RDAny &operator=(const boost::any &d) {
    RDValue::cleanup_rdvalue(m_value);
    m_value = RDValue(d);  // new boost::any(d);
    return *this;
  }

  template <class T>
  RDAny &operator=(const T &d) {
    RDValue::cleanup_rdvalue(m_value);
    boost::any *v = new boost::any(d);
    m_value = RDValue(v);
    return *this;
  }
};

////////////////////////////////////////////////////////////////
// rdany_cast
////////////////////////////////////////////////////////////////

// Const Access
template <class T>
const T rdany_cast(const RDAny &d) {
  return rdvalue_cast<T>(d.m_value);
}

// Direct access
template <class T>
T rdany_cast(RDAny &d) {
  return rdvalue_cast<T>(d.m_value);
}

template <class T>
typename boost::enable_if<boost::is_arithmetic<T>, T>::type from_rdany(
    const RDAny &arg) {
  T res;
  if (arg.m_value.getTag() == RDTypeTag::StringTag) {
    Utils::LocaleSwitcher ls;
    try {
      res = rdany_cast<T>(arg);
    } catch (const boost::bad_any_cast &exc) {
      try {
        res = boost::lexical_cast<T>(rdany_cast<std::string>(arg));
      } catch (...) {
        throw exc;
      }
    }
  } else {
    res = rdany_cast<T>(arg);
  }
  return res;
}

template <class T>
typename boost::disable_if<boost::is_arithmetic<T>, T>::type from_rdany(
    const RDAny &arg) {
  return rdany_cast<T>(arg);
}

}  // namespace RDKit
#endif
