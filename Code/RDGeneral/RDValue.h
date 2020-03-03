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
#ifndef RDKIT_RDVALUE_H
#define RDKIT_RDVALUE_H

//#define UNSAFE_RDVALUE
#ifdef UNSAFE_RDVALUE
#include "RDValue-doublemagic.h"
#else
#include "RDValue-taggedunion.h"
#endif

namespace RDKit {
//  Common Casts (POD Casts are implementation dependent)
// string casts
template <>
inline std::string rdvalue_cast<std::string>(RDValue_cast_t v) {
  if (rdvalue_is<std::string>(v)) return *v.ptrCast<std::string>();
  throw boost::bad_any_cast();
}

template <>
inline std::string &rdvalue_cast<std::string &>(RDValue_cast_t v) {
  if (rdvalue_is<std::string>(v)) return *v.ptrCast<std::string>();
  throw boost::bad_any_cast();
}

// Special Vecor Casts
template <>
inline std::vector<double> rdvalue_cast<std::vector<double>>(RDValue_cast_t v) {
  if (rdvalue_is<std::vector<double>>(v))
    return *v.ptrCast<std::vector<double>>();
  throw boost::bad_any_cast();
}

template <>
inline std::vector<double> &rdvalue_cast<std::vector<double> &>(
    RDValue_cast_t v) {
  if (rdvalue_is<std::vector<double>>(v))
    return *v.ptrCast<std::vector<double>>();
  throw boost::bad_any_cast();
}

template <>
inline std::vector<float> rdvalue_cast<std::vector<float>>(RDValue_cast_t v) {
  if (rdvalue_is<std::vector<float>>(v))
    return *v.ptrCast<std::vector<float>>();
  throw boost::bad_any_cast();
}

template <>
inline std::vector<float> &rdvalue_cast<std::vector<float> &>(
    RDValue_cast_t v) {
  if (rdvalue_is<std::vector<float>>(v))
    return *v.ptrCast<std::vector<float>>();
  throw boost::bad_any_cast();
}

template <>
inline std::vector<std::string> rdvalue_cast<std::vector<std::string>>(
    RDValue_cast_t v) {
  if (rdvalue_is<std::vector<std::string>>(v))
    return *v.ptrCast<std::vector<std::string>>();
  throw boost::bad_any_cast();
}

template <>
inline std::vector<std::string> &rdvalue_cast<std::vector<std::string> &>(
    RDValue_cast_t v) {
  if (rdvalue_is<std::vector<std::string>>(v))
    return *v.ptrCast<std::vector<std::string>>();
  throw boost::bad_any_cast();
}

template <>
inline std::vector<int> rdvalue_cast<std::vector<int>>(RDValue_cast_t v) {
  if (rdvalue_is<std::vector<int>>(v)) return *v.ptrCast<std::vector<int>>();
  throw boost::bad_any_cast();
}

template <>
inline std::vector<int> &rdvalue_cast<std::vector<int> &>(RDValue_cast_t v) {
  if (rdvalue_is<std::vector<int>>(v)) return *v.ptrCast<std::vector<int>>();
  throw boost::bad_any_cast();
}

template <>
inline std::vector<unsigned int> rdvalue_cast<std::vector<unsigned int>>(
    RDValue_cast_t v) {
  if (rdvalue_is<std::vector<unsigned int>>(v))
    return *v.ptrCast<std::vector<unsigned int>>();
  throw boost::bad_any_cast();
}

template <>
inline std::vector<unsigned int> &rdvalue_cast<std::vector<unsigned int> &>(
    RDValue_cast_t v) {
  if (rdvalue_is<std::vector<unsigned int>>(v))
    return *v.ptrCast<std::vector<unsigned int>>();
  throw boost::bad_any_cast();
}

// Get boost any
template <>
inline boost::any rdvalue_cast<boost::any>(RDValue_cast_t v) {
  if (rdvalue_is<boost::any>(v)) {
    return *v.ptrCast<boost::any>();
  }
  throw boost::bad_any_cast();
}

template <>
inline boost::any &rdvalue_cast<boost::any &>(RDValue_cast_t v) {
  if (rdvalue_is<boost::any>(v)) {
    return *v.ptrCast<boost::any>();
  }
  throw boost::bad_any_cast();
}

template <>
inline const boost::any &rdvalue_cast<const boost::any &>(RDValue_cast_t v) {
  if (rdvalue_is<boost::any>(v)) {
    return *v.ptrCast<boost::any>();
  }
  throw boost::bad_any_cast();
}

/////////////////////////////////////////////////////////////////////////////////////
// lexical casts...
template <class T>
std::string vectToString(RDValue val) {
  const std::vector<T> &tv = rdvalue_cast<std::vector<T> &>(val);
  std::ostringstream sstr;
  sstr.imbue(std::locale("C"));
  sstr << std::setprecision(17);
  sstr << "[";
  std::copy(tv.begin(), tv.end(), std::ostream_iterator<T>(sstr, ","));
  sstr << "]";
  return sstr.str();
}

inline bool rdvalue_tostring(RDValue_cast_t val, std::string &res) {
  switch (val.getTag()) {
    case RDTypeTag::StringTag:
      res = rdvalue_cast<std::string>(val);
      break;
    case RDTypeTag::IntTag:
      res = boost::lexical_cast<std::string>(rdvalue_cast<int>(val));
      break;
    case RDTypeTag::DoubleTag: {
      Utils::LocaleSwitcher ls;  // for lexical cast...
      res = boost::lexical_cast<std::string>(rdvalue_cast<double>(val));
      break;
    }
    case RDTypeTag::UnsignedIntTag:
      res = boost::lexical_cast<std::string>(rdvalue_cast<unsigned int>(val));
      break;
#ifdef RDVALUE_HASBOOL
    case RDTypeTag::BoolTag:
      res = boost::lexical_cast<std::string>(rdvalue_cast<bool>(val));
      break;
#endif
    case RDTypeTag::FloatTag: {
      Utils::LocaleSwitcher ls;  // for lexical cast...
      res = boost::lexical_cast<std::string>(rdvalue_cast<float>(val));
      break;
    }
    case RDTypeTag::VecDoubleTag: {
      // vectToString uses std::imbue for locale
      res = vectToString<double>(val);
      break;
    }
    case RDTypeTag::VecFloatTag: {
      // vectToString uses std::imbue for locale
      res = vectToString<float>(val);
      break;
    }
    case RDTypeTag::VecIntTag:
      res = vectToString<int>(val);
      break;
    case RDTypeTag::VecUnsignedIntTag:
      res = vectToString<unsigned int>(val);
      break;
    case RDTypeTag::VecStringTag:
      res = vectToString<std::string>(val);
      break;
    case RDTypeTag::AnyTag: {
      Utils::LocaleSwitcher ls;  // for lexical cast...
      try {
        res = boost::any_cast<std::string>(rdvalue_cast<boost::any &>(val));
      } catch (const boost::bad_any_cast &) {
        if (rdvalue_cast<boost::any &>(val).type() == typeid(long)) {
          res = boost::lexical_cast<std::string>(
              boost::any_cast<long>(rdvalue_cast<boost::any &>(val)));
        } else if (rdvalue_cast<boost::any &>(val).type() ==
                   typeid(unsigned long)) {
          res = boost::lexical_cast<std::string>(
              boost::any_cast<unsigned long>(rdvalue_cast<boost::any &>(val)));
        } else {
          throw;
          return false;
        }
      }
      break;
    }
    default:
      res = "";
  }
  return true;
}

// from_rdvalue -> converts string values to appropriate types
template <class T>
typename boost::enable_if<boost::is_arithmetic<T>, T>::type from_rdvalue(
    RDValue_cast_t arg) {
  T res;
  if (arg.getTag() == RDTypeTag::StringTag) {
    Utils::LocaleSwitcher ls;
    try {
      res = rdvalue_cast<T>(arg);
    } catch (const boost::bad_any_cast &exc) {
      try {
        res = boost::lexical_cast<T>(rdvalue_cast<std::string>(arg));
      } catch (...) {
        throw exc;
      }
    }
  } else {
    res = rdvalue_cast<T>(arg);
  }
  return res;
}

template <class T>
typename boost::disable_if<boost::is_arithmetic<T>, T>::type from_rdvalue(
    RDValue_cast_t arg) {
  return rdvalue_cast<T>(arg);
}
}  // namespace RDKit
#endif
