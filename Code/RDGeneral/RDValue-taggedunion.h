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
#ifndef RDKIT_RDVALUE_TAGGED_UNION_H
#define RDKIT_RDVALUE_TAGGED_UNION_H

#include <cassert>
#include "Invariant.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cstdint>
#include <RDGeneral/BoostStartInclude.h>
#include <cstdint>
#include <any>
#include <boost/utility.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include "LocaleSwitcher.h"

#define RDVALUE_HASBOOL

namespace RDKit {

// RDValue does not dynamically create POD types (kind of like
//  cdiggins::any)  However, it doesn't use RTTI type info
//  directly, it uses a companion short valued type
//  to determine what to do.
// For unregistered types, it falls back to std::any.
//  The Size of an RDAny is (sizeof(double) + sizeof(short) == 10 bytes
//   (aligned to actually 16 so hard to pass as value type)
//
//   For the sake of compatibility, errors throw std::bad_any_cast
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
//   Falls back to std::any for non registered types
//   v = boost::shared_ptr<ROMol>(new ROMol(m));
//

// RDValue does not manange memory of non-pod data
//  this must be done externally (string, Any, vector...)
// Tagged union

namespace RDTypeTag {
const short EmptyTag = 0;
const short IntTag = 1;
const short DoubleTag = 2;
const short StringTag = 3;
const short FloatTag = 4;
const short BoolTag = 5;
const short UnsignedIntTag = 6;
const short AnyTag = 7;
const short VecDoubleTag = 8;
const short VecFloatTag = 9;
const short VecIntTag = 10;
const short VecUnsignedIntTag = 11;
const short VecStringTag = 12;
template <class T>
inline short GetTag() {
  return AnyTag;
}
template <>
inline short GetTag<double>() {
  return DoubleTag;
}
template <>
inline short GetTag<float>() {
  return FloatTag;
}
template <>
inline short GetTag<int>() {
  return IntTag;
}
template <>
inline short GetTag<unsigned int>() {
  return UnsignedIntTag;
}
template <>
inline short GetTag<bool>() {
  return BoolTag;
}
template <>
inline short GetTag<std::string>() {
  return StringTag;
}
template <>
inline short GetTag<std::vector<double>>() {
  return VecDoubleTag;
}
template <>
inline short GetTag<std::vector<float>>() {
  return VecFloatTag;
}
template <>
inline short GetTag<std::vector<int>>() {
  return VecIntTag;
}
template <>
inline short GetTag<std::vector<unsigned int>>() {
  return VecUnsignedIntTag;
}
template <>
inline short GetTag<std::vector<std::string>>() {
  return VecStringTag;
}
template <>
inline short GetTag<std::any>() {
  return AnyTag;
}

namespace detail {
union Value {
  double d;
  float f;
  int i;
  unsigned u;
  bool b;
  std::string *s;
  std::any *a;
  std::vector<double> *vd;
  std::vector<float> *vf;
  std::vector<int> *vi;
  std::vector<unsigned int> *vu;
  std::vector<std::string> *vs;

  inline Value() {}
  inline Value(double v) : d(v) {}
  inline Value(float v) : f(v) {}
  inline Value(int v) : i(v) {}
  inline Value(unsigned int v) : u(v) {}
  inline Value(bool v) : b(v) {}
  inline Value(std::string *v) : s(v) {}
  inline Value(std::any *v) : a(v) {}
  inline Value(std::vector<double> *v) : vd(v) {}
  inline Value(std::vector<float> *v) : vf(v) {}
  inline Value(std::vector<int> *v) : vi(v) {}
  inline Value(std::vector<unsigned int> *v) : vu(v) {}
  inline Value(std::vector<std::string> *v) : vs(v) {}
};

template <class T>
inline T *valuePtrCast(Value value) {
  return std::any_cast<T *>(*value.a);
}
template <>
inline std::any *valuePtrCast<std::any>(Value value) {
  return value.a;
}

template <>
inline std::string *valuePtrCast<std::string>(Value value) {
  return value.s;
}
template <>
inline std::vector<double> *valuePtrCast<std::vector<double>>(Value value) {
  return value.vd;
}
template <>
inline std::vector<float> *valuePtrCast<std::vector<float>>(Value value) {
  return value.vf;
}
template <>
inline std::vector<int> *valuePtrCast<std::vector<int>>(Value value) {
  return value.vi;
}
template <>
inline std::vector<unsigned int> *valuePtrCast<std::vector<unsigned int>>(
    Value value) {
  return value.vu;
}
template <>
inline std::vector<std::string> *valuePtrCast<std::vector<std::string>>(
    Value value) {
  return value.vs;
}
}  // namespace detail
}  // namespace RDTypeTag

struct RDValue {
  RDTypeTag::detail::Value value;
  short type;
  short reserved_tag = 0;  // 16 bit alignment

  inline RDValue() : value(0.0), type(RDTypeTag::EmptyTag) {}
  // Pod Style (Direct storage)
  inline RDValue(double v) : value(v), type(RDTypeTag::DoubleTag) {}
  inline RDValue(float v) : value(v), type(RDTypeTag::FloatTag) {}
  inline RDValue(int v) : value(v), type(RDTypeTag::IntTag) {}
  inline RDValue(unsigned v) : value(v), type(RDTypeTag::UnsignedIntTag) {}
  inline RDValue(bool v) : value(v), type(RDTypeTag::BoolTag) {}

  inline RDValue(std::any *v) : value(v), type(RDTypeTag::AnyTag) {}

  // Copies passed in pointers
  inline RDValue(const std::any &v)
      : value(new std::any(v)), type(RDTypeTag::AnyTag) {}
  inline RDValue(const std::string &v)
      : value(new std::string(v)), type(RDTypeTag::StringTag) {}
  template <class T>
  inline RDValue(const T &v)
      : value(new std::any(v)), type(RDTypeTag::AnyTag) {}

  inline RDValue(const std::vector<double> &v)
      : value(new std::vector<double>(v)), type(RDTypeTag::VecDoubleTag) {}
  inline RDValue(const std::vector<float> &v)
      : value(new std::vector<float>(v)), type(RDTypeTag::VecFloatTag) {}
  inline RDValue(const std::vector<int> &v)
      : value(new std::vector<int>(v)), type(RDTypeTag::VecIntTag) {}
  inline RDValue(const std::vector<unsigned int> &v)
      : value(new std::vector<unsigned int>(v)),
        type(RDTypeTag::VecUnsignedIntTag) {}
  inline RDValue(const std::vector<std::string> &v)
      : value(new std::vector<std::string>(v)), type(RDTypeTag::VecStringTag) {}

  short getTag() const { return type; }

  // ptrCast - unsafe, use rdvalue_cast instead.
  template <class T>
  inline T *ptrCast() const {
    return RDTypeTag::detail::valuePtrCast<T>(value);
  }

  // RDValue doesn't have an explicit destructor, it must
  //  be wrapped in a container.
  // The idea is that POD types don't need to be destroyed
  //  and this allows the container optimization possibilities.
  void destroy() {
    switch (type) {
      case RDTypeTag::StringTag:
        delete value.s;
        break;
      case RDTypeTag::AnyTag:
        delete value.a;
        break;
      case RDTypeTag::VecDoubleTag:
        delete value.vd;
        break;
      case RDTypeTag::VecFloatTag:
        delete value.vf;
        break;
      case RDTypeTag::VecIntTag:
        delete value.vi;
        break;
      case RDTypeTag::VecUnsignedIntTag:
        delete value.vu;
        break;
      case RDTypeTag::VecStringTag:
        delete value.vs;
        break;
      default:
        break;
    }
    type = RDTypeTag::EmptyTag;
  }

  static  // Given a type and an RDAnyValue - delete the appropriate structure
      inline void
      cleanup_rdvalue(RDValue &rdvalue) {
    rdvalue.destroy();
  }
};

/////////////////////////////////////////////////////////////////////////////////////
// Given two RDValue::Values - copy the appropriate structure
//   RDValue doesn't have a copy constructor, the default
//   copy act's like a move for better value semantics.
//  Containers may need to copy though.
inline void copy_rdvalue(RDValue &dest, const RDValue &src) {
  if (&dest == &src) {  // don't copy over yourself
    return;
  }
  dest.destroy();
  dest.type = src.type;
  switch (src.type) {
    case RDTypeTag::StringTag:
      dest.value.s = new std::string(*src.value.s);
      break;
    case RDTypeTag::AnyTag:
      dest.value.a = new std::any(*src.value.a);
      break;
    case RDTypeTag::VecDoubleTag:
      dest.value.vd = new std::vector<double>(*src.value.vd);
      break;
    case RDTypeTag::VecFloatTag:
      dest.value.vf = new std::vector<float>(*src.value.vf);
      break;
    case RDTypeTag::VecIntTag:
      dest.value.vi = new std::vector<int>(*src.value.vi);
      break;
    case RDTypeTag::VecUnsignedIntTag:
      dest.value.vu = new std::vector<unsigned int>(*src.value.vu);
      break;
    case RDTypeTag::VecStringTag:
      dest.value.vs = new std::vector<std::string>(*src.value.vs);
      break;
    default:
      dest = src;
  }
}

#ifdef RDK_32BIT_BUILD
// avoid register pressure and spilling on 32 bit systems
typedef const RDValue &RDValue_cast_t;
#else
typedef RDValue RDValue_cast_t;
#endif

/////////////////////////////////////////////////////////////////////////////////////
// rdvalue_is<T>

template <class T>
inline bool rdvalue_is(RDValue_cast_t v) {
  const short tag =
      RDTypeTag::GetTag<typename boost::remove_reference<T>::type>();

  // If we are an Any tag, check the any type info
  //  see the template specialization below if we are
  //  looking for a boost any directly
  if (v.getTag() == RDTypeTag::AnyTag) {
    return v.value.a->type() == typeid(T);
  }

  if (v.getTag() == tag) {
    return true;
  }

  return false;
}

template <>
inline bool rdvalue_is<std::any>(RDValue_cast_t v) {
  // If we are explicitly looking for a std::any
  //  then just check the top level tag
  const short tag = RDTypeTag::GetTag<std::any>();
  return v.getTag() == tag;
}
/////////////////////////////////////////////////////////////////////////////////////
// rdvalue_cast<T>
//
//  POD types do not support reference semantics.  Other types do.
// rdvalue_cast<const std::vector<double> &>(RDValue);  // ok
// rdvalue_cast<const float &>(RDValue); // bad_any_cast

// Get stuff stored in boost any
template <class T>
inline T rdvalue_cast(RDValue_cast_t v) {
  // Disable reference and pointer casts to POD data.
  BOOST_STATIC_ASSERT(!(
      (boost::is_pointer<T>::value &&
       (boost::is_integral<typename boost::remove_pointer<T>::type>::value ||
        boost::is_floating_point<
            typename boost::remove_pointer<T>::type>::value)) ||
      (boost::is_reference<T>::value &&
       (boost::is_integral<typename boost::remove_reference<T>::type>::value ||
        boost::is_floating_point<
            typename boost::remove_reference<T>::type>::value))));

  if (rdvalue_is<std::any>(v)) {
    return std::any_cast<T>(*v.ptrCast<std::any>());
  }
  throw std::bad_any_cast();
}

// POD casts
template <>
inline double rdvalue_cast<double>(RDValue_cast_t v) {
  if (rdvalue_is<double>(v)) {
    return v.value.d;
  }
  if (rdvalue_is<float>(v)) {
    return v.value.f;
  }
  throw std::bad_any_cast();
}

template <>
inline float rdvalue_cast<float>(RDValue_cast_t v) {
  if (rdvalue_is<float>(v)) {
    return v.value.f;
  }
  if (rdvalue_is<double>(v)) {
    return boost::numeric_cast<float>(v.value.d);
  }
  throw std::bad_any_cast();
}

template <>
inline int rdvalue_cast<int>(RDValue_cast_t v) {
  if (rdvalue_is<int>(v)) {
    return v.value.i;
  }
  if (rdvalue_is<unsigned int>(v)) {
    return boost::numeric_cast<int>(v.value.u);
  }
  throw std::bad_any_cast();
}

template <>
inline std::int8_t rdvalue_cast<std::int8_t>(RDValue_cast_t v) {
  if (rdvalue_is<int>(v)) {
    return boost::numeric_cast<std::int8_t>(v.value.i);
  }
  if (rdvalue_is<unsigned int>(v)) {
    return boost::numeric_cast<std::int8_t>(v.value.u);
  }
  throw std::bad_any_cast();
}

template <>
inline std::int16_t rdvalue_cast<std::int16_t>(RDValue_cast_t v) {
  if (rdvalue_is<int>(v)) {
    return boost::numeric_cast<std::int16_t>(v.value.i);
  }
  if (rdvalue_is<unsigned int>(v)) {
    return boost::numeric_cast<std::int16_t>(v.value.u);
  }
  throw std::bad_any_cast();
}

template <>
inline std::int64_t rdvalue_cast<std::int64_t>(RDValue_cast_t v) {
  if (rdvalue_is<int>(v)) {
    return static_cast<std::int64_t>(v.value.i);
  }
  if (rdvalue_is<unsigned int>(v)) {
    return static_cast<std::int64_t>(v.value.u);
  }
  if (rdvalue_is<std::any>(v)) {
    return std::any_cast<std::int64_t>(*v.ptrCast<std::any>());
  }
  throw std::bad_any_cast();
}

template <>
inline unsigned int rdvalue_cast<unsigned int>(RDValue_cast_t v) {
  if (rdvalue_is<unsigned int>(v)) {
    return v.value.u;
  }
  if (rdvalue_is<int>(v)) {
    return boost::numeric_cast<unsigned int>(v.value.i);
  }
  throw std::bad_any_cast();
}

template <>
inline std::uint8_t rdvalue_cast<std::uint8_t>(RDValue_cast_t v) {
  if (rdvalue_is<int>(v)) {
    return boost::numeric_cast<std::uint8_t>(v.value.i);
  }
  if (rdvalue_is<unsigned int>(v)) {
    return boost::numeric_cast<std::uint8_t>(v.value.u);
  }
  throw std::bad_any_cast();
}

template <>
inline std::uint16_t rdvalue_cast<std::uint16_t>(RDValue_cast_t v) {
  if (rdvalue_is<int>(v)) {
    return boost::numeric_cast<std::uint16_t>(v.value.i);
  }
  if (rdvalue_is<unsigned int>(v)) {
    return boost::numeric_cast<std::uint16_t>(v.value.u);
  }
  throw std::bad_any_cast();
}

template <>
inline std::uint64_t rdvalue_cast<std::uint64_t>(RDValue_cast_t v) {
  if (rdvalue_is<unsigned int>(v)) {
    return static_cast<std::uint64_t>(v.value.u);
  }
  if (rdvalue_is<int>(v)) {
    return boost::numeric_cast<std::uint64_t>(v.value.i);
  }
  if (rdvalue_is<std::any>(v)) {
    return std::any_cast<std::uint64_t>(*v.ptrCast<std::any>());
  }
  throw std::bad_any_cast();
}

template <>
inline bool rdvalue_cast<bool>(RDValue_cast_t v) {
  if (rdvalue_is<bool>(v)) {
    return v.value.b;
  }
  throw std::bad_any_cast();
}

}  // namespace RDKit
#endif
