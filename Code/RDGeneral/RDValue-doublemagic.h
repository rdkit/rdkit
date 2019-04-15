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
#ifndef RDKIT_RDVALUE_PTRMAGIC_H
#define RDKIT_RDVALUE_PTRMAGIC_H

#include <cstdint>
#include <cassert>
#include <boost/any.hpp>
#include "Invariant.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/utility.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/type_traits.hpp>
#include <boost/static_assert.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <cmath>
#include "LocaleSwitcher.h"

#define RDVALUE_HASBOOL

namespace RDKit {

// Inspired by
// https://nikic.github.io/2012/02/02/Pointer-magic-for-efficient-dynamic-value-representations.html
// 16 bit storage for value types using Quiet NaN spaces in
//  doubles
// Won't work on Solaris and some other os's as mmaping maps from
// top memory down
// Example check:
//     std::string *pointer = new std::string(v);
//      assert((reinterpret_cast<boost::uint64_t>(pointer) & StringTag) == 0);

//  implementations, need a typedef at compile time to figure this out.
//  current implementation is probably little endian, need to check.

/*
  Encoding for storing other things as a double.  Use
  Quiet NaN
  Quiet NaN: // used to encode types
   F   F    F   1XXX < - X = type bits (first bit is set to one)

  seeeeeee|eeeemmmm|mmmmmmmm|mmmmmmmm|mmmmmmmm|mmmmmmmm|mmmmmmmm|mmmmmmmm
  s1111111|11111ppp|pppppppp|pppppppp|pppppppp|pppppppp|pppppppp|pppppppp
               ^- first mantissa bit 1    everything else is "payload" -^
   ^- exponent bits all 1                 and mustn't be all-zero (as it
  ^- any sign bit                         would be INF then)

  Available
  8 = 1000 MaxDouble // Not really a tag, is a sentinel
  9 = 1001 Float
  b = 1010 Int32
  a = 1011 Uint32
  C = 1100 <none>
  D = 1101 <none>
  E = 1110 <none>
  F = 1111 PtrTag (look at lower 3 bits for type)
*/

namespace RDTypeTag {
static const boost::uint64_t NaN = 0xfff7FFFFFFFFFFFF;        // signalling NaN
static const boost::uint64_t MaxDouble = 0xfff8000000000000;  //
static const boost::uint64_t DoubleTag = 0xfff8000000000000;  //
static const boost::uint64_t FloatTag = 0xfff9000000000000;   //
static const boost::uint64_t IntTag = 0xfffa000000000000;     //
static const boost::uint64_t UnsignedIntTag = 0xfffb000000000000;  //
static const boost::uint64_t BoolTag = 0xfffc000000000000;         //

// PTR Tags use the last 3 bits for typing info
static const boost::uint64_t PtrTag = 0xffff000000000000;
static const boost::uint64_t StringTag = 0xffff000000000001;          // 001
static const boost::uint64_t VecDoubleTag = 0xffff000000000002;       // 010
static const boost::uint64_t VecFloatTag = 0xffff000000000003;        // 011
static const boost::uint64_t VecIntTag = 0xffff000000000004;          // 100
static const boost::uint64_t VecUnsignedIntTag = 0xffff000000000005;  // 101
static const boost::uint64_t VecStringTag = 0xffff000000000006;       // 110
static const boost::uint64_t AnyTag = 0xffff000000000007;             // 111

// Retrieves the tag (and PtrMask) from the type
template <class T>
inline boost::uint64_t GetTag() {
  return AnyTag;
}
template <>
inline boost::uint64_t GetTag<double>() {
  return MaxDouble;
}
template <>
inline boost::uint64_t GetTag<float>() {
  return FloatTag;
}
template <>
inline boost::uint64_t GetTag<int>() {
  return IntTag;
}
template <>
inline boost::uint64_t GetTag<unsigned int>() {
  return UnsignedIntTag;
}
template <>
inline boost::uint64_t GetTag<bool>() {
  return BoolTag;
}
template <>
inline boost::uint64_t GetTag<std::string>() {
  return StringTag;
}
template <>
inline boost::uint64_t GetTag<std::vector<double>>() {
  return VecDoubleTag;
}
template <>
inline boost::uint64_t GetTag<std::vector<float>>() {
  return VecFloatTag;
}
template <>
inline boost::uint64_t GetTag<std::vector<int>>() {
  return VecIntTag;
}
template <>
inline boost::uint64_t GetTag<std::vector<unsigned int>>() {
  return VecUnsignedIntTag;
}
template <>
inline boost::uint64_t GetTag<std::vector<std::string>>() {
  return VecStringTag;
}
template <>
inline boost::uint64_t GetTag<boost::any>() {
  return AnyTag;
}
}  // namespace RDTypeTag

struct RDValue {
  // Bit Twidling for conversion from the Tag to a Pointer
  static const boost::uint64_t TagMask = 0xFFFF000000000000;
  static const boost::uint64_t PointerTagMask = 0xFFFF000000000007;
  static const boost::uint64_t ApplyMask = 0x0000FFFFFFFFFFFF;
  static const boost::uint64_t ApplyPtrMask = 0x0000FFFFFFFFFFF8;

  union {
    double doubleBits;
    boost::uint64_t otherBits;
  };

  inline RDValue() : doubleBits(0.0) {}

  inline RDValue(double number) {
    if (boost::math::isnan(number)) {
      // Store a signalling NaN for NaN's.
      //  quiet NaNs are used for other types.
      otherBits = RDTypeTag::NaN;
      assert(boost::math::isnan(doubleBits));
    } else
      doubleBits = number;
  }

  inline RDValue(float number) {
    otherBits = 0 | RDTypeTag::FloatTag;
    memcpy(((char *)&otherBits), &number, sizeof(float));
  }

  inline RDValue(int32_t number) {
    otherBits = (((boost::uint64_t)number) & ApplyMask) | RDTypeTag::IntTag;
  }

  inline RDValue(unsigned int number) {
    otherBits =
        (((boost::uint64_t)number) & ApplyMask) | RDTypeTag::UnsignedIntTag;
  }

  inline RDValue(bool number) {
    otherBits =
        (static_cast<boost::uint64_t>(number) & ApplyMask) | RDTypeTag::BoolTag;
  }

  inline RDValue(boost::any *pointer) {
    // ensure that the pointer really is only 48 bit
    assert((reinterpret_cast<boost::uint64_t>(pointer) & RDTypeTag::AnyTag) ==
           0);
    otherBits = reinterpret_cast<boost::uint64_t>(pointer) | RDTypeTag::AnyTag;
  }

  inline RDValue(const boost::any &any) {
    // ensure that the pointer really is only 48 bit
    boost::any *pointer = new boost::any(any);
    assert((reinterpret_cast<boost::uint64_t>(pointer) & RDTypeTag::AnyTag) ==
           0);
    otherBits = reinterpret_cast<boost::uint64_t>(pointer) | RDTypeTag::AnyTag;
  }

  // Unknown types are stored as boost::any
  template <class T>
  inline RDValue(const T &v) {
    boost::any *pointer = new boost::any(v);
    assert((reinterpret_cast<boost::uint64_t>(pointer) & RDTypeTag::AnyTag) ==
           0);
    otherBits = reinterpret_cast<boost::uint64_t>(pointer) | RDTypeTag::AnyTag;
  }

  inline RDValue(const std::string &v) {
    std::string *pointer = new std::string(v);
    assert((reinterpret_cast<boost::uint64_t>(pointer) &
            RDTypeTag::StringTag) == 0);
    otherBits =
        reinterpret_cast<boost::uint64_t>(pointer) | RDTypeTag::StringTag;
  }

  inline RDValue(const std::vector<double> &v) {
    std::vector<double> *pointer = new std::vector<double>(v);
    assert((reinterpret_cast<boost::uint64_t>(pointer) &
            RDTypeTag::VecDoubleTag) == 0);
    otherBits =
        reinterpret_cast<boost::uint64_t>(pointer) | RDTypeTag::VecDoubleTag;
  }

  inline RDValue(const std::vector<float> &v) {
    std::vector<float> *pointer = new std::vector<float>(v);
    assert((reinterpret_cast<boost::uint64_t>(pointer) &
            RDTypeTag::VecFloatTag) == 0);
    otherBits =
        reinterpret_cast<boost::uint64_t>(pointer) | RDTypeTag::VecFloatTag;
  }

  inline RDValue(const std::vector<int> &v) {
    std::vector<int> *pointer = new std::vector<int>(v);
    assert((reinterpret_cast<boost::uint64_t>(pointer) &
            RDTypeTag::VecIntTag) == 0);
    otherBits =
        reinterpret_cast<boost::uint64_t>(pointer) | RDTypeTag::VecIntTag;
  }

  inline RDValue(const std::vector<unsigned int> &v) {
    std::vector<unsigned int> *pointer = new std::vector<unsigned int>(v);
    assert((reinterpret_cast<boost::uint64_t>(pointer) &
            RDTypeTag::VecIntTag) == 0);
    otherBits = reinterpret_cast<boost::uint64_t>(pointer) |
                RDTypeTag::VecUnsignedIntTag;
  }

  inline RDValue(const std::vector<std::string> &v) {
    std::vector<std::string> *pointer = new std::vector<std::string>(v);
    assert((reinterpret_cast<boost::uint64_t>(pointer) &
            RDTypeTag::VecStringTag) == 0);
    otherBits =
        reinterpret_cast<boost::uint64_t>(pointer) | RDTypeTag::VecStringTag;
  }

  boost::uint64_t getTag() const {
    if (otherBits < RDTypeTag::MaxDouble ||
        (otherBits & RDTypeTag::NaN) == RDTypeTag::NaN) {
      return RDTypeTag::DoubleTag;
    }

    boost::uint64_t tag = otherBits & TagMask;
    if (tag == RDTypeTag::PtrTag) return otherBits & PointerTagMask;
    return tag;
  }

  // ptrCast - unsafe, use rdvalue_cast instead.
  template <class T>
  inline T *ptrCast() const {
    return reinterpret_cast<T *>(otherBits & ~RDTypeTag::GetTag<T>());
  }

  // RDValue doesn't have an explicit destructor, it must
  //  be wrapped in a container.
  // The idea is that POD types don't need to be destroyed
  //  and this allows the container optimization possibilities.
  inline void destroy() {
    switch (getTag()) {
      case RDTypeTag::StringTag:
        delete ptrCast<std::string>();
        break;
      case RDTypeTag::VecDoubleTag:
        delete ptrCast<std::vector<double>>();
        break;
      case RDTypeTag::VecFloatTag:
        delete ptrCast<std::vector<float>>();
        break;
      case RDTypeTag::VecIntTag:
        delete ptrCast<std::vector<int>>();
        break;
      case RDTypeTag::VecUnsignedIntTag:
        delete ptrCast<std::vector<unsigned int>>();
        break;
      case RDTypeTag::VecStringTag:
        delete ptrCast<std::vector<std::string>>();
        break;
      case RDTypeTag::AnyTag:
        delete ptrCast<boost::any>();
        break;
      default:
        break;
    }
  }

  static inline void cleanup_rdvalue(RDValue v) { v.destroy(); }
};

/////////////////////////////////////////////////////////////////////////////////////
// Given two RDValue::Values - copy the appropriate structure
//   RDValue doesn't have a copy constructor, the default
//   copy act's like a move for better value semantics.
//  Containers may need to copy though.
inline void copy_rdvalue(RDValue &dest, const RDValue &src) {
  dest.destroy();
  switch (src.getTag()) {
    case RDTypeTag::StringTag:
      dest = RDValue(*src.ptrCast<std::string>());
      break;
    case RDTypeTag::VecDoubleTag:
      dest = RDValue(*src.ptrCast<std::vector<double>>());
      break;
    case RDTypeTag::VecFloatTag:
      dest = RDValue(*src.ptrCast<std::vector<float>>());
      break;
    case RDTypeTag::VecIntTag:
      dest = RDValue(*src.ptrCast<std::vector<int>>());
      break;
    case RDTypeTag::VecUnsignedIntTag:
      dest = RDValue(*src.ptrCast<std::vector<unsigned int>>());
      break;
    case RDTypeTag::VecStringTag:
      dest = RDValue(*src.ptrCast<std::vector<std::string>>());
      break;
    case RDTypeTag::AnyTag:
      dest = RDValue(*src.ptrCast<boost::any>());
      break;
    default:
      dest = src;
  }
}

/////////////////////////////////////////////////////////////////////////////////////
// rdvalue_is<T>

template <class T>
inline bool rdvalue_is(const RDValue_cast_t) {
  const short tag =
      RDTypeTag::GetTag<typename boost::remove_reference<T>::type>();
  if (v.getTag() == tag) return true;

  // If we are an Any tag, check the any type info
  if (v.getTag() == RDTypeTag::AnyTag) {
    return v.value.a->type() == typeid(T);
  }

  return false;
}

template <>
inline bool rdvalue_is<double>(const RDValue_cast_t) {
  return v.otherBits < RDTypeTag::MaxDouble ||
         (v.otherBits & RDTypeTag::NaN) == RDTypeTag::NaN;
}

template <>
inline bool rdvalue_is<const double &>(const RDValue_cast_t) {
  return rdvalue_is<double>(v);
}

/*
template<>
inline bool rdvalue_is<bool>(const RDValue_cast_t) {
  return (v.getTag() == RDTypeTag::IntTag &&
          (static_cast<int32_t>(v.otherBits & ~RDTypeTag::IntTag) == 1 ||
           static_cast<int32_t>(v.otherBits & ~RDTypeTag::IntTag) == 0 ));
}

template<>
inline bool rdvalue_is<const bool&>(const RDValue_cast_t) {
  return rdvalue_is<bool>(v);
}
*/

/////////////////////////////////////////////////////////////////////////////////////
// rdvalue_cast<T>
//
//  POD types do not support reference semantics.  Other types do.
// rdvalue_cast<const std::vector<double> &>(RDValue);  // ok
// rdvalue_cast<const float &>(RDValue); // bad_any_cast

typedef RDValue RDValue_cast_t;
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

  if (rdvalue_is<boost::any>(v)) {
    return boost::any_cast<T>(*v.ptrCast<boost::any>());
  }
  throw boost::bad_any_cast();
}

// POD casts
template <>
inline double rdvalue_cast<double>(RDValue_cast_t v) {
  if (rdvalue_is<double>(v)) return v.doubleBits;
  throw boost::bad_any_cast();
}

template <>
inline float rdvalue_cast<float>(RDValue_cast_t v) {
  if (rdvalue_is<float>(v)) {
    float f;
    memcpy(&f, ((char *)&v.otherBits), sizeof(float));
    return f;
  }
  throw boost::bad_any_cast();
}

// n.b. with const expressions, could use ~RDTagTypes::GetTag<T>()
//  and enable_if
template <>
inline int rdvalue_cast<int>(RDValue_cast_t v) {
  if (rdvalue_is<int>(v))
    return static_cast<int32_t>(v.otherBits & ~RDTypeTag::IntTag);
  throw boost::bad_any_cast();
}
template <>
inline unsigned int rdvalue_cast<unsigned int>(RDValue_cast_t v) {
  if (rdvalue_is<unsigned int>(v))
    return static_cast<uint32_t>(v.otherBits & ~RDTypeTag::UnsignedIntTag);
  throw boost::bad_any_cast();
}

template <>
inline bool rdvalue_cast<bool>(RDValue_cast_t v) {
  if (rdvalue_is<bool>(v))
    return static_cast<bool>(v.otherBits & ~RDTypeTag::BoolTag);
  throw boost::bad_any_cast();
}

}  // namespace RDKit
#endif
