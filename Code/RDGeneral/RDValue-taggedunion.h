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
#ifndef RDKIT_RDVALUE_TAGGED_UNION_H
#define RDKIT_RDVALUE_TAGGED_UNION_H

#include <stdint.h>
#include <cassert>
#include <boost/any.hpp>
#include "Invariant.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <boost/utility.hpp>
#include <boost/lexical_cast.hpp>
#include "LocaleSwitcher.h"
#include "tags.h"

#define RDVALUE_HAS_KEY

namespace RDKit {

// RDValue does not dynamically create POD types (kind of like
//  cdiggins::any)  However, it doesn't use RTTI type info
//  directly, it uses a companion short valued type
//  to determine what to do.
// For unregistered types, it falls back to boost::any.
//  The Size of an RDAny is (sizeof(double) + sizeof(short) == 10 bytes
//   (aligned to actually 16 so hard to pass as value type)
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

// RDValue does not manange memory of non-pod data
//  this must be done externally (string, Any, vector...)
// Tagged union

namespace RDTypeTag {
  const unsigned char EmptyTag           = 0;
  const unsigned char IntTag             = 1;
  const unsigned char DoubleTag          = 2;
  const unsigned char StringTag          = 3;
  const unsigned char FloatTag           = 4;
  const unsigned char BoolTag            = 5;
  const unsigned char UnsignedIntTag     = 6;
  const unsigned char AnyTag             = 7;
  const unsigned char VecDoubleTag       = 8;
  const unsigned char VecFloatTag        = 9;
  const unsigned char VecIntTag          = 10;
  const unsigned char VecUnsignedIntTag  = 11;
  const unsigned char VecStringTag       = 12;
  template<class T> inline unsigned char GetTag() { return AnyTag; }
  template<> inline unsigned char GetTag<double>() { return DoubleTag; }
  template<> inline unsigned char GetTag<float>() { return FloatTag; }
  template<> inline unsigned char GetTag<int>() { return IntTag; }
  template<> inline unsigned char GetTag<unsigned int>() { return UnsignedIntTag; }
  template<> inline unsigned char GetTag<bool>() { return BoolTag; }
  template<> inline unsigned char GetTag<std::string>() { return StringTag; }
  template<> inline unsigned char GetTag<std::vector<double> >() { return VecDoubleTag; }
  template<> inline unsigned char GetTag<std::vector<float> >() { return VecFloatTag; }
  template<> inline unsigned char GetTag<std::vector<int> >() { return VecIntTag; }
  template<> inline unsigned char GetTag<std::vector<unsigned int> >() { return VecUnsignedIntTag; }
  template<> inline unsigned char GetTag<std::vector<std::string> >() { return VecStringTag; }
  template<> inline unsigned char GetTag<boost::any>() { return AnyTag; }

  namespace detail {
    union Value {
      double d;
      float f;
      int i;
      unsigned u;
      bool b;
      std::string *s;
      boost::any *a;
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
      inline Value(bool  v) : b(v) {}
      inline Value(std::string *v) : s(v) {}
      inline Value(boost::any *v) : a(v) {}
      inline Value(std::vector<double> *v) : vd(v) {}
      inline Value(std::vector<float> *v) : vf(v) {}
      inline Value(std::vector<int> *v) : vi(v) {}
      inline Value(std::vector<unsigned int> *v) : vu(v) {}
      inline Value(std::vector<std::string> *v) : vs(v) {}
    };
  
    template<class T>
        inline T* valuePtrCast(Value value) {
      return boost::any_cast<T*>(*value.a);
    }
    template<>
    inline boost::any *valuePtrCast<boost::any>(Value value) { return value.a; }

    template<>
    inline std::string *valuePtrCast<std::string>(Value value) { return value.s; }
    template<>
    inline std::vector<double> *valuePtrCast<std::vector<double> >(Value value) {
      return value.vd; }
    template<>
    inline std::vector<float> *valuePtrCast<std::vector<float> >(Value value) {
      return value.vf; }
    template<>
    inline std::vector<int> *valuePtrCast<std::vector<int> >(Value value) {
      return value.vi; }
    template<>
    inline std::vector<unsigned int> *valuePtrCast<std::vector<unsigned int> >(
        Value value) {
      return value.vu; }
    template<>
    inline std::vector<std::string> *valuePtrCast<std::vector<std::string> >(
        Value value) {
      return value.vs; }
  }
}

const unsigned int KEYMAX = 0x00FFFFFF;

struct RDValue {
  RDTypeTag::detail::Value value;
  short type;
  int key;

  inline RDValue(): value(0.0), type(RDTypeTag::EmptyTag), key(KEYMAX) {}
  // Pod Style (Direct storage)
 inline RDValue(double v)   : value(v), type(RDTypeTag::DoubleTag), key(KEYMAX) {}
 inline RDValue(float v)    : value(v), type(RDTypeTag::FloatTag), key(KEYMAX) {}
 inline RDValue(int v)      : value(v), type(RDTypeTag::IntTag), key(KEYMAX) {}
 inline RDValue(unsigned v) : value(v), type(RDTypeTag::UnsignedIntTag), key(KEYMAX) {}
 inline RDValue(bool v)     : value(v), type(RDTypeTag::BoolTag), key(KEYMAX) {}

 inline RDValue(boost::any *v)  : value(v),type(RDTypeTag::AnyTag), key(KEYMAX) {}

  // Copies passed in pointers
 inline RDValue(const boost::any &v)  : value(new boost::any(v)),type(RDTypeTag::AnyTag), key(KEYMAX) {}
 inline RDValue(const std::string &v) : value(new std::string(v)),type(RDTypeTag::StringTag), key(KEYMAX){};
 template <class T>
 inline RDValue(const T &v) : value(new boost::any(v)),type(RDTypeTag::AnyTag), key(KEYMAX) {}

 inline RDValue(const std::vector<double> &v) : value(new std::vector<double>(v)),
    type(RDTypeTag::VecDoubleTag), key(KEYMAX) {}
 inline RDValue(const std::vector<float> &v)  : value(new std::vector<float>(v)),
    type(RDTypeTag::VecFloatTag), key(KEYMAX) {}
 inline RDValue(const std::vector<int> &v)    : value(new std::vector<int>(v)),
    type(RDTypeTag::VecIntTag), key(KEYMAX) {}
 inline RDValue(const std::vector<unsigned int> &v) :
  value(new std::vector<unsigned int>(v)),
    type(RDTypeTag::VecUnsignedIntTag), key(KEYMAX) {}
 inline RDValue(const std::vector<std::string> &v) :
  value(new std::vector<std::string>(v)),
    type(RDTypeTag::VecStringTag), key(KEYMAX) {}

  inline unsigned char getTag() const { return type; }
  
  // we have enough space left over for a 24 bit key, might as well use it
  inline int           getKey() const { return static_cast<int>(key); }
  inline void          setKey(int k) {
    CHECK_INVARIANT(k>=0 && k < static_cast<int>((1u<<24)-1),
                    "Key index negative or too large to save");
    key=k;
  }

  // ptrCast - unsafe, use rdvalue_cast instead.
  template<class T>
  inline T* ptrCast() const {
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

  static // Given a type and an RDAnyValue - delete the appropriate structure
  inline void cleanup_rdvalue(RDValue &rdvalue) {
    rdvalue.destroy();
  }
};

/////////////////////////////////////////////////////////////////////////////////////
// Given two RDValue::Values - copy the appropriate structure
//   RDValue doesn't have a copy constructor, the default
//   copy act's like a move for better value semantics.
//  Containers may need to copy though.
inline void copy_rdvalue(RDValue &dest,
                         const RDValue &src) {
  dest.destroy();
  dest.type = src.type;
  dest.key = src.key;
  switch (src.type) {
    case RDTypeTag::StringTag:
      dest.value.s = new std::string(*src.value.s);
      break;
    case RDTypeTag::AnyTag:
      dest.value.a = new boost::any(*src.value.a);
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

/////////////////////////////////////////////////////////////////////////////////////
// rdvalue_is<T>

template<class T>
inline bool rdvalue_is(RDValue v) {
  return v.getTag() == RDTypeTag::GetTag<typename boost::remove_reference<T>::type>();
}

/////////////////////////////////////////////////////////////////////////////////////
// rdvalue_cast<T>
//
//  POD types do not support reference semantics.  Other types do.
// rdvalue_cast<const std::vector<double> &>(RDValue);  // ok
// rdvalue_cast<const float &>(RDValue); // bad_any_cast

#ifdef RDK_32BIT_BUILD
 // avoid register pressure and spilling on 32 bit systems
 typedef const RDValue& RDValue_cast_t;
#else
 typedef RDValue RDValue_cast_t;
#endif

// Get stuff stored in boost any
template<class T>
inline T rdvalue_cast(RDValue_cast_t v) {
  // Disable reference and pointer casts to POD data.
  BOOST_STATIC_ASSERT( !(
      (boost::is_pointer<T>::value && (
          boost::is_integral<typename boost::remove_pointer<T>::type>::value ||
          boost::is_floating_point<typename boost::remove_pointer<T>::type>::value)) ||
      (boost::is_reference<T>::value && (
          boost::is_integral<typename boost::remove_reference<T>::type>::value ||
          boost::is_floating_point<typename boost::remove_reference<T>::type>::value))
                         ));

  if (rdvalue_is<boost::any>(v)) {
    return boost::any_cast<T>(*v.ptrCast<boost::any>());
  }
  throw boost::bad_any_cast();
}

// POD casts
template<>
inline double rdvalue_cast<double>(RDValue_cast_t v) {
  if (rdvalue_is<double>(v)) return v.value.d;
  throw boost::bad_any_cast();
}

template<>
inline float rdvalue_cast<float>(RDValue_cast_t v) {
  if (rdvalue_is<float>(v)) return v.value.f;
  throw boost::bad_any_cast();
}

template<>
inline int rdvalue_cast<int>(RDValue_cast_t v) {
  if (rdvalue_is<int>(v)) return v.value.i;
  throw boost::bad_any_cast();
}
template<>
inline unsigned int rdvalue_cast<unsigned int>(RDValue_cast_t v) {
  if (rdvalue_is<unsigned int>(v)) return v.value.u;
  throw boost::bad_any_cast();
}

template<>
inline bool rdvalue_cast<bool>(RDValue_cast_t v) {
  if (rdvalue_is<bool>(v)) return v.value.b;
  throw boost::bad_any_cast();
}

struct KeyIntPair {
   RDValue val;
  
  KeyIntPair() : val() {}
  KeyIntPair(int k, const RDValue &v) : val(v) {
   val.setKey(k);
   TEST_ASSERT(val.getKey() == k);
  }
  KeyIntPair(const std::string &s, RDValue_cast_t v) :
   val(v) { val.setKey(tagmap.get(s));}
  
  inline void setKey(int k) { val.setKey(k); }
  inline int  getKey() const { return val.getKey(); }
  static RDTags tagmap;
};  

} // namespace rdkit
#endif

