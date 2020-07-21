//
//  Copyright (C) 2002-2008 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
#include <RDGeneral/export.h>
#ifndef _RD_STREAMOPS_H
#define _RD_STREAMOPS_H

#include "types.h"
#include "Invariant.h"
#include "RDProps.h"
#include <string>
#include <sstream>
#include <iostream>
#include <boost/cstdint.hpp>
#include <boost/predef.h>

namespace RDKit {
// this code block for handling endian problems is adapted from :
// http://stackoverflow.com/questions/105252/how-do-i-convert-between-big-endian-and-little-endian-values-in-c
enum EEndian {
  LITTLE_ENDIAN_ORDER,
  BIG_ENDIAN_ORDER,
#if defined(BOOST_ENDIAN_LITTLE_BYTE) || defined(BOOST_ENDIAN_LITTLE_WORD)
  HOST_ENDIAN_ORDER = LITTLE_ENDIAN_ORDER
#elif defined(BOOST_ENDIAN_BIG_BYTE)
  HOST_ENDIAN_ORDER = BIG_ENDIAN_ORDER
#elif defined(BOOST_ENDIAN_BIG_WORD)
#error "Cannot compile on word-swapped big-endian systems"
#else
#error "Failed to determine the system endian value"
#endif
};

// this function swap the bytes of values given it's size as a template
// parameter (could sizeof be used?).
template <class T, unsigned int size>
inline T SwapBytes(T value) {
  if (size < 2) return value;

  union {
    T value;
    char bytes[size];
  } in, out;

  in.value = value;

  for (unsigned int i = 0; i < size; ++i) {
    out.bytes[i] = in.bytes[size - 1 - i];
  }

  return out.value;
}

// Here is the function you will use. Again there is two compile-time assertion
// that use the boost libraries. You could probably comment them out, but if you
// do be cautious not to use this function for anything else than integers
// types. This function need to be called like this :
//
//     int x = someValue;
//     int i = EndianSwapBytes<HOST_ENDIAN_ORDER, BIG_ENDIAN_ORDER>(x);
//
template <EEndian from, EEndian to, class T>
inline T EndianSwapBytes(T value) {
  // A : La donnée à swapper à une taille de 2, 4 ou 8 octets
  BOOST_STATIC_ASSERT(sizeof(T) == 1 || sizeof(T) == 2 || sizeof(T) == 4 ||
                      sizeof(T) == 8);
  if (sizeof(T) == 1) return value;

  // A : La donnée à swapper est d'un type arithmetic
  // BOOST_STATIC_ASSERT(boost::is_arithmetic<T>::value);

  // Si from et to sont du même type on ne swap pas.
  if (from == to) return value;

  return SwapBytes<T, sizeof(T)>(value);
}
template <EEndian from, EEndian to>
inline char EndianSwapBytes(char value) {
  return value;
}
template <EEndian from, EEndian to>
inline unsigned char EndianSwapBytes(unsigned char value) {
  return value;
}
template <EEndian from, EEndian to>
inline signed char EndianSwapBytes(signed char value) {
  return value;
}
// --------------------------------------

//! Packs an integer and outputs it to a stream
inline void appendPackedIntToStream(std::stringstream &ss,
                                    boost::uint32_t num) {
  int nbytes, bix;
  unsigned int val, res;
  char tc;

  res = num;
  while (1) {
    if (res < (1 << 7)) {
      val = (res << 1);
      nbytes = 1;
      break;
    }
    res -= (1 << 7);
    if (res < (1 << 14)) {
      val = ((res << 2) | 1);
      nbytes = 2;
      break;
    }
    res -= (1 << 14);
    if (res < (1 << 21)) {
      val = ((res << 3) | 3);
      nbytes = 3;
      break;
    }
    res -= (1 << 21);
    if (res < (1 << 29)) {
      val = ((res << 3) | 7);
      nbytes = 4;
      break;
    } else {
      CHECK_INVARIANT(0, "ERROR: Integer too big to pack\n");
    }
  }
  // val = EndianSwapBytes<HOST_ENDIAN_ORDER,LITTLE_ENDIAN_ORDER>(val);

  for (bix = 0; bix < nbytes; bix++) {
    tc = (char)(val & 255);
    ss.write(&tc, 1);
    val >>= 8;
  }
}

//! Reads an integer from a stream in packed format and returns the result.
inline boost::uint32_t readPackedIntFromStream(std::stringstream &ss) {
  boost::uint32_t val, num;
  int shift, offset;
  char tmp;
  ss.read(&tmp, sizeof(tmp));
  if (ss.fail()) {
    throw std::runtime_error("failed to read from stream");
  }

  val = UCHAR(tmp);
  offset = 0;
  if ((val & 1) == 0) {
    shift = 1;
  } else if ((val & 3) == 1) {
    ss.read((char *)&tmp, sizeof(tmp));
    if (ss.fail()) {
      throw std::runtime_error("failed to read from stream");
    }

    val |= (UCHAR(tmp) << 8);
    shift = 2;
    offset = (1 << 7);
  } else if ((val & 7) == 3) {
    ss.read((char *)&tmp, sizeof(tmp));
    if (ss.fail()) {
      throw std::runtime_error("failed to read from stream");
    }

    val |= (UCHAR(tmp) << 8);
    ss.read((char *)&tmp, sizeof(tmp));
    if (ss.fail()) {
      throw std::runtime_error("failed to read from stream");
    }

    val |= (UCHAR(tmp) << 16);
    shift = 3;
    offset = (1 << 7) + (1 << 14);
  } else {
    ss.read((char *)&tmp, sizeof(tmp));
    if (ss.fail()) {
      throw std::runtime_error("failed to read from stream");
    }

    val |= (UCHAR(tmp) << 8);
    ss.read((char *)&tmp, sizeof(tmp));
    if (ss.fail()) {
      throw std::runtime_error("failed to read from stream");
    }

    val |= (UCHAR(tmp) << 16);
    ss.read((char *)&tmp, sizeof(tmp));
    if (ss.fail()) {
      throw std::runtime_error("failed to read from stream");
    }

    val |= (UCHAR(tmp) << 24);
    shift = 3;
    offset = (1 << 7) + (1 << 14) + (1 << 21);
  }
  num = (val >> shift) + offset;
  // num = EndianSwapBytes<LITTLE_ENDIAN_ORDER,HOST_ENDIAN_ORDER>(num);
  return num;
}

//! Reads an integer from a char * in packed format and returns the result.
//!  The argument is advanced
inline boost::uint32_t pullPackedIntFromString(const char *&text) {
  boost::uint32_t val, num;
  int shift, offset;
  char tmp;
  tmp = *text;
  text++;
  val = UCHAR(tmp);
  offset = 0;
  if ((val & 1) == 0) {
    shift = 1;
  } else if ((val & 3) == 1) {
    tmp = *text;
    text++;
    val |= (UCHAR(tmp) << 8);
    shift = 2;
    offset = (1 << 7);
  } else if ((val & 7) == 3) {
    tmp = *text;
    text++;
    val |= (UCHAR(tmp) << 8);
    tmp = *text;
    text++;
    val |= (UCHAR(tmp) << 16);
    shift = 3;
    offset = (1 << 7) + (1 << 14);
  } else {
    tmp = *text;
    text++;
    val |= (UCHAR(tmp) << 8);
    tmp = *text;
    text++;
    val |= (UCHAR(tmp) << 16);
    tmp = *text;
    text++;
    val |= (UCHAR(tmp) << 24);
    shift = 3;
    offset = (1 << 7) + (1 << 14) + (1 << 21);
  }
  num = (val >> shift) + offset;
  // num = EndianSwapBytes<LITTLE_ENDIAN_ORDER,HOST_ENDIAN_ORDER>(num);
  return num;
}

//! does a binary write of an object to a stream
template <typename T>
void streamWrite(std::ostream &ss, const T &val) {
  T tval = EndianSwapBytes<HOST_ENDIAN_ORDER, LITTLE_ENDIAN_ORDER>(val);
  ss.write((const char *)&tval, sizeof(T));
}

//! special case for string
inline void streamWrite(std::ostream &ss, const std::string &what) {
  unsigned int l = rdcast<unsigned int>(what.length());
  ss.write((const char *)&l, sizeof(l));
  ss.write(what.c_str(), sizeof(char) * l);
};

template <typename T>
void streamWriteVec(std::ostream &ss, const T &val) {
  streamWrite(ss, static_cast<boost::uint64_t>(val.size()));
  for (size_t i = 0; i < val.size(); ++i) streamWrite(ss, val[i]);
}

//! does a binary read of an object from a stream
template <typename T>
void streamRead(std::istream &ss, T &loc) {
  T tloc;
  ss.read((char *)&tloc, sizeof(T));
  if (ss.fail()) {
    throw std::runtime_error("failed to read from stream");
  }
  loc = EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(tloc);
}

//! special case for string
template <class T>
void streamRead(std::istream &ss, T &obj, int version) {
  RDUNUSED_PARAM(version);
  streamRead(ss, obj);
}

inline void streamRead(std::istream &ss, std::string &what, int version) {
  RDUNUSED_PARAM(version);
  unsigned int l;
  ss.read((char *)&l, sizeof(l));
  if (ss.fail()) {
    throw std::runtime_error("failed to read from stream");
  }
  char *buff = new char[l];
  ss.read(buff, sizeof(char) * l);
  if (ss.fail()) {
    throw std::runtime_error("failed to read from stream");
  }
  what = std::string(buff, l);
  delete[] buff;
};

template <class T>
void streamReadVec(std::istream &ss, T &val) {
  boost::uint64_t size;
  streamRead(ss, size);
  val.resize(size);

  for (size_t i = 0; i < size; ++i) streamRead(ss, val[i]);
}

inline void streamReadStringVec(std::istream &ss, std::vector<std::string> &val,
                                int version) {
  boost::uint64_t size;
  streamRead(ss, size);
  val.resize(size);

  for (size_t i = 0; i < size; ++i) streamRead(ss, val[i], version);
}

//! grabs the next line from an instream and returns it.
inline std::string getLine(std::istream *inStream) {
  std::string res;
  std::getline(*inStream, res);
  if ((res.length() > 0) && (res[res.length() - 1] == '\r')) {
    res.erase(res.length() - 1);
  }
  return res;
}
//! grabs the next line from an instream and returns it.
inline std::string getLine(std::istream &inStream) {
  return getLine(&inStream);
}

// n.b. We can't use RDTypeTag directly, they are implementation
//  specific
namespace DTags {
const unsigned char StringTag = 0;
const unsigned char IntTag = 1;
const unsigned char UnsignedIntTag = 2;
const unsigned char BoolTag = 3;
const unsigned char FloatTag = 4;
const unsigned char DoubleTag = 5;
const unsigned char VecStringTag = 6;
const unsigned char VecIntTag = 7;
const unsigned char VecUIntTag = 8;
const unsigned char VecBoolTag = 9;
const unsigned char VecFloatTag = 10;
const unsigned char VecDoubleTag = 11;

const unsigned char CustomTag = 0xFE;  // custom data
const unsigned char EndTag = 0xFF;
}  // namespace DTags

class CustomPropHandler {
 public:
  virtual ~CustomPropHandler(){};
  virtual const char *getPropName() const = 0;
  virtual bool canSerialize(const RDValue &value) const = 0;
  virtual bool read(std::istream &ss, RDValue &value) const = 0;
  virtual bool write(std::ostream &ss, const RDValue &value) const = 0;
  virtual CustomPropHandler *clone() const = 0;
};

typedef std::vector<std::shared_ptr<const CustomPropHandler>>
    CustomPropHandlerVec;

inline bool isSerializable(const Dict::Pair &pair,
                           const CustomPropHandlerVec &handlers = {}) {
  switch (pair.val.getTag()) {
    case RDTypeTag::StringTag:
    case RDTypeTag::IntTag:
    case RDTypeTag::UnsignedIntTag:
    case RDTypeTag::BoolTag:
    case RDTypeTag::FloatTag:
    case RDTypeTag::DoubleTag:

    case RDTypeTag::VecStringTag:
    case RDTypeTag::VecIntTag:
    case RDTypeTag::VecUnsignedIntTag:
    case RDTypeTag::VecFloatTag:
    case RDTypeTag::VecDoubleTag:
      return true;
    case RDTypeTag::AnyTag:
      for (auto &handler : handlers) {
        if (handler->canSerialize(pair.val)) {
          return true;
        }
      }
      return false;
    default:
      return false;
  }
}

inline bool streamWriteProp(std::ostream &ss, const Dict::Pair &pair,
                            const CustomPropHandlerVec &handlers = {}) {
  if (!isSerializable(pair, handlers)) {
    return false;
  }

  streamWrite(ss, pair.key);
  switch (pair.val.getTag()) {
    case RDTypeTag::StringTag:
      streamWrite(ss, DTags::StringTag);
      streamWrite(ss, rdvalue_cast<std::string>(pair.val));
      break;
    case RDTypeTag::IntTag:
      streamWrite(ss, DTags::IntTag);
      streamWrite(ss, rdvalue_cast<int>(pair.val));
      break;
    case RDTypeTag::UnsignedIntTag:
      streamWrite(ss, DTags::UnsignedIntTag);
      streamWrite(ss, rdvalue_cast<unsigned int>(pair.val));
      break;
    case RDTypeTag::BoolTag:
      streamWrite(ss, DTags::BoolTag);
      streamWrite(ss, rdvalue_cast<bool>(pair.val));
      break;
    case RDTypeTag::FloatTag:
      streamWrite(ss, DTags::FloatTag);
      streamWrite(ss, rdvalue_cast<float>(pair.val));
      break;
    case RDTypeTag::DoubleTag:
      streamWrite(ss, DTags::DoubleTag);
      streamWrite(ss, rdvalue_cast<double>(pair.val));
      break;

    case RDTypeTag::VecStringTag:
      streamWrite(ss, DTags::VecStringTag);
      streamWriteVec(ss, rdvalue_cast<std::vector<std::string>>(pair.val));
      break;
    case RDTypeTag::VecDoubleTag:
      streamWrite(ss, DTags::VecDoubleTag);
      streamWriteVec(ss, rdvalue_cast<std::vector<double>>(pair.val));
      break;
    case RDTypeTag::VecFloatTag:
      streamWrite(ss, DTags::VecFloatTag);
      streamWriteVec(ss, rdvalue_cast<std::vector<float>>(pair.val));
      break;
    case RDTypeTag::VecIntTag:
      streamWrite(ss, DTags::VecIntTag);
      streamWriteVec(ss, rdvalue_cast<std::vector<int>>(pair.val));
      break;
    case RDTypeTag::VecUnsignedIntTag:
      streamWrite(ss, DTags::VecUIntTag);
      streamWriteVec(ss, rdvalue_cast<std::vector<unsigned int>>(pair.val));
      break;
    default:
      for (auto &handler : handlers) {
        if (handler->canSerialize(pair.val)) {
          // The form of a custom tag is
          //  CustomTag
          //  customPropName (must be unique)
          //  custom serialization
          streamWrite(ss, DTags::CustomTag);
          streamWrite(ss, std::string(handler->getPropName()));
          handler->write(ss, pair.val);
          return true;
        }
      }

      return false;
  }
  return true;
}

inline bool streamWriteProps(std::ostream &ss, const RDProps &props,
                             bool savePrivate = false,
                             bool saveComputed = false,
                             const CustomPropHandlerVec &handlers = {}) {
  STR_VECT propsToSave = props.getPropList(savePrivate, saveComputed);
  std::set<std::string> propnames(propsToSave.begin(), propsToSave.end());

  const Dict &dict = props.getDict();
  unsigned int count = 0;
  for (Dict::DataType::const_iterator it = dict.getData().begin();
       it != dict.getData().end(); ++it) {
    if (propnames.find(it->key) != propnames.end()) {
      if (isSerializable(*it, handlers)) {
        count++;
      }
    }
  }

  streamWrite(ss, count);  // packed int?

  unsigned int writtenCount = 0;
  for (Dict::DataType::const_iterator it = dict.getData().begin();
       it != dict.getData().end(); ++it) {
    if (propnames.find(it->key) != propnames.end()) {
      if (isSerializable(*it, handlers)) {
        // note - not all properties are serializable, this may be
        //  a null op
        if (streamWriteProp(ss, *it, handlers)) {
          writtenCount++;
        }
      }
    }
  }
  POSTCONDITION(count == writtenCount,
                "Estimated property count not equal to written");
  return true;
}

template <class T>
void readRDValue(std::istream &ss, RDValue &value) {
  T v;
  streamRead(ss, v);
  value = v;
}

template <class T>
void readRDVecValue(std::istream &ss, RDValue &value) {
  std::vector<T> v;
  streamReadVec(ss, v);
  value = v;
}

inline void readRDValueString(std::istream &ss, RDValue &value) {
  std::string v;
  int version = 0;
  streamRead(ss, v, version);
  value = v;
}

inline void readRDStringVecValue(std::istream &ss, RDValue &value) {
  std::vector<std::string> v;
  int version = 0;
  streamReadStringVec(ss, v, version);
  value = v;
}

inline bool streamReadProp(std::istream &ss, Dict::Pair &pair,
                           bool &dictHasNonPOD,
                           const CustomPropHandlerVec &handlers = {}) {
  int version = 0;
  streamRead(ss, pair.key, version);

  unsigned char type;
  streamRead(ss, type);
  switch (type) {
    case DTags::IntTag:
      readRDValue<int>(ss, pair.val);
      break;
    case DTags::UnsignedIntTag:
      readRDValue<unsigned int>(ss, pair.val);
      break;
    case DTags::BoolTag:
      readRDValue<bool>(ss, pair.val);
      break;
    case DTags::FloatTag:
      readRDValue<float>(ss, pair.val);
      break;
    case DTags::DoubleTag:
      readRDValue<double>(ss, pair.val);
      break;

    case DTags::StringTag:
      readRDValueString(ss, pair.val);
      dictHasNonPOD = true;
      break;
    case DTags::VecStringTag:
      readRDStringVecValue(ss, pair.val);
      dictHasNonPOD = true;
      break;
    case DTags::VecIntTag:
      readRDVecValue<int>(ss, pair.val);
      dictHasNonPOD = true;
      break;
    case DTags::VecUIntTag:
      readRDVecValue<unsigned int>(ss, pair.val);
      dictHasNonPOD = true;
      break;
    case DTags::VecFloatTag:
      readRDVecValue<float>(ss, pair.val);
      dictHasNonPOD = true;
      break;
    case DTags::VecDoubleTag:
      readRDVecValue<double>(ss, pair.val);
      dictHasNonPOD = true;
      break;
    case DTags::CustomTag: {
      std::string propType;
      int version = 0;
      streamRead(ss, propType, version);
      for (auto &handler : handlers) {
        if (propType == handler->getPropName()) {
          handler->read(ss, pair.val);
          dictHasNonPOD = true;
          return true;
        }
      }
      return false;
    }

    default:
      return false;
  }
  return true;
}

inline unsigned int streamReadProps(std::istream &ss, RDProps &props,
                                    const CustomPropHandlerVec &handlers = {}) {
  unsigned int count;
  streamRead(ss, count);

  Dict &dict = props.getDict();
  dict.reset();  // Clear data before repopulating
  dict.getData().resize(count);
  for (unsigned index = 0; index < count; ++index) {
    CHECK_INVARIANT(streamReadProp(ss, dict.getData()[index],
                                   dict.getNonPODStatus(), handlers),
                    "Corrupted property serialization detected");
  }

  return count;
}

}  // namespace RDKit

#endif
