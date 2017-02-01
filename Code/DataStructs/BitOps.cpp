// $Id$
//
//  Copyright (C) 2003-2012 greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "BitVects.h"
#include "BitOps.h"
#include <math.h>
#include <string>
#include <iostream>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/types.h>
#include <RDGeneral/Exceptions.h>
#include <sstream>
#include <cstdlib>

#include <boost/lexical_cast.hpp>

using namespace RDKit;

int getBitId(const char*& text, int format, int size, int curr) {
  PRECONDITION(text, "no text");
  int res = -1;
  if ((format == 0) ||
      ((format == 1) && (size >= std::numeric_limits<unsigned short>::max()))) {
    int tmp;
    tmp = EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(*(int*)text);
    text += sizeof(tmp);
    res = tmp;
  } else if (format == 1) {  // version 16 and on bits sotred as short ints
    unsigned short tmp;
    tmp = EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(
        *(unsigned short*)text);
    text += sizeof(tmp);
    res = tmp;
  } else if (format == 2) {  // run length encoded format
    res = curr + RDKit::pullPackedIntFromString(text);
  }
  return res;
}

bool AllProbeBitsMatch(const std::string& probe, const std::string& ref) {
  return AllProbeBitsMatch(probe.c_str(), ref.c_str());
}

bool AllProbeBitsMatch(const char* probe, const char* ref) {
  PRECONDITION(probe, "no probe text");
  PRECONDITION(ref, "no probe text");
  int probeFormat = 0;
  int refFormat = 0;
  int version = 0;

  int probeSize =
      EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(*(int*)probe);
  probe += sizeof(probeSize);
  if (probeSize < 0) {
    version = -1 * probeSize;
    if (version == 16) {
      probeFormat = 1;
    } else if (version == 32) {
      probeFormat = 2;
    } else {
      throw("Unknown version type for the encode bit vect");
    }
    probeSize =
        EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(*(int*)probe);
    probe += sizeof(probeSize);
  }

  int refSize =
      EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(*(int*)ref);
  ref += sizeof(refSize);
  if (refSize < 0) {
    version = -1 * refSize;
    if (version == 16) {
      refFormat = 1;
    } else if (version == 32) {
      refFormat = 2;
    } else {
      throw("Unknown version type for the encode bit vect");
    }
    refSize =
        EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(*(int*)ref);
    ref += sizeof(refSize);
  }

  int nProbeOn =
      EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(*(int*)probe);
  probe += sizeof(nProbeOn);
  int nRefOn =
      EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(*(int*)ref);
  ref += sizeof(nRefOn);

  int currProbeBit = 0;
  currProbeBit = getBitId(probe, probeFormat, probeSize, currProbeBit);
  nProbeOn--;

  int currRefBit = 0;
  currRefBit = getBitId(ref, refFormat, refSize, currRefBit);
  nRefOn--;

  while (nProbeOn) {
    while (currRefBit < currProbeBit && nRefOn > 0) {
      if (refFormat == 2) currRefBit++;
      currRefBit = getBitId(ref, refFormat, refSize, currRefBit);
      nRefOn--;
    }
    if (currRefBit != currProbeBit) return false;
    if (probeFormat == 2) currProbeBit++;
    currProbeBit = getBitId(probe, probeFormat, probeSize, currProbeBit);
    nProbeOn--;
  }
  return true;
}

template <typename T1>
bool AllProbeBitsMatch(const T1& probe, const std::string& pkl) {
  const char* text = pkl.c_str();
  int format = 0;
  int nOn = 0, size, version = 0;
  size = EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(*(int*)text);
  text += sizeof(size);
  if (size < 0) {
    version = -1 * size;
    if (version == 16) {
      format = 1;
    } else if (version == 32) {
      format = 2;
    } else {
      throw("Unknown version type for the encode bit vect");
    }
    size = EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(*(int*)text);
    text += sizeof(size);
  }
  nOn = EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(*(int*)text);
  text += sizeof(nOn);

  int currBit = 0;
  currBit = getBitId(text, format, size, currBit);
  nOn--;
  std::vector<int> obl;
  probe.getOnBits(obl);

  // for(int i=0;i<probe.getNumBits();i++){
  //  if(probe.getBit(i)){
  for (std::vector<int>::const_iterator i = obl.begin(); i != obl.end(); i++) {
    while (currBit < *i && nOn > 0) {
      if (format == 2) currBit++;
      currBit = getBitId(text, format, size, currBit);
      nOn--;
    }
    if (currBit != *i) return false;
    //}
  }
  return true;
}
template bool AllProbeBitsMatch(const SparseBitVect& bv1,
                                const std::string& pkl);
template bool AllProbeBitsMatch(const ExplicitBitVect& bv1,
                                const std::string& pkl);
template <typename T1>
bool AllProbeBitsMatch(const T1& probe, const T1& ref) {
  for (unsigned int i = 0; i < probe.getNumBits(); ++i) {
    if (probe.getBit(i) && !ref.getBit(i)) {
      return false;
    }
  }
  return true;
}
template bool AllProbeBitsMatch(const SparseBitVect& bv1,
                                const SparseBitVect& bv2);
// template bool AllProbeBitsMatch(const ExplicitBitVect& bv1,const
// ExplicitBitVect &bv2);

bool AllProbeBitsMatch(const ExplicitBitVect& probe,
                       const ExplicitBitVect& ref) {
  return probe.dp_bits->is_subset_of(*(ref.dp_bits));
}

// """ -------------------------------------------------------
//
//  NumOnBitsInCommon(T1,T2)
//  Returns the number of on bits which are set in both T1 and T2.
//
// """ -------------------------------------------------------
template <typename T1, typename T2>
int NumOnBitsInCommon(const T1& bv1, const T2& bv2) {
  return static_cast<int>(OnBitsInCommon(bv1, bv2).size());
}

int NumOnBitsInCommon(const ExplicitBitVect& bv1, const ExplicitBitVect& bv2) {
  return static_cast<int>(((*bv1.dp_bits) & (*bv2.dp_bits)).count());
}

// In all these similarity metrics the notation is selected to be
//   consistent with J.W. Raymond and P. Willett, JCAMD _16_ 59-71 (2002)
// """ -------------------------------------------------------
//
//  TanimotoSimilarity(T1,T2)
//   returns the Tanamoto similarity between T1 and T2, a double.
//
//  T1 and T2 should be the same length.
//
//  C++ Notes: T1 and T2 must support operator&, getNumBits()
//  and getOnBits().
//
//  Python Notes: T1 and T2 are BitVects.
//
// """ -------------------------------------------------------
template <typename T1, typename T2>
double TanimotoSimilarity(const T1& bv1, const T2& bv2) {
  if (bv1.getNumBits() != bv2.getNumBits())
    throw ValueErrorException("BitVects must be same length");
  double x = NumOnBitsInCommon(bv1, bv2);
  double y = bv1.getNumOnBits();
  double z = bv2.getNumOnBits();
  if ((y + z - x) == 0.0)
    return 1.0;
  else
    return x / (y + z - x);
}

template <typename T1, typename T2>
double TverskySimilarity(const T1& bv1, const T2& bv2, double a, double b) {
  RANGE_CHECK(0, a, 1);
  RANGE_CHECK(0, b, 1);
  if (bv1.getNumBits() != bv2.getNumBits())
    throw ValueErrorException("BitVects must be same length");
  double x = NumOnBitsInCommon(bv1, bv2);
  double y = bv1.getNumOnBits();
  double z = bv2.getNumOnBits();
  double denom = a * y + b * z + (1 - a - b) * x;
  if (denom == 0.0)
    return 1.0;
  else
    return x / denom;
}

template <typename T1, typename T2>
double CosineSimilarity(const T1& bv1, const T2& bv2) {
  if (bv1.getNumBits() != bv2.getNumBits())
    throw ValueErrorException("BitVects must be same length");
  double x = NumOnBitsInCommon(bv1, bv2);
  double y = bv1.getNumOnBits();
  double z = bv2.getNumOnBits();

  if (y * z > 0.0) {
    return x / sqrt(y * z);
  } else {
    return 0.0;
  }
}

template <typename T1, typename T2>
double KulczynskiSimilarity(const T1& bv1, const T2& bv2) {
  if (bv1.getNumBits() != bv2.getNumBits())
    throw ValueErrorException("BitVects must be same length");
  double x = NumOnBitsInCommon(bv1, bv2);
  double y = bv1.getNumOnBits();
  double z = bv2.getNumOnBits();

  if (y * z > 0.0) {
    return x * (y + z) / (2 * y * z);
  } else {
    return 0.0;
  }
}

template <typename T1, typename T2>
double DiceSimilarity(const T1& bv1, const T2& bv2) {
  if (bv1.getNumBits() != bv2.getNumBits())
    throw ValueErrorException("BitVects must be same length");
  double x = NumOnBitsInCommon(bv1, bv2);
  double y = bv1.getNumOnBits();
  double z = bv2.getNumOnBits();

  if (y + z > 0.0) {
    return 2 * x / (y + z);
  } else {
    return 0.0;
  }
}

template <typename T1, typename T2>
double SokalSimilarity(const T1& bv1, const T2& bv2) {
  if (bv1.getNumBits() != bv2.getNumBits())
    throw ValueErrorException("BitVects must be same length");
  double x = NumOnBitsInCommon(bv1, bv2);
  double y = bv1.getNumOnBits();
  double z = bv2.getNumOnBits();

  return x / (2 * y + 2 * z - 3 * x);
}

template <typename T1, typename T2>
double McConnaugheySimilarity(const T1& bv1, const T2& bv2) {
  if (bv1.getNumBits() != bv2.getNumBits())
    throw ValueErrorException("BitVects must be same length");
  double x = NumOnBitsInCommon(bv1, bv2);
  double y = bv1.getNumOnBits();
  double z = bv2.getNumOnBits();

  if (y * z > 0.0) {
    return (x * (y + z) - (y * z)) / (y * z);
  } else {
    return 0.0;
  }
}

template <typename T>
inline T tmin(T v1, T v2) {
  if (v1 < v2) return v1;
  return v2;
}

template <typename T>
inline T tmax(T v1, T v2) {
  if (v1 > v2) return v1;
  return v2;
}

template <typename T1, typename T2>
double AsymmetricSimilarity(const T1& bv1, const T2& bv2) {
  if (bv1.getNumBits() != bv2.getNumBits())
    throw ValueErrorException("BitVects must be same length");
  double x = NumOnBitsInCommon(bv1, bv2);
  double y = bv1.getNumOnBits();
  double z = bv2.getNumOnBits();

  if (tmin(y, z) > 0) {
    return x / tmin(y, z);
  } else {
    return 0.0;
  }
}

template <typename T1, typename T2>
double BraunBlanquetSimilarity(const T1& bv1, const T2& bv2) {
  if (bv1.getNumBits() != bv2.getNumBits())
    throw ValueErrorException("BitVects must be same length");
  double x = NumOnBitsInCommon(bv1, bv2);
  double y = bv1.getNumOnBits();
  double z = bv2.getNumOnBits();

  if (tmax(y, z) > 0) {
    return x / tmax(y, z);
  } else {
    return 0.0;
  }
}

template <typename T1, typename T2>
double RusselSimilarity(const T1& bv1, const T2& bv2) {
  if (bv1.getNumBits() != bv2.getNumBits())
    throw ValueErrorException("BitVects must be same length");
  double x = NumOnBitsInCommon(bv1, bv2);
  return x / bv1.getNumBits();
}

template <typename T1, typename T2>
double RogotGoldbergSimilarity(const T1& bv1, const T2& bv2) {
  if (bv1.getNumBits() != bv2.getNumBits())
    throw ValueErrorException("BitVects must be same length");
  double x = NumOnBitsInCommon(bv1, bv2);
  double y = bv1.getNumOnBits();
  double z = bv2.getNumOnBits();
  double l = bv1.getNumBits();
  double d = l - y - z + x;
  if ((x == l) || (d == l))
    return 1.0;
  else
    return (x / (y + z) + (d) / (2 * l - y - z));
}

// """ -------------------------------------------------------
//
//  OnBitSimilarity(T1,T2)
//  Returns the percentage of possible on bits in common
//  between T1 and T2 (a double)
//
//  C++ Notes: T1 and T2 must support operator|, operator&
//  and getOnBits().
//
//  Python Notes: T1 and T2 are BitVects.
//
// """ -------------------------------------------------------
template <typename T1, typename T2>
double OnBitSimilarity(const T1& bv1, const T2& bv2) {
  if (bv1.getNumBits() != bv2.getNumBits())
    throw ValueErrorException("BitVects must be same length");

  double num = NumOnBitsInCommon(bv1, bv2);
  double denom = (bv1 | bv2).getNumOnBits();

  if (denom > 0) {
    return num / denom;
  } else {
    return 0;
  }
}

// """ -------------------------------------------------------
//
//  NumBitsInCommon(T1,T2)
//  Returns the number of bits in common (on and off)
//  between T1 and T2 (an int)
//
//  T1 and T2 should be the same length.
//
//  C++ Notes: T1 and T2 must support operator^, getNumBits().
//
//  Python Notes: T1 and T2 are BitVects.
//
// """ -------------------------------------------------------
template <typename T1, typename T2>
int NumBitsInCommon(const T1& bv1, const T2& bv2) {
  if (bv1.getNumBits() != bv2.getNumBits())
    throw ValueErrorException("BitVects must be same length");

  return bv1.getNumBits() - (bv1 ^ bv2).getNumOnBits();
}

int NumBitsInCommon(const ExplicitBitVect& bv1, const ExplicitBitVect& bv2) {
  return bv1.getNumBits() - static_cast<int>(((*bv1.dp_bits) ^ (*bv2.dp_bits)).count());
}

// """ -------------------------------------------------------
//
//  AllBitSimilarity(T1,T2)
//  Returns the percentage of bits in common (on and off)
//  between T1 and T2 (a double)
//
//  T1 and T2 should be the same length.
//
//  C++ Notes: T1 and T2 must support operator^, getNumBits()
//  and getNumOnBits().
//
//  Python Notes: T1 and T2 are BitVects.
//
// """ -------------------------------------------------------
template <typename T1, typename T2>
double AllBitSimilarity(const T1& bv1, const T2& bv2) {
  if (bv1.getNumBits() != bv2.getNumBits())
    throw ValueErrorException("BitVects must be same length");

  return double(NumBitsInCommon(bv1, bv2)) / bv1.getNumBits();
}

// """ -------------------------------------------------------
//
//  OnBitsInCommon(T1,T2)
//  Returns the on bits which are set in both T1 and T2.
//
//  T1 and T2 should be the same length.
//
//  C++ Notes: T1 and T2 must support operator&, getNumBits()
//  and getOnBits(), the return value is an IntVect.
//
//  Python Notes: T1 and T2 are BitVects, the return value
//  is a tuple of ints.
//
// """ -------------------------------------------------------
template <typename T1, typename T2>
IntVect OnBitsInCommon(const T1& bv1, const T2& bv2) {
  if (bv1.getNumBits() != bv2.getNumBits())
    throw ValueErrorException("BitVects must be same length");
  IntVect res;
  (bv1 & bv2).getOnBits(res);
  return res;
}

// """ -------------------------------------------------------
//
//  OffBitsInCommon(T1,T2)
//  Returns the off bits which are set in both T1 and T2.
//
//  T1 and T2 should be the same length.
//
//  C++ Notes: T1 and T2 must support operator|, operator~,
// getNumBits() and getOnBits(), the return value is an IntVect.
//
//  Python Notes: T1 and T2 are BitVects, the return value
//  is a tuple of ints.
//
// """ -------------------------------------------------------
template <typename T1, typename T2>
IntVect OffBitsInCommon(const T1& bv1, const T2& bv2) {
  if (bv1.getNumBits() != bv2.getNumBits())
    throw ValueErrorException("BitVects must be same length");
  IntVect res;
  (~(bv1 | bv2)).getOnBits(res);
  return res;
}

// """ -------------------------------------------------------
//
//  OnBitProjSimilarity(T1,T2)
//  Returns the projected similarity between the on bits of
//  T1 and T2.
//
//  The on bit projected similarity of T1 onto T2 is the
//  percentage of T1's on bits which are on in T2.
//
//  This type of measure may be useful for substructure-type
//  searches.
//
//  Two values are returned, the projection of T1 onto T2
//  and the projection of T2 onto T1
//
//  T1 and T2 should be the same length.
//
//  C++ Notes: T1 and T2 must support operator&, getNumBits()
//  and getNumOnBits(), the return value is an DoubleVect with
//  two elements.
//
//  Python Notes: T1 and T2 are BitVects, the return value
//  is a 2-tuple of doubles.
//
// """ -------------------------------------------------------
template <typename T1, typename T2>
DoubleVect OnBitProjSimilarity(const T1& bv1, const T2& bv2) {
  if (bv1.getNumBits() != bv2.getNumBits())
    throw ValueErrorException("BitVects must be same length");
  DoubleVect res(2, 0.0);
  double num = NumOnBitsInCommon(bv1, bv2);
  if (num) {
    res[0] = num / bv1.getNumOnBits();
    res[1] = num / bv2.getNumOnBits();
  }
  return res;
}

// """ -------------------------------------------------------
//
//  OffBitProjSimilarity(T1,T2)
//  Returns the projected similarity between the off bits of
//  T1 and T2.
//
//  The off bit projected similarity of T1 onto T2 is the
//  percentage of T1's off bits which are off in T2.
//
//  This type of measure may be useful for substructure-type
//  searches.
//
//  Two values are returned, the projection of T1 onto T2
//  and the projection of T2 onto T1
//
//  T1 and T2 should be the same length.
//
//  C++ Notes: T1 and T2 must support operator|, getNumBits()
//  and getNumOffBits(), the return value is an DoubleVect with
//  two elements.
//
//  Python Notes: T1 and T2 are BitVects, the return value
//  is a 2-tuple of doubles.
//
// """ -------------------------------------------------------
template <typename T1, typename T2>
DoubleVect OffBitProjSimilarity(const T1& bv1, const T2& bv2) {
  if (bv1.getNumBits() != bv2.getNumBits())
    throw ValueErrorException("BitVects must be same length");
  DoubleVect res(2, 0.0);
  double num = (bv1 | bv2).getNumOffBits();
  if (num) {
    res[0] = num / bv1.getNumOffBits();
    res[1] = num / bv2.getNumOffBits();
  }
  return res;
}

template <typename T1>
T1* FoldFingerprint(const T1& bv1, unsigned int factor) {
  if (factor <= 0 || factor >= bv1.getNumBits())
    throw ValueErrorException("invalid fold factor");

  int initSize = bv1.getNumBits();
  int resSize = initSize / factor;
  T1* res = new T1(resSize);

  IntVect onBits;
  bv1.getOnBits(onBits);
  for (IntVectIter iv = onBits.begin(); iv != onBits.end(); iv++) {
    int pos = (*iv) % resSize;
    res->setBit(pos);
  }
  return res;
}

template <typename T1>
std::string BitVectToText(const T1& bv1) {
  std::string res(bv1.getNumBits(), '0');
  for (unsigned int i = 0; i < bv1.getNumBits(); i++) {
    if (bv1.getBit(i)) res[i] = '1';
  }
  return res;
}

const char bin2Hex[] = {'0', '1', '2', '3', '4', '5', '6', '7',
                        '8', '9', 'a', 'b', 'c', 'd', 'e', 'f'};
template <typename T1>
std::string BitVectToFPSText(const T1& bv1) {
  unsigned int size =
      2 * (bv1.getNumBits() / 8 + (bv1.getNumBits() % 8 ? 1 : 0));
  std::string res(size, 0);
  unsigned char c = 0;
  unsigned int byte = 0;
  for (unsigned int i = 0; i < bv1.getNumBits(); i++) {
    if (bv1.getBit(i)) {
      c |= 1 << (i % 8);
    }
    if (!((i + 1) % 8)) {
      res[byte++] = bin2Hex[(c >> 4) % 16];
      res[byte++] = bin2Hex[c % 16];
      c = 0;
    }
  }
  if (byte < size) {
    res[byte++] = bin2Hex[(c >> 4) % 16];
    res[byte++] = bin2Hex[c % 16];
  }
  return res;
}

template <typename T1>
std::string BitVectToBinaryText(const T1& bv1) {
  std::string res(bv1.getNumBits() / 8 + (bv1.getNumBits() % 8 ? 1 : 0), 0);
  unsigned char c = 0;
  unsigned int byte = 0;
  for (unsigned int i = 0; i < bv1.getNumBits(); i++) {
    if (bv1.getBit(i)) {
      c |= 1 << (i % 8);
    }
    if (!((i + 1) % 8)) {
      res[byte++] = c;
      c = 0;
    }
  }
  if (bv1.getNumBits() % 8) {
    res[byte] = c;
  }
  return res;
}

template <typename T1>
void UpdateBitVectFromFPSText(T1& bv1, const std::string& fps) {
  PRECONDITION(fps.length() * 4 >= bv1.getNumBits(), "bad FPS length");
  PRECONDITION(fps.length() % 2 == 0, "bad FPS length");
  unsigned int bitIdx = 0;
  char tptr[3];
  tptr[2] = (char)0;
  for (unsigned int i = 0; i < fps.size() && bitIdx < bv1.getNumBits();
       i += 2) {
    unsigned short c = 0;
    try {
      tptr[0] = fps[i];
      tptr[1] = fps[i + 1];
      c = static_cast<unsigned short>(strtol(tptr, NULL, 16));
    } catch (...) {
      std::ostringstream errout;
      errout << "Cannot convert FPS word: " << fps.substr(i, 2) << " to int";
      std::cerr << errout.str() << std::endl;
      throw ValueErrorException(errout.str());
    }
    for (unsigned int bit = 0; bit < 8 && bitIdx < bv1.getNumBits();
         ++bit, ++bitIdx) {
      if (c & (1 << bit)) bv1.setBit(bitIdx);
    }
  }
}

template <typename T1>
void UpdateBitVectFromBinaryText(T1& bv1, const std::string& fps) {
  PRECONDITION(fps.length() * 8 >= bv1.getNumBits(), "bad FPS length");
  unsigned int bitIdx = 0;
  for (unsigned int i = 0; i < fps.size() && bitIdx < bv1.getNumBits(); i++) {
    unsigned short c = fps[i];
    for (unsigned int bit = 0; bit < 8 && bitIdx < bv1.getNumBits();
         ++bit, ++bitIdx) {
      if (c & (1 << bit)) bv1.setBit(bitIdx);
    }
  }
}

template double TanimotoSimilarity(const SparseBitVect& bv1,
                                   const SparseBitVect& bv2);
template double TverskySimilarity(const SparseBitVect& bv1,
                                  const SparseBitVect& bv2, double a, double b);
template double CosineSimilarity(const SparseBitVect& bv1,
                                 const SparseBitVect& bv2);
template double KulczynskiSimilarity(const SparseBitVect& bv1,
                                     const SparseBitVect& bv2);
template double DiceSimilarity(const SparseBitVect& bv1,
                               const SparseBitVect& bv2);
template double SokalSimilarity(const SparseBitVect& bv1,
                                const SparseBitVect& bv2);
template double McConnaugheySimilarity(const SparseBitVect& bv1,
                                       const SparseBitVect& bv2);
template double AsymmetricSimilarity(const SparseBitVect& bv1,
                                     const SparseBitVect& bv2);
template double BraunBlanquetSimilarity(const SparseBitVect& bv1,
                                        const SparseBitVect& bv2);
template double RusselSimilarity(const SparseBitVect& bv1,
                                 const SparseBitVect& bv2);
template double RogotGoldbergSimilarity(const SparseBitVect& bv1,
                                        const SparseBitVect& bv2);
template double OnBitSimilarity(const SparseBitVect& bv1,
                                const SparseBitVect& bv2);
template int NumBitsInCommon(const SparseBitVect& bv1,
                             const SparseBitVect& bv2);
template double AllBitSimilarity(const SparseBitVect& bv1,
                                 const SparseBitVect& bv2);
template int NumOnBitsInCommon(const SparseBitVect& bv1,
                               const SparseBitVect& bv2);
template IntVect OnBitsInCommon(const SparseBitVect& bv1,
                                const SparseBitVect& bv2);
template IntVect OffBitsInCommon(const SparseBitVect& bv1,
                                 const SparseBitVect& bv2);
template DoubleVect OnBitProjSimilarity(const SparseBitVect& bv1,
                                        const SparseBitVect& bv2);
template DoubleVect OffBitProjSimilarity(const SparseBitVect& bv1,
                                         const SparseBitVect& bv2);

template double TanimotoSimilarity(const ExplicitBitVect& bv1,
                                   const ExplicitBitVect& bv2);
template double TverskySimilarity(const ExplicitBitVect& bv1,
                                  const ExplicitBitVect& bv2, double a,
                                  double b);
template double CosineSimilarity(const ExplicitBitVect& bv1,
                                 const ExplicitBitVect& bv2);
template double KulczynskiSimilarity(const ExplicitBitVect& bv1,
                                     const ExplicitBitVect& bv2);
template double DiceSimilarity(const ExplicitBitVect& bv1,
                               const ExplicitBitVect& bv2);
template double SokalSimilarity(const ExplicitBitVect& bv1,
                                const ExplicitBitVect& bv2);
template double McConnaugheySimilarity(const ExplicitBitVect& bv1,
                                       const ExplicitBitVect& bv2);
template double AsymmetricSimilarity(const ExplicitBitVect& bv1,
                                     const ExplicitBitVect& bv2);
template double BraunBlanquetSimilarity(const ExplicitBitVect& bv1,
                                        const ExplicitBitVect& bv2);
template double RusselSimilarity(const ExplicitBitVect& bv1,
                                 const ExplicitBitVect& bv2);
template double RogotGoldbergSimilarity(const ExplicitBitVect& bv1,
                                        const ExplicitBitVect& bv2);
template double OnBitSimilarity(const ExplicitBitVect& bv1,
                                const ExplicitBitVect& bv2);
template int NumBitsInCommon(const ExplicitBitVect& bv1,
                             const ExplicitBitVect& bv2);
template double AllBitSimilarity(const ExplicitBitVect& bv1,
                                 const ExplicitBitVect& bv2);
template IntVect OnBitsInCommon(const ExplicitBitVect& bv1,
                                const ExplicitBitVect& bv2);
template IntVect OffBitsInCommon(const ExplicitBitVect& bv1,
                                 const ExplicitBitVect& bv2);
template DoubleVect OnBitProjSimilarity(const ExplicitBitVect& bv1,
                                        const ExplicitBitVect& bv2);
template DoubleVect OffBitProjSimilarity(const ExplicitBitVect& bv1,
                                         const ExplicitBitVect& bv2);

template SparseBitVect* FoldFingerprint(const SparseBitVect&, unsigned int);
template ExplicitBitVect* FoldFingerprint(const ExplicitBitVect&, unsigned int);

template std::string BitVectToText(const SparseBitVect&);
template std::string BitVectToText(const ExplicitBitVect&);

template std::string BitVectToFPSText(const SparseBitVect&);
template std::string BitVectToFPSText(const ExplicitBitVect&);
template void UpdateBitVectFromFPSText(SparseBitVect&, const std::string&);
template void UpdateBitVectFromFPSText(ExplicitBitVect&, const std::string&);

template std::string BitVectToBinaryText(const SparseBitVect&);
template std::string BitVectToBinaryText(const ExplicitBitVect&);
template void UpdateBitVectFromBinaryText(SparseBitVect&, const std::string&);
template void UpdateBitVectFromBinaryText(ExplicitBitVect&, const std::string&);

// from here:
// http://stackoverflow.com/questions/3849337/msvc-equivalent-to-builtin-popcount
// but corrected to get the ifdef right
#ifdef _MSC_VER
#include <intrin.h>
#ifdef _WIN64
#define BUILTIN_POPCOUNT_INSTR __popcnt64
#define BUILTIN_POPCOUNT_TYPE boost::uint64_t
#else
#define BUILTIN_POPCOUNT_INSTR __popcnt
#define BUILTIN_POPCOUNT_TYPE boost::uint32_t
#endif
#else
#define BUILTIN_POPCOUNT_INSTR __builtin_popcountll
#define BUILTIN_POPCOUNT_TYPE boost::uint64_t
#endif

// the Bitmap Tanimoto and Dice similarity code is adapted
// from Andrew Dalke's chem-fingerprints code
// http://code.google.com/p/chem-fingerprints/
namespace {
static int byte_popcounts[] = {
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};
}
unsigned int CalcBitmapPopcount(const unsigned char* afp, unsigned int nBytes) {
  PRECONDITION(afp, "no afp");
  unsigned int popcount = 0;
#ifndef USE_BUILTIN_POPCOUNT
  for (unsigned int i = 0; i < nBytes; i++) {
    popcount += byte_popcounts[afp[i]];
  }
#else
  unsigned int eidx = nBytes / sizeof(BUILTIN_POPCOUNT_TYPE);
  for (unsigned int i = 0; i < eidx; ++i) {
    popcount += static_cast<unsigned int>(BUILTIN_POPCOUNT_INSTR(((BUILTIN_POPCOUNT_TYPE*)afp)[i]));
  }
  for (unsigned int i = eidx * sizeof(BUILTIN_POPCOUNT_TYPE); i < nBytes; ++i) {
    popcount += byte_popcounts[afp[i]];
  }
#endif
  return popcount;
}
double CalcBitmapTanimoto(const unsigned char* afp, const unsigned char* bfp,
                          unsigned int nBytes) {
  PRECONDITION(afp, "no afp");
  PRECONDITION(bfp, "no bfp");
  unsigned int union_popcount = 0, intersect_popcount = 0;
#ifndef USE_BUILTIN_POPCOUNT
  for (unsigned int i = 0; i < nBytes; i++) {
    union_popcount += byte_popcounts[afp[i] | bfp[i]];
    intersect_popcount += byte_popcounts[afp[i] & bfp[i]];
  }
#else
  BUILTIN_POPCOUNT_TYPE eidx = nBytes / sizeof(BUILTIN_POPCOUNT_TYPE);
  for (BUILTIN_POPCOUNT_TYPE i = 0; i < eidx; ++i) {
    union_popcount += static_cast<unsigned int>(BUILTIN_POPCOUNT_INSTR(((BUILTIN_POPCOUNT_TYPE*)afp)[i] |
                                           ((BUILTIN_POPCOUNT_TYPE*)bfp)[i]));
    intersect_popcount += static_cast<unsigned int>(BUILTIN_POPCOUNT_INSTR(((BUILTIN_POPCOUNT_TYPE*)afp)[i] &
                                               ((BUILTIN_POPCOUNT_TYPE*)bfp)[i]));
  }
  for (BUILTIN_POPCOUNT_TYPE i = eidx * sizeof(BUILTIN_POPCOUNT_TYPE); i < nBytes; ++i) {
    union_popcount += byte_popcounts[afp[i] | bfp[i]];
    intersect_popcount += byte_popcounts[afp[i] & bfp[i]];
  }
#endif
  if (union_popcount == 0) {
    return 0.0;
  }
  return (intersect_popcount + 0.0) /
         union_popcount; /* +0.0 to coerce to double */
}

double CalcBitmapDice(const unsigned char* afp, const unsigned char* bfp,
                      unsigned int nBytes) {
  PRECONDITION(afp, "no afp");
  PRECONDITION(bfp, "no bfp");
  unsigned int intersect_popcount = 0, a_popcount = 0, b_popcount = 0;

#ifndef USE_BUILTIN_POPCOUNT
  for (unsigned int i = 0; i < nBytes; i++) {
    a_popcount += byte_popcounts[afp[i]];
    b_popcount += byte_popcounts[bfp[i]];
    intersect_popcount += byte_popcounts[afp[i] & bfp[i]];
  }
#else
  BUILTIN_POPCOUNT_TYPE eidx = nBytes / sizeof(BUILTIN_POPCOUNT_TYPE);
  for (BUILTIN_POPCOUNT_TYPE i = 0; i < eidx; ++i) {
    a_popcount += static_cast<unsigned int>(BUILTIN_POPCOUNT_INSTR(((BUILTIN_POPCOUNT_TYPE*)afp)[i]));
    b_popcount += static_cast<unsigned int>(BUILTIN_POPCOUNT_INSTR(((BUILTIN_POPCOUNT_TYPE*)bfp)[i]));
    intersect_popcount += static_cast<unsigned int>(BUILTIN_POPCOUNT_INSTR(((BUILTIN_POPCOUNT_TYPE*)afp)[i] &
                                               ((BUILTIN_POPCOUNT_TYPE*)bfp)[i]));
  }
  for (BUILTIN_POPCOUNT_TYPE i = eidx * sizeof(BUILTIN_POPCOUNT_TYPE); i < nBytes; ++i) {
    a_popcount += byte_popcounts[afp[i]];
    b_popcount += byte_popcounts[bfp[i]];
    intersect_popcount += byte_popcounts[afp[i] & bfp[i]];
  }
#endif

  if (a_popcount + b_popcount == 0) {
    return 0.0;
  }
  return (2.0 * intersect_popcount) / (a_popcount + b_popcount);
}

double CalcBitmapTversky(const unsigned char* afp, const unsigned char* bfp,
                         unsigned int nBytes, double ca, double cb) {
  PRECONDITION(afp, "no afp");
  PRECONDITION(bfp, "no bfp");
  unsigned int intersect_popcount = 0, acount = 0, bcount = 0;

#ifndef USE_BUILTIN_POPCOUNT
  for (unsigned int i = 0; i < nBytes; i++) {
    intersect_popcount += byte_popcounts[afp[i] & bfp[i]];
    acount += byte_popcounts[afp[i]];
    bcount += byte_popcounts[bfp[i]];
  }
#else
  BUILTIN_POPCOUNT_TYPE eidx = nBytes / sizeof(BUILTIN_POPCOUNT_TYPE);
  for (BUILTIN_POPCOUNT_TYPE i = 0; i < eidx; ++i) {
    intersect_popcount += static_cast<unsigned int>(BUILTIN_POPCOUNT_INSTR(((BUILTIN_POPCOUNT_TYPE*)afp)[i] &
                                               ((BUILTIN_POPCOUNT_TYPE*)bfp)[i]));
    acount += static_cast<unsigned int>(BUILTIN_POPCOUNT_INSTR(((BUILTIN_POPCOUNT_TYPE*)afp)[i]));
    bcount += static_cast<unsigned int>(BUILTIN_POPCOUNT_INSTR(((BUILTIN_POPCOUNT_TYPE*)bfp)[i]));
  }
  for (BUILTIN_POPCOUNT_TYPE i = eidx * sizeof(BUILTIN_POPCOUNT_TYPE); i < nBytes; ++i) {
    intersect_popcount += byte_popcounts[afp[i] & bfp[i]];
    acount += byte_popcounts[afp[i]];
    bcount += byte_popcounts[bfp[i]];
  }
#endif
  double denom = ca * acount + cb * bcount + (1 - ca - cb) * intersect_popcount;
  if (denom == 0.0) {
    return 0.0;
  }
  return intersect_popcount / denom;
}

bool CalcBitmapAllProbeBitsMatch(const unsigned char* probe,
                                 const unsigned char* ref,
                                 unsigned int nBytes) {
  PRECONDITION(probe, "no probe");
  PRECONDITION(ref, "no ref");

#ifndef USE_BUILTIN_POPCOUNT
  for (unsigned int i = 0; i < nBytes; i++) {
    if (byte_popcounts[probe[i] & ref[i]] != byte_popcounts[probe[i]]) {
      return false;
    }
  }
#else
  unsigned int eidx = nBytes / sizeof(BUILTIN_POPCOUNT_TYPE);
  for (unsigned int i = 0; i < eidx; ++i) {
    if (BUILTIN_POPCOUNT_INSTR(((BUILTIN_POPCOUNT_TYPE*)probe)[i] &
                           ((BUILTIN_POPCOUNT_TYPE*)ref)[i]) !=
        BUILTIN_POPCOUNT_INSTR(((BUILTIN_POPCOUNT_TYPE*)probe)[i])) {
      return false;
    }
  }
  for (unsigned int i = eidx * sizeof(BUILTIN_POPCOUNT_TYPE); i < nBytes; ++i) {
    if (byte_popcounts[probe[i] & ref[i]] != byte_popcounts[probe[i]]) {
      return false;
    }
  }
#endif
  return true;
}
