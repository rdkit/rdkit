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
#include <cmath>
#include <string>
#include <iostream>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/types.h>
#include <RDGeneral/Exceptions.h>
#include <sstream>
#include <cstdlib>
#include <algorithm>

#include <boost/lexical_cast.hpp>

#if _MSC_VER
#include <intrin.h>
#endif

// SIMD support for vectorized bitmap operations (x86/x64 only)
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
#define RDK_X86_SIMD
#ifndef _MSC_VER
#include <immintrin.h>
#endif
#endif

using namespace RDKit;

int getBitId(const char *&text, int format, int size, int curr) {
  PRECONDITION(text, "no text");
  int res = -1;
  if ((format == 0) ||
      ((format == 1) && (size >= std::numeric_limits<unsigned short>::max()))) {
    int tmp =
        EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(*(int *)text);
    text += sizeof(tmp);
    res = tmp;
  } else if (format == 1) {  // version 16 and on bits sorted as short ints
    unsigned short tmp =
        EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(
            *(unsigned short *)text);
    text += sizeof(tmp);
    res = tmp;
  } else if (format == 2) {  // run length encoded format
    res = curr + RDKit::pullPackedIntFromString(text);
  }
  return res;
}

bool AllProbeBitsMatch(const std::string &probe, const std::string &ref) {
  return AllProbeBitsMatch(probe.c_str(), ref.c_str());
}

bool AllProbeBitsMatch(const char *probe, const char *ref) {
  PRECONDITION(probe, "no probe text");
  PRECONDITION(ref, "no probe text");
  int probeFormat = 0;
  int refFormat = 0;
  int version = 0;

  int probeSize =
      EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(*(int *)probe);
  probe += sizeof(probeSize);
  if (probeSize < 0) {
    version = -1 * probeSize;
    switch (version) {
      case 16:
        probeFormat = 1;
        break;
      case 32:
        probeFormat = 2;
        break;
      default:
        throw ValueErrorException(
            "Unknown version type for the encode bit vect");
        break;
    }
    probeSize =
        EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(*(int *)probe);
    probe += sizeof(probeSize);
  }

  int refSize =
      EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(*(int *)ref);
  ref += sizeof(refSize);
  if (refSize < 0) {
    version = -1 * refSize;
    switch (version) {
      case 16:
        refFormat = 1;
        break;
      case 32:
        refFormat = 2;
        break;
      default:
        throw ValueErrorException(
            "Unknown version type for the encode bit vect");
        break;
    }
    refSize =
        EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(*(int *)ref);
    ref += sizeof(refSize);
  }

  int nProbeOn =
      EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(*(int *)probe);
  probe += sizeof(nProbeOn);
  int nRefOn =
      EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(*(int *)ref);
  ref += sizeof(nRefOn);

  int currProbeBit = 0;
  currProbeBit = getBitId(probe, probeFormat, probeSize, currProbeBit);
  nProbeOn--;

  int currRefBit = 0;
  currRefBit = getBitId(ref, refFormat, refSize, currRefBit);
  nRefOn--;

  while (nProbeOn) {
    while (currRefBit < currProbeBit && nRefOn > 0) {
      if (refFormat == 2) {
        currRefBit++;
      }
      currRefBit = getBitId(ref, refFormat, refSize, currRefBit);
      nRefOn--;
    }
    if (currRefBit != currProbeBit) {
      return false;
    }
    if (probeFormat == 2) {
      currProbeBit++;
    }
    currProbeBit = getBitId(probe, probeFormat, probeSize, currProbeBit);
    nProbeOn--;
  }
  return true;
}

template <typename T1>
bool AllProbeBitsMatch(const T1 &probe, const std::string &pkl) {
  const char *text = pkl.c_str();
  int format = 0;
  int nOn = 0, size, version = 0;
  size = EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(*(int *)text);
  text += sizeof(size);
  if (size < 0) {
    version = -1 * size;
    switch (version) {
      case 16:
        format = 1;
        break;
      case 32:
        format = 2;
        break;
      default:
        throw ValueErrorException(
            "Unknown version type for the encode bit vect");
        break;
    }
    size =
        EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(*(int *)text);
    text += sizeof(size);
  }
  nOn = EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(*(int *)text);
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
      if (format == 2) {
        currBit++;
      }
      currBit = getBitId(text, format, size, currBit);
      nOn--;
    }
    if (currBit != *i) {
      return false;
    }
    //}
  }
  return true;
}
template RDKIT_DATASTRUCTS_EXPORT bool AllProbeBitsMatch(
    const SparseBitVect &bv1, const std::string &pkl);
template RDKIT_DATASTRUCTS_EXPORT bool AllProbeBitsMatch(
    const ExplicitBitVect &bv1, const std::string &pkl);
template <typename T1>
bool AllProbeBitsMatch(const T1 &probe, const T1 &ref) {
  for (unsigned int i = 0; i < probe.getNumBits(); ++i) {
    if (probe.getBit(i) && !ref.getBit(i)) {
      return false;
    }
  }
  return true;
}
template RDKIT_DATASTRUCTS_EXPORT bool AllProbeBitsMatch(
    const SparseBitVect &bv1, const SparseBitVect &bv2);
// template bool AllProbeBitsMatch(const ExplicitBitVect& bv1,const
// ExplicitBitVect &bv2);

bool AllProbeBitsMatch(const ExplicitBitVect &probe,
                       const ExplicitBitVect &ref) {
  return probe.dp_bits->is_subset_of(*(ref.dp_bits));
}

// """ -------------------------------------------------------
//
//  NumOnBitsInCommon(T1,T2)
//  Returns the number of on bits which are set in both T1 and T2.
//
// """ -------------------------------------------------------
template <typename T1, typename T2>
int NumOnBitsInCommon(const T1 &bv1, const T2 &bv2) {
  return static_cast<int>(OnBitsInCommon(bv1, bv2).size());
}

namespace {
struct bitset_impl {
  std::vector<unsigned long> m_bits;
  std::size_t m_num_bits;
};

const bool canUseBitmapHack =
    sizeof(boost::dynamic_bitset<>) == sizeof(bitset_impl);

bool EBVToBitmap(const ExplicitBitVect &bv, const unsigned char *&fp,
                 unsigned int &nBytes) {
  if (!canUseBitmapHack) {
    return false;
  }
  const auto *p1 = (const bitset_impl *)(const void *)bv.dp_bits.get();
  // Run-time sanity check (just in case)
  if (p1->m_num_bits != bv.dp_bits->size()) {
    return false;
  }
  fp = (const unsigned char *)p1->m_bits.data();
  nBytes = (unsigned int)p1->m_num_bits / 8;
  if (p1->m_num_bits % 8) {
    ++nBytes;
  }
  return true;
}
}  // namespace

unsigned int CalcBitmapNumBitsInCommon(const unsigned char *afp,
                                       const unsigned char *bfp,
                                       unsigned int nBytes);

int NumOnBitsInCommon(const ExplicitBitVect &bv1, const ExplicitBitVect &bv2) {
  // Don't try this at home, we (hope we) know what we're doing
  const unsigned char *afp, *bfp;
  unsigned int nBytes;
  if (EBVToBitmap(bv1, afp, nBytes) && EBVToBitmap(bv2, bfp, nBytes)) {
    unsigned int result = CalcBitmapNumBitsInCommon(afp, bfp, nBytes);
    return (int)result;
  }

  return static_cast<int>(((*bv1.dp_bits) & (*bv2.dp_bits)).count());
}

// In all these similarity metrics the notation is selected to be
//   consistent with J.W. Raymond and P. Willett, JCAMD _16_ 59-71 (2002)
// """ -------------------------------------------------------
//
//  TanimotoSimilarity(T1,T2)
//   returns the Tanimoto similarity between T1 and T2, a double.
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
double TanimotoSimilarity(const T1 &bv1, const T2 &bv2) {
  if (bv1.getNumBits() != bv2.getNumBits()) {
    throw ValueErrorException("BitVects must be same length");
  }
  unsigned int total = bv1.getNumOnBits() + bv2.getNumOnBits();
  if (total == 0) {
    return 0.0;
  }
  unsigned int common = NumOnBitsInCommon(bv1, bv2);
  return (double)common / (double)(total - common);
}

template <typename T1, typename T2>
double TverskySimilarity(const T1 &bv1, const T2 &bv2, double a, double b) {
  RANGE_CHECK(0, a, 1);
  RANGE_CHECK(0, b, 1);
  if (bv1.getNumBits() != bv2.getNumBits()) {
    throw ValueErrorException("BitVects must be same length");
  }
  double x = NumOnBitsInCommon(bv1, bv2);
  auto y = bv1.getNumOnBits();
  auto z = bv2.getNumOnBits();
  if (y == 0 || z == 0) {
    return 0.0;
  }
  double denom = a * y + b * z + (1 - a - b) * x;
  if (denom == 0.0) {
    return 1.0;
  } else {
    return x / denom;
  }
}

template <typename T1, typename T2>
double CosineSimilarity(const T1 &bv1, const T2 &bv2) {
  if (bv1.getNumBits() != bv2.getNumBits()) {
    throw ValueErrorException("BitVects must be same length");
  }
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
double KulczynskiSimilarity(const T1 &bv1, const T2 &bv2) {
  if (bv1.getNumBits() != bv2.getNumBits()) {
    throw ValueErrorException("BitVects must be same length");
  }
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
double DiceSimilarity(const T1 &bv1, const T2 &bv2) {
  if (bv1.getNumBits() != bv2.getNumBits()) {
    throw ValueErrorException("BitVects must be same length");
  }
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
double SokalSimilarity(const T1 &bv1, const T2 &bv2) {
  if (bv1.getNumBits() != bv2.getNumBits()) {
    throw ValueErrorException("BitVects must be same length");
  }
  double x = NumOnBitsInCommon(bv1, bv2);
  auto y = bv1.getNumOnBits();
  auto z = bv2.getNumOnBits();
  if (y == 0 || z == 0) {
    return 0.0;
  }

  return x / (2. * y + 2. * z - 3. * x);
}

template <typename T1, typename T2>
double McConnaugheySimilarity(const T1 &bv1, const T2 &bv2) {
  if (bv1.getNumBits() != bv2.getNumBits()) {
    throw ValueErrorException("BitVects must be same length");
  }
  double x = NumOnBitsInCommon(bv1, bv2);
  double y = bv1.getNumOnBits();
  double z = bv2.getNumOnBits();

  if (y * z > 0.0) {
    return (x * (y + z) - (y * z)) / (y * z);
  } else {
    return 0.0;
  }
}

template <typename T1, typename T2>
double AsymmetricSimilarity(const T1 &bv1, const T2 &bv2) {
  if (bv1.getNumBits() != bv2.getNumBits()) {
    throw ValueErrorException("BitVects must be same length");
  }
  double x = NumOnBitsInCommon(bv1, bv2);
  double y = bv1.getNumOnBits();
  double z = bv2.getNumOnBits();

  double min = std::min(y, z);
  if (min > 0.0) {
    return x / min;
  } else {
    return 0.0;
  }
}

template <typename T1, typename T2>
double BraunBlanquetSimilarity(const T1 &bv1, const T2 &bv2) {
  if (bv1.getNumBits() != bv2.getNumBits()) {
    throw ValueErrorException("BitVects must be same length");
  }
  double x = NumOnBitsInCommon(bv1, bv2);
  double y = bv1.getNumOnBits();
  double z = bv2.getNumOnBits();

  double max = std::max(y, z);
  if (max > 0.0) {
    return x / max;
  } else {
    return 0.0;
  }
}

template <typename T1, typename T2>
double RusselSimilarity(const T1 &bv1, const T2 &bv2) {
  if (bv1.getNumBits() != bv2.getNumBits()) {
    throw ValueErrorException("BitVects must be same length");
  }

  double x = NumOnBitsInCommon(bv1, bv2);
  return x / bv1.getNumBits();
}

template <typename T1, typename T2>
double RogotGoldbergSimilarity(const T1 &bv1, const T2 &bv2) {
  if (bv1.getNumBits() != bv2.getNumBits()) {
    throw ValueErrorException("BitVects must be same length");
  }
  double x = NumOnBitsInCommon(bv1, bv2);
  auto y = bv1.getNumOnBits();
  auto z = bv2.getNumOnBits();
  if (y == 0 || z == 0) {
    return 0.0;
  }

  double l = bv1.getNumBits();
  double d = l - y - z + x;

  double denom1 = y + z;
  double denom2 = 2 * l - y - z;
  if ((x == l) || (d == l)) {
    return 1.0;
  } else if (denom1 == 0 || denom2 == 0) {
    return 0.0;
  } else {
    return (x / (y + z) + (d) / (2 * l - y - z));
  }
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
double OnBitSimilarity(const T1 &bv1, const T2 &bv2) {
  if (bv1.getNumBits() != bv2.getNumBits()) {
    throw ValueErrorException("BitVects must be same length");
  }

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
int NumBitsInCommon(const T1 &bv1, const T2 &bv2) {
  if (bv1.getNumBits() != bv2.getNumBits()) {
    throw ValueErrorException("BitVects must be same length");
  }

  return bv1.getNumBits() - (bv1 ^ bv2).getNumOnBits();
}

int NumBitsInCommon(const ExplicitBitVect &bv1, const ExplicitBitVect &bv2) {
  return bv1.getNumBits() -
         static_cast<int>(((*bv1.dp_bits) ^ (*bv2.dp_bits)).count());
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
double AllBitSimilarity(const T1 &bv1, const T2 &bv2) {
  if (bv1.getNumBits() != bv2.getNumBits()) {
    throw ValueErrorException("BitVects must be same length");
  }

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
IntVect OnBitsInCommon(const T1 &bv1, const T2 &bv2) {
  if (bv1.getNumBits() != bv2.getNumBits()) {
    throw ValueErrorException("BitVects must be same length");
  }
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
IntVect OffBitsInCommon(const T1 &bv1, const T2 &bv2) {
  if (bv1.getNumBits() != bv2.getNumBits()) {
    throw ValueErrorException("BitVects must be same length");
  }
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
DoubleVect OnBitProjSimilarity(const T1 &bv1, const T2 &bv2) {
  if (bv1.getNumBits() != bv2.getNumBits()) {
    throw ValueErrorException("BitVects must be same length");
  }
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
DoubleVect OffBitProjSimilarity(const T1 &bv1, const T2 &bv2) {
  if (bv1.getNumBits() != bv2.getNumBits()) {
    throw ValueErrorException("BitVects must be same length");
  }
  DoubleVect res(2, 0.0);
  double num = (bv1 | bv2).getNumOffBits();
  if (num) {
    res[0] = num / bv1.getNumOffBits();
    res[1] = num / bv2.getNumOffBits();
  }
  return res;
}

template <typename T1>
T1 *FoldFingerprint(const T1 &bv1, unsigned int factor) {
  if (factor <= 0 || factor >= bv1.getNumBits()) {
    throw ValueErrorException("invalid fold factor");
  }

  int initSize = bv1.getNumBits();
  int resSize = initSize / factor;
  auto *res = new T1(resSize);

  IntVect onBits;
  bv1.getOnBits(onBits);
  for (int &onBit : onBits) {
    int pos = onBit % resSize;
    res->setBit(pos);
  }
  return res;
}

template <typename T1>
std::string BitVectToText(const T1 &bv1) {
  std::string res(bv1.getNumBits(), '0');
  for (unsigned int i = 0; i < bv1.getNumBits(); i++) {
    if (bv1.getBit(i)) {
      res[i] = '1';
    }
  }
  return res;
}

const char bin2Hex[] = {'0', '1', '2', '3', '4', '5', '6', '7',
                        '8', '9', 'a', 'b', 'c', 'd', 'e', 'f'};
template <typename T1>
std::string BitVectToFPSText(const T1 &bv1) {
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
std::string BitVectToBinaryText(const T1 &bv1) {
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
void UpdateBitVectFromFPSText(T1 &bv1, const std::string &fps) {
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
      c = static_cast<unsigned short>(strtol(tptr, nullptr, 16));
    } catch (...) {
      std::ostringstream errout;
      errout << "Cannot convert FPS word: " << fps.substr(i, 2) << " to int";
      std::cerr << errout.str() << std::endl;
      throw ValueErrorException(errout.str());
    }
    for (unsigned int bit = 0; bit < 8 && bitIdx < bv1.getNumBits();
         ++bit, ++bitIdx) {
      if (c & (1 << bit)) {
        bv1.setBit(bitIdx);
      }
    }
  }
}

template <typename T1>
void UpdateBitVectFromBinaryText(T1 &bv1, const std::string &fps) {
  PRECONDITION(fps.length() * 8 >= bv1.getNumBits(), "bad FPS length");
  unsigned int bitIdx = 0;
  for (unsigned int i = 0; i < fps.size() && bitIdx < bv1.getNumBits(); i++) {
    unsigned short c = fps[i];
    for (unsigned int bit = 0; bit < 8 && bitIdx < bv1.getNumBits();
         ++bit, ++bitIdx) {
      if (c & (1 << bit)) {
        bv1.setBit(bitIdx);
      }
    }
  }
}

template RDKIT_DATASTRUCTS_EXPORT double TanimotoSimilarity(
    const SparseBitVect &bv1, const SparseBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT double TverskySimilarity(
    const SparseBitVect &bv1, const SparseBitVect &bv2, double a, double b);
template RDKIT_DATASTRUCTS_EXPORT double CosineSimilarity(
    const SparseBitVect &bv1, const SparseBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT double KulczynskiSimilarity(
    const SparseBitVect &bv1, const SparseBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT double DiceSimilarity(
    const SparseBitVect &bv1, const SparseBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT double SokalSimilarity(
    const SparseBitVect &bv1, const SparseBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT double McConnaugheySimilarity(
    const SparseBitVect &bv1, const SparseBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT double AsymmetricSimilarity(
    const SparseBitVect &bv1, const SparseBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT double BraunBlanquetSimilarity(
    const SparseBitVect &bv1, const SparseBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT double RusselSimilarity(
    const SparseBitVect &bv1, const SparseBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT double RogotGoldbergSimilarity(
    const SparseBitVect &bv1, const SparseBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT double OnBitSimilarity(
    const SparseBitVect &bv1, const SparseBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT int NumBitsInCommon(const SparseBitVect &bv1,
                                                      const SparseBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT double AllBitSimilarity(
    const SparseBitVect &bv1, const SparseBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT int NumOnBitsInCommon(
    const SparseBitVect &bv1, const SparseBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT IntVect
OnBitsInCommon(const SparseBitVect &bv1, const SparseBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT IntVect
OffBitsInCommon(const SparseBitVect &bv1, const SparseBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT DoubleVect
OnBitProjSimilarity(const SparseBitVect &bv1, const SparseBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT DoubleVect
OffBitProjSimilarity(const SparseBitVect &bv1, const SparseBitVect &bv2);

template RDKIT_DATASTRUCTS_EXPORT double TanimotoSimilarity(
    const ExplicitBitVect &bv1, const ExplicitBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT double TverskySimilarity(
    const ExplicitBitVect &bv1, const ExplicitBitVect &bv2, double a, double b);
template RDKIT_DATASTRUCTS_EXPORT double CosineSimilarity(
    const ExplicitBitVect &bv1, const ExplicitBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT double KulczynskiSimilarity(
    const ExplicitBitVect &bv1, const ExplicitBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT double DiceSimilarity(
    const ExplicitBitVect &bv1, const ExplicitBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT double SokalSimilarity(
    const ExplicitBitVect &bv1, const ExplicitBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT double McConnaugheySimilarity(
    const ExplicitBitVect &bv1, const ExplicitBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT double AsymmetricSimilarity(
    const ExplicitBitVect &bv1, const ExplicitBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT double BraunBlanquetSimilarity(
    const ExplicitBitVect &bv1, const ExplicitBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT double RusselSimilarity(
    const ExplicitBitVect &bv1, const ExplicitBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT double RogotGoldbergSimilarity(
    const ExplicitBitVect &bv1, const ExplicitBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT double OnBitSimilarity(
    const ExplicitBitVect &bv1, const ExplicitBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT int NumBitsInCommon(
    const ExplicitBitVect &bv1, const ExplicitBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT double AllBitSimilarity(
    const ExplicitBitVect &bv1, const ExplicitBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT IntVect
OnBitsInCommon(const ExplicitBitVect &bv1, const ExplicitBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT IntVect
OffBitsInCommon(const ExplicitBitVect &bv1, const ExplicitBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT DoubleVect
OnBitProjSimilarity(const ExplicitBitVect &bv1, const ExplicitBitVect &bv2);
template RDKIT_DATASTRUCTS_EXPORT DoubleVect
OffBitProjSimilarity(const ExplicitBitVect &bv1, const ExplicitBitVect &bv2);

template RDKIT_DATASTRUCTS_EXPORT SparseBitVect *FoldFingerprint(
    const SparseBitVect &, unsigned int);
template RDKIT_DATASTRUCTS_EXPORT ExplicitBitVect *FoldFingerprint(
    const ExplicitBitVect &, unsigned int);

template RDKIT_DATASTRUCTS_EXPORT std::string BitVectToText(
    const SparseBitVect &);
template RDKIT_DATASTRUCTS_EXPORT std::string BitVectToText(
    const ExplicitBitVect &);

template RDKIT_DATASTRUCTS_EXPORT std::string BitVectToFPSText(
    const SparseBitVect &);
template RDKIT_DATASTRUCTS_EXPORT std::string BitVectToFPSText(
    const ExplicitBitVect &);
template RDKIT_DATASTRUCTS_EXPORT void UpdateBitVectFromFPSText(
    SparseBitVect &, const std::string &);
template RDKIT_DATASTRUCTS_EXPORT void UpdateBitVectFromFPSText(
    ExplicitBitVect &, const std::string &);

template RDKIT_DATASTRUCTS_EXPORT std::string BitVectToBinaryText(
    const SparseBitVect &);
template RDKIT_DATASTRUCTS_EXPORT std::string BitVectToBinaryText(
    const ExplicitBitVect &);
template RDKIT_DATASTRUCTS_EXPORT void UpdateBitVectFromBinaryText(
    SparseBitVect &, const std::string &);
template RDKIT_DATASTRUCTS_EXPORT void UpdateBitVectFromBinaryText(
    ExplicitBitVect &, const std::string &);

// from here:
// http://stackoverflow.com/questions/3849337/msvc-equivalent-to-builtin-popcount
// but corrected to get the ifdef right
#ifdef _MSC_VER
#include <intrin.h>
#ifdef _WIN64
#define BUILTIN_POPCOUNT_INSTR __popcnt64
using BUILTIN_POPCOUNT_TYPE = boost::uint64_t;
#else
#define BUILTIN_POPCOUNT_INSTR __popcnt
using BUILTIN_POPCOUNT_TYPE = std::uint32_t;
#endif
#else
#define BUILTIN_POPCOUNT_INSTR __builtin_popcountll
using BUILTIN_POPCOUNT_TYPE = boost::uint64_t;
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

// ============================================================
// SIMD-accelerated bitmap operations (x86/x64)
//
// Three-tier runtime dispatch: AVX2 (256-bit) > SSSE3 (128-bit) > scalar.
// Uses the Mula vpshufb nibble-lookup popcount algorithm:
//   - Split each byte into two nibbles
//   - Use pshufb as a parallel 4-bit LUT for popcount(nibble)
//   - Accumulate per-byte counts, then horizontal-sum via psadbw
//
// GCC/Clang: __attribute__((target("avx2"))) compiles each function with
//   the correct ISA extension regardless of global -march flags.
// MSVC: intrinsics compile unconditionally; runtime CPUID gates execution.
// ============================================================
#ifdef RDK_X86_SIMD
namespace {

#ifdef _MSC_VER
#define RDK_AVX2_TARGET
#define RDK_SSSE3_TARGET
#else
#define RDK_AVX2_TARGET __attribute__((target("avx2")))
#define RDK_SSSE3_TARGET __attribute__((target("ssse3")))
#endif

// --- Runtime CPU Feature Detection ---

static bool detect_avx2() {
#ifdef _MSC_VER
  int info[4];
  __cpuid(info, 1);
  bool osxsave = (info[2] >> 27) & 1;
  if (!osxsave) return false;
  unsigned long long xcr0 = _xgetbv(0);
  if ((xcr0 & 0x6) != 0x6) return false;
  __cpuidex(info, 7, 0);
  return (info[1] >> 5) & 1;
#elif defined(__GNUC__) || defined(__clang__)
  __builtin_cpu_init();
  return __builtin_cpu_supports("avx2");
#else
  return false;
#endif
}

static bool detect_ssse3() {
#ifdef _MSC_VER
  int info[4];
  __cpuid(info, 1);
  return (info[2] >> 9) & 1;
#elif defined(__GNUC__) || defined(__clang__)
  __builtin_cpu_init();
  return __builtin_cpu_supports("ssse3");
#else
  return false;
#endif
}

static const bool g_has_avx2 = detect_avx2();
static const bool g_has_ssse3 = detect_ssse3();

// Max iterations before byte accumulator overflows (255 / 8 = 31)
static constexpr unsigned int kFlushInterval = 31;

// --- AVX2 helpers (256-bit) ---

RDK_AVX2_TARGET
static inline __m256i popcnt_bytes_avx2(__m256i v) {
  const __m256i lookup = _mm256_setr_epi8(
      0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
      0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4);
  const __m256i low_mask = _mm256_set1_epi8(0x0f);
  __m256i lo = _mm256_and_si256(v, low_mask);
  __m256i hi = _mm256_and_si256(_mm256_srli_epi16(v, 4), low_mask);
  return _mm256_add_epi8(_mm256_shuffle_epi8(lookup, lo),
                         _mm256_shuffle_epi8(lookup, hi));
}

RDK_AVX2_TARGET
static inline unsigned int hsum_bytes_avx2(__m256i acc) {
  __m256i sad = _mm256_sad_epu8(acc, _mm256_setzero_si256());
  __m128i lo = _mm256_castsi256_si128(sad);
  __m128i hi = _mm256_extracti128_si256(sad, 1);
  __m128i sum = _mm_add_epi64(lo, hi);
  sum = _mm_add_epi32(sum, _mm_srli_si128(sum, 8));
  return (unsigned int)_mm_cvtsi128_si32(sum);
}

// --- SSSE3 helpers (128-bit) ---

RDK_SSSE3_TARGET
static inline __m128i popcnt_bytes_ssse3(__m128i v) {
  const __m128i lookup =
      _mm_setr_epi8(0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4);
  const __m128i low_mask = _mm_set1_epi8(0x0f);
  __m128i lo = _mm_and_si128(v, low_mask);
  __m128i hi = _mm_and_si128(_mm_srli_epi16(v, 4), low_mask);
  return _mm_add_epi8(_mm_shuffle_epi8(lookup, lo),
                      _mm_shuffle_epi8(lookup, hi));
}

RDK_SSSE3_TARGET
static inline unsigned int hsum_bytes_ssse3(__m128i acc) {
  __m128i sad = _mm_sad_epu8(acc, _mm_setzero_si128());
  __m128i sum = _mm_add_epi32(sad, _mm_srli_si128(sad, 8));
  return (unsigned int)_mm_cvtsi128_si32(sum);
}

// === AVX2 bitmap functions ===

RDK_AVX2_TARGET
static unsigned int popcount_avx2(const unsigned char *fp, unsigned int nBytes) {
  unsigned int count = 0, i = 0, iters = 0;
  const unsigned int simd_end = (nBytes / 32) * 32;
  __m256i acc = _mm256_setzero_si256();
  for (; i < simd_end; i += 32) {
    acc = _mm256_add_epi8(
        acc, popcnt_bytes_avx2(_mm256_loadu_si256((const __m256i *)(fp + i))));
    if (++iters == kFlushInterval) {
      count += hsum_bytes_avx2(acc);
      acc = _mm256_setzero_si256();
      iters = 0;
    }
  }
  count += hsum_bytes_avx2(acc);
  for (; i < nBytes; i++) count += byte_popcounts[fp[i]];
  return count;
}

RDK_AVX2_TARGET
static unsigned int num_bits_in_common_avx2(const unsigned char *afp,
                                            const unsigned char *bfp,
                                            unsigned int nBytes) {
  unsigned int count = 0, i = 0, iters = 0;
  const unsigned int simd_end = (nBytes / 32) * 32;
  __m256i acc = _mm256_setzero_si256();
  for (; i < simd_end; i += 32) {
    __m256i a = _mm256_loadu_si256((const __m256i *)(afp + i));
    __m256i b = _mm256_loadu_si256((const __m256i *)(bfp + i));
    acc = _mm256_add_epi8(acc, popcnt_bytes_avx2(_mm256_and_si256(a, b)));
    if (++iters == kFlushInterval) {
      count += hsum_bytes_avx2(acc);
      acc = _mm256_setzero_si256();
      iters = 0;
    }
  }
  count += hsum_bytes_avx2(acc);
  for (; i < nBytes; i++) count += byte_popcounts[afp[i] & bfp[i]];
  return count;
}

RDK_AVX2_TARGET
static double tanimoto_avx2(const unsigned char *afp, const unsigned char *bfp,
                            unsigned int nBytes) {
  unsigned int union_pop = 0, inter_pop = 0, i = 0, iters = 0;
  const unsigned int simd_end = (nBytes / 32) * 32;
  __m256i u_acc = _mm256_setzero_si256();
  __m256i i_acc = _mm256_setzero_si256();
  for (; i < simd_end; i += 32) {
    __m256i a = _mm256_loadu_si256((const __m256i *)(afp + i));
    __m256i b = _mm256_loadu_si256((const __m256i *)(bfp + i));
    u_acc = _mm256_add_epi8(u_acc, popcnt_bytes_avx2(_mm256_or_si256(a, b)));
    i_acc = _mm256_add_epi8(i_acc, popcnt_bytes_avx2(_mm256_and_si256(a, b)));
    if (++iters == kFlushInterval) {
      union_pop += hsum_bytes_avx2(u_acc);
      inter_pop += hsum_bytes_avx2(i_acc);
      u_acc = _mm256_setzero_si256();
      i_acc = _mm256_setzero_si256();
      iters = 0;
    }
  }
  union_pop += hsum_bytes_avx2(u_acc);
  inter_pop += hsum_bytes_avx2(i_acc);
  for (; i < nBytes; i++) {
    union_pop += byte_popcounts[afp[i] | bfp[i]];
    inter_pop += byte_popcounts[afp[i] & bfp[i]];
  }
  if (union_pop == 0) return 0.0;
  return (inter_pop + 0.0) / union_pop;
}

RDK_AVX2_TARGET
static double dice_avx2(const unsigned char *afp, const unsigned char *bfp,
                        unsigned int nBytes) {
  unsigned int a_pop = 0, b_pop = 0, inter_pop = 0, i = 0, iters = 0;
  const unsigned int simd_end = (nBytes / 32) * 32;
  __m256i a_acc = _mm256_setzero_si256();
  __m256i b_acc = _mm256_setzero_si256();
  __m256i i_acc = _mm256_setzero_si256();
  for (; i < simd_end; i += 32) {
    __m256i a = _mm256_loadu_si256((const __m256i *)(afp + i));
    __m256i b = _mm256_loadu_si256((const __m256i *)(bfp + i));
    a_acc = _mm256_add_epi8(a_acc, popcnt_bytes_avx2(a));
    b_acc = _mm256_add_epi8(b_acc, popcnt_bytes_avx2(b));
    i_acc = _mm256_add_epi8(i_acc, popcnt_bytes_avx2(_mm256_and_si256(a, b)));
    if (++iters == kFlushInterval) {
      a_pop += hsum_bytes_avx2(a_acc);
      b_pop += hsum_bytes_avx2(b_acc);
      inter_pop += hsum_bytes_avx2(i_acc);
      a_acc = _mm256_setzero_si256();
      b_acc = _mm256_setzero_si256();
      i_acc = _mm256_setzero_si256();
      iters = 0;
    }
  }
  a_pop += hsum_bytes_avx2(a_acc);
  b_pop += hsum_bytes_avx2(b_acc);
  inter_pop += hsum_bytes_avx2(i_acc);
  for (; i < nBytes; i++) {
    a_pop += byte_popcounts[afp[i]];
    b_pop += byte_popcounts[bfp[i]];
    inter_pop += byte_popcounts[afp[i] & bfp[i]];
  }
  if (a_pop + b_pop == 0) return 0.0;
  return (2.0 * inter_pop) / (a_pop + b_pop);
}

RDK_AVX2_TARGET
static double tversky_avx2(const unsigned char *afp, const unsigned char *bfp,
                           unsigned int nBytes, double ca, double cb) {
  unsigned int a_pop = 0, b_pop = 0, inter_pop = 0, i = 0, iters = 0;
  const unsigned int simd_end = (nBytes / 32) * 32;
  __m256i a_acc = _mm256_setzero_si256();
  __m256i b_acc = _mm256_setzero_si256();
  __m256i i_acc = _mm256_setzero_si256();
  for (; i < simd_end; i += 32) {
    __m256i a = _mm256_loadu_si256((const __m256i *)(afp + i));
    __m256i b = _mm256_loadu_si256((const __m256i *)(bfp + i));
    a_acc = _mm256_add_epi8(a_acc, popcnt_bytes_avx2(a));
    b_acc = _mm256_add_epi8(b_acc, popcnt_bytes_avx2(b));
    i_acc = _mm256_add_epi8(i_acc, popcnt_bytes_avx2(_mm256_and_si256(a, b)));
    if (++iters == kFlushInterval) {
      a_pop += hsum_bytes_avx2(a_acc);
      b_pop += hsum_bytes_avx2(b_acc);
      inter_pop += hsum_bytes_avx2(i_acc);
      a_acc = _mm256_setzero_si256();
      b_acc = _mm256_setzero_si256();
      i_acc = _mm256_setzero_si256();
      iters = 0;
    }
  }
  a_pop += hsum_bytes_avx2(a_acc);
  b_pop += hsum_bytes_avx2(b_acc);
  inter_pop += hsum_bytes_avx2(i_acc);
  for (; i < nBytes; i++) {
    inter_pop += byte_popcounts[afp[i] & bfp[i]];
    a_pop += byte_popcounts[afp[i]];
    b_pop += byte_popcounts[bfp[i]];
  }
  double denom = ca * a_pop + cb * b_pop + (1 - ca - cb) * inter_pop;
  if (denom == 0.0) return 0.0;
  return inter_pop / denom;
}

RDK_AVX2_TARGET
static bool all_probe_bits_match_avx2(const unsigned char *probe,
                                      const unsigned char *ref,
                                      unsigned int nBytes) {
  unsigned int i = 0;
  const unsigned int simd_end = (nBytes / 32) * 32;
  for (; i < simd_end; i += 32) {
    __m256i p = _mm256_loadu_si256((const __m256i *)(probe + i));
    __m256i r = _mm256_loadu_si256((const __m256i *)(ref + i));
    __m256i unmatched = _mm256_andnot_si256(r, p);  // probe & ~ref
    if (!_mm256_testz_si256(unmatched, unmatched)) return false;
  }
  for (; i < nBytes; i++) {
    if (probe[i] & ~ref[i]) return false;
  }
  return true;
}

// === SSSE3 bitmap functions (128-bit, same algorithm) ===

RDK_SSSE3_TARGET
static unsigned int popcount_ssse3(const unsigned char *fp,
                                   unsigned int nBytes) {
  unsigned int count = 0, i = 0, iters = 0;
  const unsigned int simd_end = (nBytes / 16) * 16;
  __m128i acc = _mm_setzero_si128();
  for (; i < simd_end; i += 16) {
    acc = _mm_add_epi8(acc,
                       popcnt_bytes_ssse3(_mm_loadu_si128((const __m128i *)(fp + i))));
    if (++iters == kFlushInterval) {
      count += hsum_bytes_ssse3(acc);
      acc = _mm_setzero_si128();
      iters = 0;
    }
  }
  count += hsum_bytes_ssse3(acc);
  for (; i < nBytes; i++) count += byte_popcounts[fp[i]];
  return count;
}

RDK_SSSE3_TARGET
static unsigned int num_bits_in_common_ssse3(const unsigned char *afp,
                                             const unsigned char *bfp,
                                             unsigned int nBytes) {
  unsigned int count = 0, i = 0, iters = 0;
  const unsigned int simd_end = (nBytes / 16) * 16;
  __m128i acc = _mm_setzero_si128();
  for (; i < simd_end; i += 16) {
    __m128i a = _mm_loadu_si128((const __m128i *)(afp + i));
    __m128i b = _mm_loadu_si128((const __m128i *)(bfp + i));
    acc = _mm_add_epi8(acc, popcnt_bytes_ssse3(_mm_and_si128(a, b)));
    if (++iters == kFlushInterval) {
      count += hsum_bytes_ssse3(acc);
      acc = _mm_setzero_si128();
      iters = 0;
    }
  }
  count += hsum_bytes_ssse3(acc);
  for (; i < nBytes; i++) count += byte_popcounts[afp[i] & bfp[i]];
  return count;
}

RDK_SSSE3_TARGET
static double tanimoto_ssse3(const unsigned char *afp, const unsigned char *bfp,
                             unsigned int nBytes) {
  unsigned int union_pop = 0, inter_pop = 0, i = 0, iters = 0;
  const unsigned int simd_end = (nBytes / 16) * 16;
  __m128i u_acc = _mm_setzero_si128();
  __m128i i_acc = _mm_setzero_si128();
  for (; i < simd_end; i += 16) {
    __m128i a = _mm_loadu_si128((const __m128i *)(afp + i));
    __m128i b = _mm_loadu_si128((const __m128i *)(bfp + i));
    u_acc = _mm_add_epi8(u_acc, popcnt_bytes_ssse3(_mm_or_si128(a, b)));
    i_acc = _mm_add_epi8(i_acc, popcnt_bytes_ssse3(_mm_and_si128(a, b)));
    if (++iters == kFlushInterval) {
      union_pop += hsum_bytes_ssse3(u_acc);
      inter_pop += hsum_bytes_ssse3(i_acc);
      u_acc = _mm_setzero_si128();
      i_acc = _mm_setzero_si128();
      iters = 0;
    }
  }
  union_pop += hsum_bytes_ssse3(u_acc);
  inter_pop += hsum_bytes_ssse3(i_acc);
  for (; i < nBytes; i++) {
    union_pop += byte_popcounts[afp[i] | bfp[i]];
    inter_pop += byte_popcounts[afp[i] & bfp[i]];
  }
  if (union_pop == 0) return 0.0;
  return (inter_pop + 0.0) / union_pop;
}

RDK_SSSE3_TARGET
static double dice_ssse3(const unsigned char *afp, const unsigned char *bfp,
                         unsigned int nBytes) {
  unsigned int a_pop = 0, b_pop = 0, inter_pop = 0, i = 0, iters = 0;
  const unsigned int simd_end = (nBytes / 16) * 16;
  __m128i a_acc = _mm_setzero_si128();
  __m128i b_acc = _mm_setzero_si128();
  __m128i i_acc = _mm_setzero_si128();
  for (; i < simd_end; i += 16) {
    __m128i a = _mm_loadu_si128((const __m128i *)(afp + i));
    __m128i b = _mm_loadu_si128((const __m128i *)(bfp + i));
    a_acc = _mm_add_epi8(a_acc, popcnt_bytes_ssse3(a));
    b_acc = _mm_add_epi8(b_acc, popcnt_bytes_ssse3(b));
    i_acc = _mm_add_epi8(i_acc, popcnt_bytes_ssse3(_mm_and_si128(a, b)));
    if (++iters == kFlushInterval) {
      a_pop += hsum_bytes_ssse3(a_acc);
      b_pop += hsum_bytes_ssse3(b_acc);
      inter_pop += hsum_bytes_ssse3(i_acc);
      a_acc = _mm_setzero_si128();
      b_acc = _mm_setzero_si128();
      i_acc = _mm_setzero_si128();
      iters = 0;
    }
  }
  a_pop += hsum_bytes_ssse3(a_acc);
  b_pop += hsum_bytes_ssse3(b_acc);
  inter_pop += hsum_bytes_ssse3(i_acc);
  for (; i < nBytes; i++) {
    a_pop += byte_popcounts[afp[i]];
    b_pop += byte_popcounts[bfp[i]];
    inter_pop += byte_popcounts[afp[i] & bfp[i]];
  }
  if (a_pop + b_pop == 0) return 0.0;
  return (2.0 * inter_pop) / (a_pop + b_pop);
}

RDK_SSSE3_TARGET
static double tversky_ssse3(const unsigned char *afp, const unsigned char *bfp,
                            unsigned int nBytes, double ca, double cb) {
  unsigned int a_pop = 0, b_pop = 0, inter_pop = 0, i = 0, iters = 0;
  const unsigned int simd_end = (nBytes / 16) * 16;
  __m128i a_acc = _mm_setzero_si128();
  __m128i b_acc = _mm_setzero_si128();
  __m128i i_acc = _mm_setzero_si128();
  for (; i < simd_end; i += 16) {
    __m128i a = _mm_loadu_si128((const __m128i *)(afp + i));
    __m128i b = _mm_loadu_si128((const __m128i *)(bfp + i));
    a_acc = _mm_add_epi8(a_acc, popcnt_bytes_ssse3(a));
    b_acc = _mm_add_epi8(b_acc, popcnt_bytes_ssse3(b));
    i_acc = _mm_add_epi8(i_acc, popcnt_bytes_ssse3(_mm_and_si128(a, b)));
    if (++iters == kFlushInterval) {
      a_pop += hsum_bytes_ssse3(a_acc);
      b_pop += hsum_bytes_ssse3(b_acc);
      inter_pop += hsum_bytes_ssse3(i_acc);
      a_acc = _mm_setzero_si128();
      b_acc = _mm_setzero_si128();
      i_acc = _mm_setzero_si128();
      iters = 0;
    }
  }
  a_pop += hsum_bytes_ssse3(a_acc);
  b_pop += hsum_bytes_ssse3(b_acc);
  inter_pop += hsum_bytes_ssse3(i_acc);
  for (; i < nBytes; i++) {
    inter_pop += byte_popcounts[afp[i] & bfp[i]];
    a_pop += byte_popcounts[afp[i]];
    b_pop += byte_popcounts[bfp[i]];
  }
  double denom = ca * a_pop + cb * b_pop + (1 - ca - cb) * inter_pop;
  if (denom == 0.0) return 0.0;
  return inter_pop / denom;
}

RDK_SSSE3_TARGET
static bool all_probe_bits_match_ssse3(const unsigned char *probe,
                                       const unsigned char *ref,
                                       unsigned int nBytes) {
  unsigned int i = 0;
  const unsigned int simd_end = (nBytes / 16) * 16;
  const __m128i zero = _mm_setzero_si128();
  for (; i < simd_end; i += 16) {
    __m128i p = _mm_loadu_si128((const __m128i *)(probe + i));
    __m128i r = _mm_loadu_si128((const __m128i *)(ref + i));
    __m128i unmatched = _mm_andnot_si128(r, p);
    if (_mm_movemask_epi8(_mm_cmpeq_epi8(unmatched, zero)) != 0xFFFF)
      return false;
  }
  for (; i < nBytes; i++) {
    if (probe[i] & ~ref[i]) return false;
  }
  return true;
}

#undef RDK_AVX2_TARGET
#undef RDK_SSSE3_TARGET
}  // anonymous namespace
#endif  // RDK_X86_SIMD

unsigned int CalcBitmapPopcount(const unsigned char *afp, unsigned int nBytes) {
  PRECONDITION(afp, "no afp");
#ifdef RDK_X86_SIMD
  if (g_has_avx2) return popcount_avx2(afp, nBytes);
  if (g_has_ssse3) return popcount_ssse3(afp, nBytes);
#endif
  unsigned int popcount = 0;
#ifndef RDK_OPTIMIZE_POPCNT
  for (unsigned int i = 0; i < nBytes; i++) {
    popcount += byte_popcounts[afp[i]];
  }
#else
  unsigned int eidx = nBytes / sizeof(BUILTIN_POPCOUNT_TYPE);
  for (unsigned int i = 0; i < eidx; ++i) {
    popcount += static_cast<unsigned int>(
        BUILTIN_POPCOUNT_INSTR(((BUILTIN_POPCOUNT_TYPE *)afp)[i]));
  }
  for (unsigned int i = eidx * sizeof(BUILTIN_POPCOUNT_TYPE); i < nBytes; ++i) {
    popcount += byte_popcounts[afp[i]];
  }
#endif
  return popcount;
}

unsigned int CalcBitmapNumBitsInCommon(const unsigned char *afp,
                                       const unsigned char *bfp,
                                       unsigned int nBytes) {
  PRECONDITION(afp, "no afp");
  PRECONDITION(bfp, "no bfp");
#ifdef RDK_X86_SIMD
  if (g_has_avx2) return num_bits_in_common_avx2(afp, bfp, nBytes);
  if (g_has_ssse3) return num_bits_in_common_ssse3(afp, bfp, nBytes);
#endif
  unsigned int intersect_popcount = 0;
#ifndef RDK_OPTIMIZE_POPCNT
  for (unsigned int i = 0; i < nBytes; i++) {
    intersect_popcount += byte_popcounts[afp[i] & bfp[i]];
  }
#else
  BUILTIN_POPCOUNT_TYPE eidx = nBytes / sizeof(BUILTIN_POPCOUNT_TYPE);
  for (BUILTIN_POPCOUNT_TYPE i = 0; i < eidx; ++i) {
    intersect_popcount += static_cast<unsigned int>(BUILTIN_POPCOUNT_INSTR(
        ((BUILTIN_POPCOUNT_TYPE *)afp)[i] & ((BUILTIN_POPCOUNT_TYPE *)bfp)[i]));
  }
  for (BUILTIN_POPCOUNT_TYPE i = eidx * sizeof(BUILTIN_POPCOUNT_TYPE);
       i < nBytes; ++i) {
    intersect_popcount += byte_popcounts[afp[i] & bfp[i]];
  }
#endif
  return intersect_popcount;
}

double CalcBitmapTanimoto(const unsigned char *afp, const unsigned char *bfp,
                          unsigned int nBytes) {
  PRECONDITION(afp, "no afp");
  PRECONDITION(bfp, "no bfp");
#ifdef RDK_X86_SIMD
  if (g_has_avx2) return tanimoto_avx2(afp, bfp, nBytes);
  if (g_has_ssse3) return tanimoto_ssse3(afp, bfp, nBytes);
#endif
  unsigned int union_popcount = 0, intersect_popcount = 0;
#ifndef RDK_OPTIMIZE_POPCNT
  for (unsigned int i = 0; i < nBytes; i++) {
    union_popcount += byte_popcounts[afp[i] | bfp[i]];
    intersect_popcount += byte_popcounts[afp[i] & bfp[i]];
  }
#else
  BUILTIN_POPCOUNT_TYPE eidx = nBytes / sizeof(BUILTIN_POPCOUNT_TYPE);
  for (BUILTIN_POPCOUNT_TYPE i = 0; i < eidx; ++i) {
    union_popcount += static_cast<unsigned int>(BUILTIN_POPCOUNT_INSTR(
        ((BUILTIN_POPCOUNT_TYPE *)afp)[i] | ((BUILTIN_POPCOUNT_TYPE *)bfp)[i]));
    intersect_popcount += static_cast<unsigned int>(BUILTIN_POPCOUNT_INSTR(
        ((BUILTIN_POPCOUNT_TYPE *)afp)[i] & ((BUILTIN_POPCOUNT_TYPE *)bfp)[i]));
  }
  for (BUILTIN_POPCOUNT_TYPE i = eidx * sizeof(BUILTIN_POPCOUNT_TYPE);
       i < nBytes; ++i) {
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

double CalcBitmapDice(const unsigned char *afp, const unsigned char *bfp,
                      unsigned int nBytes) {
  PRECONDITION(afp, "no afp");
  PRECONDITION(bfp, "no bfp");
#ifdef RDK_X86_SIMD
  if (g_has_avx2) return dice_avx2(afp, bfp, nBytes);
  if (g_has_ssse3) return dice_ssse3(afp, bfp, nBytes);
#endif
  unsigned int intersect_popcount = 0, a_popcount = 0, b_popcount = 0;

#ifndef RDK_OPTIMIZE_POPCNT
  for (unsigned int i = 0; i < nBytes; i++) {
    a_popcount += byte_popcounts[afp[i]];
    b_popcount += byte_popcounts[bfp[i]];
    intersect_popcount += byte_popcounts[afp[i] & bfp[i]];
  }
#else
  BUILTIN_POPCOUNT_TYPE eidx = nBytes / sizeof(BUILTIN_POPCOUNT_TYPE);
  for (BUILTIN_POPCOUNT_TYPE i = 0; i < eidx; ++i) {
    a_popcount += static_cast<unsigned int>(
        BUILTIN_POPCOUNT_INSTR(((BUILTIN_POPCOUNT_TYPE *)afp)[i]));
    b_popcount += static_cast<unsigned int>(
        BUILTIN_POPCOUNT_INSTR(((BUILTIN_POPCOUNT_TYPE *)bfp)[i]));
    intersect_popcount += static_cast<unsigned int>(BUILTIN_POPCOUNT_INSTR(
        ((BUILTIN_POPCOUNT_TYPE *)afp)[i] & ((BUILTIN_POPCOUNT_TYPE *)bfp)[i]));
  }
  for (BUILTIN_POPCOUNT_TYPE i = eidx * sizeof(BUILTIN_POPCOUNT_TYPE);
       i < nBytes; ++i) {
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

double CalcBitmapTversky(const unsigned char *afp, const unsigned char *bfp,
                         unsigned int nBytes, double ca, double cb) {
  PRECONDITION(afp, "no afp");
  PRECONDITION(bfp, "no bfp");
#ifdef RDK_X86_SIMD
  if (g_has_avx2) return tversky_avx2(afp, bfp, nBytes, ca, cb);
  if (g_has_ssse3) return tversky_ssse3(afp, bfp, nBytes, ca, cb);
#endif
  unsigned int intersect_popcount = 0, acount = 0, bcount = 0;

#ifndef RDK_OPTIMIZE_POPCNT
  for (unsigned int i = 0; i < nBytes; i++) {
    intersect_popcount += byte_popcounts[afp[i] & bfp[i]];
    acount += byte_popcounts[afp[i]];
    bcount += byte_popcounts[bfp[i]];
  }
#else
  BUILTIN_POPCOUNT_TYPE eidx = nBytes / sizeof(BUILTIN_POPCOUNT_TYPE);
  for (BUILTIN_POPCOUNT_TYPE i = 0; i < eidx; ++i) {
    intersect_popcount += static_cast<unsigned int>(BUILTIN_POPCOUNT_INSTR(
        ((BUILTIN_POPCOUNT_TYPE *)afp)[i] & ((BUILTIN_POPCOUNT_TYPE *)bfp)[i]));
    acount += static_cast<unsigned int>(
        BUILTIN_POPCOUNT_INSTR(((BUILTIN_POPCOUNT_TYPE *)afp)[i]));
    bcount += static_cast<unsigned int>(
        BUILTIN_POPCOUNT_INSTR(((BUILTIN_POPCOUNT_TYPE *)bfp)[i]));
  }
  for (BUILTIN_POPCOUNT_TYPE i = eidx * sizeof(BUILTIN_POPCOUNT_TYPE);
       i < nBytes; ++i) {
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

bool CalcBitmapAllProbeBitsMatch(const unsigned char *probe,
                                 const unsigned char *ref,
                                 unsigned int nBytes) {
  PRECONDITION(probe, "no probe");
  PRECONDITION(ref, "no ref");
#ifdef RDK_X86_SIMD
  if (g_has_avx2) return all_probe_bits_match_avx2(probe, ref, nBytes);
  if (g_has_ssse3) return all_probe_bits_match_ssse3(probe, ref, nBytes);
#endif

#ifndef RDK_OPTIMIZE_POPCNT
  for (unsigned int i = 0; i < nBytes; i++) {
    if (byte_popcounts[probe[i] & ref[i]] != byte_popcounts[probe[i]]) {
      return false;
    }
  }
#else
  unsigned int eidx = nBytes / sizeof(BUILTIN_POPCOUNT_TYPE);
  for (unsigned int i = 0; i < eidx; ++i) {
    if (BUILTIN_POPCOUNT_INSTR(((BUILTIN_POPCOUNT_TYPE *)probe)[i] &
                               ((BUILTIN_POPCOUNT_TYPE *)ref)[i]) !=
        BUILTIN_POPCOUNT_INSTR(((BUILTIN_POPCOUNT_TYPE *)probe)[i])) {
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
