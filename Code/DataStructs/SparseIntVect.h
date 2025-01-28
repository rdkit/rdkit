//
//  Copyright (C) 2007-2024 Greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef __RD_SPARSE_INT_VECT_20070921__
#define __RD_SPARSE_INT_VECT_20070921__

#include <map>
#include <string>
#include <RDGeneral/Invariant.h>
#include <sstream>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/StreamOps.h>
#include <cstdint>
#include <limits>

const int ci_SPARSEINTVECT_VERSION =
    0x0001;  //!< version number to use in pickles
namespace RDKit {
//! a class for efficiently storing sparse vectors of ints
template <typename IndexType>
class SparseIntVect {
 public:
  typedef std::map<IndexType, int> StorageType;

  SparseIntVect() : d_length(0) {}

  //! initialize with a particular length
  SparseIntVect(IndexType length) : d_length(length) {}

  //! Copy constructor
  SparseIntVect(const SparseIntVect<IndexType> &other) {
    d_length = other.d_length;
    d_data.clear();
    d_data.insert(other.d_data.begin(), other.d_data.end());
  }

  //! constructor from a pickle
  SparseIntVect(const std::string &pkl) {
    initFromText(pkl.c_str(), pkl.size());
  }
  //! constructor from a pickle
  SparseIntVect(const char *pkl, const unsigned int len) {
    initFromText(pkl, len);
  }

  SparseIntVect &operator=(const SparseIntVect<IndexType> &other) {
    if (this == &other) {
      return *this;
    }
    d_length = other.d_length;
    d_data.clear();
    d_data.insert(other.d_data.begin(), other.d_data.end());
    return *this;
  }

  //! destructor (doesn't need to do anything)
  ~SparseIntVect() = default;

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wtautological-compare"
#elif (defined(__GNUC__) || defined(__GNUG__)) && \
    (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 1))
#if (__GNUC__ > 4 || __GNUC_MINOR__ > 5)
#pragma GCC diagnostic push
#endif
#pragma GCC diagnostic ignored "-Wtype-limits"
#endif
  //! return the value at an index
  int getVal(IndexType idx) const {
    if (!checkIndex(idx)) {
      throw IndexErrorException(static_cast<int>(idx));
    }
    int res = 0;
    typename StorageType::const_iterator iter = d_data.find(idx);
    if (iter != d_data.end()) {
      res = iter->second;
    }
    return res;
  }

  //! set the value at an index
  void setVal(IndexType idx, int val) {
    if (!checkIndex(idx)) {
      throw IndexErrorException(static_cast<int>(idx));
    }
    if (val != 0) {
      d_data[idx] = val;
    } else {
      d_data.erase(idx);
    }
  }
#ifdef __clang__
#pragma clang diagnostic pop
#elif (defined(__GNUC__) || defined(__GNUG__)) && \
    (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 5))
#pragma GCC diagnostic pop
#endif
  //! support indexing using []
  int operator[](IndexType idx) const { return getVal(idx); }

  //! returns the length
  IndexType getLength() const { return d_length; }

  //! returns the sum of all the elements in the vect
  //! the doAbs argument toggles summing the absolute values of the elements
  int getTotalVal(bool doAbs = false) const {
    int res = 0;
    typename StorageType::const_iterator iter;
    for (iter = d_data.begin(); iter != d_data.end(); ++iter) {
      if (!doAbs) {
        res += iter->second;
      } else {
        res += abs(iter->second);
      }
    }
    return res;
  }
  //! returns the length
  unsigned int size() const { return getLength(); }

  //! returns our nonzero elements as a map(IndexType->int)
  const StorageType &getNonzeroElements() const { return d_data; }

  //! this is a "fuzzy" intesection, the final value
  //! of each element is equal to the minimum from
  //! the two vects.
  SparseIntVect<IndexType> &operator&=(const SparseIntVect<IndexType> &other) {
    if (other.d_length != d_length) {
      throw ValueErrorException("SparseIntVect size mismatch");
    }

    typename StorageType::iterator iter = d_data.begin();
    typename StorageType::const_iterator oIter = other.d_data.begin();
    while (iter != d_data.end()) {
      // we're relying on the fact that the maps are sorted:
      while (oIter != other.d_data.end() && oIter->first < iter->first) {
        ++oIter;
      }
      if (oIter != other.d_data.end() && oIter->first == iter->first) {
        // found it:
        if (oIter->second < iter->second) {
          iter->second = oIter->second;
        }
        ++oIter;
        ++iter;
      } else {
        // not there; our value is zero, which means
        // we should remove this value:
        typename StorageType::iterator tmpIter = iter;
        ++tmpIter;
        d_data.erase(iter);
        iter = tmpIter;
      }
    }
    return *this;
  }
  const SparseIntVect<IndexType> operator&(
      const SparseIntVect<IndexType> &other) const {
    SparseIntVect<IndexType> res(*this);
    return res &= other;
  }

  //! this is a "fuzzy" union, the final value
  //! of each element is equal to the maximum from
  //! the two vects.
  SparseIntVect<IndexType> &operator|=(const SparseIntVect<IndexType> &other) {
    if (other.d_length != d_length) {
      throw ValueErrorException("SparseIntVect size mismatch");
    }

    typename StorageType::iterator iter = d_data.begin();
    typename StorageType::const_iterator oIter = other.d_data.begin();
    while (iter != d_data.end()) {
      // we're relying on the fact that the maps are sorted:
      while (oIter != other.d_data.end() && oIter->first < iter->first) {
        d_data[oIter->first] = oIter->second;
        ++oIter;
      }
      if (oIter != other.d_data.end() && oIter->first == iter->first) {
        // found it:
        if (oIter->second > iter->second) {
          iter->second = oIter->second;
        }
        ++oIter;
      }
      ++iter;
    }
    // finish up the other vect:
    while (oIter != other.d_data.end()) {
      d_data[oIter->first] = oIter->second;
      ++oIter;
    }
    return *this;
  }
  const SparseIntVect<IndexType> operator|(
      const SparseIntVect<IndexType> &other) const {
    SparseIntVect<IndexType> res(*this);
    return res |= other;
  }

  SparseIntVect<IndexType> &operator+=(const SparseIntVect<IndexType> &other) {
    if (other.d_length != d_length) {
      throw ValueErrorException("SparseIntVect size mismatch");
    }
    typename StorageType::iterator iter = d_data.begin();
    typename StorageType::const_iterator oIter = other.d_data.begin();
    while (oIter != other.d_data.end()) {
      while (iter != d_data.end() && iter->first < oIter->first) {
        ++iter;
      }
      if (iter != d_data.end() && oIter->first == iter->first) {
        // found it:
        iter->second += oIter->second;
        if (!iter->second) {
          typename StorageType::iterator tIter = iter;
          ++tIter;
          d_data.erase(iter);
          iter = tIter;
        } else {
          ++iter;
        }
      } else {
        d_data[oIter->first] = oIter->second;
      }
      ++oIter;
    }
    return *this;
  }
  const SparseIntVect<IndexType> operator+(
      const SparseIntVect<IndexType> &other) const {
    SparseIntVect<IndexType> res(*this);
    return res += other;
  }

  SparseIntVect<IndexType> &operator-=(const SparseIntVect<IndexType> &other) {
    if (other.d_length != d_length) {
      throw ValueErrorException("SparseIntVect size mismatch");
    }
    typename StorageType::iterator iter = d_data.begin();
    typename StorageType::const_iterator oIter = other.d_data.begin();
    while (oIter != other.d_data.end()) {
      while (iter != d_data.end() && iter->first < oIter->first) {
        ++iter;
      }
      if (iter != d_data.end() && oIter->first == iter->first) {
        // found it:
        iter->second -= oIter->second;
        if (!iter->second) {
          typename StorageType::iterator tIter = iter;
          ++tIter;
          d_data.erase(iter);
          iter = tIter;
        } else {
          ++iter;
        }
      } else {
        d_data[oIter->first] = -oIter->second;
      }
      ++oIter;
    }
    return *this;
  }
  const SparseIntVect<IndexType> operator-(
      const SparseIntVect<IndexType> &other) const {
    SparseIntVect<IndexType> res(*this);
    return res -= other;
  }
  SparseIntVect<IndexType> &operator*=(int v) {
    typename StorageType::iterator iter = d_data.begin();
    while (iter != d_data.end()) {
      iter->second *= v;
      ++iter;
    }
    return *this;
  }
  SparseIntVect<IndexType> &operator*(int v) {
    SparseIntVect<IndexType> res(*this);
    return res *= v;
  }
  SparseIntVect<IndexType> &operator/=(int v) {
    typename StorageType::iterator iter = d_data.begin();
    while (iter != d_data.end()) {
      iter->second /= v;
      ++iter;
    }
    return *this;
  }
  SparseIntVect<IndexType> &operator/(int v) {
    SparseIntVect<IndexType> res(*this);
    return res /= v;
  }
  SparseIntVect<IndexType> &operator+=(int v) {
    typename StorageType::iterator iter = d_data.begin();
    while (iter != d_data.end()) {
      iter->second += v;
      ++iter;
    }
    return *this;
  }
  SparseIntVect<IndexType> &operator+(int v) {
    SparseIntVect<IndexType> res(*this);
    return res += v;
  }
  SparseIntVect<IndexType> &operator-=(int v) {
    typename StorageType::iterator iter = d_data.begin();
    while (iter != d_data.end()) {
      iter->second -= v;
      ++iter;
    }
    return *this;
  }
  SparseIntVect<IndexType> &operator-(int v) {
    SparseIntVect<IndexType> res(*this);
    return res -= v;
  }

  bool operator==(const SparseIntVect<IndexType> &v2) const {
    if (d_length != v2.d_length) {
      return false;
    }
    return d_data == v2.d_data;
  }
  bool operator!=(const SparseIntVect<IndexType> &v2) const {
    return !(*this == v2);
  }

  //! returns a binary string representation (pickle)
  std::string toString() const {
    std::stringstream ss(std::ios_base::binary | std::ios_base::out |
                         std::ios_base::in);
    std::uint32_t tInt;
    tInt = ci_SPARSEINTVECT_VERSION;
    streamWrite(ss, tInt);
    tInt = sizeof(IndexType);
    streamWrite(ss, tInt);
    streamWrite(ss, d_length);
    IndexType nEntries = d_data.size();
    streamWrite(ss, nEntries);

    typename StorageType::const_iterator iter = d_data.begin();
    while (iter != d_data.end()) {
      streamWrite(ss, iter->first);
      std::int32_t tInt = iter->second;
      streamWrite(ss, tInt);
      ++iter;
    }
    return ss.str();
  }

  void fromString(const std::string &txt) {
    initFromText(txt.c_str(), txt.length());
  }

 private:
  IndexType d_length;
  StorageType d_data;

  void initFromText(const char *pkl, const unsigned int len) {
    d_data.clear();
    std::stringstream ss(std::ios_base::binary | std::ios_base::out |
                         std::ios_base::in);
    ss.write(pkl, len);

    std::uint32_t vers;
    streamRead(ss, vers);
    if (vers == 0x0001) {
      std::uint32_t tInt;
      streamRead(ss, tInt);
      if (tInt > sizeof(IndexType)) {
        throw ValueErrorException(
            "IndexType cannot accommodate index size in SparseIntVect pickle");
      }
      switch (tInt) {
        case sizeof(char):
          readVals<unsigned char>(ss);
          break;
        case sizeof(std::int32_t):
          readVals<std::uint32_t>(ss);
          break;
        case sizeof(boost::int64_t):
          readVals<boost::uint64_t>(ss);
          break;
        default:
          throw ValueErrorException("unreadable format");
      }
    } else {
      throw ValueErrorException("bad version in SparseIntVect pickle");
    }
  }
  template <typename T>
  void readVals(std::stringstream &ss) {
    PRECONDITION(sizeof(T) <= sizeof(IndexType), "invalid size");
    T tVal;
    streamRead(ss, tVal);
    d_length = tVal;
    T nEntries;
    streamRead(ss, nEntries);
    for (T i = 0; i < nEntries; ++i) {
      streamRead(ss, tVal);
      std::int32_t val;
      streamRead(ss, val);
      d_data[tVal] = val;
    }
  }
  bool checkIndex(IndexType idx) const {
    if (idx < 0 || idx > d_length ||
        (idx == d_length && d_length < std::numeric_limits<IndexType>::max())) {
      return false;
    }
    return true;
  }
};

template <typename IndexType, typename SequenceType>
void updateFromSequence(SparseIntVect<IndexType> &vect,
                        const SequenceType &seq) {
  typename SequenceType::const_iterator seqIt;
  for (seqIt = seq.begin(); seqIt != seq.end(); ++seqIt) {
    // EFF: probably not the most efficient approach
    IndexType idx = *seqIt;
    vect.setVal(idx, vect.getVal(idx) + 1);
  }
}

namespace {
template <typename IndexType>
void calcVectParams(const SparseIntVect<IndexType> &v1,
                    const SparseIntVect<IndexType> &v2, double &v1Sum,
                    double &v2Sum, double &andSum) {
  if (v1.getLength() != v2.getLength()) {
    throw ValueErrorException("SparseIntVect size mismatch");
  }
  v1Sum = v2Sum = andSum = 0.0;
  // we're doing : (v1&v2).getTotalVal(), but w/o generating
  // the other vector:
  typename SparseIntVect<IndexType>::StorageType::const_iterator iter1, iter2;
  iter1 = v1.getNonzeroElements().begin();
  if (iter1 != v1.getNonzeroElements().end()) {
    v1Sum += abs(iter1->second);
  }
  iter2 = v2.getNonzeroElements().begin();
  if (iter2 != v2.getNonzeroElements().end()) {
    v2Sum += abs(iter2->second);
  }
  while (iter1 != v1.getNonzeroElements().end()) {
    while (iter2 != v2.getNonzeroElements().end() &&
           iter2->first < iter1->first) {
      ++iter2;
      if (iter2 != v2.getNonzeroElements().end()) {
        v2Sum += abs(iter2->second);
      }
    }
    if (iter2 != v2.getNonzeroElements().end()) {
      if (iter2->first == iter1->first) {
        if (abs(iter2->second) < abs(iter1->second)) {
          andSum += abs(iter2->second);
        } else {
          andSum += abs(iter1->second);
        }
        ++iter2;
        if (iter2 != v2.getNonzeroElements().end()) {
          v2Sum += abs(iter2->second);
        }
      }
      ++iter1;
      if (iter1 != v1.getNonzeroElements().end()) {
        v1Sum += abs(iter1->second);
      }
    } else {
      break;
    }
  }
  if (iter1 != v1.getNonzeroElements().end()) {
    ++iter1;
    while (iter1 != v1.getNonzeroElements().end()) {
      v1Sum += abs(iter1->second);
      ++iter1;
    }
  }
  if (iter2 != v2.getNonzeroElements().end()) {
    ++iter2;
    while (iter2 != v2.getNonzeroElements().end()) {
      v2Sum += abs(iter2->second);
      ++iter2;
    }
  }
}
}  // namespace

template <typename IndexType>
double DiceSimilarity(const SparseIntVect<IndexType> &v1,
                      const SparseIntVect<IndexType> &v2,
                      bool returnDistance = false, double bounds = 0.0) {
  if (v1.getLength() != v2.getLength()) {
    throw ValueErrorException("SparseIntVect size mismatch");
  }
  double v1Sum = 0.0;
  double v2Sum = 0.0;
  if (!returnDistance && bounds > 0.0) {
    v1Sum = v1.getTotalVal(true);
    v2Sum = v2.getTotalVal(true);
    double denom = v1Sum + v2Sum;
    if (fabs(denom) < 1e-6) {
      // no need to worry about returnDistance here
      return 0.0;
    }
    double minV = v1Sum < v2Sum ? v1Sum : v2Sum;
    if (2. * minV / denom < bounds) {
      return 0.0;
    }
    v1Sum = 0.0;
    v2Sum = 0.0;
  }

  double numer = 0.0;

  calcVectParams(v1, v2, v1Sum, v2Sum, numer);

  double denom = v1Sum + v2Sum;
  double sim;
  if (fabs(denom) < 1e-6) {
    sim = 0.0;
  } else {
    sim = 2. * numer / denom;
  }
  if (returnDistance) {
    sim = 1. - sim;
  }
  // std::cerr<<" "<<v1Sum<<" "<<v2Sum<<" " << numer << " " << sim <<std::endl;
  return sim;
}

template <typename IndexType>
double TverskySimilarity(const SparseIntVect<IndexType> &v1,
                         const SparseIntVect<IndexType> &v2, double a, double b,
                         bool returnDistance = false, double bounds = 0.0) {
  RDUNUSED_PARAM(bounds);
  if (v1.getLength() != v2.getLength()) {
    throw ValueErrorException("SparseIntVect size mismatch");
  }
  double v1Sum = 0.0;
  double v2Sum = 0.0;
  double andSum = 0.0;

  calcVectParams(v1, v2, v1Sum, v2Sum, andSum);

  double denom = a * v1Sum + b * v2Sum + (1 - a - b) * andSum;
  double sim;

  if (fabs(denom) < 1e-6) {
    sim = 0.0;
  } else {
    sim = andSum / denom;
  }
  if (returnDistance) {
    sim = 1. - sim;
  }
  // std::cerr<<" "<<v1Sum<<" "<<v2Sum<<" " << numer << " " << sim <<std::endl;
  return sim;
}

template <typename IndexType>
double TanimotoSimilarity(const SparseIntVect<IndexType> &v1,
                          const SparseIntVect<IndexType> &v2,
                          bool returnDistance = false, double bounds = 0.0) {
  return TverskySimilarity(v1, v2, 1.0, 1.0, returnDistance, bounds);
}
}  // namespace RDKit

#endif
