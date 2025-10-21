//
//  Copyright (C) 2004-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_VECTOR_H
#define RD_VECTOR_H

#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <memory>
#include <boost/random.hpp>

static constexpr double zero_tolerance = 1.e-16;

namespace RDNumeric {

//! A class to represent vectors of numbers.
template <class TYPE>
class Vector {
 public:
  typedef std::shared_ptr<TYPE[]> DATA_SPTR;

  //! Initialize with only a size.
  constexpr explicit Vector(unsigned int N) {
    d_size = N;
    TYPE *data = new TYPE[N];
    memset(static_cast<void *>(data), 0, d_size * sizeof(TYPE));
    d_data.reset(data);
  }

  //! Initialize with a size and default value.
  constexpr Vector(unsigned int N, TYPE val) {
    d_size = N;
    d_data.reset(new TYPE[N]);
    for (auto i = 0u; i < N; i++) {
      d_data[i] = val;
    }
  }

  //! Initialize from a smart pointer.
  /*!
    <b>NOTE:</b> the data is not copied in this case
  */
  constexpr Vector(unsigned int N, DATA_SPTR data) {
    d_size = N;
    d_data = data;
  }

  //! copy constructor
  /*! We make a copy of the other vector's data.
   */
  constexpr Vector(const Vector &other) {
    d_size = other.size();
    const TYPE *otherData = other.getData();
    d_data.reset(new TYPE[d_size]);

    memcpy(static_cast<void *>(d_data.get()),
           static_cast<const void *>(otherData), d_size * sizeof(TYPE));
  }

  constexpr ~Vector() = default;

  //! return the size (dimension) of the vector
  constexpr unsigned int size() const { return d_size; }

  //! returns the value at a particular index
  TYPE getVal(unsigned int i) const {
    PRECONDITION(i < d_size, "bad index");
    return d_data[i];
  }

  //! sets the index at a particular value
  void setVal(unsigned int i, TYPE val) {
    PRECONDITION(i < d_size, "bad index");
    d_data[i] = val;
  }

  TYPE operator[](unsigned int i) const {
    PRECONDITION(i < d_size, "bad index");
    return d_data[i];
  }

  TYPE &operator[](unsigned int i) {
    PRECONDITION(i < d_size, "bad index");
    return d_data[i];
  }

  //! returns a pointer to our data array
  constexpr TYPE *getData() { return d_data.get(); }

  //! returns a const pointer to our data array
  constexpr const TYPE *getData() const {
    // return dp_data;
    return d_data.get();
  }

  //! Copy operator.
  /*! We make a copy of the other Vector's data.
   */

  Vector<TYPE> &assign(const Vector<TYPE> &other) {
    PRECONDITION(d_size == other.size(), "Size mismatch in vector copying");
    const TYPE *otherData = other.getData();
    memcpy(static_cast<void *>(d_data.get()),
           static_cast<const void *>(otherData), d_size * sizeof(TYPE));
    return *this;
  }

  //! elementwise addition, vectors must be the same size.
  Vector<TYPE> &operator+=(const Vector<TYPE> &other) {
    PRECONDITION(d_size == other.size(), "Size mismatch in vector addition");
    const TYPE *otherData = other.getData();
    TYPE *data = d_data.get();
    unsigned int i;
    for (i = 0; i < d_size; i++) {
      data[i] += otherData[i];
    }
    return *this;
  }

  //! elementwise subtraction, vectors must be the same size.
  Vector<TYPE> &operator-=(const Vector<TYPE> &other) {
    PRECONDITION(d_size == other.size(), "Size mismatch in vector subtraction");
    const TYPE *otherData = other.getData();
    TYPE *data = d_data.get();
    unsigned int i;
    for (i = 0; i < d_size; i++) {
      data[i] -= otherData[i];
    }
    return *this;
  }

  //! multiplication by a scalar
  constexpr Vector<TYPE> &operator*=(TYPE scale) {
    unsigned int i;
    for (i = 0; i < d_size; i++) {
      d_data[i] *= scale;
    }
    return *this;
  }

  //! division by a scalar
  constexpr Vector<TYPE> &operator/=(TYPE scale) {
    unsigned int i;
    for (i = 0; i < d_size; i++) {
      d_data[i] /= scale;
    }
    return *this;
  }

  //! L2 norm squared
  constexpr inline TYPE normL2Sq() const {
    TYPE res = (TYPE)0.0;
    unsigned int i;
    TYPE *data = d_data.get();
    for (i = 0; i < d_size; i++) {
      res += data[i] * data[i];
    }
    return res;
  }

  //! L2 norm
  constexpr TYPE normL2() const { return sqrt(this->normL2Sq()); }

  //! L1 norm
  constexpr TYPE normL1() const {
    TYPE res = (TYPE)0.0;
    unsigned int i;
    TYPE *data = d_data.get();
    for (i = 0; i < d_size; i++) {
      res += fabs(data[i]);
    }
    return res;
  }

  //! L-infinity norm
  constexpr TYPE normLinfinity() const {
    TYPE res = (TYPE)(-1.0);
    unsigned int i;
    TYPE *data = d_data.get();
    for (i = 0; i < d_size; i++) {
      if (fabs(data[i]) > res) {
        res = fabs(data[i]);
      }
    }
    return res;
  }

  //! \brief Gets the ID of the entry that has the largest absolute value
  //! i.e. the entry being used for the L-infinity norm
  constexpr unsigned int largestAbsValId() const {
    TYPE res = (TYPE)(-1.0);
    unsigned int i, id = d_size;
    TYPE *data = d_data.get();
    for (i = 0; i < d_size; i++) {
      if (fabs(data[i]) > res) {
        res = fabs(data[i]);
        id = i;
      }
    }
    return id;
  }

  //! \brief Gets the ID of the entry that has the largest value
  constexpr unsigned int largestValId() const {
    TYPE res = (TYPE)(-1.e8);
    unsigned int i, id = d_size;
    TYPE *data = d_data.get();
    for (i = 0; i < d_size; i++) {
      if (data[i] > res) {
        res = data[i];
        id = i;
      }
    }
    return id;
  }

  //! \brief Gets the ID of the entry that has the smallest value
  constexpr unsigned int smallestValId() const {
    TYPE res = (TYPE)(1.e8);
    unsigned int i, id = d_size;
    TYPE *data = d_data.get();
    for (i = 0; i < d_size; i++) {
      if (data[i] < res) {
        res = data[i];
        id = i;
      }
    }
    return id;
  }

  //! returns the dot product between two Vectors
  TYPE dotProduct(const Vector<TYPE> other) const {
    PRECONDITION(d_size == other.size(),
                 "Size mismatch in vector doct product");
    const TYPE *oData = other.getData();
    unsigned int i;
    TYPE res = (TYPE)(0.0);
    TYPE *data = d_data.get();
    for (i = 0; i < d_size; i++) {
      res += (data[i] * oData[i]);
    }
    return res;
  }

  //! Normalize the vector using the L2 norm
  constexpr void normalize() {
    TYPE val = this->normL2();
    if (val < zero_tolerance) {
      throw std::runtime_error("Cannot normalize a zero length vector");
    }
    (*this) /= val;
  }

  //! Set to a random unit vector
  void setToRandom(unsigned int seed = 0) {
    // we want to get our own RNG here instead of using the global
    // one.  This is related to Issue285.
    RDKit::rng_type generator(42u);
    RDKit::uniform_double dist(0, 1.0);
    RDKit::double_source_type randSource(generator, dist);
    if (seed > 0) {
      generator.seed(seed);
    } else {
      // we can't initialize using only clock(), because it's possible
      // that we'll get here fast enough that clock() will return 0
      // and generator.seed(0) is an error:
      generator.seed(clock() + 1);
    }

    unsigned int i;
    TYPE *data = d_data.get();
    for (i = 0; i < d_size; i++) {
      data[i] = randSource();
    }
    this->normalize();
  }

 private:
  unsigned int d_size;  //!< our length
  DATA_SPTR d_data;
  Vector<TYPE> &operator=(const Vector<TYPE> &other);
};

typedef Vector<double> DoubleVector;

//! returns the algebraic tanimoto similarity [defn' from JCIM 46:587-96 (2006)]
template <typename T>
constexpr double TanimotoSimilarity(const Vector<T> &v1, const Vector<T> &v2) {
  double numer = v1.dotProduct(v2);
  if (numer == 0.0) {
    return 0.0;
  }
  double denom = v1.normL2Sq() + v2.normL2Sq() - numer;
  if (denom == 0.0) {
    return 0.0;
  }
  return numer / denom;
}
}  // end of namespace RDNumeric

//! ostream operator for Vectors
template <typename TYPE>
constexpr std::ostream &operator<<(std::ostream &target,
                                   const RDNumeric::Vector<TYPE> &vec) {
  unsigned int siz = vec.size();
  target << "Size: " << siz << " [";
  unsigned int i;
  for (i = 0; i < siz; i++) {
    target << std::setw(7) << std::setprecision(3) << vec.getVal(i) << ", ";
  }
  target << "]\n";
  return target;
}

#endif
