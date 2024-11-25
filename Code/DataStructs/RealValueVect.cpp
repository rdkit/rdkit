//
//  Copyright (c) 2014-2024, Novartis Institutes for BioMedical Research and
//  other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "RealValueVect.h"
#include <RDGeneral/Invariant.h>
#include <RDGeneral/StreamOps.h>
#include "DatastructsException.h"

constexpr double VAL_TOL = 0.01;

namespace RDKit {
const int ci_REALVALUEVECTPICKLE_VERSION = 0x1;

double RealValueVect::getVal(unsigned int i) const {
  if (i >= d_length) {
    throw IndexErrorException(i);
  }
  return d_data[i];
}

void RealValueVect::setVal(unsigned int i, double val) {
  if (i >= d_length) {
    throw IndexErrorException(i);
  }
  d_data[i] = val;
}

double RealValueVect::getTotalVal() const {
  double res = 0;

  for (unsigned int i = 0; i < d_length; ++i) {
    res += d_data[i];
  }
  return res;
}

bool RealValueVect::compareVectors(const RealValueVect &other) {
  if (getLength() != other.getLength()) {
    throw ValueErrorException("Comparing vectors of different lengths");
  }
  return std::equal(this->d_data.begin(), this->d_data.end(),
                    other.d_data.begin(), [](auto v1, auto v2) {
                      return fabs((v1 - v2) / (v1 != 0 ? v1 : 1)) <= VAL_TOL;
                    });
}

double computeL1Norm(const RealValueVect &v1, const RealValueVect &v2) {
  if (v1.getLength() != v2.getLength()) {
    throw ValueErrorException("Comparing vectors of different lengths");
  }
  double res = 0.0;
  for (auto i = 0u; i < v1.getLength(); ++i) {
    res += fabs(v1.getData()[i] - v2.getData()[i]);
  }
  return res;
}

std::string RealValueVect::toString() const {
  std::stringstream ss(std::ios_base::binary | std::ios_base::out |
                       std::ios_base::in);

  boost::int32_t tVers = ci_REALVALUEVECTPICKLE_VERSION * -1;
  streamWrite(ss, tVers);
  boost::uint32_t tInt;
  tInt = d_length;
  streamWrite(ss, tInt);

#if defined(BOOST_BIG_ENDIAN)
  double *td = new double[d_length];
  for (unsigned int i = 0; i < d_length; ++i)
    td[i] = EndianSwapBytes<HOST_ENDIAN_ORDER, LITTLE_ENDIAN_ORDER>(d_data[i]);
  ss.write((const char *)td, d_length * sizeof(double));
  delete[] td;
#else
  ss.write((const char *)d_data.data(), d_length * sizeof(double));
#endif
  std::string res(ss.str());
  return res;
};

void RealValueVect::initFromText(const char *pkl, const unsigned int len) {
  std::stringstream ss(std::ios_base::binary | std::ios_base::in |
                       std::ios_base::out);
  ss.write(pkl, len);
  boost::int32_t tVers;
  streamRead(ss, tVers);
  tVers *= -1;
  if (tVers == 0x1) {
  } else {
    throw ValueErrorException("bad version in RealValueVect pickle");
  }
  boost::uint32_t tInt;
  streamRead(ss, tInt);
  d_length = tInt;
  d_data.resize(d_length);
  auto *data = d_data.data();
  ss.read((char *)data, d_length * sizeof(double));

#if defined(BOOST_BIG_ENDIAN)
  for (auto &v : d_data) {
    v = EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(v);
  }
#endif
};

RealValueVect RealValueVect::operator&(const RealValueVect &other) const {
  PRECONDITION(other.d_length == d_length, "length mismatch");
  RealValueVect ans(d_length);
  for (unsigned int i = 0; i < d_length; ++i) {
    double v1 = getVal(i);
    double v2 = other.getVal(i);
    ans.setVal(i, std::min(v1, v2));
  }
  return (ans);
};

RealValueVect RealValueVect::operator|(const RealValueVect &other) const {
  PRECONDITION(other.d_length == d_length, "length mismatch");
  RealValueVect ans(d_length);
  for (unsigned int i = 0; i < d_length; ++i) {
    double v1 = getVal(i);
    double v2 = other.getVal(i);
    ans.setVal(i, std::max(v1, v2));
  }
  return (ans);
};

RealValueVect &RealValueVect::operator+=(const RealValueVect &other) {
  PRECONDITION(other.d_length == d_length, "length mismatch");
  for (unsigned int i = 0; i < d_length; i++) {
    double v1 = getVal(i);
    double v2 = other.getVal(i);
    setVal(i, v1 + v2);
  }
  return *this;
}

RealValueVect &RealValueVect::operator-=(const RealValueVect &other) {
  PRECONDITION(other.d_length == d_length, "length mismatch");

  for (unsigned int i = 0; i < d_length; i++) {
    double v1 = getVal(i);
    double v2 = other.getVal(i);
    setVal(i, v1 - v2);
  }
  return *this;
}

RealValueVect operator+(const RealValueVect &p1, const RealValueVect &p2) {
  RealValueVect res(p1);
  res += p2;
  return res;
};

RealValueVect operator-(const RealValueVect &p1, const RealValueVect &p2) {
  RealValueVect res(p1);
  res -= p2;
  return res;
};

}  // end of namespace RDKit
