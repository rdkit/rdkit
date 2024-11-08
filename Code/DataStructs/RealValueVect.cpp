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

#define VAL_TOL 0.01

namespace RDKit {
const int ci_REALVALUEVECTPICKLE_VERSION = 0x1;

RealValueVect::RealValueVect(const RealValueVect &other) {
  d_length = other.getLength();
  const double *odata = other.getData();
  double *data = new double[d_length];
  memcpy(static_cast<void *>(data), static_cast<const void *>(odata),
         d_length * sizeof(double));
  d_data.reset(data);
}

RealValueVect &RealValueVect::operator=(const RealValueVect &other) {
  if (this == &other) {
    return *this;
  }

  d_length = other.getLength();
  const double *odata = other.getData();
  double *data = new double[d_length];
  memcpy(static_cast<void *>(data), static_cast<const void *>(odata),
         d_length * sizeof(double));
  d_data.reset(data);
  return *this;
};

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

  const double *sd1 = this->getData();
  const double *sd2 = other.getData();

  const double *send = sd1 + this->getLength();
  while (sd1 != send) {
    if (fabs((*sd1 - *sd2) / (*sd1)) > VAL_TOL) {
      return false;
    }
    sd1++;
    sd2++;
  }
  return true;
}

double computeL1Norm(const RealValueVect &v1, const RealValueVect &v2) {
  if (v1.getLength() != v2.getLength()) {
    throw ValueErrorException("Comparing vectors of different lengths");
  }

  const double *data1 = v1.getData();
  const double *data2 = v2.getData();

  double res = 0.0;

  const double *sd1 = data1;
  const double *sd2 = data2;

  const double *send = sd1 + v1.getLength();
  while (sd1 != send) {
    if (*sd1 != *sd2) {
      res += fabs((*sd1) - (*sd2));
    }
    sd1++;
    sd2++;
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
    td[i] = EndianSwapBytes<HOST_ENDIAN_ORDER, LITTLE_ENDIAN_ORDER>(
        d_data.get()[i]);
  ss.write((const char *)td, d_length * sizeof(double));
  delete[] td;
#else
  ss.write((const char *)d_data.get(), d_length * sizeof(double));
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
  double *data = new double[d_length];
  ss.read((char *)data, d_length * sizeof(double));

#if defined(BOOST_BIG_ENDIAN)
  double *td = new double[d_length];
  for (unsigned int i = 0; i < d_length; ++i)
    td[i] = EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(data[i]);
  d_data.reset(td);
  delete[] data;
#else
  d_data.reset(data);
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
