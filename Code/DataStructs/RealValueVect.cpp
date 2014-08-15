//  RealValueVect.cpp
//  Created on: Apr 7, 2014
//  Author: hahnda6
//
//  Copyright (c) 2014, Novartis Institutes for BioMedical Research Inc.
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
