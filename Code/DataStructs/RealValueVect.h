//  RealValueVect.h
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
#ifndef __RD_REAL_VALUE_VECT_20140407__
#define __RD_REAL_VALUE_VECT_20140407__

#include <boost/smart_ptr.hpp>
#include <string>
#include <cstring>

namespace RDKit {

//! a class for efficiently storing vectors of float values
//! Has additional features compared to std::vector<double>:
//! construct from and write to pickle
class RealValueVect {
 public:
  typedef boost::shared_array<double> DATA_SPTR;

  //! initialize with a particular size
  RealValueVect(unsigned int length) : d_length(length) {
    double *data = new double[d_length];
    memset(static_cast<void *>(data), 0, d_length * sizeof(double));
    d_data.reset(data);
  }

  //! initialize with a particular size and initial value
  RealValueVect(double val, unsigned int length) : d_length(length) {
    double *data = new double[d_length];
    memset(static_cast<void *>(data), 0, d_length * sizeof(double));
    d_data.reset(data);
    for (unsigned int i = 0; i < d_length; ++i) {
      d_data[i] = val;
    }
  }

  //! Copy constructor
  RealValueVect(const RealValueVect &other);
  RealValueVect &operator=(const RealValueVect &other);

  //! constructor from a pickle
  RealValueVect(const std::string &pkl) {
    initFromText(pkl.c_str(), pkl.size());
  };

  //! constructor from a pickle
  RealValueVect(const char *pkl, unsigned int len) { initFromText(pkl, len); };

  ~RealValueVect() {}

  //! return the value at an index
  double getVal(unsigned int i) const;

  //! support indexing using []
  double operator[](unsigned int idx) const { return getVal(idx); };

  //! set the value at an index
  void setVal(unsigned int i, double val);

  //! returns the sum of all the elements in the vect
  double getTotalVal() const;

  //! returns the length
  unsigned int getLength() const { return d_length; };

  //! returns the length
  unsigned int size() const { return getLength(); };

  //! return a pointer to our raw data storage
  const double *getData() const { return d_data.get(); };

  //! compares 2 vectors and returns false if different
  bool compareVectors(const RealValueVect &other);

  const DATA_SPTR &getArray() { return d_data; };

  //! support rvv3 = rvv1&rvv2
  /*!
     operator& returns the minimum value for each element.
     e.g.:
       [0,1,2,0] & [0,1,1,1] -> [0,1,1,0]

   */
  RealValueVect operator&(const RealValueVect &other) const;

  //! support rvv3 = rvv1|rvv2
  /*!

     operator| returns the maximum value for each element.
     e.g.:
       [0,1,2,0] | [0,1,1,1] -> [0,1,2,1]

   */
  RealValueVect operator|(const RealValueVect &other) const;

  //! sums up vectors
  RealValueVect &operator+=(const RealValueVect &other);

  //! subtracts vectors
  RealValueVect &operator-=(const RealValueVect &other);

  //! returns a binary string representation (pickle)
  std::string toString() const;

 private:
  unsigned int d_length;
  DATA_SPTR d_data;

  void initFromText(const char *pkl, const unsigned int len);
};  // end of declaration of class RealValueVect

//! returns L1 Norm of vectors
double computeL1Norm(const RealValueVect &v1, const RealValueVect &v2);

//! returns the sum of vectors
RealValueVect operator+(const RealValueVect &p1, const RealValueVect &p2);

//! returns the difference of vectors
RealValueVect operator-(const RealValueVect &p1, const RealValueVect &p2);

}  // end of namespace RDKit

#endif
