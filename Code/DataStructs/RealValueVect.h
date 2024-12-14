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
#include <RDGeneral/export.h>
#ifndef RD_REAL_VALUE_VECT_20140407
#define RD_REAL_VALUE_VECT_20140407

#include <boost/smart_ptr.hpp>
#include <string>
#include <cstring>

namespace RDGeom {
class UniformRealValueGrid3D;
}
namespace RDKit {

//! a class for efficiently storing vectors of double values
//! Has additional features compared to std::vector<double>:
//! construct from and write to pickle
class RDKIT_DATASTRUCTS_EXPORT RealValueVect {
 public:
  RealValueVect() = default;
  //! initialize with a particular size
  RealValueVect(unsigned int length) : d_length(length) {
    d_data.resize(d_length, 0.0);
  }

  //! initialize with a particular size and initial value
  RealValueVect(unsigned int length, double val) : d_length(length) {
    d_data.resize(d_length, val);
  }

  //! constructor from a pickle
  RealValueVect(const std::string &pkl) {
    initFromText(pkl.c_str(), pkl.size());
  };

  //! constructor from a pickle
  RealValueVect(const char *pkl, unsigned int len) { initFromText(pkl, len); };

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

  //! compares 2 vectors and returns false if different
  bool compareVectors(const RealValueVect &other) const;

  //! in-place operator&
  RealValueVect &operator&=(const RealValueVect &other);

  //! in-place operator|
  RealValueVect &operator|=(const RealValueVect &other);

  //! in-place operator+
  RealValueVect &operator+=(const RealValueVect &other);

  //! in-place operator-
  RealValueVect &operator-=(const RealValueVect &other);

  //! returns a binary string representation (pickle)
  std::string toString() const;

  void setLength(unsigned int sz) {
    d_length = sz;
    d_data.resize(sz);
  }
  void setToVal(double val) { std::fill(d_data.begin(), d_data.end(), val); }

  const std::vector<double> &getData() const { return d_data; }
  std::vector<double> &getData() { return d_data; }

 private:
  void initFromText(const char *pkl, const unsigned int len);
  unsigned int d_length = 0;
  std::vector<double> d_data;

  template <typename O>
  RealValueVect &applyBinaryOp(const RealValueVect &other, O op);
};  // end of declaration of class RealValueVect

//! returns L1 Norm of vectors
RDKIT_DATASTRUCTS_EXPORT double computeL1Norm(const RealValueVect &v1,
                                              const RealValueVect &v2);

//! returns the sum of vectors
RDKIT_DATASTRUCTS_EXPORT RealValueVect operator+(const RealValueVect &p1,
                                                 const RealValueVect &p2);

//! returns the difference of vectors
RDKIT_DATASTRUCTS_EXPORT RealValueVect operator-(const RealValueVect &p1,
                                                 const RealValueVect &p2);

//! support rvv3 = rvv1|rvv2
/*!

   operator| returns the maximum value for each element.
   e.g.:
     [0,1,2,0] | [0,1,1,1] -> [0,1,2,1]

 */
RDKIT_DATASTRUCTS_EXPORT RealValueVect operator|(const RealValueVect &p1,
                                                 const RealValueVect &p2);

//! support rvv3 = rvv1&rvv2
/*!
   operator& returns the minimum value for each element.
   e.g.:
     [0,1,2,0] & [0,1,1,1] -> [0,1,1,0]

 */
RDKIT_DATASTRUCTS_EXPORT RealValueVect operator&(const RealValueVect &p1,
                                                 const RealValueVect &p2);

}  // end of namespace RDKit

#endif
