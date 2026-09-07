//
//  Copyright (C) 2026 ETH Zurich
//  Created by: Katharina Buchthal
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#include <cstddef>
#ifndef RD_Z_MATRIX_H
#define RD_Z_MATRIX_H

#include <Numerics/SquareMatrix.h>
#include <RDGeneral/Invariant.h>
#include <iterator>
#include <optional>
#include <variant>

#include "ZMatrixUtils.h"

namespace DistGeom {

//! Class to store Z-Matrix
/*!
  Similar to a conventional ZMatrix with the following two main differences:
    * torsions are not stored as fixed values but as ranges or set of distinct
  values
    * stores information of dependend torisions, i.e., direct neighbors that are
  not the bond references with order (in Zmatrix) < than the current one
*/

template <typename T>
std::ostream &operator<<(std::ostream &os, const std::optional<T> &opt) {
  if (opt) {
    return os << opt.value();
  }
  return os << "-";
}

class RDKIT_DISTGEOMETRY_EXPORT ZMatrix {
 public:
  struct InternalCoordinateReference {
    std::optional<double> length;
    std::optional<unsigned int> bondRef;
    std::optional<double> angle;
    std::optional<unsigned int> angleRef;
    std::optional<TorsionCandidates> torsion;
    std::optional<unsigned int> torsionRef;
  };

  struct TorsionDependence {
    double offset;
    unsigned int reference;
  };

  struct ZMatrixRow {
    unsigned int atomIdx;
    const InternalCoordinateReference &internal;
    const std::optional<TorsionDependence> &torsionDependence;
  };

  class Iterator {
    // adapted from https://en.cppreference.com/cpp/iterator/iterator
   public:
    using iterator_category = std::bidirectional_iterator_tag;
    using iterator_concept = std::bidirectional_iterator_tag;
    using value_type = ZMatrixRow;
    using difference_type = std::ptrdiff_t;
    using reference = ZMatrixRow;
    using pointer = void;

    Iterator() = default;

    reference operator*() const { return (*d_matrix)[d_position]; }

    Iterator &operator++() {
      ++d_position;
      return *this;
    }
    Iterator operator++(int) {
      auto result = *this;
      ++(*this);
      return result;
    }
    Iterator &operator--() {
      --d_position;
      return *this;
    }
    Iterator operator--(int) {
      auto result = *this;
      --(*this);
      return result;
    }

    bool operator==(const Iterator &other) const {
      return d_matrix == other.d_matrix && d_position == other.d_position;
    }
    bool operator!=(const Iterator &other) const { return !(*this == other); }

   private:
    friend class ZMatrix;
    explicit Iterator(const ZMatrix *matrix, std::size_t position)
        : d_matrix(matrix), d_position(position) {}

    const ZMatrix *d_matrix = nullptr;
    std::size_t d_position = 0;
  };

  using iterator = Iterator;  // note: we only traverse on constant/finished a
                              // z-matrix, thus we do not need mutable iterators
  using const_iterator = Iterator;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  ZMatrix(unsigned int nElements)
      : d_order(),
        d_internalCoordRefs(nElements),
        d_torsionDependences(nElements) {
    d_order.reserve(nElements);
  }

  ~ZMatrix() = default;

  void addElement(unsigned int idx, std::optional<double> length = {},
                  std::optional<unsigned int> bondRef = {},
                  std::optional<double> angle = {},
                  std::optional<unsigned int> angleRef = {},
                  std::optional<TorsionCandidates> torsion = {},
                  std::optional<unsigned int> torionRef = {},
                  std::optional<TorsionDependence> torsionDep = {}) {
    PRECONDITION(idx < d_internalCoordRefs.size(), "Index out of range");
    PRECONDITION(!d_internalCoordRefs[idx], "element already set");

    d_order.emplace_back(idx);

    d_internalCoordRefs[idx] = {length,   bondRef, angle,
                                angleRef, torsion, torionRef};
    d_torsionDependences[idx] = std::move(torsionDep);
  }

  const InternalCoordinateReference &getReferences(
      const unsigned int atomIdx) const {
    PRECONDITION(atomIdx < d_internalCoordRefs.size(), "Invalid atom index");
    PRECONDITION(d_internalCoordRefs[atomIdx],
                 "Accessing an ZMatrix Element that was not set before");
    return *d_internalCoordRefs[atomIdx];
  }

  const std::optional<TorsionDependence> &getTorsionReference(
      const unsigned int atomIdx) const {
    PRECONDITION(atomIdx < d_torsionDependences.size(), "Invalid atom index");
    return d_torsionDependences[atomIdx];
  }

  void updateTorsion(const unsigned int atomIdx, TorsionRange torsion) {
    PRECONDITION(atomIdx < d_internalCoordRefs.size(), "Invalid atom index");
    d_internalCoordRefs[atomIdx]->torsion = {torsion};  // TODO check
  }

  void invertTorsionDependence(const unsigned int atomIdx) {
    PRECONDITION(atomIdx < d_torsionDependences.size(), "Invalid atom index");
    PRECONDITION(d_torsionDependences[atomIdx],
                 "Atom has no torsion dependence");
    d_torsionDependences[atomIdx]->offset =
        -d_torsionDependences[atomIdx]->offset;
  }

  ZMatrixRow operator[](unsigned int idx) const {
    PRECONDITION(idx < d_order.size(), "Invalid row index");
    auto atomIdx = d_order[idx];
    PRECONDITION(d_internalCoordRefs[atomIdx],
                 "Accessing row in matrix that was not set before");
    return {atomIdx, d_internalCoordRefs[atomIdx].value(),
            d_torsionDependences[atomIdx]};
  }

  friend std::ostream &operator<<(std::ostream &os, const ZMatrix &zmat);

  const_iterator cbegin() const { return const_iterator(this, 0); }
  const_iterator cend() const { return const_iterator(this, d_order.size()); }

  const_reverse_iterator crbegin() const {
    return const_reverse_iterator(cend());
  }
  const_reverse_iterator crend() const {
    return const_reverse_iterator(cbegin());
  }

  iterator begin() { return iterator(this, 0); }
  iterator end() { return iterator(this, d_order.size()); }
  const_iterator begin() const {
    return cbegin();
  }  // we need begin()/end() for ranges but do not want to allow iterators on a
     // mutable z matrix
  const_iterator end() const { return cend(); }

  reverse_iterator rbegin() { return reverse_iterator(end()); }
  reverse_iterator rend() { return reverse_iterator(begin()); }
  const_reverse_iterator rbegin() const { return crbegin(); }
  const_reverse_iterator rend() const { return crend(); }

 private:
  std::vector<unsigned int> d_order;
  std::vector<std::optional<InternalCoordinateReference>> d_internalCoordRefs;
  std::vector<std::optional<TorsionDependence>> d_torsionDependences;
};

typedef std::shared_ptr<ZMatrix> ZMatPtr;

inline std::ostream &operator<<(std::ostream &os, const ZMatrix &zmat) {
  os << "AtomIdx BondRefernece BondLength AngleReference BondAngle "
        "TorsionReference TorsionAngle TorsionDependence Offset\n";

  for (const auto &[atomIdx, internalReferences, torsionDependence] : zmat) {
    const auto &[bl, bondRef, ba, angleRef, torsion, torsionRef] =
        internalReferences;
    os << atomIdx << " " << bl << " " << bondRef << " " << ba << " " << angleRef
       << " " << torsionRef << " [";
    if (torsion) {
      std::visit(overloaded{[&os](const TorsionRange &r) {
                              os << r.lower << ',' << r.lower;
                            },
                            [&os](const TorsionValues &ts) {
                              for (const auto t : ts) {
                                os << dequantize(t) << "; ";
                              }
                            }},
                 *torsion);
    } else {
      os << "-,-";
    }
    os << "]";

    if (torsionDependence) {
      const auto &[offset, ref] = torsionDependence.value();
      os << " " << ref << " " << offset;
    } else {
      os << " - -";
    }
    os << "\n";
  }
  return os;
}

}  // namespace DistGeom

#endif
