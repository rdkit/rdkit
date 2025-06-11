//
//  Copyright (C) 2025 NVIDIA Corporation & Affiliates and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RD_VEC_UTILS_H
#define RD_VEC_UTILS_H
#include <vector>
#include <cstddef>

#include <RDGeneral/Invariant.h>

namespace RDKit {

//! Erases multiple indices from an iterable source. Indices must be in ascending order.
//! Uses only one copy per element in vec, and one vector resize. Does not use
//! intermediate buffers.
template <typename T, typename iterT>
void eraseMultipleIndices(std::vector<T> &vec, iterT indices, const size_t numIndices) {
  using indexT = std::remove_cv_t<std::remove_reference_t<decltype(*indices)>>;
  static_assert(std::is_integral_v<indexT>, "Indices must be integral");
  if (numIndices == 0) {
    return;
  }
  for (size_t i = 0; i < numIndices; ++i) {
    const size_t idxToRemove = *indices;
    ++indices;
    PRECONDITION((i == numIndices - 1) || idxToRemove < *indices, "Indices must be in ascending order");
    PRECONDITION(idxToRemove >= 0 && idxToRemove < vec.size(), "Index out of range");
    if (idxToRemove == vec.size() - 1) {
      break;
    }
    const size_t firstIdxToCopy = idxToRemove + 1;
    const size_t lastIdxToCopy = (i == numIndices - 1) ? vec.size() : *indices;
    const size_t destToCopy = idxToRemove - i;

    std::copy(vec.begin() + firstIdxToCopy,
              vec.begin() + lastIdxToCopy,
              vec.begin() + destToCopy);
  }
  vec.resize(vec.size() - numIndices);
}
template <typename T, typename indexT>
void eraseMultipleIndices(std::vector<T> &vec, const std::vector<indexT> indices) {
  eraseMultipleIndices(vec, indices.data(), indices.size());
}
}  // namespace RDKit

#endif // RD_VEC_UTILS_H