//
//  Copyright (C) 2026 Marc Lehner and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//! \file MatchKey.h
/// \brief the packed node match key.
/// This header is not installed. Include it only from within the DASHTree
/// library and its tests.

#ifndef RD_DASH_MATCHKEY_H
#define RD_DASH_MATCHKEY_H

#include <cstdint>

namespace RDKit {
namespace DASH {

constexpr unsigned int keyAtomTypeBits = 8;  //!< bits 0-7: feature class
constexpr unsigned int keyConAtomBits = 4;   //!< bits 8-11: attachment + 1
constexpr unsigned int keyConTypeBits = 3;   //!< bits 12-14: bond type + 1
constexpr unsigned int keyConAtomShift = keyAtomTypeBits;
constexpr unsigned int keyConTypeShift = keyAtomTypeBits + keyConAtomBits;

constexpr int maxKeyAtomType = (1 << keyAtomTypeBits) - 1;  // 255
constexpr int maxKeyConAtom = (1 << keyConAtomBits) - 2;    // 14, stored as +1
constexpr int maxKeyConType = (1 << keyConTypeBits) - 2;    // 6,  stored as +1

//! \brief a key that matches no node. Bit 15 is in no field, so no stored key
//! has it. A component out of range gives this key (the converter refuses
//! trees that need one), and withConAtom keeps it. A consumed candidate also
//! gets this key.
constexpr std::uint16_t noMatchKey = 0xFFFF;

//! \brief packs the part of a match key that does not depend on the match state
/*!
  The feature class and the bond descriptor are fixed for a given (atom,
  neighbour) pair. Only the position of the attachment atom changes as the
  match grows. With this half precomputed per bond, candidate generation is a
  bit-or.
*/
inline std::uint16_t packPartialMatchKey(int atomType, int conType) {
  if (atomType < 0 || atomType > maxKeyAtomType || conType < -1 ||
      conType > maxKeyConType) {
    return noMatchKey;
  }
  return static_cast<std::uint16_t>(atomType |
                                    ((conType + 1) << keyConTypeShift));
}

//! grafts an attachment position onto a partial key
inline std::uint16_t withConAtom(std::uint16_t partialKey, int conAtom) {
  if (conAtom < 0 || conAtom > maxKeyConAtom) {
    return noMatchKey;
  }
  return static_cast<std::uint16_t>(partialKey |
                                    ((conAtom + 1) << keyConAtomShift));
}

//! \name the components of a packed key
//! The inverse of packPartialMatchKey and withConAtom. A component that the
//! key does not have reads as -1. This is the case for the attachment of a
//! branch root, and of the heavy-atom child of a hydrogen root. withConAtom
//! did not build such a key, and packPartialMatchKey takes -1 for no bond.
//@{
inline int keyAtomType(std::uint16_t key) {
  return key & ((1 << keyAtomTypeBits) - 1);
}
inline int keyConAtom(std::uint16_t key) {
  return static_cast<int>((key >> keyConAtomShift) &
                          ((1 << keyConAtomBits) - 1)) -
         1;
}
inline int keyConType(std::uint16_t key) {
  return static_cast<int>((key >> keyConTypeShift) &
                          ((1 << keyConTypeBits) - 1)) -
         1;
}
//@}

}  // namespace DASH
}  // namespace RDKit

#endif
