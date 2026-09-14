//
//  Copyright (C) 2026 Marc Lehner and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "AtomFeatures.h"

#include <array>

#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/ROMol.h>
#include <RDGeneral/Invariant.h>

namespace RDKit {
namespace DASH {

namespace {

// The atom-feature table of the published DASH trees, in branch-index order.
// Generated from serenityff.charge.tree.atom_features.AtomFeatures.feature_list;
// the order is part of the file format and must not be changed.
// clang-format off
const AtomFeature featureTable[numAtomFeatures] = {
    {  5, 1,  0, false, 2},  // 0
    {  5, 3, -1, true , 0},  // 1
    {  5, 3,  0, false, 0},  // 2
    {  5, 3,  0, false, 2},  // 3
    {  5, 4, -1, false, 0},  // 4
    {  5, 4, -1, false, 2},  // 5
    { 35, 1,  0, false, 0},  // 6
    {  6, 1, -1, false, 0},  // 7
    {  6, 1, -1, true , 0},  // 8
    {  6, 1,  0, false, 1},  // 9
    {  6, 1,  0, false, 2},  // 10
    {  6, 1,  0, false, 3},  // 11
    {  6, 1,  0, true , 1},  // 12
    {  6, 2, -1, true , 1},  // 13
    {  6, 2,  0, false, 0},  // 14
    {  6, 2,  0, false, 1},  // 15
    {  6, 2,  0, false, 2},  // 16
    {  6, 2,  0, true , 0},  // 17
    {  6, 2,  0, true , 1},  // 18
    {  6, 3, -1, false, 0},  // 19
    {  6, 3, -1, true , 0},  // 20
    {  6, 3, -1, false, 1},  // 21
    {  6, 3, -1, true , 1},  // 22
    {  6, 3,  0, false, 0},  // 23
    {  6, 3,  0, false, 1},  // 24
    {  6, 3,  0, false, 2},  // 25
    {  6, 3,  0, true , 0},  // 26
    {  6, 3,  0, true , 1},  // 27
    {  6, 3,  0, true , 2},  // 28
    {  6, 3,  1, false, 0},  // 29
    {  6, 3,  1, true , 0},  // 30
    {  6, 4,  0, false, 0},  // 31
    {  6, 4,  0, false, 1},  // 32
    {  6, 4,  0, false, 2},  // 33
    {  6, 4,  0, false, 3},  // 34
    { 17, 1,  0, false, 0},  // 35
    {  9, 1,  0, false, 0},  // 36
    {  1, 1,  0, false, 0},  // 37
    { 53, 1,  0, false, 0},  // 38
    { 53, 3,  0, false, 0},  // 39
    { 53, 3,  0, false, 1},  // 40
    { 53, 4,  0, false, 0},  // 41
    {  7, 1, -1, false, 0},  // 42
    {  7, 1, -1, true , 0},  // 43
    {  7, 1,  0, false, 0},  // 44
    {  7, 1,  0, true , 0},  // 45
    {  7, 1,  0, true , 1},  // 46
    {  7, 1,  0, true , 2},  // 47
    {  7, 1,  1, false, 3},  // 48
    {  7, 2, -1, false, 0},  // 49
    {  7, 2, -1, false, 1},  // 50
    {  7, 2, -1, true , 0},  // 51
    {  7, 2, -1, true , 1},  // 52
    {  7, 2,  0, false, 0},  // 53
    {  7, 2,  0, false, 1},  // 54
    {  7, 2,  0, true , 0},  // 55
    {  7, 2,  0, true , 1},  // 56
    {  7, 2,  1, false, 0},  // 57
    {  7, 2,  1, false, 1},  // 58
    {  7, 2,  1, false, 2},  // 59
    {  7, 2,  1, true , 0},  // 60
    {  7, 2,  1, true , 1},  // 61
    {  7, 3,  0, false, 0},  // 62
    {  7, 3,  0, false, 1},  // 63
    {  7, 3,  0, false, 2},  // 64
    {  7, 3,  0, true , 0},  // 65
    {  7, 3,  0, true , 1},  // 66
    {  7, 3,  0, true , 2},  // 67
    {  7, 3,  1, false, 0},  // 68
    {  7, 3,  1, false, 1},  // 69
    {  7, 3,  1, true , 0},  // 70
    {  7, 3,  1, true , 1},  // 71
    {  7, 4,  1, false, 0},  // 72
    {  7, 4,  1, false, 1},  // 73
    {  7, 4,  1, false, 2},  // 74
    {  7, 4,  1, false, 3},  // 75
    {  8, 1, -1, false, 0},  // 76
    {  8, 1, -1, true , 0},  // 77
    {  8, 1,  0, false, 0},  // 78
    {  8, 1,  0, false, 1},  // 79
    {  8, 1,  0, true , 0},  // 80
    {  8, 1,  0, true , 1},  // 81
    {  8, 2,  0, false, 0},  // 82
    {  8, 2,  0, false, 1},  // 83
    {  8, 2,  0, false, 2},  // 84
    {  8, 2,  0, true , 0},  // 85
    {  8, 2,  0, true , 1},  // 86
    {  8, 2,  0, true , 2},  // 87
    {  8, 2,  1, false, 0},  // 88
    {  8, 2,  1, false, 1},  // 89
    {  8, 2,  1, true , 0},  // 90
    {  8, 3,  1, false, 0},  // 91
    {  8, 3,  1, false, 1},  // 92
    { 15, 1,  0, false, 1},  // 93
    { 15, 2,  0, false, 0},  // 94
    { 15, 2,  0, false, 1},  // 95
    { 15, 2,  0, true , 0},  // 96
    { 15, 3,  0, false, 0},  // 97
    { 15, 4,  0, false, 0},  // 98
    { 15, 4,  0, false, 1},  // 99
    { 15, 4,  0, true , 0},  // 100
    { 15, 4,  0, true , 1},  // 101
    { 15, 4,  1, false, 0},  // 102
    { 15, 5,  0, false, 0},  // 103
    { 15, 5,  0, false, 1},  // 104
    { 16, 1, -1, false, 0},  // 105
    { 16, 1, -1, true , 0},  // 106
    { 16, 1,  0, false, 0},  // 107
    { 16, 1,  0, true , 0},  // 108
    { 16, 2,  0, false, 0},  // 109
    { 16, 2,  0, true , 0},  // 110
    { 16, 2,  0, false, 1},  // 111
    { 16, 2,  1, true , 0},  // 112
    { 16, 3,  0, false, 0},  // 113
    { 16, 3,  0, true , 0},  // 114
    { 16, 3,  1, true , 0},  // 115
    { 16, 3,  1, false, 0},  // 116
    { 16, 4,  0, false, 0},  // 117
    { 16, 4,  0, false, 1},  // 118
    { 16, 4,  0, true , 0},  // 119
    { 16, 4,  1, false, 0},  // 120
    { 16, 4,  1, true , 0},  // 121
};
// clang-format on

// Bit budget of the direct-lookup table. The table covers every tuple that can
// be encoded; anything outside is rejected before the lookup, so an exotic atom
// costs one comparison rather than a search.
const unsigned int lookupAtomicNumBits = 6;  // Z <= 63 (the table tops out at I, 53)
const unsigned int lookupDegreeBits = 3;     // degree <= 7 (table: 1..5)
const unsigned int lookupChargeBits = 2;     // charge + 1 in 0..3 (table: -1..1)
const unsigned int lookupConjBits = 1;
const unsigned int lookupNumHsBits = 2;      // <= 3
const unsigned int lookupBits = lookupAtomicNumBits + lookupDegreeBits +
                                lookupChargeBits + lookupConjBits +
                                lookupNumHsBits;  // 14 -> 16 kB

inline unsigned int lookupKey(unsigned int atomicNum, unsigned int degree,
                              int formalCharge, bool conjugated,
                              unsigned int numHs) {
  unsigned int k = atomicNum;
  k = (k << lookupDegreeBits) | degree;
  k = (k << lookupChargeBits) | static_cast<unsigned int>(formalCharge + 1);
  k = (k << lookupConjBits) | (conjugated ? 1u : 0u);
  k = (k << lookupNumHsBits) | numHs;
  return k;
}

//! the direct-lookup table, built once
class FeatureLookup {
 public:
  FeatureLookup() {
    d_index.fill(-1);
    for (int i = 0; i < static_cast<int>(numAtomFeatures); ++i) {
      const AtomFeature &f = featureTable[i];
      d_index[lookupKey(f.atomicNum, f.degree, f.formalCharge, f.conjugated,
                        f.numHs)] = static_cast<std::int8_t>(i);
    }
  }
  std::int8_t at(unsigned int key) const { return d_index[key]; }

 private:
  std::array<std::int8_t, (1u << lookupBits)> d_index;
};

const FeatureLookup &featureLookup() {
  static const FeatureLookup lookup;
  return lookup;
}

}  // namespace

const AtomFeature *getAtomFeatureTable() { return featureTable; }

int atomFeatureIndex(unsigned int atomicNum, unsigned int degree,
                     int formalCharge, bool conjugated, unsigned int numHs) {
  if (atomicNum >= (1u << lookupAtomicNumBits) ||
      degree >= (1u << lookupDegreeBits) || formalCharge < -1 ||
      formalCharge > 2 || numHs >= (1u << lookupNumHsBits)) {
    return -1;
  }
  return featureLookup().at(
      lookupKey(atomicNum, degree, formalCharge, conjugated, numHs));
}

int atomFeatureIndex(const Atom *atom) {
  PRECONDITION(atom, "bad atom pointer");
  bool conjugated = false;
  for (const auto bond : atom->getOwningMol().atomBonds(atom)) {
    if (bond->getIsConjugated()) {
      conjugated = true;
      break;
    }
  }
  return atomFeatureIndex(atom->getAtomicNum(), atom->getDegree(),
                          atom->getFormalCharge(), conjugated,
                          atom->getTotalNumHs(true));
}

int dashBondType(const Bond *bond) {
  if (!bond) {
    return -1;
  }
  if (bond->getIsConjugated()) {
    return 4;
  }
  return static_cast<int>(bond->getBondType());
}

}  // namespace DASH
}  // namespace RDKit
