//
//  Copyright (C) 2026 Marc Lehner and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "DASHTreeImpl.h"

#include <string>

#include <GraphMol/Atom.h>
#include <GraphMol/MolOps.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/Invariant.h>

namespace RDKit {
namespace DASH {

namespace {

// Bit budget of the direct-lookup table. The table covers every tuple that can
// be encoded; anything outside is rejected before the lookup, so an exotic atom
// costs one comparison rather than a search.
constexpr unsigned int lookupAtomicNumBits =
    6;  // Z <= 63 (the published trees top out at I, 53)
constexpr unsigned int lookupDegreeBits = 3;  // degree <= 7 (published: 1..5)
constexpr unsigned int lookupChargeBits =
    2;  // charge + 1 in 0..3 (published: -1..1)
constexpr unsigned int lookupConjBits = 1;
constexpr unsigned int lookupNumHsBits = 2;  // <= 3
constexpr unsigned int lookupBits = lookupAtomicNumBits + lookupDegreeBits +
                                    lookupChargeBits + lookupConjBits +
                                    lookupNumHsBits;  // 14 -> 32 kB of int16

inline bool fitsLookup(unsigned int atomicNum, unsigned int degree,
                       int formalCharge, unsigned int numHs) {
  return atomicNum < (1u << lookupAtomicNumBits) &&
         degree < (1u << lookupDegreeBits) && formalCharge >= -1 &&
         formalCharge <= 2 && numHs < (1u << lookupNumHsBits);
}

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

std::string describe(const AtomFeature &f) {
  return "Z " + std::to_string(f.atomicNum) + ", degree " +
         std::to_string(f.degree) + ", charge " +
         std::to_string(f.formalCharge) + ", " +
         (f.conjugated ? "conjugated, " : "not conjugated, ") +
         std::to_string(f.numHs) + " H";
}

}  // namespace

AtomFeatureLookup::AtomFeatureLookup(const AtomFeature *features,
                                     unsigned int numFeatures)
    : d_features(features),
      d_numFeatures(numFeatures),
      d_index(1u << lookupBits, -1) {
  PRECONDITION(features || !numFeatures, "bad feature table");
  for (unsigned int i = 0; i < numFeatures; ++i) {
    const AtomFeature &f = features[i];
    if (!fitsLookup(f.atomicNum, f.degree, f.formalCharge, f.numHs)) {
      throw ValueErrorException("atom feature " + std::to_string(i) + " (" +
                                describe(f) +
                                ") is outside the range this build can "
                                "classify");
    }
    std::int16_t &slot = d_index[lookupKey(
        f.atomicNum, f.degree, f.formalCharge, f.conjugated != 0, f.numHs)];
    if (slot >= 0) {
      throw ValueErrorException("atom features " + std::to_string(slot) +
                                " and " + std::to_string(i) +
                                " are the same tuple (" + describe(f) + ")");
    }
    slot = static_cast<std::int16_t>(i);
  }
}

int AtomFeatureLookup::index(unsigned int atomicNum, unsigned int degree,
                             int formalCharge, bool conjugated,
                             unsigned int numHs) const {
  if (d_index.empty() || !fitsLookup(atomicNum, degree, formalCharge, numHs)) {
    return -1;
  }
  return d_index[lookupKey(atomicNum, degree, formalCharge, conjugated, numHs)];
}

int AtomFeatureLookup::index(const Atom *atom) const {
  PRECONDITION(atom, "bad atom pointer");
  return index(atom->getAtomicNum(), atom->getDegree(), atom->getFormalCharge(),
               MolOps::atomHasConjugatedBond(atom), atom->getTotalNumHs(true));
}

}  // namespace DASH
}  // namespace RDKit
