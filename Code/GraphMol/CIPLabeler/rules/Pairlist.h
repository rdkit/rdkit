//
//
//  Copyright (C) 2020 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once

#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

#include "../Descriptor.h"

namespace RDKit {
namespace CIPLabeler {

/**
 * Implementation of a descriptor list that allows descriptors to be added and
 * ignored. The list maintains an integer value throughout which stores the
 * pairing of descriptors and allows easy comparison between descriptor lists in
 * that higher priority descriptor pairing will always have a higher integer
 * value. The integer value can be access via the {@link #getPairing()} method.
 *
 * @see Descriptor
 */
class PairList {
 public:
  static Descriptor ref(Descriptor descriptor) {
    switch (descriptor) {
      case Descriptor::R:
      case Descriptor::M:
      case Descriptor::seqCis:
        return Descriptor::R;
      case Descriptor::S:
      case Descriptor::P:
      case Descriptor::seqTrans:
        return Descriptor::S;
      default:
        return Descriptor::NONE;
    }
  }

  PairList() = default;

  PairList(Descriptor ref) { add(ref); }

  /**
   * Creates a new list from a provided head and tail. The head and tail
   * ignored descriptors are first transferred and then their descriptors. In
   * either list, descriptors that are ignored by the other will be not be
   * added to the new instance.
   *
   * @param head the head of the list (prefix)
   * @param tail the tail of the list (suffix)
   */
  PairList(const PairList &head, const PairList &tail) {
    // add descriptors to the new instance (ignored descriptors not added)
    addAll(head.d_descriptors);
    addAll(tail.d_descriptors);
  }

  Descriptor getRefDescriptor() const { return ref(d_descriptors[0]); }

  /**
   * Adds a descriptor to the descriptor list. If the provided descriptor is
   * present in the ignore set the descriptor will not be added.
   *
   * @param descriptor the descriptor to add.
   * @return whether the descriptor was added to the list
   */

  bool add(Descriptor descriptor) {
    switch (descriptor) {
      case Descriptor::R:
      case Descriptor::S:
      case Descriptor::M:
      case Descriptor::P:
      case Descriptor::seqTrans:
      case Descriptor::seqCis:
        addAndPair(descriptor);
        return true;
      default:
        return false;
    }
  }

  /**
   * Adds multiple descriptors to the descriptor list. If the descriptor is
   * present in the ignore set it will not be added to the list.
   *
   * @param descriptors a collection of descriptors to be added
   */
  template <typename T>
  void addAll(const T &descriptors) {
    for (const auto &descriptor : descriptors) {
      add(descriptor);
    }
  }

  /**
   * Access a positive integer that represents the like/unlike pairings of
   * this descriptor list. The like/unlike is represented by set bits in an
   * integer value and means larger integer values indicates a higher
   * descriptor pairing preference.
   *
   * @return an integer representing the descriptor pairings
   */
  std::uint32_t getPairing() const { return d_pairing; }

  int compareTo(const PairList &that) const {
    if (d_descriptors.size() != that.d_descriptors.size()) {
      throw std::runtime_error("Descriptor lists should be the same length!");
    }
    Descriptor thisRef = d_descriptors[0];
    Descriptor thatRef = that.d_descriptors[0];
    for (auto i = 1u; i < d_descriptors.size(); ++i) {
      if (thisRef == d_descriptors[i] && thatRef != that.d_descriptors[i]) {
        return +1;
      }
      if (thisRef != d_descriptors[i] && thatRef == that.d_descriptors[i]) {
        return -1;
      }
    }
    return 0;
  }

  bool operator<(const PairList &that) const { return compareTo(that) == -1; }

  std::string toString() const {
    // handles cases that would break the toString method
    if (d_descriptors.empty() || d_descriptors[0] == Descriptor::NONE) {
      return "";
    }

    std::stringstream ss;
    auto basis = d_descriptors[0];
    ss << to_string(basis) << ':';

    basis = ref(basis);

    // build like (l) / unlike (u) descriptor pairing
    for (auto it = d_descriptors.begin() + 1; it != d_descriptors.end(); ++it) {
      ss << (basis == ref(*it) ? "l" : "u");
    }

    return ss.str();
  }

 private:
  std::vector<Descriptor> d_descriptors;

  std::uint32_t d_pairing = 0;

  /**
   * Adds the descriptor to the descriptor list and stores the pair in an set
   * bit (32-bit integer).
   *
   * @param descriptor the descriptor to add an pair
   * @return whether the descriptor was added
   */
  void addAndPair(Descriptor descriptor) {
    // if this isn't the first descriptor - check the pairing
    if (!d_descriptors.empty() && d_descriptors[0] == descriptor) {
      // set the bit to indicate a pair
      d_pairing |= 0x1 << (31 - d_descriptors.size());
    }
    d_descriptors.push_back(ref(descriptor));
  }
};

}  // namespace CIPLabeler
}  // namespace RDKit
