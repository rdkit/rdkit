//
//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
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

#include "../Descriptor.hpp"

namespace RDKit {
namespace NewCIPLabelling {

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

private:
  std::vector<Descriptor> descriptors;
  std::uint32_t pairing = 0;

public:
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
    addAll(head.descriptors);
    addAll(tail.descriptors);
  }

  Descriptor getRefDescriptor() const { return ref(descriptors[0]); }

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
  template <typename T> void addAll(const T &descriptors) {
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
  std::uint32_t getPairing() const { return pairing; }

  /**
   * Appends multiple descriptor lists. If more then one list is provided the
   * head (this list) is duplicate across the multiple tails (provided). If
   * the contents of this list is 'RRSS' and we invoke append with two lists
   * 'SRS' and 'RSR'. Two new lists will be returned with their contents
   * 'RRSSSRS' and 'RRSSRSR' respectively.
   * <br>
   * Empty descriptor lists are not appended, if all descriptor lists are
   * empty then 'this' list is the single returned list
   *
   * @param lists multiple descriptor lists to be appended to this list.
   * @return modified list of descriptors based on the provided input lists
   */
  template <typename T> std::vector<PairList> append(const T &lists) const {
    auto created = std::vector<PairList>();
    created.reserve(lists.size());

    for (const auto &list : lists) {
      // tail isn't empty  - create a new list with this list as the head
      if (!list.descriptors.empty()) {
        created.emplace_back(*this, list);
      }
    }

    // no modifications - make sure we maintain this descriptor list
    if (created.empty()) {
      created.push_back(*this);
    }

    return created;
  }

  int compareTo(const PairList &that) const {
    if (descriptors.size() != that.descriptors.size()) {
      throw std::runtime_error("Descriptor lists should be the same length!");
    }
    Descriptor thisRef = this->descriptors[0];
    Descriptor thatRef = that.descriptors[0];
    for (auto i = 1u; i < this->descriptors.size(); ++i) {
      if (thisRef == this->descriptors[i] && thatRef != that.descriptors[i]) {
        return +1;
      }
      if (thisRef != this->descriptors[i] && thatRef == that.descriptors[i]) {
        return -1;
      }
    }
    return 0;
  }

  bool operator<(const PairList &that) const {
    return this->compareTo(that) == -1;
  }

  /**
   * Clear the descriptor list and resets the pair value. The ignore list is
   * not cleared.
   */
  void clear() {
    pairing = 0;
    descriptors.clear();
  }

  std::string toString() const {
    // handles cases that would break the toString method
    if (descriptors.empty() || descriptors[0] == Descriptor::NONE) {
      return "";
    }

    std::stringstream ss;
    auto basis = descriptors[0];
    ss << to_string(basis) << ':';

    basis = ref(basis);

    // build like (l) / unlike (u) descriptor pairing
    for (auto it = descriptors.begin() + 1; it != descriptors.end(); ++it) {
      ss << (basis == ref(*it) ? "l" : "u");
    }

    return ss.str();
  }

private:
  /**
   * Adds the descriptor to the descriptor list and stores the pair in an set
   * bit (32-bit integer).
   *
   * @param descriptor the descriptor to add an pair
   * @return whether the descriptor was added
   */
  void addAndPair(Descriptor descriptor) {
    // if this isn't the first descriptor - check the pairing
    if (!descriptors.empty() && descriptors[0] == descriptor) {
      // set the bit to indicate a pair
      pairing |= 0x1 << (31 - descriptors.size());
    }
    descriptors.push_back(ref(descriptor));
  }
};

} // namespace NewCIPLabelling
} // namespace RDKit
