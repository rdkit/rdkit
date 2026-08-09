//
//  Copyright (C) 2004-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_RINGINFO_H
#define RD_RINGINFO_H

#include <map>
#include <cstdint>
#include <iterator>
#include <span>
#include <stdexcept>
#include <vector>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/dynamic_bitset.hpp>
#include <boost/shared_ptr.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <RingDecomposerLib.h>

namespace RDKit {
//! A class to store information about a molecule's rings
/*!

 */
typedef enum {
  FIND_RING_TYPE_FAST,
  FIND_RING_TYPE_SSSR,
  FIND_RING_TYPE_SYMM_SSSR,
  FIND_RING_TYPE_OTHER_OR_UNKNOWN
} FIND_RING_TYPE;

class RDKIT_GRAPHMOL_EXPORT RingInfo {
  friend class MolPickler;

 public:
  typedef std::vector<int> MemberType;
  typedef std::vector<MemberType> DataType;
  typedef std::vector<int> INT_VECT;
  typedef std::vector<INT_VECT> VECT_INT_VECT;

  //! A read-only, non-owning view of a collection of rings.
  /*!
    Each element is a `std::span<const int>` containing the atom or bond
    indices for one ring. The view supports indexed access and forward
    iteration without allocating or copying the underlying ring data.

    A `RingsView` does not extend the lifetime of its backing vectors. Views,
    iterators, and spans obtained from a view may be invalidated by an
    operation that modifies the originating `RingInfo`, and must not outlive
    it.

    The `values` vector used to construct a view stores all ring members
    contiguously. The `begins` vector stores the offset of each ring and a
    final end offset, so a view containing N rings has N + 1 offsets.
  */
  class RingsView {
   public:
    //! A forward iterator whose value is the span for one ring.
    class const_iterator {
     public:
      using iterator_category = std::forward_iterator_tag;
      using iterator_concept = std::forward_iterator_tag;
      using value_type = std::span<const int>;
      using difference_type = std::ptrdiff_t;
      using reference = value_type;

      const_iterator() = default;
      //! Constructs an iterator at `idx` in the supplied ring storage.
      /*!
        `values` and `begins` must describe the same valid ring collection and
        must remain alive and unchanged while the iterator is used. `idx` may
        equal the number of rings to construct the end iterator.
      */
      const_iterator(const int *values, const uint32_t *begin)
          : dp_values(values), dp_begin(begin) {}
      reference operator*() const {
        const auto offset = dp_begin[0];
        const auto size = dp_begin[1] - offset;
        return {size ? dp_values + offset : nullptr, size};
      }
      const_iterator &operator++() {
        ++dp_begin;
        return *this;
      }
      const_iterator operator++(int) {
        auto res = *this;
        ++*this;
        return res;
      }
      friend bool operator==(const const_iterator &lhs,
                             const const_iterator &rhs) {
        return lhs.dp_begin == rhs.dp_begin;
      }

     private:
      const int *dp_values{nullptr};
      const uint32_t *dp_begin{nullptr};
    };

    //! Constructs an empty view with no backing storage.
    RingsView() = default;
    //! Constructs a view over flattened ring values and their offsets.
    /*!
      `values` and `begins` must remain alive and unchanged while the view is
      used. `begins` must contain a leading zero, monotonically increasing
      offsets into `values`, and a final offset equal to `values->size()`.
    */
    RingsView(const std::vector<int> *values,
              const std::vector<uint32_t> *begins)
        : dp_values(values ? values->data() : nullptr),
          dp_begins(begins && !begins->empty() ? begins->data() : nullptr),
          d_size(begins && !begins->empty() ? begins->size() - 1 : 0) {}
    size_t size() const { return d_size; }
    bool empty() const { return d_size == 0; }
    std::span<const int> operator[](size_t idx) const {
      const auto begin = dp_begins[idx];
      const auto size = dp_begins[idx + 1] - begin;
      return {size ? dp_values + begin : nullptr, size};
    }
    std::span<const int> at(size_t idx) const {
      if (idx >= size()) {
        throw std::out_of_range("RingInfo::RingsView index out of range");
      }
      return (*this)[idx];
    }
    std::span<const int> front() const { return (*this)[0]; }
    std::span<const int> back() const { return (*this)[size() - 1]; }
    const_iterator begin() const { return {dp_values, dp_begins}; }
    const_iterator end() const {
      return {dp_values, dp_begins ? dp_begins + d_size : nullptr};
    }
    const_iterator cbegin() const { return begin(); }
    const_iterator cend() const { return end(); }

   private:
    const int *dp_values{nullptr};
    const uint32_t *dp_begins{nullptr};
    size_t d_size{0};
  };

  RingInfo() {}
  RingInfo(const RingInfo &other) = default;
  RingInfo &operator=(const RingInfo &other) = default;
  RingInfo(RingInfo &&other) noexcept = default;
  RingInfo &operator=(RingInfo &&other) noexcept = default;
  //! checks to see if we've been properly initialized
  bool isInitialized() const { return df_init; }
  //! does initialization
  void initialize(
      RDKit::FIND_RING_TYPE ringType = FIND_RING_TYPE_OTHER_OR_UNKNOWN);
  RDKit::FIND_RING_TYPE getRingType() const { return df_find_type_type; };
  //! blows out all current data and de-initializes
  void reset(bool resetRingFamilies = true);

  bool isFindFastOrBetter() const {
    return df_init && (df_find_type_type == FIND_RING_TYPE_FAST ||
                       df_find_type_type == FIND_RING_TYPE_SSSR ||
                       df_find_type_type == FIND_RING_TYPE_SYMM_SSSR);
  }

  bool isSssrOrBetter() const {
    return df_init && (df_find_type_type == FIND_RING_TYPE_SSSR ||
                       df_find_type_type == FIND_RING_TYPE_SYMM_SSSR);
  }

  bool isSymmSssr() const {
    return df_init && df_find_type_type == FIND_RING_TYPE_SYMM_SSSR;
  }

  //! adds a ring to our data
  /*!
    \param atomIndices the integer indices of the atoms involved in the ring
    \param bondIndices the integer indices of the bonds involved in the ring,
      this must be the same size as \c atomIndices.

    \return the number of rings

    <b>Notes:</b>
      - the object must be initialized before calling this

  */
  unsigned int addRing(std::span<const int> atomIndices,
                       std::span<const int> bondIndices);
  //! Adds multiple rings while building the membership index only once.
  unsigned int addRings(const VECT_INT_VECT &atomRings,
                        const VECT_INT_VECT &bondRings);
  //! Copies selected rings from another RingInfo and rebuilds membership once.
  unsigned int addRings(const RingInfo &source,
                        const std::vector<unsigned int> &ringIndices);

  //! \name Atom information
  //! @{

  //! returns a vector with sizes of the rings that atom with index \c idx is
  //! in.
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  INT_VECT atomRingSizes(unsigned int idx) const;
  //! returns whether or not the atom with index \c idx is in a \c size - ring.
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  bool isAtomInRingOfSize(unsigned int idx, unsigned int size) const;
  //! returns the number of rings atom \c idx is involved in
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  unsigned int numAtomRings(unsigned int idx) const;
  //! returns the size of the smallest ring atom \c idx is involved in
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  unsigned int minAtomRingSize(unsigned int idx) const;

  //! returns a read-only view of the atom indices in each ring
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  RingsView atomRings() const { return {&d_atomsInRings, &d_atomRingBegins}; }

  //! returns a read-only span of the rings containing atom idx, or an empty
  //! span if the atom is not in any ring.
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  std::span<const int> atomMembers(unsigned int idx) const;

  //! returns whether or not atoms with indices \c idx1 and \c idx2 belong to
  //! the same ring.
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  bool areAtomsInSameRing(unsigned int idx1, unsigned int idx2) const {
    return areAtomsInSameRingOfSize(idx1, idx2, 0);
  }

  //! returns whether or not atoms with indices \c idx1 and \c idx2 belong to
  //! the same ring of size \c size.
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  bool areAtomsInSameRingOfSize(unsigned int idx1, unsigned int idx2,
                                unsigned int size) const;

  //! @}

  //! \name Bond information
  //! @{

  //! returns a vector with sizes of the rings that bond with index \c idx is
  //! in.
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  INT_VECT bondRingSizes(unsigned int idx) const;
  //! returns whether or not the bond with index \c idx is in a \c size - ring.
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  bool isBondInRingOfSize(unsigned int idx, unsigned int size) const;
  //! returns the number of rings bond \c idx is involved in
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  unsigned int numBondRings(unsigned int idx) const;
  //! returns the size of the smallest ring bond \c idx is involved in
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  unsigned int minBondRingSize(unsigned int idx) const;

  //! returns the total number of rings
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
      - if the RDKit has been built with URF support, this returns the number
        of ring families.
  */
  unsigned int numRings() const;

  //! returns a read-only view of the bond indices in each ring
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  RingsView bondRings() const { return {&d_bondsInRings, &d_bondRingBegins}; }

  //! returns a read-only span of the rings containing bond idx, or an empty
  //! span if the bond is not in any ring.
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  std::span<const int> bondMembers(unsigned int idx) const;

  //! returns whether or not bonds with indices \c idx1 and \c idx2 belong to
  //! the same ring.
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  bool areBondsInSameRing(unsigned int idx1, unsigned int idx2) const {
    return areBondsInSameRingOfSize(idx1, idx2, 0);
  }

  //! returns whether or not bonds with indices \c idx1 and \c idx2 belong to
  //! the same ring of size \c size.
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  bool areBondsInSameRingOfSize(unsigned int idx1, unsigned int idx2,
                                unsigned int size) const;

  //! returns whether ring with index \c ringIdx is fused with other rings.
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  bool isRingFused(unsigned int ringIdx);

  //! returns whether rings with indices \c ring1Idx and \c ring2Idx have
  //! at least one bond in common.
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  bool areRingsFused(unsigned int ring1Idx, unsigned int ring2Idx);

  //! returns the number of bonds shared with other rings in ring with index
  //! \c ringIdx.
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  unsigned int numFusedBonds(unsigned int ringIdx);

  //! returns the number of rings which have at least one bond
  //! in common with ring with index \c ringIdx.
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  unsigned int numFusedRingNeighbors(unsigned int ringIdx);

  //! returns the indices of rings which have at least one bond
  //! in common with ring with index \c ringIdx.
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  std::vector<unsigned int> fusedRingNeighbors(unsigned int ringIdx);

  //! adds a ring family to our data
  /*!
    \param atomIndices the integer indices of the atoms involved in the
                       ring family
    \param bondIndices the integer indices of the bonds involved in the
                       ring family,
      this must be the same size as \c atomIndices.

    \return the number of ring families

    <b>Notes:</b>
      - the object must be initialized before calling this

  */
  unsigned int addRingFamily(const INT_VECT &atomIndices,
                             const INT_VECT &bondIndices);
  //! returns the total number of ring families
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  unsigned int numRingFamilies() const;

  //! returns the total number of relevant cycles
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  unsigned int numRelevantCycles() const;

  //! returns our atom ring family vectors
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  const VECT_INT_VECT &atomRingFamilies() const { return d_atomRingFamilies; }
  VECT_INT_VECT atomRelevantCycles() const;

  //! returns our bond ring family vectors
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  const VECT_INT_VECT &bondRingFamilies() const { return d_bondRingFamilies; }

  //! check if the ring families have been initialized
  bool areRingFamiliesInitialized() const { return dp_urfData != nullptr; }

  //! reset ring family information
  void resetRingFamilies();

  //! @}

  //! pre-allocates some memory to save time later
  void preallocate(unsigned int numAtoms, unsigned int numBonds);

 private:
  std::span<const int> atomMembersUnchecked(unsigned int idx) const {
    if (d_atomMembershipBegins.empty() ||
        idx >= d_atomMembershipBegins.size() - 1) {
      return {};
    }
    const auto begin = d_atomMembershipBegins[idx];
    const auto size = d_atomMembershipBegins[idx + 1] - begin;
    return {size ? d_atomMemberships.data() + begin : nullptr, size};
  }
  std::span<const int> bondMembersUnchecked(unsigned int idx) const {
    if (d_bondMembershipBegins.empty() ||
        idx >= d_bondMembershipBegins.size() - 1) {
      return {};
    }
    const auto begin = d_bondMembershipBegins[idx];
    const auto size = d_bondMembershipBegins[idx + 1] - begin;
    return {size ? d_bondMemberships.data() + begin : nullptr, size};
  }
  void appendRingUnchecked(std::span<const int> atomIndices,
                           std::span<const int> bondIndices);
  unsigned int finishRingUpdates();
  void initFusedRings();
  void rebuildMemberships();
  void invalidateFusedRings();
  bool df_init{false};
  FIND_RING_TYPE df_find_type_type{FIND_RING_TYPE_OTHER_OR_UNKNOWN};
  std::vector<uint32_t> d_atomRingBegins{0};
  std::vector<uint32_t> d_atomMembershipBegins{0};
  std::vector<int> d_atomMemberships;
  std::vector<uint32_t> d_bondRingBegins{0};
  std::vector<uint32_t> d_bondMembershipBegins{0};
  std::vector<int> d_bondMemberships;
  std::vector<int> d_atomsInRings;
  std::vector<int> d_bondsInRings;
  VECT_INT_VECT d_atomRingFamilies, d_bondRingFamilies;
  std::vector<bool> d_fusedRings;
  std::vector<unsigned int> d_numFusedBonds;
  bool d_fusedRingsInitialized{false};

 public:
  boost::shared_ptr<RDL_data> dp_urfData;
};

}  // namespace RDKit

#endif
