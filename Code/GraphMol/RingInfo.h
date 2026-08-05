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
#include <stdexcept>
#include <compare>
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

  //! A non-owning view of one ring or one atom/bond membership list.
  //! The view is invalidated by non-const changes to this RingInfo.
  class IntView {
   public:
    using value_type = int;
    using const_iterator = const int *;

    IntView() = default;
    IntView(const int *data, size_t size) : dp_data(data), d_size(size) {}
    const_iterator begin() const { return dp_data; }
    const_iterator end() const { return d_size ? dp_data + d_size : dp_data; }
    const_iterator cbegin() const { return begin(); }
    const_iterator cend() const { return end(); }
    size_t size() const { return d_size; }
    bool empty() const { return d_size == 0; }
    const int &operator[](size_t idx) const { return dp_data[idx]; }
    const int &at(size_t idx) const {
      if (idx >= d_size) {
        throw std::out_of_range("RingInfo::IntView index out of range");
      }
      return dp_data[idx];
    }
    const int &front() const { return dp_data[0]; }
    const int &back() const { return dp_data[d_size - 1]; }

   private:
    const int *dp_data{nullptr};
    size_t d_size{0};
  };

  //! A non-owning indexed view of all ordinary atom or bond rings.
  //! Iteration is read-only and views must not be retained across mutation.
  class RingsView {
   public:
    class const_iterator {
     public:
      using iterator_category = std::forward_iterator_tag;
      using iterator_concept = std::forward_iterator_tag;
      using value_type = IntView;
      using difference_type = std::ptrdiff_t;
      using reference = IntView;

      reference operator*() const { return valueAt(d_idx); }
      const IntView *operator->() const {
        d_cachedValue = valueAt(d_idx);
        return &d_cachedValue;
      }
      reference operator[](difference_type offset) const {
        return valueAt(d_idx + offset);
      }
      const_iterator &operator++() {
        ++d_idx;
        return *this;
      }
      const_iterator operator++(int) {
        auto res = *this;
        ++*this;
        return res;
      }
      const_iterator &operator--() {
        --d_idx;
        return *this;
      }
      const_iterator operator--(int) {
        auto res = *this;
        --*this;
        return res;
      }
      const_iterator &operator+=(difference_type offset) {
        d_idx += offset;
        return *this;
      }
      const_iterator &operator-=(difference_type offset) {
        d_idx -= offset;
        return *this;
      }
      friend const_iterator operator+(const_iterator it,
                                      difference_type offset) {
        return it += offset;
      }
      friend const_iterator operator+(difference_type offset,
                                      const_iterator it) {
        return it += offset;
      }
      friend const_iterator operator-(const_iterator it,
                                      difference_type offset) {
        return it -= offset;
      }
      friend difference_type operator-(const const_iterator &lhs,
                                       const const_iterator &rhs) {
        return static_cast<difference_type>(lhs.d_idx) -
               static_cast<difference_type>(rhs.d_idx);
      }
      friend bool operator==(const const_iterator &lhs,
                             const const_iterator &rhs) {
        return lhs.dp_values == rhs.dp_values &&
               lhs.dp_begins == rhs.dp_begins && lhs.d_idx == rhs.d_idx;
      }
      friend auto operator<=>(const const_iterator &lhs,
                              const const_iterator &rhs) {
        return lhs.d_idx <=> rhs.d_idx;
      }

     private:
      friend class RingsView;
      const_iterator(const std::vector<int> *values,
                     const std::vector<uint32_t> *begins, size_t idx)
          : dp_values(values), dp_begins(begins), d_idx(idx) {}
      IntView valueAt(size_t idx) const {
        const auto begin = (*dp_begins)[idx];
        const auto size = (*dp_begins)[idx + 1] - begin;
        return {size ? dp_values->data() + begin : nullptr, size};
      }
      const std::vector<int> *dp_values{nullptr};
      const std::vector<uint32_t> *dp_begins{nullptr};
      size_t d_idx{0};
      mutable IntView d_cachedValue;
    };

    RingsView() = default;
    RingsView(const std::vector<int> *values,
              const std::vector<uint32_t> *begins)
        : dp_values(values), dp_begins(begins) {}
    size_t size() const {
      return dp_begins && !dp_begins->empty() ? dp_begins->size() - 1 : 0;
    }
    bool empty() const { return size() == 0; }
    IntView operator[](size_t idx) const {
      const auto begin = (*dp_begins)[idx];
      const auto size = (*dp_begins)[idx + 1] - begin;
      return {size ? dp_values->data() + begin : nullptr, size};
    }
    IntView at(size_t idx) const {
      if (idx >= size()) {
        throw std::out_of_range("RingInfo::RingsView index out of range");
      }
      return (*this)[idx];
    }
    IntView front() const { return (*this)[0]; }
    IntView back() const { return (*this)[size() - 1]; }
    const_iterator begin() const { return {dp_values, dp_begins, 0}; }
    const_iterator end() const { return {dp_values, dp_begins, size()}; }
    const_iterator cbegin() const { return begin(); }
    const_iterator cend() const { return end(); }

   private:
    const std::vector<int> *dp_values{nullptr};
    const std::vector<uint32_t> *dp_begins{nullptr};
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
  unsigned int addRing(const INT_VECT &atomIndices,
                       const INT_VECT &bondIndices);
  unsigned int addRing(IntView atomIndices, IntView bondIndices);
  //! Adds multiple rings while building the membership index only once.
  unsigned int addRings(const VECT_INT_VECT &atomRings,
                        const VECT_INT_VECT &bondRings);

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

  //! returns our \c atom-rings vectors, i.e. a vector of int vectors
  //! reporting the atom indices which are part of each ring
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  RingsView atomRings() const { return {&d_atomsInRings, &d_atomRingBegins}; }
  //! Materializes ordinary atom rings for APIs which require owning vectors.
  VECT_INT_VECT atomRingsAsVectors() const;

  //! returns our \c atom-members vector for atom idx (i.e.,
  //! a vector of ints reporting the ring indices that
  //! atom idx is member of), or an empty vector if the atom is
  //! not in any ring.
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  IntView atomMembers(unsigned int idx) const;

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

  //! returns our \c bond-rings vectors, i.e. a vector of int vectors
  //! reporting the bond indices which are part of each ring
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  RingsView bondRings() const { return {&d_bondsInRings, &d_bondRingBegins}; }
  //! Materializes ordinary bond rings for APIs which require owning vectors.
  VECT_INT_VECT bondRingsAsVectors() const;

  //! returns our \c bond-members vector for bond idx (i.e.,
  //! a vector of ints reporting the ring indices that
  //! bond idx is member of), or an empty vector if the bond is
  //! not in any ring.
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  IntView bondMembers(unsigned int idx) const;

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
  void initFusedRings();
  void rebuildMemberships();
  void invalidateFusedRings();
  bool df_init{false};
  FIND_RING_TYPE df_find_type_type{FIND_RING_TYPE_OTHER_OR_UNKNOWN};
  std::vector<uint32_t> d_atomRingBegins{0}, d_bondRingBegins{0};
  std::vector<int> d_atomsInRings, d_bondsInRings;
  std::vector<uint32_t> d_atomMembershipBegins{0}, d_bondMembershipBegins{0};
  std::vector<int> d_atomMemberships, d_bondMemberships;
  VECT_INT_VECT d_atomRingFamilies, d_bondRingFamilies;
  std::vector<bool> d_fusedRings;
  std::vector<unsigned int> d_numFusedBonds;
  bool d_fusedRingsInitialized{false};

 public:
  boost::shared_ptr<RDL_data> dp_urfData;
};

}  // namespace RDKit

#endif
