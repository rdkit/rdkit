//
//  Copyright (C) 2004-2022 Greg Landrum and other RDKit contributors
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
#include <vector>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/dynamic_bitset.hpp>
#ifdef RDK_USE_URF
#include <boost/shared_ptr.hpp>
#endif
#include <RDGeneral/BoostEndInclude.h>
#ifdef RDK_USE_URF
#include <RingDecomposerLib.h>
#endif

namespace RDKit {
//! A class to store information about a molecule's rings
/*!

 */
class RDKIT_GRAPHMOL_EXPORT RingInfo {
  friend class MolPickler;

 public:
  typedef std::vector<int> MemberType;
  typedef std::vector<MemberType> DataType;
  typedef std::vector<int> INT_VECT;
  typedef std::vector<INT_VECT> VECT_INT_VECT;

  RingInfo() {}
  RingInfo(const RingInfo &other) = default;
  RingInfo &operator=(const RingInfo &other) = default;
  RingInfo(RingInfo &&other) noexcept = default;
  RingInfo &operator=(RingInfo &&other) noexcept = default;
  //! checks to see if we've been properly initialized
  bool isInitialized() const { return df_init; }
  //! does initialization
  void initialize();

  //! blows out all current data and de-initializes
  void reset();

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
  const VECT_INT_VECT &atomRings() const { return d_atomRings; }

  //! returns our \c atom-members vector for atom idx (i.e.,
  //! a vector of ints reporting the ring indices that
  //! atom idx is member of), or an empty vector if the atom is
  //! not in any ring.
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  const INT_VECT &atomMembers(unsigned int idx) const;

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
  const VECT_INT_VECT &bondRings() const { return d_bondRings; }

  //! returns our \c bond-members vector for bond idx (i.e.,
  //! a vector of ints reporting the ring indices that
  //! bond idx is member of), or an empty vector if the bond is
  //! not in any ring.
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  const INT_VECT &bondMembers(unsigned int idx) const;

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

#ifdef RDK_USE_URF
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

  //! returns our bond ring family vectors
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  const VECT_INT_VECT &bondRingFamilies() const { return d_bondRingFamilies; }

  //! check if the ring families have been initialized
  bool areRingFamiliesInitialized() const { return dp_urfData != nullptr; }
#endif

  //! @}

 private:
  //! pre-allocates some memory to save time later
  void preallocate(unsigned int numAtoms, unsigned int numBonds);
  void initFusedRings();
  bool df_init{false};
  DataType d_atomMembers, d_bondMembers;
  VECT_INT_VECT d_atomRings, d_bondRings;
  VECT_INT_VECT d_atomRingFamilies, d_bondRingFamilies;
  std::vector<boost::dynamic_bitset<>> d_fusedRings;
  std::vector<unsigned int> d_numFusedBonds;

#ifdef RDK_USE_URF
 public:
  boost::shared_ptr<RDL_data> dp_urfData;
#endif
};
}  // namespace RDKit

#endif
