//
//  Copyright (C) 2004-2019 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef _RD_RINGINFO_H
#define _RD_RINGINFO_H

#include <map>
#include <vector>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/shared_ptr.hpp>
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

  RingInfo()  {};
  RingInfo(const RingInfo &other)
      : df_init(other.df_init),
        d_atomMembers(other.d_atomMembers),
        d_bondMembers(other.d_bondMembers),
        d_atomRings(other.d_atomRings),
        d_bondRings(other.d_bondRings),
        d_atomRingFamilies(other.d_atomRingFamilies),
        d_bondRingFamilies(other.d_bondRingFamilies)
#ifdef RDK_USE_URF
        ,dp_urfData(other.dp_urfData)
#endif
            {};

  //! checks to see if we've been properly initialized
  bool isInitialized() const { return df_init; };
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
  //@{

  //! returns a vector with sizes of the rings that atom with index \c idx is in.
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

  //! returns our \c atom-rings vectors
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  const VECT_INT_VECT &atomRings() const { return d_atomRings; };

  //@}

  //! \name Bond information
  //@{

  //! returns a vector with sizes of the rings that bond with index \c idx is in.
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

  //! returns our \c bond-rings vectors
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  const VECT_INT_VECT &bondRings() const { return d_bondRings; };

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
  const VECT_INT_VECT &atomRingFamilies() const { return d_atomRingFamilies; };

  //! returns our bond ring family vectors
  /*!
    <b>Notes:</b>
      - the object must be initialized before calling this
  */
  const VECT_INT_VECT &bondRingFamilies() const { return d_bondRingFamilies; };

  //! check if the ring families have been initialized
  bool areRingFamiliesInitialized() const { return dp_urfData != nullptr; }
#endif

  //@}

 private:
  //! pre-allocates some memory to save time later
  void preallocate(unsigned int numAtoms, unsigned int numBonds);

  bool df_init{false};
  DataType d_atomMembers, d_bondMembers;
  VECT_INT_VECT d_atomRings, d_bondRings;
  VECT_INT_VECT d_atomRingFamilies, d_bondRingFamilies;

#ifdef RDK_USE_URF
 public:
  boost::shared_ptr<RDL_data> dp_urfData;
#endif
};
}  // namespace RDKit

#endif
