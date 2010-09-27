//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_RINGINFO_H
#define _RD_RINGINFO_H

#include <map>
#include <vector>

namespace RDKit {
  //! A class to store information about a molecule's rings
  /*!

   */
  class RingInfo {
    friend class MolPickler;
  public:
    typedef std::vector<int> MemberType;
    typedef std::vector<MemberType > DataType;
    typedef std::vector<int> INT_VECT;
    typedef std::vector< INT_VECT > VECT_INT_VECT;

    RingInfo() : df_init(false) {};
    RingInfo(const RingInfo &other) : df_init(other.df_init),
				      d_atomMembers(other.d_atomMembers),
				      d_bondMembers(other.d_bondMembers),
				      d_atomRings(other.d_atomRings),
				      d_bondRings(other.d_bondRings) {};
    
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
    unsigned int addRing(const INT_VECT &atomIndices,const INT_VECT &bondIndices);

  
    //! \name Atom information
    //@{

    //! returns whether or not the atom with index \c idx is in a \c size - ring.
    /*!
      <b>Notes:</b>
        - the object must be initialized before calling this
    */
    bool isAtomInRingOfSize(unsigned int idx,unsigned int size) const;
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

    //! returns whether or not the bond with index \c idx is in a \c size - ring.
    /*!
      <b>Notes:</b>
        - the object must be initialized before calling this
    */
    bool isBondInRingOfSize(unsigned int idx,unsigned int size) const;
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
    */
    unsigned int numRings() const;

    //! returns our \c bond-rings vectors
    /*!
      <b>Notes:</b>
        - the object must be initialized before calling this
    */
    const VECT_INT_VECT &bondRings() const { return d_bondRings; };

    //@}
    
  private:
    //! pre-allocates some memory to save time later
    void preallocate(unsigned int numAtoms,unsigned int numBonds);

    bool df_init;
    DataType d_atomMembers,d_bondMembers;
    VECT_INT_VECT d_atomRings,d_bondRings;
  };
}

#endif
