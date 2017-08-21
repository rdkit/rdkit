//
//  Copyright (C) 2001-2017 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_BOND_H
#define _RD_BOND_H

// std stuff
#include <iostream>

// Ours
#include <RDGeneral/Invariant.h>
#include <Query/QueryObjects.h>
#include <RDGeneral/types.h>
#include <RDGeneral/RDProps.h>
#include <GraphMol/details.h>
#include <boost/foreach.hpp>

namespace RDKit {
class ROMol;
class RWMol;
class Atom;
typedef boost::shared_ptr<Atom> ATOM_SPTR;

//! class for representing a bond
/*!

  <b>Notes:</b>
    - many of the methods of Atom require that the Atom be associated
      with a molecule (an ROMol).
    - each Bond maintains a Dict of \c properties:
        - Each \c property is keyed by name and can store an
          arbitrary type.
        - \c Properties can be marked as \c calculated, in which case
          they will be cleared when the \c clearComputedProps() method
          is called.
        - Because they have no impact upon chemistry, all \c property
          operations are \c const, this allows extra flexibility for
          clients who need to store extra data on Bond objects.

*/
class Bond : public RDProps {
  friend class RWMol;
  friend class ROMol;

 public:
  typedef boost::shared_ptr<Bond> BOND_SPTR;
  // FIX: grn...
  typedef Queries::Query<int, Bond const *, true> QUERYBOND_QUERY;

  //! the type of Bond
  typedef enum {
    UNSPECIFIED = 0,
    SINGLE,
    DOUBLE,
    TRIPLE,
    QUADRUPLE,
    QUINTUPLE,
    HEXTUPLE,
    ONEANDAHALF,
    TWOANDAHALF,
    THREEANDAHALF,
    FOURANDAHALF,
    FIVEANDAHALF,
    AROMATIC,
    IONIC,
    HYDROGEN,
    THREECENTER,
    DATIVEONE,  //!< one-electron dative (e.g. from a C in a Cp ring to a metal)
    DATIVE,     //!< standard two-electron dative
    DATIVEL,    //!< standard two-electron dative
    DATIVER,    //!< standard two-electron dative
    OTHER,
    ZERO  //!< Zero-order bond (from
    // http://pubs.acs.org/doi/abs/10.1021/ci200488k)
  } BondType;

  //! the bond's direction (for chirality)
  typedef enum {
    NONE = 0,    //!< no special style
    BEGINWEDGE,  //!< wedged: narrow at begin
    BEGINDASH,   //!< dashed: narrow at begin
    // FIX: this may not really be adequate
    ENDDOWNRIGHT,  //!< for cis/trans
    ENDUPRIGHT,    //!<  ditto
    EITHERDOUBLE,  //!< a "crossed" double bond
    UNKNOWN,       //!< intentionally unspecified stereochemistry
  } BondDir;

  //! the nature of the bond's stereochem (for cis/trans)
  typedef enum {     // stereochemistry of double bonds
    STEREONONE = 0,  // no special style
    STEREOANY,       // intentionally unspecified
    // -- Put any true specifications about this point so
    // that we can do comparisons like if(bond->getStereo()>Bond::STEREOANY)
    STEREOZ,     // Z double bond
    STEREOE,     // E double bond
    STEREOCIS,   // cis double bond
    STEREOTRANS  // trans double bond
  } BondStereo;

  Bond();
  //! construct with a particular BondType
  explicit Bond(BondType bT);
  Bond(const Bond &other);
  virtual ~Bond();
  Bond &operator=(const Bond &other);

  //! returns a copy
  /*!
    <b>Note:</b> the caller is responsible for <tt>delete</tt>ing
     the returned pointer.
  */
  virtual Bond *copy() const;

  //! returns our \c bondType
  BondType getBondType() const { return static_cast<BondType>(d_bondType); };
  //! sets our \c bondType
  void setBondType(BondType bT) { d_bondType = bT; };
  //! \brief returns our \c bondType as a double
  //!   (e.g. SINGLE->1.0, AROMATIC->1.5, etc.)
  double getBondTypeAsDouble() const;

  //! returns our contribution to the explicit valence of an Atom
  /*!
    <b>Notes:</b>
      - requires an owning molecule
  */
  double getValenceContrib(const Atom *at) const;
  // \overload
  double getValenceContrib(ATOM_SPTR at) const;

  //! sets our \c isAromatic flag
  void setIsAromatic(bool what) { df_isAromatic = what; };
  //! returns the status of our \c isAromatic flag
  bool getIsAromatic() const { return df_isAromatic; };

  //! sets our \c isConjugated flag
  void setIsConjugated(bool what) { df_isConjugated = what; };
  //! returns the status of our \c isConjugated flag
  bool getIsConjugated() const { return df_isConjugated; };

  //! returns a reference to the ROMol that owns this Bond
  ROMol &getOwningMol() const {
    PRECONDITION(dp_mol, "no owner");
    return *dp_mol;
  };
  //! sets our owning molecule
  void setOwningMol(ROMol *other);
  //! sets our owning molecule
  void setOwningMol(ROMol &other) { setOwningMol(&other); };

  //! returns our index within the ROMol
  /*!
    <b>Notes:</b>
      - this makes no sense if we do not have an owning molecule

  */
  unsigned int getIdx() const { return d_index; };
  //! sets our index within the ROMol
  /*!
    <b>Notes:</b>
      - this makes no sense if we do not have an owning molecule
      - the index should be <tt>< this->getOwningMol()->getNumBonds()</tt>
  */
  void setIdx(unsigned int index) { d_index = index; };

  //! returns the index of our begin Atom
  /*!
    <b>Notes:</b>
      - this makes no sense if we do not have an owning molecule
  */
  unsigned int getBeginAtomIdx() const { return d_beginAtomIdx; };

  //! returns the index of our end Atom
  /*!
    <b>Notes:</b>
      - this makes no sense if we do not have an owning molecule
  */
  unsigned int getEndAtomIdx() const { return d_endAtomIdx; };

  //! given the index of one Atom, returns the index of the other
  /*!
    <b>Notes:</b>
      - this makes no sense if we do not have an owning molecule
  */
  unsigned int getOtherAtomIdx(unsigned int thisIdx) const;

  //! sets the index of our begin Atom
  /*!
    <b>Notes:</b>
      - requires an owning molecule
  */
  void setBeginAtomIdx(unsigned int what);
  //! sets the index of our end Atom
  /*!
    <b>Notes:</b>
      - requires an owning molecule
  */
  void setEndAtomIdx(unsigned int what);

  //! sets our begin Atom
  /*!
    <b>Notes:</b>
      - requires an owning molecule
  */
  void setBeginAtom(Atom *at);
  //! \overload
  void setBeginAtom(ATOM_SPTR at);
  //! sets our end Atom
  /*!
    <b>Notes:</b>
      - requires an owning molecule
  */
  void setEndAtom(Atom *at);
  //! \overload
  void setEndAtom(ATOM_SPTR at);

  //! returns a pointer to our begin Atom
  /*!
    <b>Notes:</b>
      - requires an owning molecule
  */
  Atom *getBeginAtom() const;
  //! returns a pointer to our end Atom
  /*!
    <b>Notes:</b>
      - requires an owning molecule
  */
  Atom *getEndAtom() const;
  //! returns a pointer to the other Atom
  /*!
    <b>Notes:</b>
      - requires an owning molecule
  */
  Atom *getOtherAtom(Atom const *what) const;

  // ------------------------------------
  // Please see the note in Atom.h for some explanation
  // of these methods
  // ------------------------------------

  // This method can be used to distinguish query bonds from standard bonds
  virtual bool hasQuery() const { return false; };

  // FIX: the const crap here is all mucked up.
  //! NOT CALLABLE
  virtual void setQuery(QUERYBOND_QUERY *what);
  //! NOT CALLABLE
  virtual QUERYBOND_QUERY *getQuery() const;

  //! NOT CALLABLE
  virtual void expandQuery(
      QUERYBOND_QUERY *what,
      Queries::CompositeQueryType how = Queries::COMPOSITE_AND,
      bool maintainOrder = true);

  //! returns whether or not we match the argument
  /*!
      <b>Notes:</b>
        - for Bond objects, "match" means that either one of the Bonds
          has \c bondType Bond::UNSPECIFIED or both Bonds have the
          same \c bondType.
  */
  virtual bool Match(Bond const *what) const;
  //! \overload
  virtual bool Match(const Bond::BOND_SPTR what) const;

  //! sets our direction
  void setBondDir(BondDir what) { d_dirTag = what; };
  //! returns our direction
  BondDir getBondDir() const { return static_cast<BondDir>(d_dirTag); };

  //! sets our stereo code
  /*!
      STEREONONE, STEREOANY, STEREOE and STEREOZ can be set without
      neighboring atoms specified in getStereoAtoms since they are
      defined by the topology of the molecular graph. In order to set
      STEREOCIS or STEREOTRANS the neighboring atoms must be set first
      (using setStereoBonds()) to know what atoms are being considered.

      <b>Notes:</b>
        - MolOps::findPotentialStereoBonds can be used to set
          getStereoAtoms before setting CIS/TRANS
  */
  void setStereo(BondStereo what) {
    PRECONDITION(what <= STEREOE || getStereoAtoms().size() == 2,
                 "Stereo atoms should be specified before specifying CIS/TRANS "
                 "bond stereochemistry")
    d_stereo = what;
  };
  //! returns our stereo code
  BondStereo getStereo() const { return static_cast<BondStereo>(d_stereo); };

  //! sets the atoms to be considered as reference points for bond stereo
  /*!
      These do not necessarily need to be the highest 'ranking' atoms
      like CIP stereo requires. They can be any arbitrary atoms
      neighboring the begin and end atoms of this bond
      respectively. STEREOCIS or STEREOTRANS is then set relative to
      only these atoms.

      If CIP rankings are desired, use
      MolOps::findPotentialStereoBonds, but this is a more costly
      function as it takes the whole molecule topology into account.
  */
  void setStereoAtoms(unsigned int bgnIdx, unsigned int endIdx);

  //! returns the indices of our stereo atoms
  const INT_VECT &getStereoAtoms() const {
    if (!dp_stereoAtoms) {
      const_cast<Bond *>(this)->dp_stereoAtoms = new INT_VECT();
    }
    return *dp_stereoAtoms;
  };
  //! \overload
  INT_VECT &getStereoAtoms() {
    if (!dp_stereoAtoms) dp_stereoAtoms = new INT_VECT();
    return *dp_stereoAtoms;
  };

  //! calculates any of our lazy \c properties
  /*!
    <b>Notes:</b>
      - requires an owning molecule
  */
  void updatePropertyCache(bool strict = true) { (void)strict; }

 protected:
  //! sets our owning molecule
  // void setOwningMol(ROMol *other);
  //! sets our owning molecule
  // void setOwningMol(ROMol &other) {setOwningMol(&other);};
  bool df_isAromatic;
  bool df_isConjugated;
  boost::uint8_t d_bondType;
  boost::uint8_t d_dirTag;
  boost::uint8_t d_stereo;
  atomindex_t d_index;
  atomindex_t d_beginAtomIdx, d_endAtomIdx;
  ROMol *dp_mol;
  INT_VECT *dp_stereoAtoms;

  void initBond();
};
};

//! allows Bond objects to be dumped to streams
extern std::ostream &operator<<(std::ostream &target, const RDKit::Bond &b);

#endif
