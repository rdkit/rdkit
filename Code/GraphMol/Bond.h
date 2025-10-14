//
//  Copyright (C) 2001-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_BOND_H
#define RD_BOND_H

// std stuff
#include <utility>

// Ours
#include <RDGeneral/Invariant.h>
#include <Query/QueryObjects.h>
#include <RDGeneral/types.h>
#include <RDGeneral/RDProps.h>
#include <GraphMol/details.h>
#include <GraphMol/RDMol.h>
#include <GraphMol/rdmol_throw.h>

#include "GraphMolEnums.h"

namespace RDKit {
class ROMol;
class RWMol;
class Atom;
class ConstRDMolBond;

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
class RDKIT_GRAPHMOL_EXPORT Bond {
  friend class ROMol;
  friend class RWMol;
  friend class RDMol;

  // When a bond is created on its own, it owns dp_dataMol to store its data
  // and dp_owningMol is nullptr. Calling setOwningMol will set dp_owningMol,
  // but as before, this does not actually insert the Bond's data into the
  // owning RDMol, so dp_dataMol will continue keeping the data until
  // ROMol::addBond copies the data and frees dp_dataMol, setting dp_dataMol
  // to be equal to dp_owningMol, to indicate that dp_dataMol is no longer
  // owned by this Bond, and setting d_index accordingly.
  //
  // The 3 states are:
  // 1) dp_owningMol == nullptr: dp_dataMol owned by this & no owning mol
  // 2) else if dp_dataMol != dp_owningMol: dp_dataMol owned by this
  // 3) dp_dataMol == dp_owningMol: dp_owningMol owns this
  //
  // This is necessary for compatibility with previous workflows that would
  // call setOwningMol and separately ROMol::addBond.
  RDMol *dp_dataMol = nullptr;
  RDMol *dp_owningMol = nullptr;
  uint32_t d_index = 0;

 protected:
  Bond(RDMol *mol, uint32_t bondIndex)
      : dp_dataMol(mol), dp_owningMol(mol), d_index(bondIndex) {}

 public:
  // FIX: grn...
  typedef Queries::Query<int, Bond const *, true> QUERYBOND_QUERY;
  //typedef Queries::Query<int, ConstRDMolBond, true> QUERYBOND_QUERY;

  //! the type of Bond
  enum BondType {
    UNSPECIFIED = BondEnums::BondType::UNSPECIFIED,
    SINGLE = BondEnums::BondType::SINGLE,
    DOUBLE = BondEnums::BondType::DOUBLE,
    TRIPLE = BondEnums::BondType::TRIPLE,
    QUADRUPLE = BondEnums::BondType::QUADRUPLE,
    QUINTUPLE = BondEnums::BondType::QUINTUPLE,
    HEXTUPLE = BondEnums::BondType::HEXTUPLE,
    ONEANDAHALF = BondEnums::BondType::ONEANDAHALF,
    TWOANDAHALF = BondEnums::BondType::TWOANDAHALF,
    THREEANDAHALF = BondEnums::BondType::THREEANDAHALF,
    FOURANDAHALF = BondEnums::BondType::FOURANDAHALF,
    FIVEANDAHALF = BondEnums::BondType::FIVEANDAHALF,
    AROMATIC = BondEnums::BondType::AROMATIC,
    IONIC = BondEnums::BondType::IONIC,
    HYDROGEN = BondEnums::BondType::HYDROGEN,
    THREECENTER = BondEnums::BondType::THREECENTER,
    //!< one-electron dative (e.g. from a C in a Cp ring to a metal)
    DATIVEONE = BondEnums::BondType::DATIVEONE,
    DATIVE = BondEnums::BondType::DATIVE,    //!< standard two-electron dative
    DATIVEL = BondEnums::BondType::DATIVEL,  //!< standard two-electron dative
    DATIVER = BondEnums::BondType::DATIVER,  //!< standard two-electron dative
    OTHER = BondEnums::BondType::OTHER,
    ZERO = BondEnums::BondType::ZERO  //!< Zero-order bond (from
    // http://pubs.acs.org/doi/abs/10.1021/ci200488k)
  };

  //! the bond's direction (for chirality)
  enum BondDir {
    NONE = BondEnums::BondDir::NONE,              //!< no special style
    BEGINWEDGE = BondEnums::BondDir::BEGINWEDGE,  //!< wedged: narrow at begin
    BEGINDASH = BondEnums::BondDir::BEGINDASH,    //!< dashed: narrow at begin
    // FIX: this may not really be adequate
    ENDDOWNRIGHT = BondEnums::BondDir::ENDDOWNRIGHT,  //!< for cis/trans
    ENDUPRIGHT = BondEnums::BondDir::ENDUPRIGHT,      //!<  ditto
    //!< a "crossed" double bond
    EITHERDOUBLE = BondEnums::BondDir::EITHERDOUBLE,
    //!< intentionally unspecified stereochemistry
    UNKNOWN = BondEnums::BondDir::UNKNOWN,
  };

  //! the nature of the bond's stereochem (for cis/trans)
  enum BondStereo {  // stereochemistry of double bonds
    STEREONONE = BondEnums::BondStereo::STEREONONE,         // no special style
    STEREOANY = BondEnums::BondStereo::STEREOANY,           // intentionally unspecified
    // -- Put any true specifications about this point so
    // that we can do comparisons like if(bond->getStereo()>Bond::STEREOANY)
    STEREOZ = BondEnums::BondStereo::STEREOZ,               // Z double bond
    STEREOE = BondEnums::BondStereo::STEREOE,               // E double bond
    STEREOCIS = BondEnums::BondStereo::STEREOCIS,           // cis double bond
    STEREOTRANS = BondEnums::BondStereo::STEREOTRANS,       // trans double bond
    STEREOATROPCW = BondEnums::BondStereo::STEREOATROPCW,   // cis double bond
    STEREOATROPCCW = BondEnums::BondStereo::STEREOATROPCCW  // trans double bond
  };

  Bond();
  //! construct with a particular BondType
  explicit Bond(BondType bT);
  Bond(const Bond &other);

  virtual ~Bond();
  Bond &operator=(const Bond &other);

  // NOTE: the move methods are somewhat fraught for bonds associated with
  // molecules since the molecule will still be pointing to the original object
  Bond(Bond &&o) noexcept;
  Bond &operator=(Bond &&o) noexcept;

  //! returns a copy
  /*!
    <b>Note:</b> the caller is responsible for <tt>delete</tt>ing
     the returned pointer.
  */
  virtual Bond *copy() const;

  //! returns whether or not this instance belongs to a molecule
  bool hasOwningMol() const { return dp_owningMol != nullptr; }

  //! returns a reference to the ROMol that owns this instance
  ROMol &getOwningMol() const;

  //! returns a reference to the RDMol that owns this instance
  const RDMol &getRDMol() const {
    PRECONDITION(dp_owningMol, "no owner");
    return *dp_owningMol;
  }
  //! overload
  RDMol &getRDMol() {
    PRECONDITION(dp_owningMol, "no owner");
    return *dp_owningMol;
  }

 protected:
  //! returns a reference to the RDMol that contains the data for this bond,
  //! for internal use only
  const RDMol &getDataRDMol() const {
    return *dp_dataMol;
  }
  //! overload
  RDMol &getDataRDMol() {
    return *dp_dataMol;
  }
 public:

  //! sets our owning molecule
  void setOwningMol(ROMol *other);
  //! \overload
  void setOwningMol(ROMol &other) { setOwningMol(&other); }
  //! \overload
  void setOwningMol(RDMol *other);

  // inverts the chirality of an atropisomer
  bool invertChirality();

  //! returns our index within the ROMol
  unsigned int getIdx() const { return d_index; }
  //! sets our index within the ROMol
  /*!
    <b>Notes:</b>
      - this makes no sense if we do not have an owning molecule
      - the index should be <tt>< this->getOwningMol()->getNumBonds()</tt>
  */
  void setIdx(unsigned int index) {d_index = index;}

  //! returns our \c bondType
  BondType getBondType() const;
  //! sets our \c bondType
  void setBondType(BondType bT);
  //! \brief returns our \c bondType as a double
  //!   (e.g. SINGLE->1.0, AROMATIC->1.5, etc.)
  double getBondTypeAsDouble() const;

  //! returns our contribution to the explicit valence of an Atom
  /*!
    <b>Notes:</b>
      - requires an owning molecule
  */
  virtual double getValenceContrib(const Atom *at) const;

  //! sets our \c isAromatic flag
  void setIsAromatic(bool what);
  //! returns the status of our \c isAromatic flag
  bool getIsAromatic() const;

  //! sets our \c isConjugated flag
  void setIsConjugated(bool what);
  //! returns the status of our \c isConjugated flag
  bool getIsConjugated() const;

  //! returns the index of our begin Atom
  /*!
    <b>Notes:</b>
      - this makes no sense if we do not have an owning molecule
  */
  unsigned int getBeginAtomIdx() const;

  //! returns the index of our end Atom
  /*!
    <b>Notes:</b>
      - this makes no sense if we do not have an owning molecule
  */
  unsigned int getEndAtomIdx() const;

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
  //! sets our end Atom
  /*!
    <b>Notes:</b>
      - requires an owning molecule
  */
  void setEndAtom(Atom *at);

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
  virtual bool hasQuery() const { return false; }

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

  //! sets our direction
  void setBondDir(BondDir what);
  //! returns our direction
  BondDir getBondDir() const;

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
  void setStereo(BondStereo what);
  //! returns our stereo code
  BondStereo getStereo() const;

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
  const INT_VECT &getStereoAtoms() const;
  //! \overload
  INT_VECT &getStereoAtoms();

  //! calculates any of our lazy \c properties
  /*!
    <b>Notes:</b>
      - requires an owning molecule
  */
  void updatePropertyCache(bool strict = true);

  //! Flags that can be used by to store information on bonds.
  //!   These are not serialized and should be treated as temporary values.
  //!   No guarantees are made about preserving these flags across library
  //!   calls.
  void setFlags(std::uint64_t flags);
  std::uint64_t getFlags() const;
  std::uint64_t &getFlags();

  //! returns a list with the names of our \c properties
  STR_VECT getPropList(bool includePrivate = true,
                       bool includeComputed = true) const;

  //! sets a \c property value
  /*!
    \param key the name under which the \c property should be stored.
    If a \c property is already stored under this name, it will be
    replaced.
    \param val the value to be stored
    \param computed (optional) allows the \c property to be flagged
    \c computed.

    <b>Notes:</b>
    - This must be allowed on a const Bond for backward compatibility,
    but DO NOT call this while this bond might be being modified in another
    thread, NOR while any properties on the same molecule might be being read or
    written.
  */

  //! \overload
  template <typename T>
  void setProp(const std::string_view &key, T val,
               bool computed = false) const {
    // Continue support for mixed-type cases like:
    // bond0->setProp("a", 1.7); bond1->setProp("a", 123);
    // If there is a type mismatch, this option will convert to RDValue
    dp_dataMol->setSingleBondProp(PropToken(key), getIdx(), val, computed, true);
  }
  void setProp(const std::string_view &key, const char *val,
               bool computed = false) const {
    std::string sv(val);
    // Continue support for mixed-type cases like:
    // bond0->setProp("a", 1.7); bond1->setProp("a", 123);
    // If there is a type mismatch, this option will convert to RDValue
    dp_dataMol->setSingleBondProp(PropToken(key), getIdx(), sv, computed, true);
  }

  //! allows retrieval of a particular property value
  /*!
    \param key the name under which the \c property should be stored.
    If a \c property is already stored under this name, it will be
    replaced.
    \param res a reference to the storage location for the value.

    <b>Notes:</b>
    - if no \c property with name \c key exists, a KeyErrorException will be
    thrown.
    - the \c boost::lexical_cast machinery is used to attempt type
    conversions.
    If this fails, a \c boost::bad_lexical_cast exception will be thrown.
  */
  //! \overload
  template <typename T>
  void getProp(const std::string_view &key, T &res) const {
    PropToken token(key);
    if constexpr (std::is_same_v<T, STR_VECT>) {
      if (token == detail::computedPropNameToken) {
        dp_dataMol->getComputedPropList(res, RDMol::Scope::BOND, getIdx());
        return;
      }
    }
    bool found = dp_dataMol->getBondPropIfPresent(token, getIdx(), res);
    if (!found) {
      throw KeyErrorException(key);
    }
  }

  //! \overload
  template <typename T>
  T getProp(const std::string_view &key) const {
    PropToken token(key);
    if constexpr (std::is_same_v<T, STR_VECT>) {
      if (token == detail::computedPropNameToken) {
        STR_VECT res;
        dp_dataMol->getComputedPropList(res, RDMol::Scope::BOND, getIdx());
        return res;
      }
    }
    return dp_dataMol->getBondProp<T>(token, getIdx());
  }

  //! returns whether or not we have a \c property with name \c key
  //!  and assigns the value if we do
  //! \overload
  template <typename T>
  bool getPropIfPresent(const std::string_view &key, T &res) const {
    PropToken token(key);
    if constexpr (std::is_same_v<T, STR_VECT>) {
      if (token == detail::computedPropNameToken) {
        dp_dataMol->getComputedPropList(res, RDMol::Scope::BOND, getIdx());
        return true;
      }
    }
    return dp_dataMol->getBondPropIfPresent(token, getIdx(), res);
  }

  //! \overload
  bool hasProp(const std::string_view &key) const;

  //! clears the value of a \c property
  /*!
    <b>Notes:</b>
    - if no \c property with name \c key exists, a KeyErrorException
    will be thrown.
    - if the \c property is marked as \c computed, it will also be removed
    from our list of \c computedProperties

    <b>Notes:</b>
    - This must be allowed on a const Bond for backward compatibility,
    but DO NOT call this while this bond might be being modified in another
    thread, NOR while any properties on the same molecule might be being read or
    written.
  */
  //! \overload
  void clearProp(const std::string_view &key) const;

  //! clears all of our \c computed \c properties
  /*!
    <b>Notes:</b>
    - This must be allowed on a const Bond for backward compatibility,
    but DO NOT call this while this bond might be being modified in another
    thread, NOR while any properties on the same molecule might be being read or
    written.
  */
  void clearComputedProps() const;

  //! update the properties from another
  /*
    \param source    Source to update the properties from
    \param preserve  Existing If true keep existing data, else override from
    the source
  */
  void updateProps(const Bond& source, bool preserveExisting = false);
  void updateProps(const RDProps& source, bool preserveExisting = false);

  //! Clears all properties of this Bond, but leaves everything else
  void clear();
  [[noreturn]] Dict &getDict() {
    raiseNonImplementedFunction("GetDict");
  }
  [[noreturn]] const Dict &getDict() const {
    raiseNonImplementedFunction("GetDict");
  }
 private:
  void initFromOther(const Bond& other, bool preserveProps = false);
};

inline bool isDative(const Bond::BondType bt) {
  return bt == Bond::BondType::DATIVE || bt == Bond::BondType::DATIVEL ||
         bt == Bond::BondType::DATIVER || bt == Bond::BondType::DATIVEONE;
}

inline bool isDative(const Bond &bond) {
  auto bt = bond.getBondType();
  return isDative(bt);
}

inline bool canSetDoubleBondStereo(const Bond &bond) {
  auto bondType = bond.getBondType();
  return (bondType == Bond::SINGLE || bondType == Bond::AROMATIC ||
          isDative(bond));
}

inline bool canHaveDirection(const Bond &bond) {
  auto bondType = bond.getBondType();
  return (bondType == Bond::SINGLE || bondType == Bond::AROMATIC);
}

//! returns twice the \c bondType
//! (e.g. SINGLE->2, AROMATIC->3, etc.)
RDKIT_GRAPHMOL_EXPORT extern uint8_t getTwiceBondType(RDKit::Bond::BondType b);
RDKIT_GRAPHMOL_EXPORT extern uint8_t getTwiceBondType(const RDKit::Bond &b);

};  // namespace RDKit

//! allows Bond objects to be dumped to streams
RDKIT_GRAPHMOL_EXPORT extern std::ostream &operator<<(std::ostream &target,
                                                      const RDKit::Bond &b);

#endif
