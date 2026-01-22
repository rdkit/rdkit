//
//  Copyright (C) 2001-2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file Atom.h

  \brief Defines the Atom class and associated typedefs

*/
#include <RDGeneral/export.h>
#ifndef _RD_ATOM_H
#define _RD_ATOM_H

// ours
#include <RDGeneral/Invariant.h>
#include <Query/QueryObjects.h>
#include <RDGeneral/types.h>
#include <RDGeneral/RDProps.h>
#include <GraphMol/details.h>
#include <GraphMol/RDMol.h>
#include <GraphMol/rdmol_throw.h>
#include "GraphMolEnums.h"

namespace RDKit {
class Atom;
}
//! allows Atom objects to be dumped to streams
RDKIT_GRAPHMOL_EXPORT std::ostream &operator<<(std::ostream &target,
                                               const RDKit::Atom &at);

namespace RDKit {
class ROMol;
class RWMol;
class AtomMonomerInfo;
class ConstRDMolAtom;

//! The class for representing atoms
/*!

  <b>Notes:</b>
    - many of the methods of Atom require that the Atom be associated
      with a molecule (an ROMol).
    - each Atom maintains a Dict of \c properties:
        - Each \c property is keyed by name and can store an
          arbitrary type.
        - \c Properties can be marked as \c calculated, in which case
          they will be cleared when the \c clearComputedProps() method
          is called.
        - Because they have no impact upon chemistry, all \c property
          operations are \c const, this allows extra flexibility for
          clients who need to store extra data on Atom objects.
    - Atom objects are lazy about computing their explicit and implicit valence
      values.  These will not be computed until their values are requested.

  <b>Chirality:</b>

  The chirality of an Atom is determined by two things:
    - its \c chiralTag
    - the input order of its bonds (see note below for handling of
      implicit Hs)

  For tetrahedral coordination, the \c chiralTag tells you what
  direction you have to rotate to get from bond 2 to bond 3 while looking
  down bond 1. This is pretty much identical to the SMILES representation of
  chirality.

  NOTE: if an atom has an implicit H, the bond to that H is considered to be
  at the *end* of the list of other bonds.

*/
class RDKIT_GRAPHMOL_EXPORT Atom {
  friend class RDMol;
  friend class ROMol;
  friend class RWMol;
  friend class MolPickler;
  friend std::ostream &(::operator<<)(std::ostream &target, const Atom &at);

  // When an atom is created on its own, it owns dp_dataMol to store its data
  // and dp_owningMol is nullptr. Calling setOwningMol will set dp_owningMol,
  // but as before, this does not actually insert the Atom's data into the
  // owning RDMol, so dp_dataMol will continue keeping the data until
  // ROMol::addAtom copies the data and frees dp_dataMol, setting dp_dataMol
  // to be equal to dp_owningMol, to indicate that dp_dataMol is no longer
  // owned by this Atom, and setting d_index accordingly.
  //
  // The 3 states are:
  // 1) dp_owningMol == nullptr: dp_dataMol owned by this & no owning mol
  // 2) else if dp_dataMol != dp_owningMol: dp_dataMol owned by this
  // 3) dp_dataMol == dp_owningMol: dp_owningMol owns this
  //
  // This is necessary for compatibility with previous workflows that would
  // call setOwningMol and separately ROMol::addAtom.
  RDMol *dp_dataMol = nullptr;
  RDMol *dp_owningMol = nullptr;
  atomindex_t d_index = 0;

 protected:
  Atom(RDMol *mol, atomindex_t atomIndex)
      : dp_dataMol(mol), dp_owningMol(mol), d_index(atomIndex) {}

 public:
  // FIX: grn...
  typedef Queries::Query<int, Atom const *, true> QUERYATOM_QUERY;
  // typedef Queries::Query<int, ConstRDMolAtom, true> QUERYATOM_QUERY;

  //! store hybridization
  enum HybridizationType {
    //!< hybridization that hasn't been specified
    UNSPECIFIED = AtomEnums::HybridizationType::UNSPECIFIED,
    S = AtomEnums::HybridizationType::S,
    SP = AtomEnums::HybridizationType::SP,
    SP2 = AtomEnums::HybridizationType::SP2,
    SP3 = AtomEnums::HybridizationType::SP3,
    SP2D = AtomEnums::HybridizationType::SP2D,
    SP3D = AtomEnums::HybridizationType::SP3D,
    SP3D2 = AtomEnums::HybridizationType::SP3D2,
    OTHER = AtomEnums::HybridizationType::OTHER  //!< unrecognized hybridization
  };

  //! store type of chirality
  enum ChiralType {
    //!< chirality that hasn't been specified
    CHI_UNSPECIFIED = AtomEnums::ChiralType::CHI_UNSPECIFIED,
    //!< tetrahedral: clockwise rotation (SMILES \@\@)
    CHI_TETRAHEDRAL_CW = AtomEnums::ChiralType::CHI_TETRAHEDRAL_CW,
    //!< tetrahedral: counter-clockwise rotation (SMILES \@)
    CHI_TETRAHEDRAL_CCW = AtomEnums::ChiralType::CHI_TETRAHEDRAL_CCW,
    //!< some unrecognized type of chirality
    CHI_OTHER = AtomEnums::ChiralType::CHI_OTHER,
    //!< tetrahedral, use permutation flag
    CHI_TETRAHEDRAL = AtomEnums::ChiralType::CHI_TETRAHEDRAL,
    //!< allene, use permutation flag
    CHI_ALLENE = AtomEnums::ChiralType::CHI_ALLENE,
    //!< square planar, use permutation flag
    CHI_SQUAREPLANAR = AtomEnums::ChiralType::CHI_SQUAREPLANAR,
    //!< trigonal bipyramidal, use permutation flag
    CHI_TRIGONALBIPYRAMIDAL = AtomEnums::ChiralType::CHI_TRIGONALBIPYRAMIDAL,
    //!< octahedral, use permutation flag
    CHI_OCTAHEDRAL = AtomEnums::ChiralType::CHI_OCTAHEDRAL
  };

  // this can be done with a `using` definition, but that breaks the SWIG
  // bindings
  enum class ValenceType : std::uint8_t {
    IMPLICIT = 0,
    EXPLICIT
  };

  Atom();
  //! construct an Atom with a particular atomic number
  explicit Atom(unsigned int num);
  //! construct an Atom with a particular symbol (looked up in the
  /// PeriodicTable)
  explicit Atom(const std::string &what);
  Atom(const Atom &other);
  Atom &operator=(const Atom &other);
  // NOTE: the move methods are somewhat fraught for atoms associated with
  // molecules since the molecule will still be pointing to the original object
  Atom(Atom &&other);
  Atom &operator=(Atom &&other);

 private:
  void initFromOther(const Atom &other, bool preserveProps = false);

 public:
  virtual ~Atom();

  //! makes a copy of this Atom and returns a pointer to it.
  /*!
    <b>Note:</b> the caller is responsible for <tt>delete</tt>ing the result
  */
  virtual Atom *copy() const;

  //! makes a copy of this Atom and returns a pointer to it.
  //! FIXME: Should this be private or deleted?
  /*!
    <b>Note:</b> the caller is responsible for <tt>delete</tt>ing the result
  */
  // virtual Atom *copy() const;

  //! returns our atomic number
  int getAtomicNum() const;
  //! sets our atomic number
  void setAtomicNum(int newNum);

  //! returns our symbol (determined by our atomic number)
  std::string getSymbol() const;

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
  //! returns a reference to the RDMol that contains the data for this atom,
  //! for internal use only
  const RDMol &getDataRDMol() const { return *dp_dataMol; }
  //! overload
  RDMol &getDataRDMol() { return *dp_dataMol; }
  friend RDKIT_GRAPHMOL_EXPORT int getAtomRLabel(const Atom *atom);
  friend std::string getAtomStringProp(const Atom *atom,
                                       const PropToken &token);
  friend void setAtomStringProp(Atom *atom, const PropToken &token,
                                const std::string &value);

 public:
  //! returns our index within the ROMol
  unsigned int getIdx() const { return d_index; }
  //! sets our index within the ROMol
  //! FIXME: Should this be private or deleted?
  /*!
    <b>Notes:</b>
      - this makes no sense if we do not have an owning molecule
      - the index should be <tt>< this->getOwningMol()->getNumAtoms()</tt>
  */
  void setIdx(unsigned int index) { d_index = index; }
  //! overload
  //! FIXME: Should this be private or deleted?
  template <class U>
  void setIdx(const U index) {
    setIdx(rdcast<unsigned int>(index));
  }

  //! returns the explicit degree of the Atom (number of bonded
  //!   neighbors in the graph)
  unsigned int getDegree() const;

  //! returns the total degree of the Atom (number of bonded
  //!   neighbors + number of Hs)
  unsigned int getTotalDegree() const;

  //! \brief returns the total number of Hs (implicit and explicit) that
  //! this Atom is bound to
  unsigned int getTotalNumHs(bool includeNeighbors = false) const;

  //! \brief returns the total valence (implicit and explicit)
  //! for an atom
  unsigned int getTotalValence() const;

  //! returns the number of implicit Hs this Atom is bound to
  unsigned int getNumImplicitHs() const;

  //! returns the valence (explicit or implicit) of this atom
  unsigned int getValence(ValenceType which) const;

  //! returns the explicit valence (including Hs) of this atom
  [[deprecated("please use getValence(ValenceType::EXPLICIT)")]]
  int getExplicitValence() const;

  //! returns the implicit valence for this Atom
  [[deprecated("please use getValence(ValenceType::IMPLICIT)")]]
  int getImplicitValence() const;

  //! returns whether the atom has a valency violation or not
  bool hasValenceViolation() const;

  //! returns the number of radical electrons for this Atom
  unsigned int getNumRadicalElectrons() const;
  void setNumRadicalElectrons(unsigned int num);

  //! returns the formal charge of this atom
  int getFormalCharge() const;
  //! set's the formal charge of this atom
  void setFormalCharge(int what);

  //! \brief sets our \c noImplicit flag, indicating whether or not
  //!  we are allowed to have implicit Hs
  void setNoImplicit(bool what);
  //! returns the \c noImplicit flag
  bool getNoImplicit() const;

  //! sets our number of explicit Hs
  void setNumExplicitHs(unsigned int what);
  //! returns our number of explicit Hs
  unsigned int getNumExplicitHs() const;

  //! sets our \c isAromatic flag, indicating whether or not we are aromatic
  void setIsAromatic(bool what);
  //! returns our \c isAromatic flag
  bool getIsAromatic() const;

  //! returns our mass
  double getMass() const;

  //! sets our isotope number
  void setIsotope(unsigned int what);
  //! returns our isotope number
  unsigned int getIsotope() const;

  //! sets our \c chiralTag
  void setChiralTag(ChiralType what);
  //! inverts our \c chiralTag, returns whether or not a change was made
  bool invertChirality();
  //! returns our \c chiralTag
  ChiralType getChiralTag() const;

  //! sets our hybridization
  void setHybridization(HybridizationType what);
  //! returns our hybridization
  HybridizationType getHybridization() const;

  // ------------------------------------
  // Some words of explanation before getting down into
  // the query stuff.
  // These query functions are really only here so that they
  //  can have real functionality in subclasses (like QueryAtoms).
  // Since pretty much it's gonna be a mistake to call any of these
  //  (ever), we're saddling them all with a precondition which
  //  is guaranteed to fail.  I'd like to have them be pure virtual,
  //  but that doesn't work since we need to be able to instantiate
  //  Atoms.
  // ------------------------------------

  // This method can be used to distinguish query atoms from standard atoms:
  virtual bool hasQuery() const { return false; }

  virtual std::string getQueryType() const { return ""; }

  //! NOT CALLABLE
  virtual void setQuery(QUERYATOM_QUERY *what);

  //! NOT CALLABLE
  virtual QUERYATOM_QUERY *getQuery() const;
  //! NOT CALLABLE
  virtual void expandQuery(
      QUERYATOM_QUERY *what,
      Queries::CompositeQueryType how = Queries::COMPOSITE_AND,
      bool maintainOrder = true);

  //! returns whether or not we match the argument
  /*!
      <b>Notes:</b>
        The general rule is that if a property on this atom has a non-default
     value,
        the property on the other atom must have the same value.
        The exception to this is H counts, which are ignored. These turns out to
     be
          impossible to handle generally, so rather than having odd and
     hard-to-explain
          exceptions, we ignore them entirely.

        Here are the rules for atom-atom matching:
        | This    | Other   | Match | Reason
        | CCO     | CCO     | Yes   |
        | CCO     | CC[O-]  | Yes   |
        | CC[O-]  | CCO     | No    | Charge
        | CC[O-]  | CC[O-]  | Yes   |
        | CC[OH]  | CC[O-]  | Yes   |
        | CC[OH]  | CCOC    | Yes   |
        | CCO     | CCOC    | Yes   |
        | CCC     | CCC     | Yes   |
        | CCC     | CC[14C] | Yes   |
        | CC[14C] | CCC     | No    | Isotope
        | CC[14C] | CC[14C] | Yes   |
        | C       | OCO     | Yes   |
        | [CH]    | OCO     | Yes   |
        | [CH2]   | OCO     | Yes   |
        | [CH3]   | OCO     | No    | Radical
        | C       | O[CH2]O | Yes   |
        | [CH2]   | O[CH2]O | Yes   |
  */
  virtual bool Match(Atom const *what) const;

  //! returns the perturbation order for a list of integers
  /*!

    This value is associated with chirality.

    \param probe a list of bond indices.  This must be the same
      length as our number of incoming bonds (our degree).

    \return the number of swaps required to convert the ordering
      of the probe list to match the order of our incoming bonds:
      e.g. if our incoming bond order is: <tt>[0,1,2,3]</tt>
      \verbatim
      getPerturbationOrder([1,0,2,3]) = 1
      getPerturbationOrder([1,2,3,0]) = 3
      getPerturbationOrder([1,2,0,3]) = 2
      \endverbatim

    See the class documentation for a more detailed description
    of our representation of chirality.

    <b>Notes:</b>
      - requires an owning molecule

  */
  int getPerturbationOrder(const INT_LIST &probe) const;

  //! calculates any of our lazy \c properties
  /*!
    <b>Notes:</b>
      - requires an owning molecule
      - the current lazy \c properties are implicit and explicit valence
  */
  void updatePropertyCache(bool strict = true);

  bool needsUpdatePropertyCache() const;
  void clearPropertyCache();

  //! calculates and returns our explicit valence
  /*!
    <b>Notes:</b>
      - requires an owning molecule
  */
  int calcExplicitValence(bool strict = true);

  //! calculates and returns our implicit valence
  /*!
    <b>Notes:</b>
      - requires an owning molecule
  */
  int calcImplicitValence(bool strict = true);

  AtomMonomerInfo *getMonomerInfo();
  const AtomMonomerInfo *getMonomerInfo() const;
  //! takes ownership of the pointer
  void setMonomerInfo(AtomMonomerInfo *info);

  //! Set the atom map Number of the atom
  void setAtomMapNum(int mapno, bool strict = true);
  //! Gets the atom map Number of the atom, if no atom map exists, 0 is
  //! returned.
  int getAtomMapNum() const;

  //! Flags that can be used by to store information on atoms.
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
    - This must be allowed on a const Atom for backward compatibility,
    but DO NOT call this while this atom might be being modified in another
    thread, NOR while any properties on the same molecule might be being read or
    written.
  */
  //! \overload
  template <typename T>
  void setProp(const std::string_view &key, T val,
               bool computed = false) const {
    // Continue support for mixed-type cases like:
    // atom0->setProp("a", 1.7); atom1->setProp("a", 123);
    // If there is a type mismatch, this option will convert to RDValue
    dp_dataMol->setSingleAtomProp(PropToken(key), getIdx(), val, computed,
                                  true);
  }
  void setProp(const std::string_view &key, const char *val,
               bool computed = false) const {
    std::string sv(val);
    // Continue support for mixed-type cases like:
    // atom0->setProp("a", 1.7); atom1->setProp("a", 123);
    // If there is a type mismatch, this option will convert to RDValue
    dp_dataMol->setSingleAtomProp(PropToken(key), getIdx(), sv, computed, true);
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
        dp_dataMol->getComputedPropList(res, Properties::Scope::ATOM, getIdx());
        return;
      }
    }
    bool found = dp_dataMol->getAtomPropIfPresent(token, getIdx(), res);
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
        dp_dataMol->getComputedPropList(res, Properties::Scope::ATOM, getIdx());
        return res;
      }
    }
    return dp_dataMol->getAtomProp<T>(token, getIdx());
  }

  //! returns whether or not we have a \c property with name \c key
  //!  and assigns the value if we do
  //! \overload
  template <typename T>
  bool getPropIfPresent(const std::string_view &key, T &res) const {
    PropToken token(key);
    if constexpr (std::is_same_v<T, STR_VECT>) {
      if (token == detail::computedPropNameToken) {
        dp_dataMol->getComputedPropList(res, Properties::Scope::ATOM, getIdx());
        return true;
      }
    }
    return dp_dataMol->getAtomPropIfPresent(token, getIdx(), res);
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
    - This must be allowed on a const Atom for backward compatibility,
    but DO NOT call this while this atom might be being modified in another
    thread, NOR while any properties on the same molecule might be being read or
    written.
  */
  //! \overload
  void clearProp(const std::string_view &key) const;

  //! clears all of our \c computed \c properties
  /*!
    <b>Notes:</b>
    - This must be allowed on a const Atom for backward compatibility,
    but DO NOT call this while this atom might be being modified in another
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
  void updateProps(const Atom &source, bool preserveExisting = false);
  void updateProps(const RDProps &source, bool preserveExisting = false);

  //! Clears all properties of this Atom, but leaves everything else
  void clear();
  [[noreturn]] Dict &getDict() { raiseNonImplementedFunction("GetDict"); }
  [[noreturn]] const Dict &getDict() const {
    raiseNonImplementedFunction("GetDict");
  }

 protected:
  //! sets our owning molecule
  void setOwningMol(ROMol *other);
  //! \overload
  void setOwningMol(ROMol &other) { setOwningMol(&other); }
  //! \overload
  void setOwningMol(RDMol *other);

  int getExplicitValencePrivate() const;
  int getImplicitValencePrivate() const;
};

//! Set the atom's MDL integer RLabel
///  Setting to 0 clears the rlabel.  Rlabel must be in the range [0..99]
RDKIT_GRAPHMOL_EXPORT void setAtomRLabel(Atom *atom, int rlabel);
RDKIT_GRAPHMOL_EXPORT int getAtomRLabel(const Atom *atom);

//! Set the atom's MDL atom alias
///  Setting to an empty string clears the alias
RDKIT_GRAPHMOL_EXPORT void setAtomAlias(Atom *atom, const std::string &alias);
RDKIT_GRAPHMOL_EXPORT std::string getAtomAlias(const Atom *atom);

//! Set the atom's MDL atom value
///  Setting to an empty string clears the value
///  This is where recursive smarts get stored in MolBlock Queries
RDKIT_GRAPHMOL_EXPORT void setAtomValue(Atom *atom, const std::string &value);
RDKIT_GRAPHMOL_EXPORT std::string getAtomValue(const Atom *atom);

//! Sets the supplemental label that will follow the atom when writing
///  smiles strings.
RDKIT_GRAPHMOL_EXPORT void setSupplementalSmilesLabel(Atom *atom,
                                                      const std::string &label);
RDKIT_GRAPHMOL_EXPORT std::string getSupplementalSmilesLabel(const Atom *atom);

//! returns true if the atom is aromatic or has an aromatic bond
RDKIT_GRAPHMOL_EXPORT bool isAromaticAtom(const Atom &atom);
//! overload
RDKIT_GRAPHMOL_EXPORT bool isAromaticAtom(ConstRDMolAtom atom);

//! returns the number of pi electrons on the atom
RDKIT_GRAPHMOL_EXPORT unsigned int numPiElectrons(const Atom &atom);
};  // namespace RDKit

//! allows Atom objects to be dumped to streams
RDKIT_GRAPHMOL_EXPORT std::ostream &operator<<(std::ostream &target,
                                               const RDKit::Atom &at);
#endif
