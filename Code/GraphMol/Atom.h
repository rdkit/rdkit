//
//  Copyright (C) 2001-2021 Greg Landrum and other RDKit contributors
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

// Std stuff
#include <iostream>

// ours
#include <RDGeneral/Invariant.h>
#include <Query/QueryObjects.h>
#include <RDGeneral/types.h>
#include <RDGeneral/RDProps.h>
#include <GraphMol/details.h>

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
class RDKIT_GRAPHMOL_EXPORT Atom : public RDProps {
  friend class MolPickler;  //!< the pickler needs access to our privates
  friend class ROMol;
  friend class RWMol;
  friend std::ostream &(::operator<<)(std::ostream &target, const Atom &at);

 public:
  // FIX: grn...
  typedef Queries::Query<int, Atom const *, true> QUERYATOM_QUERY;

  //! store hybridization
  typedef enum {
    UNSPECIFIED = 0,  //!< hybridization that hasn't been specified
    S,
    SP,
    SP2,
    SP3,
    SP2D,
    SP3D,
    SP3D2,
    OTHER  //!< unrecognized hybridization
  } HybridizationType;

  //! store type of chirality
  typedef enum {
    CHI_UNSPECIFIED = 0,  //!< chirality that hasn't been specified
    CHI_TETRAHEDRAL_CW,   //!< tetrahedral: clockwise rotation (SMILES \@\@)
    CHI_TETRAHEDRAL_CCW,  //!< tetrahedral: counter-clockwise rotation (SMILES
                          //\@)
    CHI_OTHER,            //!< some unrecognized type of chirality
    CHI_TETRAHEDRAL,      //!< tetrahedral, use permutation flag
    CHI_ALLENE,           //!< allene, use permutation flag
    CHI_SQUAREPLANAR,     //!< square planar, use permutation flag
    CHI_TRIGONALBIPYRAMIDAL,  //!< trigonal bipyramidal, use permutation flag
    CHI_OCTAHEDRAL            //!< octahedral, use permutation flag
  } ChiralType;

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
  Atom(Atom &&other) = default;
  Atom &operator=(Atom &&other) = default;

  virtual ~Atom();

  //! makes a copy of this Atom and returns a pointer to it.
  /*!
    <b>Note:</b> the caller is responsible for <tt>delete</tt>ing the result
  */
  virtual Atom *copy() const;

  //! returns our atomic number
  int getAtomicNum() const { return d_atomicNum; }
  //! sets our atomic number
  void setAtomicNum(int newNum) { d_atomicNum = newNum; }

  //! returns our symbol (determined by our atomic number)
  std::string getSymbol() const;

  //! returns whether or not this instance belongs to a molecule
  bool hasOwningMol() const { return dp_mol != nullptr; }

  //! returns a reference to the ROMol that owns this instance
  ROMol &getOwningMol() const {
    PRECONDITION(dp_mol, "no owner");
    return *dp_mol;
  }

  //! returns our index within the ROMol
  unsigned int getIdx() const { return d_index; }
  //! sets our index within the ROMol
  /*!
    <b>Notes:</b>
      - this makes no sense if we do not have an owning molecule
      - the index should be <tt>< this->getOwningMol()->getNumAtoms()</tt>
  */
  void setIdx(unsigned int index) { d_index = index; }
  //! overload
  template <class U>
  void setIdx(const U index) {
    setIdx(rdcast<unsigned int>(index));
  }
  //! returns the explicit degree of the Atom (number of bonded
  //!   neighbors in the graph)
  /*!
    <b>Notes:</b>
      - requires an owning molecule
  */
  unsigned int getDegree() const;

  //! returns the total degree of the Atom (number of bonded
  //!   neighbors + number of Hs)
  /*!
    <b>Notes:</b>
      - requires an owning molecule
  */
  unsigned int getTotalDegree() const;

  //! \brief returns the total number of Hs (implicit and explicit) that
  //! this Atom is bound to
  /*!
    <b>Notes:</b>
      - requires an owning molecule
  */
  unsigned int getTotalNumHs(bool includeNeighbors = false) const;

  //! \brief returns the total valence (implicit and explicit)
  //! for an atom
  /*!
    <b>Notes:</b>
      - requires an owning molecule
  */
  unsigned int getTotalValence() const;

  //! returns the number of implicit Hs this Atom is bound to
  /*!
    <b>Notes:</b>
      - requires an owning molecule
  */
  unsigned int getNumImplicitHs() const;

  //! returns the explicit valence (including Hs) of this atom
  int getExplicitValence() const;

  //! returns the implicit valence for this Atom
  /*!
    <b>Notes:</b>
      - requires an owning molecule
  */
  int getImplicitValence() const;

  //! returns whether the atom has a valency violation or not
  /*!
    <b>Notes:</b>
      - requires an owning molecule
  */
  bool hasValenceViolation() const;

  //! returns the number of radical electrons for this Atom
  /*!
    <b>Notes:</b>
      - requires an owning molecule
  */
  unsigned int getNumRadicalElectrons() const { return d_numRadicalElectrons; }
  void setNumRadicalElectrons(unsigned int num) { d_numRadicalElectrons = num; }

  //! returns the formal charge of this atom
  int getFormalCharge() const { return d_formalCharge; }
  //! set's the formal charge of this atom
  void setFormalCharge(int what) { d_formalCharge = what; }

  //! \brief sets our \c noImplicit flag, indicating whether or not
  //!  we are allowed to have implicit Hs
  void setNoImplicit(bool what) { df_noImplicit = what; }
  //! returns the \c noImplicit flag
  bool getNoImplicit() const { return df_noImplicit; }

  //! sets our number of explicit Hs
  void setNumExplicitHs(unsigned int what) { d_numExplicitHs = what; }
  //! returns our number of explicit Hs
  unsigned int getNumExplicitHs() const { return d_numExplicitHs; }

  //! sets our \c isAromatic flag, indicating whether or not we are aromatic
  void setIsAromatic(bool what) { df_isAromatic = what; }
  //! returns our \c isAromatic flag
  bool getIsAromatic() const { return df_isAromatic; }

  //! returns our mass
  double getMass() const;

  //! sets our isotope number
  void setIsotope(unsigned int what);
  //! returns our isotope number
  unsigned int getIsotope() const { return d_isotope; }

  //! sets our \c chiralTag
  void setChiralTag(ChiralType what) { d_chiralTag = what; }
  //! inverts our \c chiralTag, returns whether or not a change was made
  bool invertChirality();
  //! returns our \c chiralTag
  ChiralType getChiralTag() const {
    return static_cast<ChiralType>(d_chiralTag);
  }

  //! sets our hybridization
  void setHybridization(HybridizationType what) { d_hybrid = what; }
  //! returns our hybridization
  HybridizationType getHybridization() const {
    return static_cast<HybridizationType>(d_hybrid);
  }

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

  AtomMonomerInfo *getMonomerInfo() { return dp_monomerInfo; }
  const AtomMonomerInfo *getMonomerInfo() const { return dp_monomerInfo; }
  //! takes ownership of the pointer
  void setMonomerInfo(AtomMonomerInfo *info) { dp_monomerInfo = info; }

  //! Set the atom map Number of the atom
  void setAtomMapNum(int mapno, bool strict = true) {
    PRECONDITION(
        !strict || (mapno >= 0 && mapno < 1000),
        "atom map number out of range [0..1000], use strict=false to override");
    if (mapno) {
      setProp(common_properties::molAtomMapNumber, mapno);
    } else if (hasProp(common_properties::molAtomMapNumber)) {
      clearProp(common_properties::molAtomMapNumber);
    }
  }
  //! Gets the atom map Number of the atom, if no atom map exists, 0 is
  //! returned.
  int getAtomMapNum() const {
    int mapno = 0;
    getPropIfPresent(common_properties::molAtomMapNumber, mapno);
    return mapno;
  }

 protected:
  //! sets our owning molecule
  void setOwningMol(ROMol *other);
  //! sets our owning molecule
  void setOwningMol(ROMol &other) { setOwningMol(&other); }

  bool df_isAromatic;
  bool df_noImplicit;
  std::uint8_t d_numExplicitHs;
  std::int8_t d_formalCharge;
  std::uint8_t d_atomicNum;
  // NOTE that these cannot be signed, they are calculated using
  // a lazy scheme and are initialized to -1 to indicate that the
  // calculation has not yet been done.
  std::int8_t d_implicitValence, d_explicitValence;
  std::uint8_t d_numRadicalElectrons;
  std::uint8_t d_chiralTag;
  std::uint8_t d_hybrid;

  std::uint16_t d_isotope;
  atomindex_t d_index;

  ROMol *dp_mol;
  AtomMonomerInfo *dp_monomerInfo;
  void initAtom();
  void initFromOther(const Atom &other);
};

//! Set the atom's MDL integer RLabel
///  Setting to 0 clears the rlabel.  Rlabel must be in the range [0..99]
RDKIT_GRAPHMOL_EXPORT void setAtomRLabel(Atom *atm, int rlabel);
RDKIT_GRAPHMOL_EXPORT int getAtomRLabel(const Atom *atm);

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

//! returns true if the atom is to the left of C
RDKIT_GRAPHMOL_EXPORT bool isEarlyAtom(int atomicNum);
//! returns true if the atom is aromatic or has an aromatic bond
RDKIT_GRAPHMOL_EXPORT bool isAromaticAtom(const Atom &atom);
//! returns the number of pi electrons on the atom
RDKIT_GRAPHMOL_EXPORT unsigned int numPiElectrons(const Atom &atom);
};  // namespace RDKit

//! allows Atom objects to be dumped to streams
RDKIT_GRAPHMOL_EXPORT std::ostream &operator<<(std::ostream &target,
                                               const RDKit::Atom &at);
#endif
