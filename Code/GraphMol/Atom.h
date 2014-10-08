//
//  Copyright (C) 2001-2010 Greg Landrum and Rational Discovery LLC
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
#ifndef _RD_ATOM_H
#define _RD_ATOM_H

// Std stuff
#include <iostream>
#include <boost/foreach.hpp>

// ours
#include <Query/QueryObjects.h>
#include <RDGeneral/types.h>
#include <RDGeneral/Dict.h>

namespace RDKit{
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
  class Atom {
    friend class MolPickler; //!< the pickler needs access to our privates
    friend class ROMol;
    friend class RWMol;
  public:

    typedef boost::shared_ptr<Atom>    ATOM_SPTR;
    typedef boost::shared_ptr<const Atom>    C_ATOM_SPTR;
    // FIX: grn...
    typedef Queries::Query<int,Atom const *,true> QUERYATOM_QUERY;

    //! store hybridization
    typedef enum {
      UNSPECIFIED=0, //!< hybridization that hasn't been specified
      S,
      SP,
      SP2,
      SP3,
      SP3D,
      SP3D2,
      OTHER          //!< unrecognized hybridization
    } HybridizationType;

    //! store type of chirality
    typedef enum {
      CHI_UNSPECIFIED=0,  //!< chirality that hasn't been specified
      CHI_TETRAHEDRAL_CW, //!< tetrahedral: clockwise rotation (SMILES \@\@)
      CHI_TETRAHEDRAL_CCW,//!< tetrahedral: counter-clockwise rotation (SMILES \@)
      CHI_OTHER           //!< some unrecognized type of chirality
    } ChiralType;

    Atom();
    //! construct an Atom with a particular atomic number
    explicit Atom(unsigned int num);
    //! construct an Atom with a particular symbol (looked up in the PeriodicTable)
    explicit Atom(std::string what);
    Atom(const Atom & other);
    virtual ~Atom();

    //! makes a copy of this Atom and returns a pointer to it.
    /*!
      <b>Note:</b> the caller is responsible for <tt>delete</tt>ing the result
    */
    virtual Atom *copy() const;

    //! returns our atomic number
    int getAtomicNum() const { return d_atomicNum; };
    //! sets our atomic number
    void setAtomicNum(int newNum) { d_atomicNum = newNum; };

    //! returns our symbol (determined by our atomic number)
    std::string getSymbol() const;

    //! returns a reference to the ROMol that owns this Atom
    ROMol &getOwningMol() const { return *dp_mol; };
  
    //! returns our index within the ROMol
    unsigned int getIdx() const {return d_index;};
    //! sets our index within the ROMol
    /*!
      <b>Notes:</b>
        - this makes no sense if we do not have an owning molecule
	- the index should be <tt>< this->getOwningMol()->getNumAtoms()</tt>
    */
    void setIdx(unsigned int index) {d_index=index;};

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
    unsigned int getTotalNumHs(bool includeNeighbors=false) const;

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

    //! returns the number of radical electrons for this Atom
    /*!
      <b>Notes:</b>
        - requires an owning molecule
    */
    unsigned int getNumRadicalElectrons() const { return d_numRadicalElectrons; };
    void setNumRadicalElectrons(unsigned int num) { d_numRadicalElectrons=num; };


    //! returns the formal charge of this atom
    int getFormalCharge() const { return d_formalCharge; };
    //! set's the formal charge of this atom
    void setFormalCharge(int what) { d_formalCharge = what;} ;

    //! \brief sets our \c noImplicit flag, indicating whether or not
    //!  we are allowed to have implicit Hs
    void setNoImplicit( bool what ) { df_noImplicit = what; };
    //! returns the \c noImplicit flag
    bool getNoImplicit() const { return df_noImplicit; };
    
    //! sets our number of explict Hs
    void setNumExplicitHs(unsigned int what) { d_numExplicitHs = what; };
    //! returns our number of explict Hs
    unsigned int getNumExplicitHs() const { return d_numExplicitHs; };
  
    //! sets our \c isAromatic flag, indicating whether or not we are aromatic
    void setIsAromatic( bool what ) { df_isAromatic = what; };
    //! returns our \c isAromatic flag
    bool getIsAromatic() const { return df_isAromatic; };

    //! sets our mass (should no longer be used)
    void setMass( double what) { d_mass = what; };
    //! returns our mass
    double getMass() const {return d_mass; };

    //! sets our isotope number
    void setIsotope(unsigned int what);
    //! returns our isotope number
    unsigned int getIsotope() const {return d_isotope; };

    //! sets our \c dativeFlag
    // intended to be used only in construction.
    // NOTE: the dative flag is not currently used anywhere
    void setDativeFlag(int what) {
      d_dativeFlag = what;
    };
    //! returns our \c dativeFlag
    // NOTE: the dative flag is not currently used anywhere
    int getDativeFlag() const {
      return d_dativeFlag;
    };
    bool hasDativeFlag(int what) const {
      return d_dativeFlag==what;
    };
    //! clears our \c dativeFlag
    // NOTE: the dative flag is not currently used anywhere
    void clearDativeFlag(){ d_dativeFlag = 0; };

    //! sets our \c chiralTag
    void setChiralTag(ChiralType what) { d_chiralTag = what; };
    //! inverts our \c chiralTag
    void invertChirality();
    //! returns our \c chiralTag
    ChiralType getChiralTag() const { return d_chiralTag; };

    //! sets our hybridization
    void setHybridization(HybridizationType what) { d_hybrid = what; };
    //! returns our hybridization
    HybridizationType getHybridization() const { return d_hybrid; };

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
    virtual bool hasQuery() const { return false; };
   
    //! NOT CALLABLE
    virtual void setQuery(QUERYATOM_QUERY *what);

    //! NOT CALLABLE
    virtual QUERYATOM_QUERY *getQuery() const;
    //! NOT CALLABLE
    virtual void expandQuery(QUERYATOM_QUERY *what,
			     Queries::CompositeQueryType how=Queries::COMPOSITE_AND,
			     bool maintainOrder=true);

    //! returns whether or not we match the argument
    /*!
        <b>Notes:</b>
          The general rule is that if a property on this atom has a non-default value,
          the property on the other atom must have the same value.
          The exception to this is H counts, which are ignored. These turns out to be
            impossible to handle generally, so rather than having odd and hard-to-explain
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
    //! \overload
    virtual inline bool Match(const ATOM_SPTR &what) const {
      return Match(what.get());
    };
  

    // ------------------------------------
    //  Local Property Dict functionality
    //  all setProp functions are const because they
    //     are not meant to change the atom chemically 
    // ------------------------------------
    //! returns a list with the names of our \c properties
    STR_VECT getPropList() const {
      return dp_props->keys();
    }

    //! sets a \c property value
    /*!
       \param key the name under which the \c property should be stored.
           If a \c property is already stored under this name, it will be
	   replaced.
       \param val the value to be stored
       \param computed (optional) allows the \c property to be flagged
           \c computed.
     */
    template <typename T>
    void setProp(const char *key, T val, bool computed=false) const{
    
      //if(!dp_props) dp_props = new Dict();
      std::string what(key);
      setProp(what,val, computed);
    }

    //! \overload
    template <typename T>
    void setProp(const std::string key, T val, bool computed=false) const {
      if (computed) {
	STR_VECT compLst;
	if(hasProp(detail::computedPropName)) getProp(detail::computedPropName, compLst);
	if (std::find(compLst.begin(), compLst.end(), key) == compLst.end()) {
	  compLst.push_back(key);
	  dp_props->setVal(detail::computedPropName, compLst);
	}
      }
      //setProp(key.c_str(),val);
      dp_props->setVal(key, val);
    }

    //! allows retrieval of a particular property value
    /*!

       \param key the name under which the \c property should be stored.
           If a \c property is already stored under this name, it will be
	   replaced.
       \param res a reference to the storage location for the value.

       <b>Notes:</b>
         - if no \c property with name \c key exists, a KeyErrorException will be thrown.
	 - the \c boost::lexical_cast machinery is used to attempt type conversions.
	   If this fails, a \c boost::bad_lexical_cast exception will be thrown.

    */
    template <typename T>
    void getProp(const char *key,T &res) const {
      dp_props->getVal(key,res);
    }
    //! \overload
    template <typename T>
    void getProp(const std::string key,T &res) const {
      dp_props->getVal(key,res);
    }
    //! \overload
    template <typename T>
    T getProp(const char *key) const {
      return dp_props->getVal<T>(key);
    }
    //! \overload
    template <typename T>
    T getProp(const std::string key) const {
      return dp_props->getVal<T>(key);
    }

    //! returns whether or not we have a \c property with name \c key
    bool hasProp(const char *key) const {
      if(!dp_props) return false;
      return dp_props->hasVal(key);
    };
    //! \overload
    bool hasProp(const std::string key) const {
      if(!dp_props) return false;
      return dp_props->hasVal(key);
    };

    
    //! clears the value of a \c property
    /*!
       <b>Notes:</b>
         - if no \c property with name \c key exists, a KeyErrorException
	   will be thrown.
	 - if the \c property is marked as \c computed, it will also be removed
	   from our list of \c computedProperties
    */
    void clearProp(const char *key) const {
      std::string what(key);
      clearProp(what);
    };
    //! \overload
    void clearProp(const std::string key) const {
      if(hasProp(detail::computedPropName)){
	STR_VECT compLst;
	getProp(detail::computedPropName, compLst);
	STR_VECT_I svi = std::find(compLst.begin(), compLst.end(), key);
	if (svi != compLst.end()) {
	  compLst.erase(svi);
	  dp_props->setVal(detail::computedPropName, compLst);
	}
      }    
      dp_props->clearVal(key);
    };

    //! clears all of our \c computed \c properties
    void clearComputedProps() const {
      if(!hasProp(detail::computedPropName)) return;
      STR_VECT compLst;
      getProp(detail::computedPropName, compLst);
      BOOST_FOREACH(const std::string &sv,compLst){
	dp_props->clearVal(sv);
      }
      compLst.clear();
      dp_props->setVal(detail::computedPropName, compLst);
    }

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
    int getPerturbationOrder(INT_LIST probe) const;

    //! calculates any of our lazy \c properties
    /*!
      <b>Notes:</b>
        - requires an owning molecule
	- the current lazy \c properties are implicit and explicit valence
    */
    void updatePropertyCache(bool strict=true);

    //! calculates and returns our explicit valence
    /*!
      <b>Notes:</b>
        - requires an owning molecule
    */
    int calcExplicitValence(bool strict=true);

    //! calculates and returns our implicit valence
    /*!
      <b>Notes:</b>
        - requires an owning molecule
    */
    int calcImplicitValence(bool strict=true);

    AtomMonomerInfo *getMonomerInfo() { return dp_monomerInfo; };
    const AtomMonomerInfo *getMonomerInfo() const { return dp_monomerInfo; };
    //! takes ownership of the pointer
    void setMonomerInfo(AtomMonomerInfo *info) { dp_monomerInfo=info; };
    
  protected:
    //! sets our owning molecule
    void setOwningMol(ROMol *other);
    //! sets our owning molecule
    void setOwningMol(ROMol &other) {setOwningMol(&other);};

    bool df_isAromatic; 
    bool df_noImplicit;
    int d_dativeFlag;
    unsigned int d_numExplicitHs;
    int d_formalCharge;
    unsigned int d_atomicNum;
    unsigned int d_index;
    // NOTE that these cannot be signed ints, they are calculated using
    // a lazy scheme and are initialized to -1 to indicate that the
    // calculation has not yet been done.
    int d_implicitValence, d_explicitValence;
    unsigned int d_numRadicalElectrons;
    ChiralType d_chiralTag;
    HybridizationType d_hybrid;
    double d_mass;
    unsigned int d_isotope;
    ROMol *dp_mol;
    Dict *dp_props;
    AtomMonomerInfo *dp_monomerInfo;
    void initAtom();
  };

};
//! allows Atom objects to be dumped to streams
std::ostream & operator<<(std::ostream& target, const RDKit::Atom &at);

#endif
