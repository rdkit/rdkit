//
//  Copyright (C) 2001-2010 Greg Landrum and Rational Discovery LLC
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
// FIX: grn...
#include <Query/QueryObjects.h>
#include <RDGeneral/types.h>

namespace RDKit{
  class ROMol;
  class RWMol;
  class Atom;
  typedef boost::shared_ptr<Atom>    ATOM_SPTR;

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
  class Bond {
    friend class RWMol;
    friend class ROMol;
  public:
    typedef boost::shared_ptr<Bond>    BOND_SPTR;
    // FIX: grn...
    typedef Queries::Query<int,Bond const *,true> QUERYBOND_QUERY;

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
      DATIVEONE,   //!< one-electron dative (e.g. from a C in a Cp ring to a metal)
      DATIVE,      //!< standard two-electron dative
      DATIVEL,     //!< standard two-electron dative
      DATIVER,     //!< standard two-electron dative
      OTHER,       
      ZERO         //!< Zero-order bond (from http://pubs.acs.org/doi/abs/10.1021/ci200488k)
    } BondType;

    //! the bond's direction (for chirality)
    typedef enum {   
      NONE=0,         //!< no special style
      BEGINWEDGE,     //!< wedged: narrow at begin
      BEGINDASH,      //!< dashed: narrow at begin
      // FIX: this may not really be adequate
      ENDDOWNRIGHT,   //!< for cis/trans
      ENDUPRIGHT,     //!<  ditto
      EITHERDOUBLE,   //!< a "crossed" double bond 
      UNKNOWN,        //!< intentionally unspecified stereochemistry
    } BondDir;
    
    //! the nature of the bond's stereochem (for cis/trans)
    typedef enum {   // stereochemistry of double bonds
      STEREONONE=0,  // no special style
      STEREOANY,     // intentionally unspecified
      // -- Put any true specifications about this point so
      // that we can do comparisons like if(bond->getStereo()>Bond::STEREOANY) 
      STEREOZ,       // Z double bond
      STEREOE,       // E double bond
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
    BondType getBondType() const { return d_bondType; };
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
    void setIsAromatic( bool what ) { df_isAromatic = what; };
    //! returns the status of our \c isAromatic flag
    bool getIsAromatic() const { return df_isAromatic; };

    //! sets our \c isConjugated flag
    void setIsConjugated(bool what) {df_isConjugated = what;};
    //! returns the status of our \c isConjugated flag
    bool getIsConjugated() const {return df_isConjugated;};
  
    //! returns a reference to the ROMol that owns this Bond
    ROMol &getOwningMol() const { return *dp_mol; };
    //! sets our owning molecule
    void setOwningMol(ROMol *other);
    //! sets our owning molecule
    void setOwningMol(ROMol &other) {setOwningMol(&other);};

    //! returns our index within the ROMol
    /*!
      <b>Notes:</b>
        - this makes no sense if we do not have an owning molecule

    */
    unsigned int getIdx() const {return d_index;};
    //! sets our index within the ROMol
    /*!
      <b>Notes:</b>
        - this makes no sense if we do not have an owning molecule
	- the index should be <tt>< this->getOwningMol()->getNumBonds()</tt>
    */
    void setIdx(unsigned int index) {d_index=index;};

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
    virtual void expandQuery(QUERYBOND_QUERY *what,
                             Queries::CompositeQueryType how=Queries::COMPOSITE_AND,
                             bool maintainOrder=true);

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
    BondDir getBondDir() const { return d_dirTag; };
  
    //! sets our stereo code
    void setStereo(BondStereo what) { d_stereo = what; };
    //! returns our stereo code
    BondStereo getStereo() const { return d_stereo; };

    //! returns the indices of our stereo atoms
    const INT_VECT &getStereoAtoms() const { return d_stereoAtoms; };
    //! \overload
    INT_VECT &getStereoAtoms() { return d_stereoAtoms; };
  
    // ------------------------------------
    //  Local Property Dict functionality
    //  FIX: at some point this stuff should go in a mixin class
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
    void setProp(const char *key,T val, bool computed=false) const{
      //if(!dp_props) dp_props = new Dict();
      std::string what(key);
      setProp(what,val, computed);
    }
    //! \overload
    template <typename T>
    void setProp(const std::string key,T val, bool computed=false ) const{
      //setProp(key.c_str(),val);
      if (computed) {
	STR_VECT compLst;
	if(hasProp(detail::computedPropName)) getProp(detail::computedPropName, compLst);
	if (std::find(compLst.begin(), compLst.end(), key) == compLst.end()) {
	  compLst.push_back(key);
	  dp_props->setVal(detail::computedPropName, compLst);
	}
      }
      dp_props->setVal(key,val);
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
      PRECONDITION(dp_props,"getProp called on empty property dict");
      dp_props->getVal(key,res);
    }
    //! \overload
    template <typename T>
    void getProp(const std::string key,T &res) const {
      PRECONDITION(dp_props,"getProp called on empty property dict");
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
      STR_VECT_CI svi;
      for (svi = compLst.begin(); svi != compLst.end(); svi++) {
	dp_props->clearVal(*svi);
      }
      compLst.clear();
      dp_props->setVal(detail::computedPropName, compLst);
    }

    //! calculates any of our lazy \c properties
    /*!
      <b>Notes:</b>
        - requires an owning molecule
    */
    void updatePropertyCache(bool strict=true) { (void)strict; }

  protected:
    //! sets our owning molecule
    //void setOwningMol(ROMol *other);
    //! sets our owning molecule
    //void setOwningMol(ROMol &other) {setOwningMol(&other);};
    BondType d_bondType;
    ROMol *dp_mol;
    bool df_isAromatic;
    bool df_isConjugated;
    unsigned int d_index;
    unsigned int d_beginAtomIdx,d_endAtomIdx;
    BondDir d_dirTag;
    BondStereo d_stereo;
    INT_VECT d_stereoAtoms;
    //Atom::ATOM_SPTR dsp_beginAtom,dsp_endAtom;
    Dict *dp_props;

    void initBond();
  };

};

//! allows Bond objects to be dumped to streams
extern std::ostream & operator<<(std::ostream& target, const RDKit::Bond &b);

#endif
