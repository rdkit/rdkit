//
//  Copyright (C) 2003-2010 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file ROMol.h

  \brief Defines the primary molecule class \c ROMol as well as associated typedefs

*/  

#ifndef __RD_ROMOL_H__
#define __RD_ROMOL_H__

/// Std stuff
#include <utility>
#include <map>

// boost stuff
#include <boost/graph/adjacency_list.hpp>
#include <boost/smart_ptr.hpp>

// our stuff
#include <RDGeneral/types.h>
#include "Atom.h"
#include "Bond.h"

#include "Conformer.h"

namespace RDKit{
  class Atom;
  class Bond;
  typedef boost::shared_ptr<Atom>    ATOM_SPTR;
  typedef boost::shared_ptr<Bond>    BOND_SPTR;

  //! This is the BGL type used to store the topology:
  typedef boost::adjacency_list< boost::vecS,
                                 boost::vecS,
                                 boost::undirectedS,
                                 ATOM_SPTR,
                                 BOND_SPTR> MolGraph; 
  class MolPickler;
  class RWMol;
  class QueryAtom;
  class QueryBond;
  class RingInfo;

  template <class T1,class T2>
  class AtomIterator_;
  class BondIterator_;
  class ConstBondIterator_;

  template <class T1,class T2>
  class AromaticAtomIterator_;
  template <class T1,class T2>
  class HeteroatomIterator_;
  template <class T1,class T2>
  class QueryAtomIterator_;
  template <class T1,class T2>
  class MatchingAtomIterator_;




  extern const int ci_RIGHTMOST_ATOM;
  extern const int ci_LEADING_BOND;
  extern const int ci_ATOM_HOLDER;


  //! ROMol is a molecule class that is intended to have a fixed topology
  /*!
    This is the primary class for most molecule operations.

    If you need to be manipulating the molecule (e.g. adding or deleting
    atoms or bonds, use an RWMol instead.

    <b>Notes:</b>
      - each ROMol maintains a Dict of \c properties:  
          - Each \c property is keyed by name and can store an
            arbitrary type.
          - \c Properties can be marked as \c calculated, in which case
            they will be cleared when the \c clearComputedProps() method
            is called.
          - Because they have no impact upon chemistry, all \c property
            operations are \c const, this allows extra flexibility for
            clients who need to store extra data on ROMol objects.

      - each ROMol has collections of \c bookmarks for Atoms and Bonds:
          - the Atom bookmarks and Bond bookmarks are stored separately
            from each other
          - each \c bookmark, an integer, can map to more than one
            Atom or Bond
          - these are currently used in molecule construction, but
            could also be useful for reaction mapping and the like

      - information about rings (SSSR and the like) is stored in the
        molecule's RingInfo pointer.
    
   */
  
  class ROMol {
  public:
    friend class MolPickler;
    friend class RWMol;
    
    //! \cond TYPEDEFS

    //! \name typedefs
    //@{
    typedef MolGraph::vertex_descriptor vertex_descriptor;
    typedef MolGraph::edge_descriptor edge_descriptor;

    typedef MolGraph::edge_iterator EDGE_ITER;
    typedef MolGraph::out_edge_iterator OEDGE_ITER;
    typedef MolGraph::vertex_iterator VERTEX_ITER;
    typedef MolGraph::adjacency_iterator ADJ_ITER;
    typedef std::pair<EDGE_ITER,EDGE_ITER> BOND_ITER_PAIR;
    typedef std::pair<OEDGE_ITER,OEDGE_ITER> OBOND_ITER_PAIR;
    typedef std::pair<VERTEX_ITER,VERTEX_ITER> ATOM_ITER_PAIR;
    typedef std::pair<ADJ_ITER,ADJ_ITER> ADJ_ITER_PAIR;

    typedef std::vector<ATOM_SPTR> ATOM_SPTR_VECT;
    typedef ATOM_SPTR_VECT::iterator ATOM_SPTR_VECT_I;
    typedef ATOM_SPTR_VECT::const_iterator ATOM_SPTR_VECT_CI;
    typedef std::vector<BOND_SPTR> BOND_SPTR_VECT;
    typedef BOND_SPTR_VECT::iterator BOND_SPTR_VECT_I;
    typedef BOND_SPTR_VECT::const_iterator BOND_SPTR_VECT_CI;
  
    typedef std::vector<Atom *> ATOM_PTR_VECT;
    typedef ATOM_PTR_VECT::iterator ATOM_PTR_VECT_I;
    typedef ATOM_PTR_VECT::const_iterator ATOM_PTR_VECT_CI;
    typedef std::vector<Bond *> BOND_PTR_VECT;
    typedef BOND_PTR_VECT::iterator BOND_PTR_VECT_I;
    typedef BOND_PTR_VECT::const_iterator BOND_PTR_VECT_CI;

    typedef std::list<Atom *> ATOM_PTR_LIST;
    typedef ATOM_PTR_LIST::iterator ATOM_PTR_LIST_I;
    typedef ATOM_PTR_LIST::const_iterator ATOM_PTR_LIST_CI;
    typedef std::list<Bond *> BOND_PTR_LIST;
    typedef BOND_PTR_LIST::iterator BOND_PTR_LIST_I;
    typedef BOND_PTR_LIST::const_iterator BOND_PTR_LIST_CI;

    // list of conformations
    typedef std::list<CONFORMER_SPTR> CONF_SPTR_LIST;
    typedef CONF_SPTR_LIST::iterator CONF_SPTR_LIST_I;
    typedef CONF_SPTR_LIST::const_iterator CONF_SPTR_LIST_CI;
    typedef std::pair<CONF_SPTR_LIST_I, CONF_SPTR_LIST_I> CONFS_I_PAIR;

    // ROFIX: these will need to be readonly somehow?
    typedef std::map<int,ATOM_PTR_LIST> ATOM_BOOKMARK_MAP;
    typedef std::map<int,BOND_PTR_LIST> BOND_BOOKMARK_MAP;

    typedef class AtomIterator_<Atom,ROMol> AtomIterator;
    typedef class AtomIterator_<const Atom,const ROMol> ConstAtomIterator;
    typedef class BondIterator_ BondIterator;
    typedef class ConstBondIterator_ ConstBondIterator;
    typedef class AromaticAtomIterator_<Atom,ROMol> AromaticAtomIterator;
    typedef class AromaticAtomIterator_<const Atom,const ROMol> ConstAromaticAtomIterator;
    typedef class HeteroatomIterator_<Atom,ROMol> HeteroatomIterator;
    typedef class HeteroatomIterator_<const Atom,const ROMol> ConstHeteroatomIterator;
    typedef class QueryAtomIterator_<Atom,ROMol> QueryAtomIterator;
    typedef class QueryAtomIterator_<const Atom,const ROMol> ConstQueryAtomIterator;
    typedef class MatchingAtomIterator_<Atom,ROMol> MatchingAtomIterator;
    typedef class MatchingAtomIterator_<const Atom,const ROMol> ConstMatchingAtomIterator;


    typedef CONF_SPTR_LIST_I ConformerIterator;
    typedef  CONF_SPTR_LIST_CI ConstConformerIterator;

    //@}
    //! \endcond

    ROMol() { initMol(); }

    //! copy constructor with a twist
    /*!
      \param other     the molecule to be copied
      \param quickCopy (optional) if this is true, the resulting ROMol will not
           copy any of the properties or bookmarks and conformers from \c other.  This can
           make the copy substantially faster (thus the name).
    */
    ROMol(const ROMol &other,bool quickCopy=false) {dp_props=0;dp_ringInfo=0;initFromOther(other,quickCopy);};
    //! construct a molecule from a pickle string
    ROMol(const std::string &binStr);

    virtual ~ROMol() { destroy(); };

    //! \name Atoms
    //@{

    //! returns our number of atoms
    unsigned int getNumAtoms(bool onlyExplicit=1) const;
    //! returns our number of heavy atoms (atomic number > 1)
    unsigned int getNumHeavyAtoms() const;
    //! returns a pointer to a particular Atom
    Atom *getAtomWithIdx(unsigned int idx);
    //! \overload
    const Atom *getAtomWithIdx(unsigned int idx) const;
    //! returns the degree (number of neighbors) of an Atom in the graph
    unsigned int getAtomDegree(const Atom *at) const;
    //! \overload
    unsigned int getAtomDegree(ATOM_SPTR at) const;
    //@}

    //! \name Bonds
    //@{

    //! returns our number of Bonds
    unsigned int getNumBonds(bool onlyHeavy=1) const; 
    //! returns a pointer to a particular Bond
    Bond *getBondWithIdx(unsigned int idx);
    //! \overload
    const Bond * getBondWithIdx(unsigned int idx) const;
    //! returns a pointer to the bond between two atoms, Null on failure
    Bond *getBondBetweenAtoms(unsigned int idx1,unsigned int idx2);
    //! \overload
    const Bond *getBondBetweenAtoms(unsigned int idx1,unsigned int idx2) const;
    //@}


    //! \name Bookmarks
    //@{

    //! associates an Atom pointer with a bookmark
    void setAtomBookmark(ATOM_SPTR at,int mark) {d_atomBookmarks[mark].push_back(at.get());};
    //! \overload
    void setAtomBookmark(Atom *at,int mark) {d_atomBookmarks[mark].push_back(at);};
    //! associates an Atom pointer with a bookmark
    void replaceAtomBookmark(ATOM_SPTR at,int mark) {
      d_atomBookmarks[mark].clear();
      d_atomBookmarks[mark].push_back(at.get());
    };
    //! \overload
    void replaceAtomBookmark(Atom *at,int mark) {
      d_atomBookmarks[mark].clear();
      d_atomBookmarks[mark].push_back(at);
    };
    //! returns the first Atom associated with the \c bookmark provided
    Atom *getAtomWithBookmark(int mark);
    //! returns all Atoms associated with the \c bookmark provided
    ATOM_PTR_LIST &getAllAtomsWithBookmark(int mark);
    //! removes a \c bookmark from our collection
    void clearAtomBookmark(const int mark);
    //! removes a particular Atom from the list associated with the \c bookmark
    void clearAtomBookmark(const int mark,const Atom *atom);
    //! \overload
    void clearAtomBookmark(const int mark,ATOM_SPTR atom) {clearAtomBookmark(mark,atom.get());};
    //! blows out all atomic \c bookmarks
    void clearAllAtomBookmarks() { d_atomBookmarks.clear(); };
    //! queries whether or not any atoms are associated with a \c bookmark
    bool hasAtomBookmark(int mark) const {return d_atomBookmarks.count(mark);};
    //! returns a pointer to all of our atom \c bookmarks
    ATOM_BOOKMARK_MAP *getAtomBookmarks() { return &d_atomBookmarks; };

    //! associates a Bond pointer with a bookmark
    void setBondBookmark(BOND_SPTR bond,int mark) {d_bondBookmarks[mark].push_back(bond.get());};
    //! \overload
    void setBondBookmark(Bond *bond,int mark) {d_bondBookmarks[mark].push_back(bond);};
    //! returns the first Bond associated with the \c bookmark provided
    Bond *getBondWithBookmark(int mark);
    //! returns all bonds associated with the \c bookmark provided
    BOND_PTR_LIST &getAllBondsWithBookmark(int mark);
    //! removes a \c bookmark from our collection
    void clearBondBookmark(int mark);
    //! removes a particular Bond from the list associated with the \c bookmark
    void clearBondBookmark(int mark,const Bond *bond);
    //! \overload
    void clearBondBookmark(int mark,BOND_SPTR bond) {clearBondBookmark(mark,bond.get());};
    //! blows out all bond \c bookmarks
    void clearAllBondBookmarks() { d_bondBookmarks.clear(); };
    //! queries whether or not any bonds are associated with a \c bookmark
    bool hasBondBookmark(int mark) const {return d_bondBookmarks.count(mark);};
    //! returns a pointer to all of our bond \c bookmarks
    BOND_BOOKMARK_MAP *getBondBookmarks() { return &d_bondBookmarks; };

    //@}


    //! \name Conformers
    //@{

    //! return the conformer with a specified ID
    //! if the ID is negative the first conformation will be returned
    const Conformer &getConformer(int id=-1) const;
    
    //! return the conformer with a specified ID
    //! if the ID is negative the first conformation will be returned
    Conformer &getConformer(int id=-1);

    //! Delete the conformation with the specified ID
    void removeConformer(unsigned int id);
    
    //! Clear all the conformations on the molecule
    void clearConformers() {d_confs.clear();}

    //! Add a new conformation to the molecule
    /*!
      \param conf - conformation to be added to the molecule, this molecule takes ownership 
                    of the conformer
      \param assignId - a unique ID will be assigned to the the conformation if true
                        otherwise it is assumed that the conformation already has an (unique) ID set
    */
    unsigned int addConformer(Conformer * conf, bool assignId=false);

    inline unsigned int getNumConformers() const {
      return d_confs.size();
    }

    //@}


    //! \name Topology
    //@{

    //! returns a pointer to our RingInfo structure
    //! <b>Note:</b> the client should not delete this.
    RingInfo *getRingInfo() const { return dp_ringInfo; };

    //! provides access to all neighbors around an Atom
    /*!
      \param at the atom whose neighbors we are looking for

      <b>Usage</b>
      \code
        ... molPtr is a const ROMol & ...
        ... atomPtr is a const Atom * ...
        ROMol::ADJ_ITER nbrIdx,endNbrs;
        boost::tie(nbrIdx,endNbrs) = molPtr.getAtomNeighbors(atomPtr);
        while(nbrIdx!=endNbrs){
          const ATOM_SPTR at=molPtr[*nbrIdx];
          ... do something with the Atom ...
          ++nbrIdx;
        }
      \endcode

    */
    ADJ_ITER_PAIR getAtomNeighbors(Atom const *at) const;
    //! \overload
    ADJ_ITER_PAIR getAtomNeighbors(ATOM_SPTR at) const;

    //! provides access to all Bond objects connected to an Atom
    /*!
      \param at the atom whose neighbors we are looking for

      <b>Usage</b>
      \code
        ... molPtr is a const ROMol * ...
        ... atomPtr is a const Atom * ...
        ROMol::OEDGE_ITER beg,end;
        boost::tie(beg,end) = molPtr->getAtomBonds(atomPtr);
        while(beg!=end){
          const BOND_SPTR bond=(*molPtr)[*beg];
          ... do something with the Bond ...
          ++beg;
        }
      \endcode
      or, if you need a non-const Bond *:
      \code
        ... molPtr is a ROMol * ...
        ... atomPtr is a const Atom * ...
        ROMol::OEDGE_ITER beg,end;
        boost::tie(beg,end) = molPtr->getAtomBonds(atomPtr);
        while(beg!=end){
          BOND_SPTR bond=(*molPtr)[*beg];
          ... do something with the Bond ...
          ++beg;
        }
      \endcode
      
      
    */
    OBOND_ITER_PAIR getAtomBonds(Atom const *at) const;

    //! returns an iterator pair for looping over all Atoms
    /*!

      <b>Usage</b>
      \code

        ROMol::VERTEX_ITER atBegin,atEnd;
        boost::tie(atBegin,atEnd) = mol.getVertices();  
        while(atBegin!=atEnd){
          ATOM_SPTR at2=mol[*atBegin];
          ... do something with the Atom ...
          ++atBegin;
        }
      \endcode
    */
    ATOM_ITER_PAIR getVertices();
    //! returns an iterator pair for looping over all Bonds
    /*!

      <b>Usage</b>
      \code

        ROMol::EDGE_ITER firstB,lastB;
        boost::tie(firstB,lastB) = mol.getEdges();
        while(firstB!=lastB){
          BOND_SPTR bond = mol[*firstB];
          ... do something with the Bond ...
          ++firstB;
        }
      \endcode
    */
    BOND_ITER_PAIR getEdges();
    //! \overload
    ATOM_ITER_PAIR getVertices() const;
    //! \overload
    BOND_ITER_PAIR getEdges() const;

    //! brief returns a pointer to our underlying BGL object
    /*!
        This can be useful if you need to call other BGL algorithms:

        Here's an example:
        \code
           ... mol is a const ROMol ...
           ... mapping is an INT_VECT ...
           mapping.resize(mol.getNumAtoms());
           const MolGraph &G_p = mol.getTopology();
           int res = boost::connected_components(G_p,&mapping[0]);
        \endcode
     */
    MolGraph const &getTopology() const { return d_graph; };
    //@}


    //! \name Iterators
    //@{

    //! get an AtomIterator pointing at our first Atom
    AtomIterator beginAtoms();
    //! \overload
    ConstAtomIterator beginAtoms() const;
    //! get an AtomIterator pointing at the end of our Atoms
    AtomIterator endAtoms();
    //! \overload
    ConstAtomIterator endAtoms() const;
    //! get a BondIterator pointing at our first Bond
    BondIterator beginBonds();
    //! \overload
    ConstBondIterator beginBonds() const;
    //! get a BondIterator pointing at the end of our Bonds
    BondIterator endBonds();
    //! \overload
    ConstBondIterator endBonds() const;
  
    //! get an AtomIterator pointing at our first aromatic Atom
    AromaticAtomIterator beginAromaticAtoms();
    //! \overload
    ConstAromaticAtomIterator beginAromaticAtoms() const;
    //! get an AtomIterator pointing at the end of our Atoms
    AromaticAtomIterator endAromaticAtoms();
    //! \overload
    ConstAromaticAtomIterator endAromaticAtoms() const;

    //! get an AtomIterator pointing at our first hetero Atom
    HeteroatomIterator beginHeteros();
    //! \overload
    ConstHeteroatomIterator beginHeteros() const;
    //! get an AtomIterator pointing at the end of our Atoms
    HeteroatomIterator endHeteros();
    //! \overload
    ConstHeteroatomIterator endHeteros() const;

    //! get an AtomIterator pointing at our first Atom that matches \c query
    QueryAtomIterator beginQueryAtoms(QueryAtom const *query);
    //! \overload
    ConstQueryAtomIterator beginQueryAtoms(QueryAtom const *) const;
    //! get an AtomIterator pointing at the end of our Atoms
    QueryAtomIterator endQueryAtoms();
    //! \overload
    ConstQueryAtomIterator endQueryAtoms() const;

    //! get an AtomIterator pointing at our first Atom that matches \c query
    MatchingAtomIterator beginMatchingAtoms(bool (*query)(Atom *));
    //! \overload
    ConstMatchingAtomIterator beginMatchingAtoms(bool (*query)(const Atom *)) const;
    //! get an AtomIterator pointing at the end of our Atoms
    MatchingAtomIterator endMatchingAtoms();
    //! \overload
    ConstMatchingAtomIterator endMatchingAtoms() const;

    inline ConformerIterator beginConformers() {
      return d_confs.begin();
    }

    inline ConformerIterator endConformers() {
      return d_confs.end();
    }

    inline ConstConformerIterator beginConformers() const {
      return d_confs.begin();
    }

    inline ConstConformerIterator endConformers() const {
      return d_confs.end();
    }

    //@}

    //! \name Properties
    //@{

    //! returns a list with the names of our \c properties
    STR_VECT getPropList(bool includePrivate=true,
                         bool includeComputed=true) const {
      const STR_VECT &tmp=dp_props->keys();
      STR_VECT res,computed;
      if(!includeComputed && hasProp(detail::computedPropName)){
        getProp(detail::computedPropName,computed);
        computed.push_back(detail::computedPropName);
      }
      
      STR_VECT::const_iterator pos = tmp.begin();
      while(pos!=tmp.end()){
        if((includePrivate || (*pos)[0]!='_') &&
	 std::find(computed.begin(),computed.end(),*pos)==computed.end()){
          res.push_back(*pos);
	}
        pos++;
      }
      return res;
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
    void setProp(const char *key, T val, bool computed=false) const {
      std::string what(key);
      setProp(what,val, computed);
    }
    //! \overload
    template <typename T>
    void setProp(const std::string key, T val, bool computed=false) const {
      if (computed) {
        STR_VECT compLst;
        dp_props->getVal(detail::computedPropName,compLst);
        if (std::find(compLst.begin(), compLst.end(), key) == compLst.end()) {
          compLst.push_back(key);
          dp_props->setVal(detail::computedPropName, compLst);
        }
      }
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
    void getProp(const char *key, T &res) const {
      dp_props->getVal(key, res);
    }
    //! \overload
    template <typename T>
    void getProp(const std::string key, T &res) const {
      dp_props->getVal(key, res);
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
      if (!dp_props) return false;
      return dp_props->hasVal(key);
    }
    //! \overload
    bool hasProp(const std::string key) const {
      if (!dp_props) return false;
      return dp_props->hasVal(key);
      //return hasProp(key.c_str());
    }

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
      STR_VECT compLst;
      getProp(detail::computedPropName, compLst);
      STR_VECT_I svi = std::find(compLst.begin(), compLst.end(), key);
      if (svi != compLst.end()) {
        compLst.erase(svi);
        dp_props->setVal(detail::computedPropName, compLst);
      }
    
      dp_props->clearVal(key);
    };

    //! clears all of our \c computed \c properties
    void clearComputedProps(bool includeRings=true) const;
    //! calculates any of our lazy \c properties
    /*!
      <b>Notes:</b>
         - this calls \c updatePropertyCache() on each of our Atoms and Bonds
    */
    void updatePropertyCache(bool strict=true);

    //@}


    //! \name Misc
    //@{
    //! sends some debugging info to a stream
    void debugMol(std::ostream& str) const;
    //@}


    ATOM_SPTR operator[](const vertex_descriptor &v) { return d_graph[v]; };
    const ATOM_SPTR operator[](const vertex_descriptor &v)  const { return d_graph[v]; };
    
    BOND_SPTR operator[](const edge_descriptor &e) { return d_graph[e]; };
    const BOND_SPTR operator[](const edge_descriptor &e)  const { return d_graph[e]; };
    
  private:
    MolGraph d_graph;
    ATOM_BOOKMARK_MAP d_atomBookmarks;
    BOND_BOOKMARK_MAP d_bondBookmarks;
    Dict *dp_props;
    RingInfo *dp_ringInfo;
    CONF_SPTR_LIST d_confs;
    ROMol &operator=(const ROMol &); // disable assignment, RWMol's support assignment 

#ifdef WIN32
  protected:
#endif
    void initMol();
    virtual void destroy();
    //! adds an Atom to our collection
    /*!
      \param atom          pointer to the Atom to add
      \param updateLabel   (optional) if this is true, the new Atom will be
                           our \c activeAtom
      \param takeOwnership (optional) if this is true, we take ownership of \c atom
                           instead of copying it.

      \return the new number of atoms
    */
    unsigned int addAtom(Atom *atom,bool updateLabel=true,bool takeOwnership=false);
    //! adds an Atom to our collection
    /*!
      \param atom          pointer to the Atom to add
      \param updateLabel   (optional) if this is true, the new Atom will be
                           our \c activeAtom


      \return the new number of atoms

      <b>Note:</b> since this is using a smart pointer, we don't need to worry about
      issues of ownership.

    */
    unsigned int addAtom(ATOM_SPTR,bool updateLabel=true);
    //! adds a Bond to our collection
    /*!
      \param bond          pointer to the Bond to add
      \param takeOwnership (optional) if this is true, we take ownership of \c bond
                           instead of copying it.

      \return the new number of bonds
    */
    unsigned int addBond(Bond *bond,bool takeOwnership=false);
    //! adds a Bond to our collection
    /*!
      \param bond          pointer to the Bond to add

      \return the new number of bonds

      <b>Note:</b> since this is using a smart pointer, we don't need to worry about
      issues of ownership.
    */
    unsigned int addBond(BOND_SPTR bsp);


    //! initializes from the contents of another molecule
    /*!
      \param other     the molecule to be copied
      \param quickCopy (optional) if this is true, we will not 
           copy any of the properties or bookmarks and conformers from \c other.  This can
           make the copy substantially faster (thus the name).
    */
    void initFromOther(const ROMol &other,bool quickCopy);

  };

  typedef std::vector<ROMol> MOL_VECT;
  typedef boost::shared_ptr<ROMol>    ROMOL_SPTR;
  typedef std::vector<ROMol *> MOL_PTR_VECT;
  typedef std::vector<ROMOL_SPTR> MOL_SPTR_VECT;

  typedef MOL_PTR_VECT::const_iterator MOL_PTR_VECT_CI;
  typedef MOL_PTR_VECT::iterator MOL_PTR_VECT_I;

}; // end of RDKit namespace
#endif
