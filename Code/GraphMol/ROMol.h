//
//  Copyright (C) 2003-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#ifndef __RD_ROMOL_H__
#define __RD_ROMOL_H__

/// Std stuff
#include <utility>
#include <map>

// boost stuff
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
// our stuff
#include "AtomProps.h"
#include "BondProps.h"

#include "Conformer.h"

namespace RDKit{
  //! This is the BGL type used to store the topology:
  typedef boost::adjacency_list< boost::vecS,
				 boost::vecS,
				 boost::undirectedS,
				 AtomProperty,
				 BondProperty> MolGraph; 
  class MolPickler;
  class RWMol;
  class Atom;
  class Bond;
  class QueryAtom;
  class QueryBond;
  class RingInfo;

  template <class T1,class T2,class T3>
  class AtomIterator_;
  template <class T1,class T2,class T3>
  class AromaticAtomIterator_;
  template <class T1,class T2,class T3>
  class HeteroatomIterator_;
  template <class T1,class T2,class T3>
  class QueryAtomIterator_;
  class BondIterator_;
  class ConstBondIterator_;

  typedef boost::shared_ptr<Atom>    ATOM_SPTR;
  typedef boost::shared_ptr<Bond>    BOND_SPTR;

  extern const int ci_RIGHTMOST_ATOM;
  extern const int ci_LEADING_BOND;
  extern const int ci_ATOM_HOLDER;


  //! ROMol is a molecule class that is intended to have a fixed topology
  /*!

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
    

    typedef boost::property_map<MolGraph,vertex_atom_t>  GRAPH_MOL_ATOM_PMAP;
    typedef boost::property_map<MolGraph,edge_bond_t>  GRAPH_MOL_BOND_PMAP;
    typedef boost::graph_traits<MolGraph> GRAPH_MOL_TRAITS;
    typedef GRAPH_MOL_TRAITS::edge_iterator EDGE_ITER;
    typedef GRAPH_MOL_TRAITS::out_edge_iterator OEDGE_ITER;
    typedef GRAPH_MOL_TRAITS::vertex_iterator VERTEX_ITER;
    typedef GRAPH_MOL_TRAITS::adjacency_iterator ADJ_ITER;
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
    typedef Atom * GRAPH_NODE_TYPE;
    typedef Bond * GRAPH_EDGE_TYPE;
    typedef Atom const * GRAPH_NODE_CONST_TYPE;
    typedef Bond const * GRAPH_EDGE_CONST_TYPE;
    typedef std::map<int,ATOM_PTR_LIST> ATOM_BOOKMARK_MAP;
    typedef std::map<int,BOND_PTR_LIST> BOND_BOOKMARK_MAP;

    typedef class AtomIterator_<Atom,ROMol,GRAPH_MOL_ATOM_PMAP::type> AtomIterator;
    typedef class AtomIterator_<const Atom,const ROMol,GRAPH_MOL_ATOM_PMAP::const_type> ConstAtomIterator;
    typedef class AromaticAtomIterator_<Atom,ROMol,GRAPH_MOL_ATOM_PMAP::type> AromaticAtomIterator;
    typedef class AromaticAtomIterator_<const Atom,const ROMol,GRAPH_MOL_ATOM_PMAP::const_type> ConstAromaticAtomIterator;
    typedef class HeteroatomIterator_<Atom,ROMol,GRAPH_MOL_ATOM_PMAP::type> HeteroatomIterator;
    typedef class HeteroatomIterator_<const Atom,const ROMol,GRAPH_MOL_ATOM_PMAP::const_type> ConstHeteroatomIterator;
    typedef class QueryAtomIterator_<Atom,ROMol,GRAPH_MOL_ATOM_PMAP::type> QueryAtomIterator;
    typedef class QueryAtomIterator_<const Atom,const ROMol,GRAPH_MOL_ATOM_PMAP::const_type> ConstQueryAtomIterator;
    typedef class BondIterator_ BondIterator;
    typedef class ConstBondIterator_ ConstBondIterator;

    typedef CONF_SPTR_LIST_I ConformerIterator;
    typedef  CONF_SPTR_LIST_CI ConstConformerIterator;

    ROMol() { initMol(); }

    //! copy constructor with a twist
    /*!
      \param other     the molecule to be copied
      \param quickCopy (optional) if this is true, the resulting ROMol will not
           copy any of the properties or bookmarks and conformers from \c other.  This can
	   make the copy substantially faster (thus the name).
    */
    ROMol(const ROMol &other,bool quickCopy=false);
    //! construct a molecule from a pickle string
    ROMol(const std::string &binStr);

    virtual ~ROMol() { destroy(); };
  

    // --------------------------------------------
    //
    //  Atoms
    //
    // --------------------------------------------

    //! returns our number of Atoms
    unsigned int getNumAtoms(bool onlyHeavy=1) const;
    //! returns a pointer to a particular Atom
    GRAPH_NODE_TYPE getAtomWithIdx(unsigned int idx);
    //! \overload
    GRAPH_NODE_CONST_TYPE getAtomWithIdx(unsigned int idx) const;
    //! returns the degree (number of neighbors) of an Atom in the graph
    unsigned int getAtomDegree(const Atom *at) const;
    //! \overload
    unsigned int getAtomDegree(ATOM_SPTR at) const;

    // --------------------------------------------
    //
    //  Bonds
    //
    // --------------------------------------------

    //! returns our number of Bonds
    unsigned int getNumBonds(bool onlyHeavy=1) const; 
    //! returns a pointer to a particular Bond
    GRAPH_EDGE_TYPE getBondWithIdx(unsigned int idx);
    //! \overload
    GRAPH_EDGE_CONST_TYPE getBondWithIdx(unsigned int idx) const;
    //! returns a pointer to the bond between two atoms, Null on failure
    GRAPH_EDGE_TYPE getBondBetweenAtoms(unsigned int idx1,unsigned int idx2);
    //! \overload
    GRAPH_EDGE_CONST_TYPE getBondBetweenAtoms(unsigned int idx1,unsigned int idx2) const;


    // --------------------------------------------
    //
    //  Bookmarks
    //
    // --------------------------------------------

    //! associates an Atom pointer with a bookmark
    void setAtomBookmark(ATOM_SPTR at,int mark) {d_atomBookmarks[mark].push_back(at.get());};
    //! \overload
    void setAtomBookmark(Atom *at,int mark) {d_atomBookmarks[mark].push_back(at);};
    //! returns the first Atom associated with the \c bookmark provided
    GRAPH_NODE_TYPE  getAtomWithBookmark(int mark);
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
    GRAPH_EDGE_TYPE getBondWithBookmark(int mark);
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
    bool hasBondBookmark(int mark) {return d_bondBookmarks.count(mark);};
    //! returns a pointer to all of our bond \c bookmarks
    BOND_BOOKMARK_MAP *getBondBookmarks() { return &d_bondBookmarks; };


    // --------------------------------------------
    //
    //  Collections
    //
    // --------------------------------------------

    //! provides access to all neighbors around an Atom
    /*!
      \param at the atom whose neighbors we are looking for

      <b>Usage</b>
      \verbatim
        ... molPtr is a const ROMol * ...
        ... atomPtr is a const Atom * ...
        ROMol::ADJ_ITER nbrIdx,endNbrs;
	boost::tie(nbrIdx,endNbrs) = molPtr->getAtomNeighbors(atomPtr);
	while(nbrIdx!=endNbrs){
	  const Atom *at=molPtr->getAtomWithIdx(*nbrIdx);
	  ... do something with the Atom ...
	  nbrIdx++;
	}
      \endverbatim

      <b>Notes:</b>
        - technically, we're probably suppposed to be using the atom pmap here
	  (accessible using ROMol::getAtomPMap()), but that's not actually required.
	  
    */
    ADJ_ITER_PAIR getAtomNeighbors(Atom const *at) const;
    //! \overload
    ADJ_ITER_PAIR getAtomNeighbors(ATOM_SPTR at) const;

    //! provides access to all Bond objects connected to an Atom
    /*!
      \param at the atom whose neighbors we are looking for

      <b>Usage</b>
      \verbatim
        ... molPtr is a const ROMol * ...
        ... atomPtr is a const Atom * ...
        ROMol::OEDGE_ITER beg,end;
        ROMol::GRAPH_MOL_BOND_PMAP::const_type pMap = molPtr->getBondPMap();
	boost::tie(beg,end) = molPtr->getAtomNeighbors(atomPtr);
	while(beg!=end){
	  const Bond *bond=pMap[*beg];
	  ... do something with the Bond ...
	  beg++;
	}
      \endverbatim
      or, if you need a non-const Bond *:
      \verbatim
        ... molPtr is a ROMol * ...
        ... atomPtr is a const Atom * ...
        ROMol::OEDGE_ITER beg,end;
        ROMol::GRAPH_MOL_BOND_PMAP::type pMap = molPtr->getBondPMap();
	boost::tie(beg,end) = molPtr->getAtomNeighbors(atomPtr);
	while(beg!=end){
	  Bond *bond=pMap[*beg];
	  ... do something with the Bond ...
	  beg++;
	}
      \endverbatim
      
      
    */
    OBOND_ITER_PAIR getAtomBonds(Atom const *at) const;
    //! returns the atom PMap
    GRAPH_MOL_ATOM_PMAP::type getAtomPMap();
    //! returns the bond PMap (required to use Bond iterators)
    GRAPH_MOL_BOND_PMAP::type getBondPMap();
    //! \overload
    GRAPH_MOL_ATOM_PMAP::const_type getAtomPMap() const;
    //! \overload
    GRAPH_MOL_BOND_PMAP::const_type getBondPMap() const;

    //! returns an iterator pair for looping over all Atoms
    /*!

      <b>Usage</b>
      \verbatim
        ... molPtr is an ROMol * ...
	ROMol::GRAPH_MOL_ATOM_PMAP::type atomMap = molPtr->getAtomPMap();
	ROMol::VERTEX_ITER atBegin,atEnd;
	while(atBegin!=atEnd){
   	  Atom *at2=atomMap[*atBegin];
	  ... do something with the Atom ...
          atBegin++;
        }
      \endverbatim
    */
    ATOM_ITER_PAIR getVertices();
    //! returns an iterator pair for looping over all Bonds
    /*!

      <b>Usage</b>
      \verbatim
        ... molPtr is a ROMol * ...
        ROMol::EDGE_ITER firstB,lastB;
	boost::tie(firstB,lastB) = mol.getEdges();
	ROMol::GRAPH_MOL_BOND_PMAP::type bondMap = mol.getBondPMap();
	while(firstB!=lastB){
   	  Bond *bond = bondMap[*firstB];
	  ... do something with the Bond ...
	  firstB++;
	}
      \endverbatim
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
	\verbatim
	   ... mol is a const ROMol ...
	   ... mapping is an INT_VECT ...
	   mapping.resize(mol.getNumAtoms());
           const MolGraph *G_p = mol.getTopology();
           int res = boost::connected_components(*G_p,&mapping[0]);
	\endverbatim
     */
    MolGraph const *getTopology() const { return &d_graph; };

    //! sends some debugging info to a stream
    void debugMol(std::ostream& str) const;

    // --------------------------------------------
    //
    //  Iterators
    //
    // --------------------------------------------

    //! get an AtomIterator pointing at our first Atom
    AtomIterator beginAtoms();
    //! \overload
    ConstAtomIterator beginAtoms() const;
    //! get an AtomIterator pointing at the end of our Atoms
    AtomIterator endAtoms();
    //! \overload
    ConstAtomIterator endAtoms() const;
  
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
    //! gte an AtomIterator pointing at the end of our Atoms
    QueryAtomIterator endQueryAtoms();
    //! \overload
    ConstQueryAtomIterator endQueryAtoms() const;

    //! get a BondIterator pointing at our first Bond
    BondIterator beginBonds();
    //! \overload
    ConstBondIterator beginBonds() const;
    //! get a BondIterator pointing at the end of our Bonds
    BondIterator endBonds();
    //! \overload
    ConstBondIterator endBonds() const;
  
    //---------------------------------------------------
    //
    //    Properties
    //
    //---------------------------------------------------

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
    void setProp(const char *key, T val, bool computed=false) const {
      std::string what(key);
      setProp(what,val, computed);
    }
    //! \overload
    template <typename T>
    void setProp(const std::string key, T val, bool computed=false) const {
      if (computed) {
	STR_VECT compLst;
	getProp("computedProps", compLst);
	if (std::find(compLst.begin(), compLst.end(), key) == compLst.end()) {
	  compLst.push_back(key);
	  dp_props->setVal("computedProps", compLst);
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
      //getProp(key.c_str(), res);
      dp_props->getVal(key, res);
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
      getProp("computedProps", compLst);
      STR_VECT_I svi = std::find(compLst.begin(), compLst.end(), key);
      if (svi != compLst.end()) {
	compLst.erase(svi);
	dp_props->setVal("computedProps", compLst);
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

    //! returns a pointer to our RingInfo structure
    RingInfo *getRingInfo() const { return dp_ringInfo; };

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
      \param conf - conformation to be added to the molecule
      \param assignId - a unique ID will be assigned to the the conformation if true
                        otherwise it is assumed that the conformation already has an (unique) ID set
    */
    unsigned int addConformer(Conformer * conf, bool assignId=false);

    //! return a reference to the list of conformers
    //CONF_SPTR_LIST &getConformers() {
    //  return d_confs;
    //}

    // return a const reference to the list of conformers
    //const CONF_SPTR_LIST &getConformers() const {
    //  return d_confs;
    //}

    inline unsigned int getNumConformers() const {
      return d_confs.size();
    }

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

  private:
    MolGraph d_graph;
    ATOM_BOOKMARK_MAP d_atomBookmarks;
    BOND_BOOKMARK_MAP d_bondBookmarks;
    Dict *dp_props;
    RingInfo *dp_ringInfo;
    CONF_SPTR_LIST d_confs;
    ROMol &operator=(const ROMol &); // disable assignment

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

  };

  typedef std::vector<ROMol> MOL_VECT;
  typedef boost::shared_ptr<ROMol>    ROMOL_SPTR;
  typedef std::vector<ROMol *> MOL_PTR_VECT;
  typedef std::vector<ROMOL_SPTR> MOL_SPTR_VECT;

  typedef MOL_PTR_VECT::const_iterator MOL_PTR_VECT_CI;
  typedef MOL_PTR_VECT::iterator MOL_PTR_VECT_I;

}; // end of RDKit namespace
#endif
