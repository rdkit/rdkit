// $Id$
//
//  Copyright (C) 2003-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#include <iostream>

// our stuff
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include "ROMol.h"
#include "Atom.h"
#include "QueryAtom.h"
#include "Bond.h"
#include "QueryBond.h"
#include "MolPickler.h"

#if 0
#include "AtomIterators.h"
#include "BondIterators.h"
#endif

#include "Conformer.h"

namespace RDKit{
  class QueryAtom;
  class QueryBond;

  const int ci_RIGHTMOST_ATOM=-0xBADBEEF;
  const int ci_LEADING_BOND=-0xBADBEEF+1;
  const int ci_ATOM_HOLDER=-0xDEADD06;

  void ROMol::destroy(){

#if 0
    // blow out atoms and edges:
    MolGraph::vertex_iterator begVerts,endVerts;
    boost::tie(begVerts,endVerts)= boost::vertices(d_graph);
    while(begVerts != endVerts){
      MolGraph::vertex_descriptor u=*begVerts;
      ++begVerts;
      // first remove all connections to the atom:
      boost::clear_vertex(u,d_graph);
      // now remove the atom itself:
      boost::remove_vertex(u,d_graph);
    }
#endif
    
    d_atomBookmarks.clear();
    d_bondBookmarks.clear();
    d_graph.clear();
    if (dp_props) {
      delete dp_props;
      dp_props=0;
    }
    if (dp_ringInfo) {
      delete dp_ringInfo;
      dp_ringInfo=0;
    }
  };

  ROMol::ROMol(const std::string &pickle) {
    initMol();
    MolPickler::molFromPickle(pickle,*this);
  }
  
  void ROMol::initFromOther(const ROMol &other,bool quickCopy){
    // copy over the atoms
    const MolGraph &oGraph=other.d_graph;
    ROMol::ATOM_ITER_PAIR atItP = other.getVertices();
    while(atItP.first != atItP.second ){
      addAtom(oGraph[*atItP.first]->copy(),false,true);
      ++atItP.first;
    }

    // and the bonds:
    ROMol::BOND_ITER_PAIR bondItP = other.getEdges();
    while(bondItP.first != bondItP.second ){
      addBond(oGraph[*(bondItP.first++)]->copy(),true);
    }

    // ring information
    if(dp_ringInfo) delete dp_ringInfo;
    if(other.dp_ringInfo){
      dp_ringInfo = new RingInfo(*(other.dp_ringInfo));
    } else {
      dp_ringInfo = new RingInfo();
    }

    if(dp_props) delete dp_props;
    dp_props=0;
    if(!quickCopy){
      // copy conformations
      ConstConformerIterator ci;
    
      for (ci = other.beginConformers(); ci != other.endConformers(); ci++) {
        Conformer *conf = new Conformer(*(*ci));
        this->addConformer(conf);
      }

      if (other.dp_props) {
        dp_props = new Dict(*other.dp_props);
      }
  
      // Bookmarks should be copied as well:
      ATOM_BOOKMARK_MAP::const_iterator abmI;
      for(abmI=other.d_atomBookmarks.begin();abmI!=other.d_atomBookmarks.end();abmI++){
        ATOM_PTR_LIST::const_iterator aplI;
        for(aplI=abmI->second.begin();aplI!=abmI->second.end();aplI++){
          int idx=(*aplI)->getIdx();
          int first=abmI->first;
          Atom *at=getAtomWithIdx(idx);
          setAtomBookmark(at,first);
        }
      }
      BOND_BOOKMARK_MAP::const_iterator bbmI;
      for(bbmI=other.d_bondBookmarks.begin();bbmI!=other.d_bondBookmarks.end();bbmI++){
        BOND_PTR_LIST::const_iterator bplI;
        for(bplI=bbmI->second.begin();bplI!=bbmI->second.end();bplI++){
          setBondBookmark(getBondWithIdx((*bplI)->getIdx()),bbmI->first);
        }
      }
    }
    if(!dp_props){
      dp_props = new Dict();
      STR_VECT computed;
      dp_props->setVal("computedProps", computed);
    }
  }

  void ROMol::initMol() {
    dp_props = new Dict();
    dp_ringInfo = new RingInfo();
    // ok every molecule contains a property entry called "computedProps" which provides
    //  list of property keys that correspond to value that have been computed
    // this can used to blow out all computed properties while leaving the rest along
    // initialize this list to an empty vector of strings
    STR_VECT computed;
    dp_props->setVal("computedProps", computed);
  }
  
  unsigned int ROMol::getAtomDegree(const Atom *at) const {
    return (boost::out_degree(at->getIdx(),d_graph));
  };
  unsigned int ROMol::getAtomDegree(Atom::ATOM_SPTR at) const {
    return getAtomDegree(at.get());
  };
  unsigned int ROMol::getNumAtoms(bool onlyHeavy) const {
    int res = boost::num_vertices(d_graph);
    if (!onlyHeavy) {
      // if we are interested in hydrogens as well add them up from 
      // each heavy atom
      ConstAtomIterator ai;
      for (ai = beginAtoms(); ai != endAtoms(); ai++) {
        res += (*ai)->getTotalNumHs();
      }
    }
    return res;
  };


  Atom *ROMol::getAtomWithIdx(unsigned int idx)
  {
    PRECONDITION(getNumAtoms()>0,"no atoms");
    RANGE_CHECK(0,idx,getNumAtoms()-1);
  
    MolGraph::vertex_descriptor vd = boost::vertex(idx,d_graph);
    Atom *res = d_graph[vd].get();
    POSTCONDITION(res,"");
    return res;
  }

  const Atom *ROMol::getAtomWithIdx(unsigned int idx) const 
  {
    PRECONDITION(getNumAtoms()>0,"no atoms");
    RANGE_CHECK(0,idx,getNumAtoms()-1);
  
    MolGraph::vertex_descriptor vd = boost::vertex(idx,d_graph);
    const Atom *res = d_graph[vd].get();

    POSTCONDITION(res,"");
    return res;
  }

  // returns the first inserted atom with the given bookmark
  Atom * ROMol::getAtomWithBookmark(const int mark) {
    PRECONDITION(d_atomBookmarks.count(mark) != 0,"atom bookmark not found");
    PRECONDITION(d_atomBookmarks[mark].begin() != d_atomBookmarks[mark].end(),"atom bookmark not found");

    return *(d_atomBookmarks[mark].begin());
  };

  // returns all atoms with the given bookmark
  ROMol::ATOM_PTR_LIST &ROMol::getAllAtomsWithBookmark(const int mark) {
    PRECONDITION(d_atomBookmarks.count(mark) != 0,"atom bookmark not found");
    return d_atomBookmarks[mark];
  };

  // returns the first inserted bond with the given bookmark
  Bond *ROMol::getBondWithBookmark(const int mark) {
    PRECONDITION(d_bondBookmarks.count(mark) != 0,"bond bookmark not found");
    PRECONDITION(d_bondBookmarks[mark].begin() != d_bondBookmarks[mark].end(),"bond bookmark not found");
    return *(d_bondBookmarks[mark].begin());
  };

  // returns all bonds with the given bookmark
  ROMol::BOND_PTR_LIST &ROMol::getAllBondsWithBookmark(const int mark) {
    PRECONDITION(d_bondBookmarks.count(mark) != 0,"bond bookmark not found");
    return d_bondBookmarks[mark];
  };


  void ROMol::clearAtomBookmark(const int mark){
    d_atomBookmarks.erase(mark);
  }

  void ROMol::clearAtomBookmark(const int mark,const Atom *atom){
    if(d_atomBookmarks.count(mark) != 0){
      ATOM_PTR_LIST *entry=&d_atomBookmarks[mark];
      ATOM_PTR_LIST::iterator i;
      unsigned int tgtIdx=atom->getIdx();
      for(i=entry->begin();i!=entry->end();i++){
        if((*i)->getIdx()==tgtIdx){
          entry->erase(i);
          break;
        }
      }
      if(entry->begin() == entry->end()){
        d_atomBookmarks.erase(mark);
      }
    }

  }

  void ROMol::clearBondBookmark(const int mark){
    d_bondBookmarks.erase(mark);
  }
  void ROMol::clearBondBookmark(const int mark,const Bond *bond){
    if(d_bondBookmarks.count(mark) != 0){
      BOND_PTR_LIST *entry=&d_bondBookmarks[mark];
      BOND_PTR_LIST::iterator i;    
      unsigned int tgtIdx=bond->getIdx();
      for(i=entry->begin();i!=entry->end();i++){
        if((*i)->getIdx()==tgtIdx){
          entry->erase(i);
          break;
        }
      }
      if(entry->begin() == entry->end()){
        d_bondBookmarks.erase(mark);
      }
    }
  }

  unsigned int ROMol::getNumBonds(bool onlyHeavy) const {
    // By default resturn the bonds that connect only the heavy atoms
    // hydrogen connecting bonds are ignores
    int res = boost::num_edges(d_graph);
    if (!onlyHeavy) {
      // If we need hydrogen connecting bonds add them up
      ConstAtomIterator ai;
      for (ai = beginAtoms(); ai != endAtoms(); ai++) {
        res += (*ai)->getTotalNumHs();
      }
    }
    return res;
  }


  Bond *ROMol::getBondWithIdx(unsigned int idx){
    PRECONDITION(getNumBonds()>0,"no bonds");
    RANGE_CHECK(0,idx,getNumBonds()-1);

    BOND_ITER_PAIR bIter = getEdges();
    for(unsigned int i=0;i<idx;i++) ++bIter.first;
    Bond *res = d_graph[*(bIter.first)].get();

    POSTCONDITION(res!=0,"Invalid bond requested");
    return res;
  }
  

  const Bond * ROMol::getBondWithIdx(unsigned int idx) const {
    PRECONDITION(getNumBonds()>0,"no bonds");
    RANGE_CHECK(0,idx,getNumBonds()-1);

    BOND_ITER_PAIR bIter = getEdges();
    for(unsigned int i=0;i<idx;i++) ++bIter.first;
    const Bond *res = d_graph[*(bIter.first)].get();

    POSTCONDITION(res!=0,"Invalid bond requested");
    return res;
  }
  

  Bond * ROMol::getBondBetweenAtoms(unsigned int idx1,unsigned int idx2){
    RANGE_CHECK(0,idx1,getNumAtoms()-1);
    RANGE_CHECK(0,idx2,getNumAtoms()-1);
    Bond *res=0;

    MolGraph::edge_descriptor edge;
    bool found;
    boost::tie(edge,found) = boost::edge(boost::vertex(idx1,d_graph),
                                         boost::vertex(idx2,d_graph),
                                         d_graph);
    if(found){
      res = d_graph[edge].get();
    }
    return res;
  }

  const Bond * ROMol::getBondBetweenAtoms(unsigned int idx1,unsigned int idx2) const{
    RANGE_CHECK(0,idx1,getNumAtoms()-1);
    RANGE_CHECK(0,idx2,getNumAtoms()-1);
    const Bond *res=0;

    MolGraph::edge_descriptor edge;
    bool found;
    boost::tie(edge,found) = boost::edge(boost::vertex(idx1,d_graph),
                                         boost::vertex(idx2,d_graph),
                                         d_graph);
    if(found){
      res = d_graph[edge].get();
    }
    return res;
  }
  
  

  ROMol::ADJ_ITER_PAIR ROMol::getAtomNeighbors(Atom const *at) const {
    return boost::adjacent_vertices(at->getIdx(),d_graph);
  };
  ROMol::ADJ_ITER_PAIR ROMol::getAtomNeighbors(Atom::ATOM_SPTR at) const {
    return boost::adjacent_vertices(at->getIdx(),d_graph);
  };
  ROMol::OBOND_ITER_PAIR ROMol::getAtomBonds(Atom const *at) const {
    return boost::out_edges(at->getIdx(),d_graph);
  }

  ROMol::ATOM_ITER_PAIR ROMol::getVertices() { return boost::vertices(d_graph);}
  ROMol::BOND_ITER_PAIR ROMol::getEdges() {return boost::edges(d_graph);}
  ROMol::ATOM_ITER_PAIR ROMol::getVertices() const { return boost::vertices(d_graph);}
  ROMol::BOND_ITER_PAIR ROMol::getEdges() const {return boost::edges(d_graph);}


  unsigned int ROMol::addAtom(Atom *atom_pin,bool updateLabel,bool takeOwnership){
    PRECONDITION(atom_pin,"null atom passed in");
    Atom *atom_p;
    if(!takeOwnership) atom_p = atom_pin->copy();
    else atom_p = atom_pin;

    atom_p->setOwningMol(this);
    MolGraph::vertex_descriptor which=boost::add_vertex(d_graph);
    d_graph[which].reset(atom_p);
    atom_p->setIdx(which);
    if(updateLabel){
      if(hasAtomBookmark(ci_RIGHTMOST_ATOM)) clearAtomBookmark(ci_RIGHTMOST_ATOM);
      setAtomBookmark(atom_p,ci_RIGHTMOST_ATOM);
    }
    ConformerIterator cfi;
    // EFF: I think this is going to cause a lot of re-allocation when
    // building molecules from SMILES
    for (cfi = this->beginConformers(); cfi != this->endConformers(); cfi++) {
      (*cfi)->setAtomPos(which, RDGeom::Point3D(0.0, 0.0, 0.0));
    }
    return which;
  };
  unsigned int ROMol::addAtom(Atom::ATOM_SPTR atom_sp,bool updateLabel)
  {
    return addAtom(atom_sp.get(),updateLabel,false);
  }
  unsigned int ROMol::addBond(Bond *bond_pin,bool takeOwnership){
    PRECONDITION(bond_pin,"null bond passed in");
    RANGE_CHECK(0,bond_pin->getBeginAtomIdx(),getNumAtoms()-1);
    RANGE_CHECK(0,bond_pin->getEndAtomIdx(),getNumAtoms()-1);
    PRECONDITION(bond_pin->getBeginAtomIdx()!=bond_pin->getEndAtomIdx(),
                 "attempt to add self-bond");
    PRECONDITION(!getBondBetweenAtoms(bond_pin->getBeginAtomIdx(),
                                      bond_pin->getEndAtomIdx()),"bond already exists");
    Bond *bond_p;
    if(!takeOwnership) bond_p = bond_pin->copy();
    else bond_p = bond_pin;

    bond_p->setOwningMol(this);
    bool ok;
    MolGraph::edge_descriptor which;
    boost::tie(which,ok) = boost::add_edge(bond_p->getBeginAtomIdx(),bond_p->getEndAtomIdx(),d_graph);
    CHECK_INVARIANT(ok,"bond could not be added");
    d_graph[which].reset(bond_p);
    int res = boost::num_edges(d_graph);
    bond_p->setIdx(res-1);
    return res;
  }
  unsigned int ROMol::addBond(Bond::BOND_SPTR bsp){
    return addBond(bsp.get());
  }


  void ROMol::debugMol(std::ostream &str) const{
    ATOM_ITER_PAIR atItP = getVertices();
    BOND_ITER_PAIR bondItP = getEdges();

    str << "Atoms:" << std::endl;
    while(atItP.first != atItP.second ){
      str << "\t" << *d_graph[*(atItP.first++)].get() << std::endl;
    }

    str << "Bonds:" << std::endl;
    while(bondItP.first != bondItP.second ){
      str << "\t" << *d_graph[*(bondItP.first++)].get() << std::endl;
    }
  }


  // --------------------------------------------
  //
  //  Iterators
  //
  // --------------------------------------------
  ROMol::AtomIterator ROMol::beginAtoms() {
    return AtomIterator(this);
  }
  ROMol::ConstAtomIterator ROMol::beginAtoms() const {
    return ConstAtomIterator(this);
  }
  ROMol::AtomIterator ROMol::endAtoms(){
    return AtomIterator(this,getNumAtoms());
  }
  ROMol::ConstAtomIterator ROMol::endAtoms() const{
    return ConstAtomIterator(this,getNumAtoms());
  }

  ROMol::AromaticAtomIterator ROMol::beginAromaticAtoms(){
    return AromaticAtomIterator(this);
  }
  ROMol::ConstAromaticAtomIterator ROMol::beginAromaticAtoms() const{
    return ConstAromaticAtomIterator(this);
  }
  ROMol::AromaticAtomIterator ROMol::endAromaticAtoms(){
    return AromaticAtomIterator(this,getNumAtoms());
  }
  ROMol::ConstAromaticAtomIterator ROMol::endAromaticAtoms() const{
    return ConstAromaticAtomIterator(this,getNumAtoms());
  }

  ROMol::HeteroatomIterator ROMol::beginHeteros(){
    return HeteroatomIterator(this);
  }
  ROMol::ConstHeteroatomIterator ROMol::beginHeteros() const{
    return ConstHeteroatomIterator(this);
  }
  ROMol::HeteroatomIterator ROMol::endHeteros(){
    return HeteroatomIterator(this,getNumAtoms());
  }
  ROMol::ConstHeteroatomIterator ROMol::endHeteros() const{
    return ConstHeteroatomIterator(this,getNumAtoms());
  }

  ROMol::QueryAtomIterator ROMol::beginQueryAtoms(QueryAtom const *what) {
    return QueryAtomIterator(this,what);
  }
  ROMol::ConstQueryAtomIterator ROMol::beginQueryAtoms(QueryAtom const *what) const {
    return ConstQueryAtomIterator(this,what);
  }
  ROMol::QueryAtomIterator ROMol::endQueryAtoms(){
    return QueryAtomIterator(this,getNumAtoms());
  }
  ROMol::ConstQueryAtomIterator ROMol::endQueryAtoms() const{
    return ConstQueryAtomIterator(this,getNumAtoms());
  }

  ROMol::BondIterator ROMol::beginBonds(){
    return BondIterator(this);
  }
  ROMol::ConstBondIterator ROMol::beginBonds() const {
    return ConstBondIterator(this);
  }
  ROMol::BondIterator ROMol::endBonds(){
    EDGE_ITER beg,end;
    boost::tie(beg,end)=getEdges();
    return BondIterator(this,end);
  }
  ROMol::ConstBondIterator ROMol::endBonds() const{
    EDGE_ITER beg,end;
    boost::tie(beg,end)=getEdges();
    return ConstBondIterator(this,end);
  }


  void ROMol::clearComputedProps(bool includeRings) const {
    // the SSSR information:
    if(includeRings) this->dp_ringInfo->reset();

    STR_VECT compLst;
    getProp("computedProps", compLst);
    STR_VECT_CI svi;
    for (svi = compLst.begin(); svi != compLst.end(); svi++) {
      dp_props->clearVal(*svi);
    }
    compLst.clear();
    dp_props->setVal("computedProps", compLst);
    for(ConstAtomIterator atomIt=this->beginAtoms();
        atomIt!=this->endAtoms();
        atomIt++){
      (*atomIt)->clearComputedProps();
    }
    for(ConstBondIterator bondIt=this->beginBonds();
        bondIt!=this->endBonds();
        bondIt++){
      (*bondIt)->clearComputedProps();
    }
  }

  void ROMol::updatePropertyCache(bool strict) {
    for(AtomIterator atomIt=this->beginAtoms();
        atomIt!=this->endAtoms();
        atomIt++){
      (*atomIt)->updatePropertyCache(strict);
    }
    for(BondIterator bondIt=this->beginBonds();
        bondIt!=this->endBonds();
        bondIt++){
      (*bondIt)->updatePropertyCache(strict);
    }
  }

  const Conformer &ROMol::getConformer(int id) const {
    // make sure we have more than one conformation
    if (d_confs.size() == 0) {
      throw ConformerException("No conformations available on the molecule");
    }
    
    ConstConformerIterator ci;
    if (id < 0) {
      return *(d_confs.front());
    }
    unsigned int cid = (unsigned int)id;
    for (ci = this->beginConformers(); ci != this->endConformers(); ci++) {
      if ((*ci)->getId() == cid) {
        return *(*ci);
      }
    }
    // we did not find a coformation with the specified ID
    std::string mesg = "Can't find conformation with ID: ";
    mesg += id;
    throw ConformerException(mesg);
  } 

  Conformer &ROMol::getConformer(int id) {
    // make sure we have more than one conformation
    if (d_confs.size() == 0) {
      throw ConformerException("No conformations available on the molecule");
    }
    
    if (id < 0) {
      return *(d_confs.front());
    }
    ConformerIterator ci;
    unsigned int cid = (unsigned int)id;
    for (ci = this->beginConformers(); ci != this->endConformers(); ci++) {
      if ((*ci)->getId() == cid) {
        return *(*ci);
      }
    }
    // we did not find a coformation with the specified ID
    std::string mesg = "Can't find conformation with ID: ";
    mesg += id;
    throw ConformerException(mesg);
  }

  void ROMol::removeConformer(unsigned int id) {
    CONF_SPTR_LIST_I ci;
    for (ci = d_confs.begin(); ci != d_confs.end(); ci++) {
      if ((*ci)->getId() == id) {
        d_confs.erase(ci);
        return;
      }
    }
  }

  unsigned int ROMol::addConformer(Conformer * conf, bool assignId) {
    PRECONDITION(conf->getNumAtoms() == this->getNumAtoms(), "Number of atom mismatch")
      if (assignId) {
        int maxId = -1;
        CONF_SPTR_LIST_I ci;
        for (ci = d_confs.begin(); ci != d_confs.end(); ci++) {
          if ((int)(*ci)->getId() > maxId) {
            maxId = (int)(*ci)->getId();
          }
        }
        maxId++;
        conf->setId((unsigned int)maxId);
      }
    conf->setOwningMol(this);
    CONFORMER_SPTR nConf(conf);
    d_confs.push_back(nConf);
    return conf->getId();
  }
} // end o' namespace
