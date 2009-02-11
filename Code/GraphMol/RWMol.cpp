// $Id$
//
//  Copyright (C) 2003-2009 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

// our stuff
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include "RWMol.h"
#include "Atom.h"
#include "Bond.h"
#include "BondIterators.h"



namespace RDKit{
  void RWMol::destroy() {
    ROMol::destroy();
    d_partialBonds.clear();
    d_partialBonds.resize(0);
  };

  void RWMol::insertMol(const ROMol &other)
  {

    VERTEX_ITER firstA,lastA;
    boost::tie(firstA,lastA) = boost::vertices(other.d_graph);
    while(firstA!=lastA){
      Atom *newAt = other.d_graph[*firstA]->copy();
      addAtom(newAt,false,true);
      // set a bookmark for this atom's original index so we can update bonds
      setAtomBookmark(newAt,ci_ATOM_HOLDER + other.d_graph[*firstA]->getIdx());
      ++firstA;
    }

    EDGE_ITER firstB,lastB;
    boost::tie(firstB,lastB) = boost::edges(other.d_graph);
    while(firstB != lastB){
      Bond *bond_p = other.d_graph[*firstB]->copy();
      unsigned int idx1,idx2;
      idx1 = getAtomWithBookmark(ci_ATOM_HOLDER + bond_p->getBeginAtomIdx())->getIdx();
      idx2 = getAtomWithBookmark(ci_ATOM_HOLDER + bond_p->getEndAtomIdx())->getIdx();
      bond_p->setOwningMol(this);
      bond_p->setBeginAtomIdx(idx1);
      bond_p->setEndAtomIdx(idx2);
      addBond(bond_p,true);
      ++firstB;
    }
    // blow out those bookmarks we set
    boost::tie(firstA,lastA) = boost::vertices(other.d_graph);
    while(firstA!=lastA){
      clearAtomBookmark(ci_ATOM_HOLDER + *firstA);
      ++firstA;
    }
  }

  unsigned int RWMol::addAtom(bool updateLabel){
    Atom *atom_p = new Atom();
    atom_p->setOwningMol(this);
    MolGraph::vertex_descriptor which = boost::add_vertex(d_graph);
    d_graph[which].reset(atom_p);
    atom_p->setIdx(which);
    if(updateLabel){
      if(hasAtomBookmark(ci_RIGHTMOST_ATOM)) clearAtomBookmark(ci_RIGHTMOST_ATOM);
      setAtomBookmark(atom_p,ci_RIGHTMOST_ATOM);
    }
    return which;
  }
  
  void RWMol::replaceAtom(unsigned int idx,Atom *atom_pin,bool updateLabel){
    PRECONDITION(atom_pin,"bad atom passed to replaceAtom");
    RANGE_CHECK(0,idx,getNumAtoms()-1);
    Atom *atom_p = atom_pin->copy();
    atom_p->setOwningMol(this);
    atom_p->setIdx(idx);
    MolGraph::vertex_descriptor vd = boost::vertex(idx,d_graph);
    d_graph[vd].reset(atom_p);
    // FIX: do something about bookmarks
  };

  Atom *RWMol::getActiveAtom() {
    if(hasAtomBookmark(ci_RIGHTMOST_ATOM)) return getAtomWithBookmark(ci_RIGHTMOST_ATOM);
    else return getLastAtom();
  };

  void RWMol::setActiveAtom(Atom *at) {
    if(hasAtomBookmark(ci_RIGHTMOST_ATOM)) clearAtomBookmark(ci_RIGHTMOST_ATOM);
    setAtomBookmark(at,ci_RIGHTMOST_ATOM);
  };
  void RWMol::setActiveAtom(unsigned int idx) {
    setActiveAtom(getAtomWithIdx(idx));
  };


  void RWMol::removeAtom(unsigned int idx) {
    Atom *oatom = getAtomWithIdx(idx);
    // remove any bookmarks which point to this atom:
    ATOM_BOOKMARK_MAP *marks = getAtomBookmarks();
    ATOM_BOOKMARK_MAP::iterator markI=marks->begin();
    while(markI != marks->end()){
      ATOM_PTR_LIST &atoms=markI->second;
      // we need to copy the iterator then increment it, because the
      // deletion we're going to do in clearAtomBookmark will invalidate
      // it.
      ATOM_BOOKMARK_MAP::iterator tmpI=markI;
      markI++;
      if(std::find(atoms.begin(),atoms.end(),oatom)!=atoms.end()){
        clearAtomBookmark(tmpI->first,oatom);
      }
    }

    // remove bonds attached to the atom
    std::vector<std::pair<unsigned int,unsigned int> > nbrs;
    ADJ_ITER b1,b2;
    boost::tie(b1,b2)=getAtomNeighbors(oatom);
    while(b1!=b2){
      nbrs.push_back(std::make_pair(oatom->getIdx(),*b1));
      ++b1;
    }
    for(unsigned int i=0;i<nbrs.size();++i){
      removeBond(nbrs[i].first,nbrs[i].second);
    }

    // loop over all atoms with higher indices and update their indices
    for(unsigned int i=idx+1;i<getNumAtoms();i++){
      Atom *atom = getAtomWithIdx(i);
      atom->setIdx(i-1);
    }

    // do the same with the coordinates in the conformations
    CONF_SPTR_LIST_I ci;
    for (ci = d_confs.begin(); ci != d_confs.end(); ci++) {
      RDGeom::POINT3D_VECT &positions = (*ci)->getPositions();
      RDGeom::POINT3D_VECT_I pi = positions.begin();
      for (unsigned int i = 0; i < getNumAtoms()-1;i++) {
        pi++;
        if (i >= idx) {
          positions[i] = positions[i+1];
        }
      }
      positions.erase(pi);
    }
    // now deal with bonds:
    //   their end indices may need to be decremented and their
    //   indices will need to be handled
    BondIterator bondIt;
    unsigned int nBonds=0;
    for(bondIt=beginBonds();bondIt!=endBonds();bondIt++){
      Bond *bond = *bondIt;
      unsigned int tmpIdx = bond->getBeginAtomIdx();
      if( tmpIdx > idx) bond->setBeginAtomIdx(tmpIdx-1);
      tmpIdx = bond->getEndAtomIdx();
      if( tmpIdx > idx) bond->setEndAtomIdx(tmpIdx-1);
      bond->setIdx(nBonds++);
    }

    // remove all connections to the atom:
    MolGraph::vertex_descriptor vd = boost::vertex(idx,d_graph);
    boost::clear_vertex(vd,d_graph);
    // finally remove the vertex itself
    boost::remove_vertex(vd,d_graph);
  }

  void RWMol::removeAtom(Atom *atom) {
    removeAtom(atom->getIdx());
  }

  unsigned int RWMol::addBond(unsigned int atomIdx1,unsigned int atomIdx2,
                              Bond::BondType bondType){
    RANGE_CHECK(0,atomIdx1,getNumAtoms()-1);
    RANGE_CHECK(0,atomIdx2,getNumAtoms()-1);
    PRECONDITION(atomIdx1!=atomIdx2,"attempt to add self-bond");
    PRECONDITION(!getBondBetweenAtoms(atomIdx1,atomIdx2),"bond already exists");
    Bond *b = new Bond(bondType);
    b->setOwningMol(this);
    if(bondType==Bond::AROMATIC){
      b->setIsAromatic(1);
      //
      // assume that aromatic bonds connect aromatic atoms
      //   This is relevant for file formats like MOL, where there
      //   is no such thing as an aromatic atom, but bonds can be
      //   marked aromatic.
      //
      getAtomWithIdx(atomIdx1)->setIsAromatic(1);
      getAtomWithIdx(atomIdx2)->setIsAromatic(1);
    }
    bool ok;
    MolGraph::edge_descriptor which;
    boost::tie(which,ok) = boost::add_edge(atomIdx1,atomIdx2,d_graph);
    d_graph[which].reset(b);
    unsigned int res = boost::num_edges(d_graph);
    b->setIdx(res-1);
    b->setBeginAtomIdx(atomIdx1);
    b->setEndAtomIdx(atomIdx2);
    
    return res;
  }

  unsigned int RWMol::addBond(Atom *atom1,Atom *atom2,
                              Bond::BondType bondType){
    PRECONDITION(atom1&&atom2,"NULL atom passed in");
    return addBond(atom1->getIdx(),atom2->getIdx(),bondType);
  }

  unsigned int RWMol::addBond(Atom::ATOM_SPTR atom1,Atom::ATOM_SPTR atom2,
                              Bond::BondType bondType){
    return addBond(atom1->getIdx(),atom2->getIdx(),bondType);
  }

  void RWMol::removeBond(unsigned int aid1, unsigned int aid2) {
    RANGE_CHECK(0,aid1,getNumAtoms()-1);
    RANGE_CHECK(0,aid2,getNumAtoms()-1);
    Bond *bnd = getBondBetweenAtoms(aid1, aid2);
    if(!bnd) return;
    
    // remove any bookmarks which point to this bond:
    BOND_BOOKMARK_MAP *marks = getBondBookmarks();
    BOND_BOOKMARK_MAP::iterator markI=marks->begin();
    while(markI != marks->end()){
      BOND_PTR_LIST &bonds=markI->second;
      // we need to copy the iterator then increment it, because the
      // deletion we're going to do in clearBondBookmark will invalidate
      // it.
      BOND_BOOKMARK_MAP::iterator tmpI=markI;
      markI++;
      if(std::find(bonds.begin(),bonds.end(),bnd)!=bonds.end()){
        clearBondBookmark(tmpI->first,bnd);
      }
    }
    MolGraph::vertex_descriptor vd1 = boost::vertex(aid1,d_graph);
    MolGraph::vertex_descriptor vd2 = boost::vertex(aid2,d_graph);
    boost::remove_edge(vd1, vd2, d_graph);
  }

  Bond *RWMol::createPartialBond(unsigned int atomIdx1,Bond::BondType bondType){
    RANGE_CHECK(0,atomIdx1,getNumAtoms()-1);

    Bond *b = new Bond(bondType);
    b->setOwningMol(this);
    b->setBeginAtomIdx(atomIdx1);
    return b;
  }

  unsigned int RWMol::finishPartialBond(unsigned int atomIdx2,int bondBookmark,
                                        Bond::BondType bondType){
    PRECONDITION(hasBondBookmark(bondBookmark),"no such partial bond");
    RANGE_CHECK(0,atomIdx2,getNumAtoms()-1);

    Bond *bsp = getBondWithBookmark(bondBookmark);
    if(bondType==Bond::UNSPECIFIED){
      bondType = bsp->getBondType();
    }
    return addBond(bsp->getBeginAtomIdx(),atomIdx2,bondType);
  }

} // end o' namespace
