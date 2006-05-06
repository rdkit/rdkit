// $Id: AtomIterators.cpp 4961 2006-02-18 00:14:47Z glandrum $
//
//  Copyright (C) 2002-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include "AtomIterators.h"
#include "RDKitBase.h"
#include "RDKitQueries.h"

namespace RDKit{
  template <class Atom_, class Mol_, class PMAP_>
  AtomIterator_<Atom_, Mol_, PMAP_>::AtomIterator_(Mol_ * mol) {
    _mol = mol;
    _pMap = mol->getAtomPMap();
    _pos = 0;
    _max = mol->getNumAtoms();
  };
  template <class Atom_, class Mol_, class PMAP_>
  AtomIterator_<Atom_, Mol_, PMAP_>::AtomIterator_(Mol_ * mol,int pos) {
    _mol = mol;
    _pMap = mol->getAtomPMap();
    _pos = pos;
    _max = mol->getNumAtoms();
  };
  template <class Atom_, class Mol_, class PMAP_>
  AtomIterator_<Atom_, Mol_, PMAP_>::AtomIterator_(const AtomIterator_<Atom_, Mol_, PMAP_> &other){
    _mol = other._mol;
    _pMap = other._pMap;
    _pos = other._pos;
    _max = other._max;
  }
  template <class Atom_, class Mol_, class PMAP_>
  AtomIterator_<Atom_, Mol_, PMAP_> &AtomIterator_<Atom_, Mol_, PMAP_>::operator=(const AtomIterator_<Atom_, Mol_, PMAP_> &other){
    _mol = other._mol;
    _pMap = other._pMap;
    _pos = other._pos;
    _max = other._max;
    return *this;
  }


  template <class Atom_, class Mol_, class PMAP_>
  AtomIterator_<Atom_, Mol_, PMAP_> &AtomIterator_<Atom_, Mol_, PMAP_>::operator+=(int val){
    _pos += val;
    if(_pos<0 || _pos>_max) _pos = _max;
    return *this;
  }
  template <class Atom_, class Mol_, class PMAP_>
  AtomIterator_<Atom_, Mol_, PMAP_> &AtomIterator_<Atom_, Mol_, PMAP_>::operator-=(int val){
    _pos -= val;
    if(_pos<0 || _pos>_max) _pos = _max;
    return *this;
  }
  template <class Atom_, class Mol_, class PMAP_>
  AtomIterator_<Atom_, Mol_, PMAP_> AtomIterator_<Atom_, Mol_, PMAP_>::operator+(int val){
    AtomIterator_<Atom_, Mol_, PMAP_> res(*this);
    res += val;
    // += takes care of the pre/post conditions for us, so we're safe to return
    return res;
  }
  template <class Atom_, class Mol_, class PMAP_>
  AtomIterator_<Atom_, Mol_, PMAP_> AtomIterator_<Atom_, Mol_, PMAP_>::operator-(int val){
    AtomIterator_<Atom_, Mol_, PMAP_> res(*this);
    // -= takes care of the pre/post conditions for us, so we're safe to return
    res -= val;
    return res;
  }

  // iterator subtraction
  template <class Atom_, class Mol_, class PMAP_>
  int AtomIterator_<Atom_, Mol_, PMAP_>::operator-(AtomIterator_<Atom_, Mol_, PMAP_> &other){
    PRECONDITION(_mol==other._mol,"bad operator- call");
    return _pos - other._pos;
  }

  // dereference 
  template <class Atom_, class Mol_, class PMAP_>
  Atom_ * AtomIterator_<Atom_, Mol_, PMAP_>::operator*() {
    RANGE_CHECK(0,_pos,_max-1);
    return _pMap[_pos];
  }
  // random access
  template <class Atom_, class Mol_, class PMAP_>
  Atom_ * AtomIterator_<Atom_, Mol_, PMAP_>::operator[](const int which){
    RANGE_CHECK(0,which,_max-1);
    return _pMap[which];
  }

  template <class Atom_, class Mol_, class PMAP_>
  bool AtomIterator_<Atom_, Mol_, PMAP_>::operator==(const AtomIterator_<Atom_, Mol_, PMAP_> &other){
    return _mol==other._mol && _pos==other._pos;
  }
  template <class Atom_, class Mol_, class PMAP_>
  bool AtomIterator_<Atom_, Mol_, PMAP_>::operator!=(const AtomIterator_<Atom_, Mol_, PMAP_> &other){
    return _mol!=other._mol || _pos!=other._pos;
  }
  template <class Atom_, class Mol_, class PMAP_>
  bool AtomIterator_<Atom_, Mol_, PMAP_>::operator<(const AtomIterator_<Atom_, Mol_, PMAP_> &other){
    return _mol==other._mol && _pos<other._pos;
  }
  template <class Atom_, class Mol_, class PMAP_>
  bool AtomIterator_<Atom_, Mol_, PMAP_>::operator<=(const AtomIterator_<Atom_, Mol_, PMAP_> &other){
    return _mol==other._mol && _pos<=other._pos;
  }
  template <class Atom_, class Mol_, class PMAP_>
  bool AtomIterator_<Atom_, Mol_, PMAP_>::operator>(const AtomIterator_<Atom_, Mol_, PMAP_> &other){
    return _mol==other._mol && _pos>other._pos;
  }
  template <class Atom_, class Mol_, class PMAP_>
  bool AtomIterator_<Atom_, Mol_, PMAP_>::operator>=(const AtomIterator_<Atom_, Mol_, PMAP_> &other){
    return _mol==other._mol && _pos>=other._pos;
  }


  // pre-increment
  template <class Atom_, class Mol_, class PMAP_>
  AtomIterator_<Atom_, Mol_, PMAP_> &AtomIterator_<Atom_, Mol_, PMAP_>::operator++() {
    RANGE_CHECK(0,_pos,_max-1);
    _pos++;
    return *this;
  }
  template <class Atom_, class Mol_, class PMAP_>
  AtomIterator_<Atom_, Mol_, PMAP_> AtomIterator_<Atom_, Mol_, PMAP_>::operator++(int) {
    RANGE_CHECK(0,_pos,_max-1);
    AtomIterator_<Atom_,Mol_,PMAP_> res(*this);
    _pos++;
    return res;
  }
  // pre-decrement
  template <class Atom_, class Mol_, class PMAP_>
  AtomIterator_<Atom_, Mol_, PMAP_> &AtomIterator_<Atom_, Mol_, PMAP_>::operator--() {
    _pos--;
    RANGE_CHECK(0,_pos,_max-1);
    return *this;
  }
  template <class Atom_, class Mol_, class PMAP_>
  AtomIterator_<Atom_, Mol_, PMAP_> AtomIterator_<Atom_, Mol_, PMAP_>::operator--(int) {
    RANGE_CHECK(0,_pos,_max-1);
    AtomIterator_<Atom_,Mol_,PMAP_> res(*this);
    if(_pos-1 < 0) _pos = _max;
    else _pos--;
    return res;
  }



  //-----------------------------------------
  //
  //  HeteroatomIterator
  //
  //-----------------------------------------
  template <class Atom_, class Mol_, class PMAP_>
  HeteroatomIterator_<Atom_,Mol_,PMAP_>::HeteroatomIterator_(Mol_ * mol){
    _mol = mol;
    _qA = new QueryAtom(6);
    _qA->getQuery()->setNegation(true);
    _pMap = mol->getAtomPMap();
    _end = mol->getNumAtoms();
    _pos = _findNext(0);
  };
  template <class Atom_, class Mol_, class PMAP_>
  HeteroatomIterator_<Atom_,Mol_,PMAP_>::HeteroatomIterator_(Mol_ * mol,int pos){
    _mol = mol;
    _qA = new QueryAtom(6);
    _qA->getQuery()->setNegation(true);
    _pMap = mol->getAtomPMap();
    _end =  mol->getNumAtoms();
    _pos = pos;
  };

  template <class Atom_, class Mol_, class PMAP_>
  HeteroatomIterator_<Atom_,Mol_,PMAP_>::~HeteroatomIterator_() {
    if (_qA) delete _qA;
  }

  template <class Atom_, class Mol_, class PMAP_>
  HeteroatomIterator_<Atom_,Mol_,PMAP_>::HeteroatomIterator_(const ThisType &other){
    _mol = other._mol;
    _end = other._end;
    _pMap = other._pMap;
    _pos = other._pos;
    //_qA = new QueryAtom(6);
    //_qA->getQuery()->setNegation(true);
    _qA = static_cast<QueryAtom *>(other._qA->copy());
  }

  template <class Atom_, class Mol_, class PMAP_>
  HeteroatomIterator_<Atom_,Mol_,PMAP_> &HeteroatomIterator_<Atom_,Mol_,PMAP_>::operator=(const ThisType &other){
    _mol = other._mol;
    _end = other._end;
    _pMap = other._pMap;
    _pos = other._pos;
    _qA = static_cast<QueryAtom *>(other._qA->copy());
    return *this;
  }


  template <class Atom_, class Mol_, class PMAP_>
  bool HeteroatomIterator_<Atom_,Mol_,PMAP_>::operator==(const ThisType &other){
    return _mol==other._mol && _pos==other._pos;
  }
  template <class Atom_, class Mol_, class PMAP_>
  bool HeteroatomIterator_<Atom_,Mol_,PMAP_>::operator!=(const ThisType &other){
    return _mol!=other._mol || _pos!=other._pos;
  }

  template <class Atom_, class Mol_, class PMAP_>
  Atom_ * HeteroatomIterator_<Atom_,Mol_,PMAP_>::operator*() {
    return _pMap[_pos];
  }
  // pre-increment
  template <class Atom_, class Mol_, class PMAP_>
  HeteroatomIterator_<Atom_,Mol_,PMAP_> &HeteroatomIterator_<Atom_,Mol_,PMAP_>::operator++() {
    _pos = _findNext(_pos+1);
    RANGE_CHECK(0,_pos,_end-1);
    return *this;
  }
  template <class Atom_, class Mol_, class PMAP_>
  HeteroatomIterator_<Atom_,Mol_,PMAP_> HeteroatomIterator_<Atom_,Mol_,PMAP_>::operator++(int) {
    HeteroatomIterator_<Atom_,Mol_,PMAP_> res(*this); 
    _pos = _findNext(_pos+1);
    return res;
  }
  // pre-decrement
  template <class Atom_, class Mol_, class PMAP_>
  HeteroatomIterator_<Atom_,Mol_,PMAP_> &HeteroatomIterator_<Atom_,Mol_,PMAP_>::operator--() {
    _pos = _findPrev(_pos-1);
    RANGE_CHECK(0,_pos,_end-1);
    return *this;
  }
  template <class Atom_, class Mol_, class PMAP_>
  HeteroatomIterator_<Atom_,Mol_,PMAP_>
  HeteroatomIterator_<Atom_,Mol_,PMAP_>::operator--(int) {
    HeteroatomIterator_<Atom_,Mol_,PMAP_> res(*this); 
    _pos = _findPrev(_pos-1);
    return res;
  }
  template <class Atom_, class Mol_, class PMAP_>
  int HeteroatomIterator_<Atom_,Mol_,PMAP_>::_findNext(int from){
    while(from<_end){
      if(_qA->Match(_pMap[from])) break;
      else from++;
    } 
    return from;
  }
  
  template <class Atom_, class Mol_, class PMAP_>
  int HeteroatomIterator_<Atom_,Mol_,PMAP_>::_findPrev(int from){
    while(from>0){
      if(_qA->Match(_pMap[from])) break;
      else from--;
    }
    if(from<0) from = _end;
    return from;
  }


  //-----------------------------------------
  //
  //  AromaticAtomIterator
  //
  //-----------------------------------------
  template <class Atom_, class Mol_, class PMAP_>
  AromaticAtomIterator_<Atom_,Mol_,PMAP_>::AromaticAtomIterator_(Mol_ * mol){
    _mol = mol;
    _pMap = mol->getAtomPMap();
    _end = mol->getNumAtoms();
    _pos = _findNext(0);
  };
  template <class Atom_, class Mol_, class PMAP_>
  AromaticAtomIterator_<Atom_,Mol_,PMAP_>::AromaticAtomIterator_(Mol_ * mol,int pos){
    _mol = mol;
    _pMap = mol->getAtomPMap();
    _end =  mol->getNumAtoms();
    _pos = pos;
  };

  template <class Atom_, class Mol_, class PMAP_>
  AromaticAtomIterator_<Atom_,Mol_,PMAP_>::~AromaticAtomIterator_() {
  }

  template <class Atom_, class Mol_, class PMAP_>
  AromaticAtomIterator_<Atom_,Mol_,PMAP_>::AromaticAtomIterator_(const ThisType &other){
    _mol = other._mol;
    _end = other._end;
    _pMap = other._pMap;
    _pos = other._pos;
  }

  template <class Atom_, class Mol_, class PMAP_>
  AromaticAtomIterator_<Atom_,Mol_,PMAP_> &AromaticAtomIterator_<Atom_,Mol_,PMAP_>::operator=(const ThisType &other){
    _mol = other._mol;
    _end = other._end;
    _pMap = other._pMap;
    _pos = other._pos;
    return *this;
  }


  template <class Atom_, class Mol_, class PMAP_>
  bool AromaticAtomIterator_<Atom_,Mol_,PMAP_>::operator==(const ThisType &other){
    return _mol==other._mol && _pos==other._pos;
  }
  template <class Atom_, class Mol_, class PMAP_>
  bool AromaticAtomIterator_<Atom_,Mol_,PMAP_>::operator!=(const ThisType &other){
    return _mol!=other._mol || _pos!=other._pos;
  }

  template <class Atom_, class Mol_, class PMAP_>
  Atom_ * AromaticAtomIterator_<Atom_,Mol_,PMAP_>::operator*() {
    return _pMap[_pos];
  }
  // pre-increment
  template <class Atom_, class Mol_, class PMAP_>
  AromaticAtomIterator_<Atom_,Mol_,PMAP_> &AromaticAtomIterator_<Atom_,Mol_,PMAP_>::operator++() {
    _pos = _findNext(_pos+1);
    RANGE_CHECK(0,_pos,_end-1);
    return *this;
  }
  template <class Atom_, class Mol_, class PMAP_>
  AromaticAtomIterator_<Atom_,Mol_,PMAP_> AromaticAtomIterator_<Atom_,Mol_,PMAP_>::operator++(int) {
    AromaticAtomIterator_<Atom_,Mol_,PMAP_> res(*this); 
    _pos = _findNext(_pos+1);
    return res;
  }
  // pre-decrement
  template <class Atom_, class Mol_, class PMAP_>
  AromaticAtomIterator_<Atom_,Mol_,PMAP_> &AromaticAtomIterator_<Atom_,Mol_,PMAP_>::operator--() {
    _pos = _findPrev(_pos-1);
    RANGE_CHECK(0,_pos,_end-1);
    return *this;
  }
  template <class Atom_, class Mol_, class PMAP_>
  AromaticAtomIterator_<Atom_,Mol_,PMAP_>
  AromaticAtomIterator_<Atom_,Mol_,PMAP_>::operator--(int) {
    AromaticAtomIterator_<Atom_,Mol_,PMAP_> res(*this); 
    _pos = _findPrev(_pos-1);
    return res;
  }
  template <class Atom_, class Mol_, class PMAP_>
  int AromaticAtomIterator_<Atom_,Mol_,PMAP_>::_findNext(int from){
    while(from<_end){
      if(_pMap[from]->getIsAromatic()) break;
      else from++;
    } 
    return from;
  }
  
  template <class Atom_, class Mol_, class PMAP_>
  int AromaticAtomIterator_<Atom_,Mol_,PMAP_>::_findPrev(int from){
    while(from>0){
      if(_pMap[from]->getIsAromatic()) break;
      else from--;
    }
    if(from<0) from = _end;
    return from;
  }






  
  //-----------------------------------------
  //
  //  QueryAtomIterator
  //
  //-----------------------------------------
  template <class Atom_, class Mol_, class PMAP_>
  QueryAtomIterator_<Atom_,Mol_,PMAP_>::QueryAtomIterator_(Mol_ * mol,QueryAtom const *what){
    _mol = mol;
    _qA = static_cast<QueryAtom *>(what->copy());
    _pMap = mol->getAtomPMap();
    _end = mol->getNumAtoms();
    _pos = _findNext(0);
  };
  template <class Atom_, class Mol_, class PMAP_>
  QueryAtomIterator_<Atom_,Mol_,PMAP_>::QueryAtomIterator_(Mol_ * mol,int pos){
    _mol = mol;
    _qA = NULL;
    _pMap = mol->getAtomPMap();
    _end = mol->getNumAtoms();
    _pos = pos;
  };
  template <class Atom_, class Mol_, class PMAP_>
  QueryAtomIterator_<Atom_,Mol_,PMAP_>::~QueryAtomIterator_() {
    if(_qA) delete _qA;
  }
  template <class Atom_, class Mol_, class PMAP_>
  QueryAtomIterator_<Atom_,Mol_,PMAP_>::QueryAtomIterator_(const QueryAtomIterator_<Atom_,Mol_,PMAP_> &other){
    _mol = other._mol;
    _pMap = other._pMap;
    _pos = other._pos;
    _end = other._end;
    _qA = static_cast<QueryAtom *>(other._qA->copy());
  }

  template <class Atom_, class Mol_, class PMAP_>
  QueryAtomIterator_<Atom_,Mol_,PMAP_> &
  QueryAtomIterator_<Atom_,Mol_,PMAP_>::operator =(const QueryAtomIterator_<Atom_,Mol_,PMAP_> &other){
    _mol = other._mol;
    _pMap = other._pMap;
    _pos = other._pos;
    _end = other._end;
    _qA = static_cast<QueryAtom *>(other._qA->copy());
    return *this;
  }
  template <class Atom_, class Mol_, class PMAP_>
  bool QueryAtomIterator_<Atom_,Mol_,PMAP_>::operator==(const QueryAtomIterator_<Atom_,Mol_,PMAP_> &other){
    return _mol==other._mol && _pos==other._pos;
  }
  template <class Atom_, class Mol_, class PMAP_>
  bool QueryAtomIterator_<Atom_,Mol_,PMAP_>::operator!=(const QueryAtomIterator_<Atom_,Mol_,PMAP_> &other){
    return _mol!=other._mol || _pos!=other._pos;
  }

  template <class Atom_, class Mol_, class PMAP_>
  Atom_ * QueryAtomIterator_<Atom_,Mol_,PMAP_>::operator*() {
    return _pMap[_pos];
  }
  // pre-increment
  template <class Atom_, class Mol_, class PMAP_>
  QueryAtomIterator_<Atom_,Mol_,PMAP_> &QueryAtomIterator_<Atom_,Mol_,PMAP_>::operator++() {
    _pos = _findNext(_pos+1);
    RANGE_CHECK(0,_pos,_end-1);
    return *this;
  }
  template <class Atom_, class Mol_, class PMAP_>
  QueryAtomIterator_<Atom_,Mol_,PMAP_> QueryAtomIterator_<Atom_,Mol_,PMAP_>::operator++(int) {
    QueryAtomIterator_ res(*this); 
    _pos = _findNext(_pos+1);
    return res;
  }
  // pre-decrement
  template <class Atom_, class Mol_, class PMAP_>
  QueryAtomIterator_<Atom_,Mol_,PMAP_> &QueryAtomIterator_<Atom_,Mol_,PMAP_>::operator--() {
    _pos = _findPrev(_pos-1);
    RANGE_CHECK(0,_pos,_end-1);
    return *this;
  }
  template <class Atom_, class Mol_, class PMAP_>
  QueryAtomIterator_<Atom_,Mol_,PMAP_> QueryAtomIterator_<Atom_,Mol_,PMAP_>::operator--(int) {
    QueryAtomIterator_<Atom_,Mol_,PMAP_> res(*this); 
    _pos = _findPrev(_pos-1);
    return res;
  }
  template <class Atom_, class Mol_, class PMAP_>
  int QueryAtomIterator_<Atom_,Mol_,PMAP_>::_findNext(int from){
    PRECONDITION(_qA!=NULL,"no query set");
    while(from<_end){
      if(_qA->Match(_pMap[from])) break;
      else from++;
    }
    return from;
  }
  
  template <class Atom_, class Mol_, class PMAP_>
  int QueryAtomIterator_<Atom_,Mol_,PMAP_>::_findPrev(int from){
    PRECONDITION(_qA!=NULL,"no query set");
    while(from>0){
      if(_qA->Match(_pMap[from])) break;
      else from--;
    }
    if(from<0) from = _end;
    return from;
  }

  template <class T>
  void forceIt(T &arg){
    ROMol *m;
    T arg2(m,2);
    *arg;
    arg++;
    ++arg;
    arg--;
    --arg;
    arg==arg2;
    if(arg != arg2){
      arg!=++arg;
      arg=++arg;
    }
  };

  template class AtomIterator_<Atom,ROMol,ROMol::GRAPH_MOL_ATOM_PMAP::type>;
  template class AtomIterator_<const Atom,const ROMol,ROMol::GRAPH_MOL_ATOM_PMAP::const_type>;
  template class AromaticAtomIterator_<Atom,ROMol,ROMol::GRAPH_MOL_ATOM_PMAP::type>;
  template class AromaticAtomIterator_<const Atom,const ROMol,ROMol::GRAPH_MOL_ATOM_PMAP::const_type>;
  template class HeteroatomIterator_<Atom,ROMol,ROMol::GRAPH_MOL_ATOM_PMAP::type>;
  template class HeteroatomIterator_<const Atom,const ROMol,ROMol::GRAPH_MOL_ATOM_PMAP::const_type>;
  template class QueryAtomIterator_<Atom,ROMol,ROMol::GRAPH_MOL_ATOM_PMAP::type>;
  template class QueryAtomIterator_<const Atom,const ROMol,ROMol::GRAPH_MOL_ATOM_PMAP::const_type>;


}; // end o' namespace


