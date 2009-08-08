// $Id$
//
//  Copyright (C) 2002-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include "AtomIterators.h"
#include "RDKitBase.h"
#include "RDKitQueries.h"

namespace RDKit{
  template <class Atom_, class Mol_>
  AtomIterator_<Atom_, Mol_>::AtomIterator_(Mol_ * mol) {
    _mol = mol;
    _pos = 0;
    _max = mol->getNumAtoms();
  };
  template <class Atom_, class Mol_>
  AtomIterator_<Atom_, Mol_>::AtomIterator_(Mol_ * mol,int pos) {
    _mol = mol;
    _pos = pos;
    _max = mol->getNumAtoms();
  };
  template <class Atom_, class Mol_>
  AtomIterator_<Atom_, Mol_>::AtomIterator_(const AtomIterator_<Atom_, Mol_> &other){
    _mol = other._mol;
    _pos = other._pos;
    _max = other._max;
  }
  template <class Atom_, class Mol_>
  AtomIterator_<Atom_, Mol_> &AtomIterator_<Atom_, Mol_>::operator=(const AtomIterator_<Atom_, Mol_> &other){
    _mol = other._mol;
    _pos = other._pos;
    _max = other._max;
    return *this;
  }


  template <class Atom_, class Mol_>
  AtomIterator_<Atom_, Mol_> &AtomIterator_<Atom_, Mol_>::operator+=(int val){
    _pos += val;
    if(_pos<0 || _pos>_max) _pos = _max;
    return *this;
  }
  template <class Atom_, class Mol_>
  AtomIterator_<Atom_, Mol_> &AtomIterator_<Atom_, Mol_>::operator-=(int val){
    _pos -= val;
    if(_pos<0 || _pos>_max) _pos = _max;
    return *this;
  }
  template <class Atom_, class Mol_>
  AtomIterator_<Atom_, Mol_> AtomIterator_<Atom_, Mol_>::operator+(int val){
    AtomIterator_<Atom_, Mol_> res(*this);
    res += val;
    // += takes care of the pre/post conditions for us, so we're safe to return
    return res;
  }
  template <class Atom_, class Mol_>
  AtomIterator_<Atom_, Mol_> AtomIterator_<Atom_, Mol_>::operator-(int val){
    AtomIterator_<Atom_, Mol_> res(*this);
    // -= takes care of the pre/post conditions for us, so we're safe to return
    res -= val;
    return res;
  }

  // iterator subtraction
  template <class Atom_, class Mol_>
  int AtomIterator_<Atom_, Mol_>::operator-(AtomIterator_<Atom_, Mol_> &other){
    PRECONDITION(_mol==other._mol,"bad operator- call");
    return _pos - other._pos;
  }

  // dereference 
  template <class Atom_, class Mol_>
  Atom_ * AtomIterator_<Atom_, Mol_>::operator*() {
    RANGE_CHECK(0,_pos,_max-1);
    return (*_mol)[_pos].get();
  }
  // random access
  template <class Atom_, class Mol_>
  Atom_ * AtomIterator_<Atom_, Mol_>::operator[](const int which){
    RANGE_CHECK(0,which,_max-1);
    return (*_mol)[which].get();
  }

  template <class Atom_, class Mol_>
  bool AtomIterator_<Atom_, Mol_>::operator==(const AtomIterator_<Atom_, Mol_> &other){
    return _mol==other._mol && _pos==other._pos;
  }
  template <class Atom_, class Mol_>
  bool AtomIterator_<Atom_, Mol_>::operator!=(const AtomIterator_<Atom_, Mol_> &other){
    return _mol!=other._mol || _pos!=other._pos;
  }
  template <class Atom_, class Mol_>
  bool AtomIterator_<Atom_, Mol_>::operator<(const AtomIterator_<Atom_, Mol_> &other){
    return _mol==other._mol && _pos<other._pos;
  }
  template <class Atom_, class Mol_>
  bool AtomIterator_<Atom_, Mol_>::operator<=(const AtomIterator_<Atom_, Mol_> &other){
    return _mol==other._mol && _pos<=other._pos;
  }
  template <class Atom_, class Mol_>
  bool AtomIterator_<Atom_, Mol_>::operator>(const AtomIterator_<Atom_, Mol_> &other){
    return _mol==other._mol && _pos>other._pos;
  }
  template <class Atom_, class Mol_>
  bool AtomIterator_<Atom_, Mol_>::operator>=(const AtomIterator_<Atom_, Mol_> &other){
    return _mol==other._mol && _pos>=other._pos;
  }


  // pre-increment
  template <class Atom_, class Mol_>
  AtomIterator_<Atom_, Mol_> &AtomIterator_<Atom_, Mol_>::operator++() {
    RANGE_CHECK(0,_pos,_max-1);
    _pos++;
    return *this;
  }
  template <class Atom_, class Mol_>
  AtomIterator_<Atom_, Mol_> AtomIterator_<Atom_, Mol_>::operator++(int) {
    RANGE_CHECK(0,_pos,_max-1);
    AtomIterator_<Atom_,Mol_> res(*this);
    _pos++;
    return res;
  }
  // pre-decrement
  template <class Atom_, class Mol_>
  AtomIterator_<Atom_, Mol_> &AtomIterator_<Atom_, Mol_>::operator--() {
    _pos--;
    RANGE_CHECK(0,_pos,_max-1);
    return *this;
  }
  template <class Atom_, class Mol_>
  AtomIterator_<Atom_, Mol_> AtomIterator_<Atom_, Mol_>::operator--(int) {
    RANGE_CHECK(0,_pos,_max-1);
    AtomIterator_<Atom_,Mol_> res(*this);
    if(_pos-1 < 0) _pos = _max;
    else _pos--;
    return res;
  }


  template class AtomIterator_<Atom,ROMol>;
  template class AtomIterator_<const Atom,const ROMol>;

}; // end o' namespace


