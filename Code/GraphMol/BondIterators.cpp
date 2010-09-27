// $Id$
//
//  Copyright (C) 2002-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "BondIterators.h"

namespace RDKit{

  BondIterator_::BondIterator_(ROMol *mol) {
    _mol = mol;
    boost::tie(_beg,_end)=mol->getEdges();
    _pos = _beg;
  };
  BondIterator_::BondIterator_(ROMol *mol,ROMol::EDGE_ITER pos) {
    _mol = mol;
    boost::tie(_beg,_end)=mol->getEdges();

    _pos = pos;
  };
  BondIterator_::BondIterator_(const BondIterator_ &other){
    _mol = other._mol;
    _pos = other._pos;
    _beg = other._beg;
    _end = other._end;
  }

  BondIterator_ &BondIterator_::operator=(const BondIterator_ &other){
    _mol = other._mol;
    _pos = other._pos;
    _beg = other._beg;
    _end = other._end;
    return *this;
  }

  bool BondIterator_::operator==(const BondIterator_ &other){
    return _mol==other._mol && _pos==other._pos;
  }
  bool BondIterator_::operator!=(const BondIterator_ &other){
    return _mol!=other._mol || _pos!=other._pos;
  }

  Bond *BondIterator_::operator*() {
    return (*_mol)[*_pos].get();
  }
  // pre-increment
  BondIterator_ &BondIterator_::operator++() {
    PRECONDITION(_pos!=_end,"bad initial position")
    _pos++;
    return *this;
  }
  BondIterator_ BondIterator_::operator++(int) {
    PRECONDITION(_pos!=_end,"bad initial position")
    BondIterator_ res(*this);
    _pos++;
    return res;
  }
  // pre-decrement
  BondIterator_ &BondIterator_::operator--() {
    if(_pos == _beg) _pos = _end;
    else _pos--;
    return *this;
  }
  BondIterator_ BondIterator_::operator--(int) {
    BondIterator_ res(*this);
    if(_pos == _beg) _pos = _end;
    else _pos--;
    return res;
  }


  // CONST
  ConstBondIterator_::ConstBondIterator_(ROMol const *mol) {
    _mol = mol;
    boost::tie(_beg,_end)=mol->getEdges();
    _pos = _beg;
  };
  ConstBondIterator_::ConstBondIterator_(ROMol const *mol,ROMol::EDGE_ITER pos) {
    _mol = mol;
    boost::tie(_beg,_end)=mol->getEdges();

    _pos = pos;
  };
  ConstBondIterator_::ConstBondIterator_(const ConstBondIterator_ &other){
    _mol = other._mol;
    _pos = other._pos;
    _beg = other._beg;
    _end = other._end;
  }
  ConstBondIterator_ &ConstBondIterator_::operator=(const ConstBondIterator_ &other){
    _mol = other._mol;
    _pos = other._pos;
    _beg = other._beg;
    _end = other._end;
    return *this;
  }
  bool ConstBondIterator_::operator==(const ConstBondIterator_ &other){
    return _mol==other._mol && _pos==other._pos;
  }
  bool ConstBondIterator_::operator!=(const ConstBondIterator_ &other){
    return _mol!=other._mol || _pos!=other._pos;
  }

  Bond const *ConstBondIterator_::operator*() {
    return (*_mol)[*_pos].get();
  }
  // pre-increment
  ConstBondIterator_ &ConstBondIterator_::operator++() {
    PRECONDITION(_pos!=_end,"bad initial position")
    _pos++;
    return *this;
  }
  ConstBondIterator_ ConstBondIterator_::operator++(int) {
    ConstBondIterator_ res(*this);
    PRECONDITION(_pos!=_end,"bad initial position")
    _pos++;
    return res;
  }
  // pre-decrement
  ConstBondIterator_ &ConstBondIterator_::operator--() {
    if(_pos == _beg) _pos = _end;
    else _pos--;
    return *this;
  }
  ConstBondIterator_ ConstBondIterator_::operator--(int) {
    ConstBondIterator_ res(*this);
    if(_pos == _beg) _pos = _end;
    else _pos--;
    return res;
  }


}
