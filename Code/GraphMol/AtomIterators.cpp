//
//  Copyright (C) 2002-2017 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "AtomIterators.h"
#include "RDKitBase.h"
#include "RDKitQueries.h"

namespace RDKit {
template <class Atom_, class Mol_>
AtomIterator_<Atom_, Mol_>::AtomIterator_(Mol_ *mol) {
  _mol = mol;
  _pos = 0;
  _max = mol->getNumAtoms();
};
template <class Atom_, class Mol_>
AtomIterator_<Atom_, Mol_>::AtomIterator_(Mol_ *mol, int pos) {
  _mol = mol;
  _pos = pos;
  _max = mol->getNumAtoms();
};
template <class Atom_, class Mol_>
AtomIterator_<Atom_, Mol_>::AtomIterator_(
    const AtomIterator_<Atom_, Mol_> &other) {
  _mol = other._mol;
  _pos = other._pos;
  _max = other._max;
}
template <class Atom_, class Mol_>
AtomIterator_<Atom_, Mol_> &AtomIterator_<Atom_, Mol_>::operator=(
    const AtomIterator_<Atom_, Mol_> &other) {
  _mol = other._mol;
  _pos = other._pos;
  _max = other._max;
  return *this;
}

template <class Atom_, class Mol_>
AtomIterator_<Atom_, Mol_> &AtomIterator_<Atom_, Mol_>::operator+=(int val) {
  _pos += val;
  if (_pos < 0 || _pos > _max) {
    _pos = _max;
  }
  return *this;
}
template <class Atom_, class Mol_>
AtomIterator_<Atom_, Mol_> &AtomIterator_<Atom_, Mol_>::operator-=(int val) {
  _pos -= val;
  if (_pos < 0 || _pos > _max) {
    _pos = _max;
  }
  return *this;
}
template <class Atom_, class Mol_>
AtomIterator_<Atom_, Mol_> AtomIterator_<Atom_, Mol_>::operator+(
    int val) const {
  AtomIterator_<Atom_, Mol_> res(*this);
  res += val;
  // += takes care of the pre/post conditions for us, so we're safe to return
  return res;
}
template <class Atom_, class Mol_>
AtomIterator_<Atom_, Mol_> AtomIterator_<Atom_, Mol_>::operator-(
    int val) const {
  AtomIterator_<Atom_, Mol_> res(*this);
  // -= takes care of the pre/post conditions for us, so we're safe to return
  res -= val;
  return res;
}

// iterator subtraction
template <class Atom_, class Mol_>
int AtomIterator_<Atom_, Mol_>::operator-(
    AtomIterator_<Atom_, Mol_> &other) const {
  PRECONDITION(_mol == other._mol, "bad operator- call");
  return _pos - other._pos;
}

// dereference
template <class Atom_, class Mol_>
Atom_ *AtomIterator_<Atom_, Mol_>::operator*() const {
  PRECONDITION(_mol != nullptr, "no molecule");
  RANGE_CHECK(0, _pos, _max - 1);
  return (*_mol)[_pos];
}
// random access
template <class Atom_, class Mol_>
Atom_ *AtomIterator_<Atom_, Mol_>::operator[](const int which) const {
  PRECONDITION(_mol != nullptr, "no molecule");
  RANGE_CHECK(0, which, _max - 1);
  return (*_mol)[which];
}

template <class Atom_, class Mol_>
bool AtomIterator_<Atom_, Mol_>::operator==(
    const AtomIterator_<Atom_, Mol_> &other) const {
  return _mol == other._mol && _pos == other._pos;
}
template <class Atom_, class Mol_>
bool AtomIterator_<Atom_, Mol_>::operator!=(
    const AtomIterator_<Atom_, Mol_> &other) const {
  return _mol != other._mol || _pos != other._pos;
}
template <class Atom_, class Mol_>
bool AtomIterator_<Atom_, Mol_>::operator<(
    const AtomIterator_<Atom_, Mol_> &other) const {
  return _mol == other._mol && _pos < other._pos;
}
template <class Atom_, class Mol_>
bool AtomIterator_<Atom_, Mol_>::operator<=(
    const AtomIterator_<Atom_, Mol_> &other) const {
  return _mol == other._mol && _pos <= other._pos;
}
template <class Atom_, class Mol_>
bool AtomIterator_<Atom_, Mol_>::operator>(
    const AtomIterator_<Atom_, Mol_> &other) const {
  return _mol == other._mol && _pos > other._pos;
}
template <class Atom_, class Mol_>
bool AtomIterator_<Atom_, Mol_>::operator>=(
    const AtomIterator_<Atom_, Mol_> &other) const {
  return _mol == other._mol && _pos >= other._pos;
}

// pre-increment
template <class Atom_, class Mol_>
AtomIterator_<Atom_, Mol_> &AtomIterator_<Atom_, Mol_>::operator++() {
  _pos++;
  return *this;
}
template <class Atom_, class Mol_>
AtomIterator_<Atom_, Mol_> AtomIterator_<Atom_, Mol_>::operator++(int) {
  AtomIterator_<Atom_, Mol_> res(*this);
  _pos++;
  return res;
}
// pre-decrement
template <class Atom_, class Mol_>
AtomIterator_<Atom_, Mol_> &AtomIterator_<Atom_, Mol_>::operator--() {
  _pos--;
  return *this;
}
template <class Atom_, class Mol_>
AtomIterator_<Atom_, Mol_> AtomIterator_<Atom_, Mol_>::operator--(int) {
  AtomIterator_<Atom_, Mol_> res(*this);
  if (_pos - 1 < 0) {
    _pos = _max;
  } else {
    _pos--;
  }
  return res;
}

//-----------------------------------------
//
//  HeteroatomIterator
//
//-----------------------------------------
template <class Atom_, class Mol_>
HeteroatomIterator_<Atom_, Mol_>::HeteroatomIterator_(Mol_ *mol) {
  _mol = mol;
  _qA = new QueryAtom(6);
  _qA->getQuery()->setNegation(true);
  _end = mol->getNumAtoms();
  _pos = _findNext(0);
};
template <class Atom_, class Mol_>
HeteroatomIterator_<Atom_, Mol_>::HeteroatomIterator_(Mol_ *mol, int pos) {
  _mol = mol;
  _qA = new QueryAtom(6);
  _end = mol->getNumAtoms();
  _pos = pos;
};

template <class Atom_, class Mol_>
HeteroatomIterator_<Atom_, Mol_>::~HeteroatomIterator_() {
  delete _qA;
  _qA = nullptr;
}

template <class Atom_, class Mol_>
HeteroatomIterator_<Atom_, Mol_>::HeteroatomIterator_(const ThisType &other) {
  _mol = other._mol;
  _end = other._end;
  _pos = other._pos;
  _qA = static_cast<QueryAtom *>(other._qA->copy());
}

template <class Atom_, class Mol_>
HeteroatomIterator_<Atom_, Mol_> &HeteroatomIterator_<Atom_, Mol_>::operator=(
    const ThisType &other) {
  _mol = other._mol;
  _end = other._end;
  _pos = other._pos;
  _qA = static_cast<QueryAtom *>(other._qA->copy());
  return *this;
}

template <class Atom_, class Mol_>
bool HeteroatomIterator_<Atom_, Mol_>::operator==(const ThisType &other) const {
  return _mol == other._mol && _pos == other._pos;
}
template <class Atom_, class Mol_>
bool HeteroatomIterator_<Atom_, Mol_>::operator!=(const ThisType &other) const {
  return _mol != other._mol || _pos != other._pos;
}

template <class Atom_, class Mol_>
Atom_ *HeteroatomIterator_<Atom_, Mol_>::operator*() const {
  PRECONDITION(_mol != nullptr, "no molecule");
  return (*_mol)[_pos];
}
// pre-increment
template <class Atom_, class Mol_>
HeteroatomIterator_<Atom_, Mol_> &
HeteroatomIterator_<Atom_, Mol_>::operator++() {
  _pos = _findNext(_pos + 1);
  return *this;
}
template <class Atom_, class Mol_>
HeteroatomIterator_<Atom_, Mol_> HeteroatomIterator_<Atom_, Mol_>::operator++(
    int) {
  HeteroatomIterator_<Atom_, Mol_> res(*this);
  _pos = _findNext(_pos + 1);
  return res;
}
// pre-decrement
template <class Atom_, class Mol_>
HeteroatomIterator_<Atom_, Mol_> &
HeteroatomIterator_<Atom_, Mol_>::operator--() {
  _pos = _findPrev(_pos - 1);
  return *this;
}
template <class Atom_, class Mol_>
HeteroatomIterator_<Atom_, Mol_> HeteroatomIterator_<Atom_, Mol_>::operator--(
    int) {
  HeteroatomIterator_<Atom_, Mol_> res(*this);
  _pos = _findPrev(_pos - 1);
  return res;
}
template <class Atom_, class Mol_>
int HeteroatomIterator_<Atom_, Mol_>::_findNext(int from) {
  while (from < _end) {
    if (_qA->Match((*_mol)[from])) {
      break;
    } else {
      from++;
    }
  }
  return from;
}

template <class Atom_, class Mol_>
int HeteroatomIterator_<Atom_, Mol_>::_findPrev(int from) {
  while (from > 0) {
    if (_qA->Match((*_mol)[from])) {
      break;
    } else {
      from--;
    }
  }
  if (from < 0) {
    from = _end;
  }
  return from;
}

//-----------------------------------------
//
//  AromaticAtomIterator
//
//-----------------------------------------
template <class Atom_, class Mol_>
AromaticAtomIterator_<Atom_, Mol_>::AromaticAtomIterator_(Mol_ *mol) {
  _mol = mol;
  _end = mol->getNumAtoms();
  _pos = _findNext(0);
};
template <class Atom_, class Mol_>
AromaticAtomIterator_<Atom_, Mol_>::AromaticAtomIterator_(Mol_ *mol, int pos) {
  _mol = mol;
  _end = mol->getNumAtoms();
  _pos = pos;
};

template <class Atom_, class Mol_>
AromaticAtomIterator_<Atom_, Mol_>::~AromaticAtomIterator_() {}

template <class Atom_, class Mol_>
AromaticAtomIterator_<Atom_, Mol_>::AromaticAtomIterator_(
    const ThisType &other) {
  _mol = other._mol;
  _end = other._end;
  _pos = other._pos;
}

template <class Atom_, class Mol_>
AromaticAtomIterator_<Atom_, Mol_> &
AromaticAtomIterator_<Atom_, Mol_>::operator=(const ThisType &other) {
  _mol = other._mol;
  _end = other._end;
  _pos = other._pos;
  return *this;
}

template <class Atom_, class Mol_>
bool AromaticAtomIterator_<Atom_, Mol_>::operator==(
    const ThisType &other) const {
  return _mol == other._mol && _pos == other._pos;
}
template <class Atom_, class Mol_>
bool AromaticAtomIterator_<Atom_, Mol_>::operator!=(
    const ThisType &other) const {
  return _mol != other._mol || _pos != other._pos;
}

template <class Atom_, class Mol_>
Atom_ *AromaticAtomIterator_<Atom_, Mol_>::operator*() const {
  PRECONDITION(_mol != nullptr, "no molecule");
  return (*_mol)[_pos];
}
// pre-increment
template <class Atom_, class Mol_>
AromaticAtomIterator_<Atom_, Mol_> &
AromaticAtomIterator_<Atom_, Mol_>::operator++() {
  _pos = _findNext(_pos + 1);
  return *this;
}
template <class Atom_, class Mol_>
AromaticAtomIterator_<Atom_, Mol_>
AromaticAtomIterator_<Atom_, Mol_>::operator++(int) {
  AromaticAtomIterator_<Atom_, Mol_> res(*this);
  _pos = _findNext(_pos + 1);
  return res;
}
// pre-decrement
template <class Atom_, class Mol_>
AromaticAtomIterator_<Atom_, Mol_> &
AromaticAtomIterator_<Atom_, Mol_>::operator--() {
  _pos = _findPrev(_pos - 1);
  return *this;
}
template <class Atom_, class Mol_>
AromaticAtomIterator_<Atom_, Mol_>
AromaticAtomIterator_<Atom_, Mol_>::operator--(int) {
  AromaticAtomIterator_<Atom_, Mol_> res(*this);
  _pos = _findPrev(_pos - 1);
  return res;
}
template <class Atom_, class Mol_>
int AromaticAtomIterator_<Atom_, Mol_>::_findNext(int from) {
  while (from < _end) {
    if ((*_mol)[from]->getIsAromatic()) {
      break;
    } else {
      from++;
    }
  }
  return from;
}

template <class Atom_, class Mol_>
int AromaticAtomIterator_<Atom_, Mol_>::_findPrev(int from) {
  while (from > 0) {
    if ((*_mol)[from]->getIsAromatic()) {
      break;
    } else {
      from--;
    }
  }
  if (from < 0) {
    from = _end;
  }
  return from;
}

//-----------------------------------------
//
//  QueryAtomIterator
//
//-----------------------------------------
template <class Atom_, class Mol_>
QueryAtomIterator_<Atom_, Mol_>::QueryAtomIterator_(Mol_ *mol,
                                                    QueryAtom const *what) {
  PRECONDITION(what, "bad query atom");
  _mol = mol;
  _qA = static_cast<QueryAtom *>(what->copy());
  _end = mol->getNumAtoms();
  _pos = _findNext(0);
};
template <class Atom_, class Mol_>
QueryAtomIterator_<Atom_, Mol_>::QueryAtomIterator_(Mol_ *mol, int pos) {
  _mol = mol;
  _qA = nullptr;
  _end = mol->getNumAtoms();
  _pos = pos;
};
template <class Atom_, class Mol_>
QueryAtomIterator_<Atom_, Mol_>::~QueryAtomIterator_() {
  delete _qA;
  _qA = nullptr;
}
template <class Atom_, class Mol_>
QueryAtomIterator_<Atom_, Mol_>::QueryAtomIterator_(
    const QueryAtomIterator_<Atom_, Mol_> &other) {
  _mol = other._mol;
  _pos = other._pos;
  _end = other._end;
  if (other._qA) {
    _qA = static_cast<QueryAtom *>(other._qA->copy());
  } else {
    _qA = nullptr;
  }
}

template <class Atom_, class Mol_>
QueryAtomIterator_<Atom_, Mol_> &QueryAtomIterator_<Atom_, Mol_>::operator=(
    const QueryAtomIterator_<Atom_, Mol_> &other) {
  if (this != &other) {
    _mol = other._mol;
    _pos = other._pos;
    _end = other._end;
    delete _qA;
    if (other._qA) {
      _qA = static_cast<QueryAtom *>(other._qA->copy());
    } else {
      _qA = nullptr;
    }
  }
  return *this;
}
template <class Atom_, class Mol_>
bool QueryAtomIterator_<Atom_, Mol_>::operator==(
    const QueryAtomIterator_<Atom_, Mol_> &other) const {
  return _mol == other._mol && _pos == other._pos;
}
template <class Atom_, class Mol_>
bool QueryAtomIterator_<Atom_, Mol_>::operator!=(
    const QueryAtomIterator_<Atom_, Mol_> &other) const {
  return _mol != other._mol || _pos != other._pos;
}

template <class Atom_, class Mol_>
Atom_ *QueryAtomIterator_<Atom_, Mol_>::operator*() const {
  PRECONDITION(_mol != nullptr, "no molecule");
  return (*_mol)[_pos];
}
// pre-increment
template <class Atom_, class Mol_>
QueryAtomIterator_<Atom_, Mol_> &QueryAtomIterator_<Atom_, Mol_>::operator++() {
  _pos = _findNext(_pos + 1);
  return *this;
}
template <class Atom_, class Mol_>
QueryAtomIterator_<Atom_, Mol_> QueryAtomIterator_<Atom_, Mol_>::operator++(
    int) {
  QueryAtomIterator_ res(*this);
  _pos = _findNext(_pos + 1);
  return res;
}
// pre-decrement
template <class Atom_, class Mol_>
QueryAtomIterator_<Atom_, Mol_> &QueryAtomIterator_<Atom_, Mol_>::operator--() {
  _pos = _findPrev(_pos - 1);
  return *this;
}
template <class Atom_, class Mol_>
QueryAtomIterator_<Atom_, Mol_> QueryAtomIterator_<Atom_, Mol_>::operator--(
    int) {
  QueryAtomIterator_<Atom_, Mol_> res(*this);
  _pos = _findPrev(_pos - 1);
  return res;
}
template <class Atom_, class Mol_>
int QueryAtomIterator_<Atom_, Mol_>::_findNext(int from) {
  PRECONDITION(_mol != nullptr, "no molecule");
  PRECONDITION(_qA != nullptr, "no query set");
  while (from < _end) {
    if (_qA->Match((*_mol)[from])) {
      break;
    } else {
      from++;
    }
  }
  return from;
}

template <class Atom_, class Mol_>
int QueryAtomIterator_<Atom_, Mol_>::_findPrev(int from) {
  PRECONDITION(_mol != nullptr, "no molecule");
  PRECONDITION(_qA != nullptr, "no query set");
  while (from > 0) {
    if (_qA->Match((*_mol)[from])) {
      break;
    } else {
      from--;
    }
  }
  if (from < 0) {
    from = _end;
  }
  return from;
}

//-----------------------------------------
//
//  MatchingAtomIterator
//
//-----------------------------------------
template <class Atom_, class Mol_>
MatchingAtomIterator_<Atom_, Mol_>::MatchingAtomIterator_(Mol_ *mol,
                                                          bool (*fn)(Atom_ *)) {
  PRECONDITION(fn, "bad query function");
  _mol = mol;
  _qF = fn;
  _end = mol->getNumAtoms();
  _pos = _findNext(0);
};
template <class Atom_, class Mol_>
MatchingAtomIterator_<Atom_, Mol_>::MatchingAtomIterator_(Mol_ *mol, int pos) {
  _mol = mol;
  _qF = nullptr;
  _end = mol->getNumAtoms();
  _pos = pos;
};

template <class Atom_, class Mol_>
MatchingAtomIterator_<Atom_, Mol_>::~MatchingAtomIterator_() {}

template <class Atom_, class Mol_>
MatchingAtomIterator_<Atom_, Mol_>::MatchingAtomIterator_(
    const MatchingAtomIterator_<Atom_, Mol_> &other) {
  _mol = other._mol;
  _pos = other._pos;
  _end = other._end;
  _qF = other._qF;
}

template <class Atom_, class Mol_>
MatchingAtomIterator_<Atom_, Mol_> &
MatchingAtomIterator_<Atom_, Mol_>::operator=(
    const MatchingAtomIterator_<Atom_, Mol_> &other) {
  if (this != &other) {
    _mol = other._mol;
    _pos = other._pos;
    _end = other._end;
    _qF = other._qF;
  }
  return *this;
}
template <class Atom_, class Mol_>
bool MatchingAtomIterator_<Atom_, Mol_>::operator==(
    const MatchingAtomIterator_<Atom_, Mol_> &other) const {
  return _mol == other._mol && _pos == other._pos;
}
template <class Atom_, class Mol_>
bool MatchingAtomIterator_<Atom_, Mol_>::operator!=(
    const MatchingAtomIterator_<Atom_, Mol_> &other) const {
  return _mol != other._mol || _pos != other._pos;
}

template <class Atom_, class Mol_>
Atom_ *MatchingAtomIterator_<Atom_, Mol_>::operator*() const {
  PRECONDITION(_mol != nullptr, "no molecule");
  return (*_mol)[_pos];
}
// pre-increment
template <class Atom_, class Mol_>
MatchingAtomIterator_<Atom_, Mol_> &
MatchingAtomIterator_<Atom_, Mol_>::operator++() {
  _pos = _findNext(_pos + 1);
  return *this;
}
template <class Atom_, class Mol_>
MatchingAtomIterator_<Atom_, Mol_>
MatchingAtomIterator_<Atom_, Mol_>::operator++(int) {
  MatchingAtomIterator_ res(*this);
  _pos = _findNext(_pos + 1);
  return res;
}
// pre-decrement
template <class Atom_, class Mol_>
MatchingAtomIterator_<Atom_, Mol_> &
MatchingAtomIterator_<Atom_, Mol_>::operator--() {
  _pos = _findPrev(_pos - 1);
  return *this;
}
template <class Atom_, class Mol_>
MatchingAtomIterator_<Atom_, Mol_>
MatchingAtomIterator_<Atom_, Mol_>::operator--(int) {
  MatchingAtomIterator_<Atom_, Mol_> res(*this);
  _pos = _findPrev(_pos - 1);
  return res;
}
template <class Atom_, class Mol_>
int MatchingAtomIterator_<Atom_, Mol_>::_findNext(int from) {
  PRECONDITION(_mol != nullptr, "no molecule");
  PRECONDITION(_qF != nullptr, "no query set");
  while (from < _end) {
    if (_qF((*_mol)[from])) {
      break;
    } else {
      ++from;
    }
  }
  return from;
}

template <class Atom_, class Mol_>
int MatchingAtomIterator_<Atom_, Mol_>::_findPrev(int from) {
  PRECONDITION(_mol != nullptr, "no molecule");
  PRECONDITION(_qF != nullptr, "no query set");
  while (from > 0) {
    if (_qF((*_mol)[from])) {
      break;
    } else {
      --from;
    }
  }
  if (from < 0) {
    from = _end;
  }
  return from;
}

template class AtomIterator_<Atom, ROMol>;
template class AtomIterator_<const Atom, const ROMol>;

template class AromaticAtomIterator_<Atom, ROMol>;
template class AromaticAtomIterator_<const Atom, const ROMol>;
template class HeteroatomIterator_<Atom, ROMol>;
template class HeteroatomIterator_<const Atom, const ROMol>;
template class QueryAtomIterator_<Atom, ROMol>;
template class QueryAtomIterator_<const Atom, const ROMol>;
template class MatchingAtomIterator_<Atom, ROMol>;
template class MatchingAtomIterator_<const Atom, const ROMol>;

};  // namespace RDKit
