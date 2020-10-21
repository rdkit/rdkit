//
//  Copyright (C) 2003-2006 greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _SEQS_HPP_
#define _SEQS_HPP_

#include <GraphMol/RDKitBase.h>
#include <iostream>
#include <RDBoost/python.h>
namespace python = boost::python;

namespace RDKit {

class AtomCountFunctor {
 private:
  const ROMol &_mol;

 public:
  AtomCountFunctor(const ROMol &mol) : _mol(mol){};
  unsigned int operator()() const { return _mol.getNumAtoms(); };
};
class BondCountFunctor {
 private:
  const ROMol &_mol;

 public:
  BondCountFunctor(const ROMol &mol) : _mol(mol){};
  unsigned int operator()() const { return _mol.getNumBonds(); };
};
class ConformerCountFunctor {
 private:
  const ROMol &_mol;

 public:
  ConformerCountFunctor(const ROMol &mol) : _mol(mol){};
  unsigned int operator()() const { return _mol.getNumConformers(); };
};

// Note: T1 should be some iterator type,
//       T2 is the value when we dereference T
template <class T1, class T2, class T3>
class ReadOnlySeq : public RDChild {
 private:
  T1 _start, _end, _pos;
  int _size;
  T3 _lenFunc;
  size_t _origLen;
  ROMol *_mol;

 public:
  ~ReadOnlySeq() {}
  ReadOnlySeq(ROMol *mol, T1 start, T1 end, T3 lenFunc)
      : RDChild(mol),
        _start(start),
        _end(end),
        _pos(start),
        _size(-1),
        _lenFunc(lenFunc),
        _origLen(lenFunc()),
        _mol(mol) {}
  ReadOnlySeq(const ReadOnlySeq<T1, T2, T3> &other)
      : RDChild(other._mol),
        _start(other._start),
        _end(other._end),
        _pos(other._pos),
        _size(other._size),
        _lenFunc(other._lenFunc),
        _origLen(other._origLen),
        _mol(other._mol) {}
  void throwIfParentMolDeleted() const {
    if (!hasParent()) {
      PyErr_SetString(PyExc_RuntimeError, "Underlying ROMol object deleted");
      throw python::error_already_set();
    }
  }
  void throwIfSequenceModified() const {
    if (_lenFunc() != _origLen) {
      PyErr_SetString(PyExc_RuntimeError, "Sequence modified during iteration");
      throw python::error_already_set();
    }
  }
  void throwEndOfSequenceHit() const {
    PyErr_SetString(PyExc_StopIteration, "End of sequence hit");
    throw python::error_already_set();
  }
  void throwIndexOutOfBounds() const {
    PyErr_SetString(PyExc_IndexError, "IndexOutOfBounds");
    throw python::error_already_set();
  }

  void reset() {
    // std::cerr << "**** reset ****" << this<< std::endl;
    _pos = _start;
  }
  ReadOnlySeq<T1, T2, T3> *__iter__() {
    // std::cerr << "**** ITER ****" << this << std::endl;
    reset();
    // std::cerr << "  finish ****" << this << std::endl;
    return this;
  }
  T2 next() {
    // std::cerr << "\tnext: " << _pos._pos << " " << _end._pos << std::endl;
    throwIfParentMolDeleted();
    if (_pos == _end) {
      throwEndOfSequenceHit();
    }
    throwIfSequenceModified();
    T2 res = *_pos;
    ++_pos;
    return res;
  }
  T2 get_item(int which) {
    // std::cerr << "get_item: " <<which<< std::endl;
    throwIfParentMolDeleted();
    if (which >= len()) {
      throwIndexOutOfBounds();
    }
    throwIfSequenceModified();
    T1 it = _start;
    for (int i = 0; i < which; i++) {
      ++it;
    }
    return *it;
  }
  int len() {
    // std::cerr << "len " << std::endl;
    if (_size < 0) {
      throwIfParentMolDeleted();
      _size = 0;
      for (T1 tmp = _start; tmp != _end; ++tmp, ++_size)
        ;
    }
    // std::cerr << "\tret" << std::endl;
    return _size;
  }
};

typedef ReadOnlySeq<ROMol::AtomIterator, Atom *, AtomCountFunctor> AtomIterSeq;
typedef ReadOnlySeq<ROMol::QueryAtomIterator, Atom *, AtomCountFunctor>
    QueryAtomIterSeq;
typedef ReadOnlySeq<ROMol::BondIterator, Bond *, BondCountFunctor> BondIterSeq;
typedef ReadOnlySeq<ROMol::ConformerIterator, CONFORMER_SPTR &, ConformerCountFunctor> ConformerIterSeq;
}
#endif
