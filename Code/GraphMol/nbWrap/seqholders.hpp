
//
//  Copyright (C) 2026 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryAtom.h>

namespace RDKit {

using AtomsIterator =
    CXXAtomIterator<MolGraph, Atom *, MolGraph::vertex_iterator, true>;
using AtomNeighborsIterator =
    CXXAtomIterator<MolGraph, Atom *, MolGraph::adjacency_iterator, true, true>;
using BondsIterator =
    CXXBondIterator<MolGraph, Bond *, MolGraph::edge_iterator, true>;
using AtomBondsIterator =
    CXXBondIterator<MolGraph, Bond *, MolGraph::out_edge_iterator, true>;

template <typename IterT = AtomsIterator>
struct AtomSeqHolder {
  IterT iter;
  ROMol *mol;
  AtomSeqHolder(ROMol &mol) : iter(mol.checkedAtoms()), mol(&mol) {};
  AtomSeqHolder(ROMol &mol, const IterT &iter) : iter(iter), mol(&mol) {};
  size_t size() const { return iter.size(); }
  typename IterT::CXXAtomIter begin() { return iter.begin(); }
  typename IterT::CXXAtomIter end() { return iter.end(); }
  Atom *operator[](int idx) {
    if (idx < 0 || static_cast<size_t>(idx) >= iter.size()) {
      throw IndexErrorException(idx);
    }
    return mol->getAtomWithIdx(idx);
  }
};

template <typename IterT = BondsIterator>
struct BondSeqHolder {
  IterT iter;
  ROMol *mol;
  BondSeqHolder(ROMol &mol) : iter(mol.checkedBonds()), mol(&mol) {};
  BondSeqHolder(ROMol &mol, const IterT &iter) : iter(iter), mol(&mol) {};
  size_t size() const { return iter.size(); }
  typename IterT::CXXBondIter begin() { return iter.begin(); }
  typename IterT::CXXBondIter end() { return iter.end(); }
  Bond *operator[](int idx) {
    if (idx < 0 || static_cast<size_t>(idx) >= iter.size()) {
      throw IndexErrorException(idx);
    }
    return mol->getBondWithIdx(idx);
  }
};

class AtomCountFunctor {
 private:
  const ROMol &_mol;

 public:
  AtomCountFunctor(const ROMol &mol) : _mol(mol) {}
  unsigned int operator()() const { return _mol.getNumAtoms(); }
};

class ConformerCountFunctor {
 private:
  const ROMol &_mol;

 public:
  ConformerCountFunctor(const ROMol &mol) : _mol(mol) {}
  unsigned int operator()() const { return _mol.getNumConformers(); }
};

// Note: T1 should be some iterator type,
//       T2 is the value when we dereference T
template <class T1, class T2, class T3, class Mapper = std::function<T2(T1)>>
class ReadOnlySeq {
 private:
  T1 _start, _end, _pos;
  int _size;
  T3 _lenFunc;
  size_t _origLen;
  Mapper _mapper;
  const ROMol &_mol;

 public:
  ~ReadOnlySeq() {
    std::cerr << "\n**** delete ReadOnlySeq ****" << this << std::endl;
  }
  // ~ReadOnlySeq() = default;
  ReadOnlySeq(
      const ROMol &mol, T1 start, T1 end, T3 lenFunc,
      Mapper mapper = [](T1 it) { return *it; })
      : _start(start),
        _end(end),
        _pos(start),
        _size(-1),
        _lenFunc(lenFunc),
        _origLen(lenFunc()),
        _mapper(mapper),
        _mol(mol) {
    std::cerr << "\n**** new ReadOnlySeq ****" << this << std::endl;
  }
  ReadOnlySeq(const ReadOnlySeq<T1, T2, T3, Mapper> &other)
      : _start(other._start),
        _end(other._end),
        _pos(other._pos),
        _size(other._size),
        _lenFunc(other._lenFunc),
        _origLen(other._origLen),
        _mapper(other._mapper),
        _mol(other._mol) {
    std::cerr << "\n**** new ReadOnlySeq2 ****" << this << std::endl;
  }
  T1 begin() { return _start; }
  T1 end() { return _end; }
  void reset() {
    // std::cerr << "**** reset ****" << this<< std::endl;
    _pos = _start;
  }
  ReadOnlySeq<T1, T2, T3> *__iter__() {
    std::cerr << "**** ITER ****" << this << std::endl;
    reset();
    std::cerr << "  finish ****" << this << std::endl;
    return this;
  }
  T2 next() {
    if (_pos == _end) {
      throw nb::stop_iteration();
    }
    std::cerr << "\tnext: " << std::endl;
    if (_lenFunc() != _origLen) {
      throw std::runtime_error("Sequence modified during iteration");
    }
    T2 res = _mapper(_pos);
    ++_pos;
    return res;
  }
  T2 get_item(int which) {
    std::cerr << "get_item: " << which << std::endl;
    if (which >= len()) {
      throw nb::index_error("End of sequence hit");
    }
    if (_lenFunc() != _origLen) {
      throw std::runtime_error("Sequence modified during iteration");
    }
    T1 it = _start;
    for (int i = 0; i < which; i++) {
      ++it;
    }
    return _mapper(it);
  }
  int len() {
    std::cerr << "len " << std::endl;
    if (_size < 0) {
      _size = 0;
      for (T1 tmp = _start; tmp != _end; tmp++) {
        std::cerr << "\tincr: " << std::endl;
        _size++;
      }
    }
    std::cerr << "\tret" << std::endl;
    return _size;
  }
};

class ConformerIterSeq {
 private:
  const ROMol &_mol;
  std::vector<Conformer *> _confs;

 public:
  ConformerIterSeq(const ROMol &mol) : _mol(mol) {
    _confs.reserve(mol.getNumConformers());
    std::transform(mol.beginConformers(), mol.endConformers(),
                   std::back_inserter(_confs),
                   [](const CONFORMER_SPTR &conf) { return conf.get(); });
  };
  std::vector<Conformer *>::iterator begin() { return _confs.begin(); }
  std::vector<Conformer *>::iterator end() { return _confs.end(); }
  Conformer *operator[](int idx) {
    if(_mol.getNumConformers() != _confs.size()) {
      throw std::runtime_error("Sequence modified during iteration");
    }
    if (idx < 0 || static_cast<size_t>(idx) >= _confs.size()) {
      throw IndexErrorException(idx);
    }
    return _confs[idx];
  }
  size_t size() { return _confs.size(); }
};

class QueryAtomIterSeq {
 private:
  const ROMol &_mol;
  const QueryAtom *_qa;
  int _size = -1;
  bool _ownsQA = false;

 public:
  QueryAtomIterSeq(const ROMol &mol, const QueryAtom *qa, bool ownsQA = false)
      : _mol(mol), _qa(qa), _ownsQA(ownsQA) {};
  QueryAtomIterSeq(const QueryAtomIterSeq &other)
      : _mol(other._mol),
        _qa(static_cast<const QueryAtom *>(other._qa->copy())),
        _ownsQA(true) {};
  QueryAtomIterSeq &operator=(const QueryAtomIterSeq &other) = delete;
  ~QueryAtomIterSeq() {
    if (_ownsQA) {
      delete _qa;
    }
  }
  ROMol::ConstQueryAtomIterator begin() { return _mol.beginQueryAtoms(_qa); }
  ROMol::ConstQueryAtomIterator end() { return _mol.endQueryAtoms(); }
  size_t size() {
    if (_size < 0) {
      _size = 0;
      for (auto it = begin(); it != end(); ++it) {
        _size++;
      }
    }
    return static_cast<size_t>(_size);
  }
  const Atom *operator[](int idx) {
    if (idx < 0) {
      throw IndexErrorException(idx);
    }
    auto it = begin();
    for (int i = 0; i < idx; i++) {
      ++it;
      if (it == end()) {
        throw IndexErrorException(idx);
      }
    }
    return *it;
  }
};
}  // namespace RDKit