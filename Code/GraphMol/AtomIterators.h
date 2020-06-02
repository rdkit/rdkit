//
//  Copyright (C) 2002-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file AtomIterators.h

  \brief various tools for iterating over a molecule's Atoms.

   <b>WARNING:</b> If you go changing the molecule underneath one of
   these iterators you will be sad...
*/
#include <RDGeneral/export.h>
#ifndef __RD_ATOM_ITERATORS_H__
#define __RD_ATOM_ITERATORS_H__

#ifdef _MSC_VER
#pragma warning(disable : 4661)  // no suitable definition provided for explicit
                                 // template instantiation request
#endif

namespace RDKit {
class QueryAtom;

//! A general random access iterator
template <class Atom_, class Mol_>
class RDKIT_GRAPHMOL_EXPORT AtomIterator_ {
 public:
  typedef AtomIterator_<Atom_, Mol_> ThisType;
  AtomIterator_() : _mol(nullptr){};
  AtomIterator_(Mol_ *mol);
  AtomIterator_(Mol_ *mol, int pos);
  AtomIterator_(const ThisType &other);
  AtomIterator_ &operator=(const ThisType &other);
  AtomIterator_ &operator+=(int val);
  AtomIterator_ &operator-=(int val);
  AtomIterator_ operator+(int val) const;
  AtomIterator_ operator-(int val) const;

  // iterator subtraction
  int operator-(ThisType &other) const;

  // dereference
  Atom_ *operator*() const;
  // random access
  Atom_ *operator[](const int which) const;
  bool operator==(const ThisType &other) const;
  bool operator!=(const ThisType &other) const;
  bool operator<(const ThisType &other) const;
  bool operator<=(const ThisType &other) const;
  bool operator>(const ThisType &other) const;
  bool operator>=(const ThisType &other) const;

  // pre-increment
  ThisType &operator++();
  ThisType operator++(int);

  // pre-decrement
  ThisType &operator--();
  ThisType operator--(int);

 private:
  int _pos{-1};
  int _max{-1};
  Mol_ *_mol;
};

//! Iterate over heteroatoms, this is bidirectional
template <class Atom_, class Mol_>
class RDKIT_GRAPHMOL_EXPORT HeteroatomIterator_ {
 public:
  typedef HeteroatomIterator_<Atom_, Mol_> ThisType;
  HeteroatomIterator_() : _mol(nullptr){};
  HeteroatomIterator_(Mol_ *mol);
  HeteroatomIterator_(Mol_ *mol, int pos);
  ~HeteroatomIterator_();
  HeteroatomIterator_(const ThisType &other);
  HeteroatomIterator_ &operator=(const ThisType &other);
  bool operator==(const ThisType &other) const;
  bool operator!=(const ThisType &other) const;

  Atom_ *operator*() const;

  // pre-increment
  ThisType &operator++();
  ThisType operator++(int);

  // pre-decrement
  ThisType &operator--();
  ThisType operator--(int);

 private:
  int _end{-1};
  int _pos{-1};
  Mol_ *_mol;
  // FIX: somehow changing the following to a pointer make the regression test
  // pass
  // QueryAtom _qA;
  QueryAtom *_qA;

  int _findNext(int from);
  int _findPrev(int from);
};

//! Iterate over aromatic atoms, this is bidirectional
template <class Atom_, class Mol_>
class RDKIT_GRAPHMOL_EXPORT AromaticAtomIterator_ {
 public:
  typedef AromaticAtomIterator_<Atom_, Mol_> ThisType;
  AromaticAtomIterator_() : _mol(nullptr){};
  AromaticAtomIterator_(Mol_ *mol);
  AromaticAtomIterator_(Mol_ *mol, int pos);
  ~AromaticAtomIterator_();
  AromaticAtomIterator_(const ThisType &other);
  AromaticAtomIterator_ &operator=(const ThisType &other);
  bool operator==(const ThisType &other) const;
  bool operator!=(const ThisType &other) const;

  Atom_ *operator*() const;

  // pre-increment
  ThisType &operator++();
  ThisType operator++(int);

  // pre-decrement
  ThisType &operator--();
  ThisType operator--(int);

 private:
  int _end{-1};
  int _pos{-1};
  Mol_ *_mol;

  int _findNext(int from);
  int _findPrev(int from);
};

//! Iterate over atoms matching a query. This is bidirectional.
template <class Atom_, class Mol_>
class RDKIT_GRAPHMOL_EXPORT QueryAtomIterator_ {
 public:
  typedef QueryAtomIterator_<Atom_, Mol_> ThisType;
  QueryAtomIterator_() : _mol(nullptr){};
  QueryAtomIterator_(Mol_ *mol, QueryAtom const *what);
  QueryAtomIterator_(Mol_ *mol, int pos);
  ~QueryAtomIterator_();
  QueryAtomIterator_(const ThisType &other);
  QueryAtomIterator_ &operator=(const ThisType &other);
  bool operator==(const ThisType &other) const;
  bool operator!=(const ThisType &other) const;

  Atom_ *operator*() const;

  // pre-increment
  ThisType &operator++();
  ThisType operator++(int);

  // pre-decrement
  ThisType &operator--();
  ThisType operator--(int);

 private:
  int _end{-1};
  int _pos{-1};
  Mol_ *_mol;
  QueryAtom *_qA{nullptr};

  int _findNext(int from);
  int _findPrev(int from);
};

//! Iterate over atoms matching a query function. This is bidirectional.
template <class Atom_, class Mol_>
class RDKIT_GRAPHMOL_EXPORT MatchingAtomIterator_ {
 public:
  typedef MatchingAtomIterator_<Atom_, Mol_> ThisType;
  MatchingAtomIterator_() : _mol(nullptr), _qF(nullptr){};
  MatchingAtomIterator_(Mol_ *mol, bool (*fn)(Atom_ *));
  MatchingAtomIterator_(Mol_ *mol, int pos);
  ~MatchingAtomIterator_();
  MatchingAtomIterator_(const ThisType &other);
  MatchingAtomIterator_ &operator=(const ThisType &other);
  bool operator==(const ThisType &other) const;
  bool operator!=(const ThisType &other) const;

  Atom_ *operator*() const;

  // pre-increment
  ThisType &operator++();
  ThisType operator++(int);

  // pre-decrement
  ThisType &operator--();
  ThisType operator--(int);

 private:
  int _end{-1};
  int _pos{-1};
  Mol_ *_mol;
  bool (*_qF)(Atom_ *);

  int _findNext(int from);
  int _findPrev(int from);
};

}  // namespace RDKit

#endif
