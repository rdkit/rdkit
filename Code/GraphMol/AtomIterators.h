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
#ifndef __RD_ATOM_ITERATORS_H__
#define __RD_ATOM_ITERATORS_H__

#ifdef _MSC_VER
#pragma warning (disable: 4661)  // no suitable definition provided for explicit template instantiation request
#endif

namespace RDKit{
  class QueryAtom;
  
  //! A general random access iterator
  template <class Atom_, class Mol_>
  class AtomIterator_ {
  public:
    typedef AtomIterator_<Atom_,Mol_> ThisType;
    AtomIterator_() : _pos(0),_max(-1),_mol(0) {}; 
    AtomIterator_(Mol_ * mol);
    AtomIterator_(Mol_ * mol,int pos);
    AtomIterator_(const ThisType &other);
    AtomIterator_ &operator=(const ThisType &other);
    AtomIterator_ &operator+=(int val);
    AtomIterator_ &operator-=(int val);
    AtomIterator_ operator+(int val);
    AtomIterator_ operator-(int val);

    // iterator subtraction
    int operator-(ThisType &other);

    // dereference 
    Atom_ * operator*();
    // random access
    Atom_ * operator[](const int which);
    bool operator==(const ThisType &other);
    bool operator!=(const ThisType &other);
    bool operator<(const ThisType &other);
    bool operator<=(const ThisType &other);
    bool operator>(const ThisType &other);
    bool operator>=(const ThisType &other);

    // pre-increment
    ThisType &operator++();
    ThisType operator++(int);

    // pre-decrement
    ThisType &operator--();
    ThisType operator--(int);
  
    //private:
    int _pos,_max;
    Mol_ * _mol;
  };


  //! Iterate over heteroatoms, this is bidirectional
  template <class Atom_, class Mol_>
  class HeteroatomIterator_ {
  public:
    typedef HeteroatomIterator_<Atom_,Mol_> ThisType;
    HeteroatomIterator_() : _mol(0) {}; 
    HeteroatomIterator_(Mol_ * mol);
    HeteroatomIterator_(Mol_ * mol,int pos);
    ~HeteroatomIterator_();
    HeteroatomIterator_(const ThisType &other);
    HeteroatomIterator_ &operator=(const ThisType &other);
    bool operator==(const ThisType &other);
    bool operator!=(const ThisType &other);

    Atom_ * operator*();

    // pre-increment
    ThisType &operator++();
    ThisType operator++(int);

    // pre-decrement
    ThisType &operator--();
    ThisType operator--(int);
  private:
    int _end,_pos;
    Mol_ * _mol;
    //FIX: somehow changing the following to a pointer make the regression test pass
    // QueryAtom _qA;
    QueryAtom *_qA;

    int _findNext(int from);  
    int _findPrev(int from);
  };

  //! Iterate over aromatic atoms, this is bidirectional
  template <class Atom_, class Mol_>
  class AromaticAtomIterator_ {
  public:
    typedef AromaticAtomIterator_<Atom_,Mol_> ThisType;
    AromaticAtomIterator_() : _mol(0) {}; 
    AromaticAtomIterator_(Mol_ * mol);
    AromaticAtomIterator_(Mol_ * mol,int pos);
    ~AromaticAtomIterator_();
    AromaticAtomIterator_(const ThisType &other);
    AromaticAtomIterator_ &operator=(const ThisType &other);
    bool operator==(const ThisType &other);
    bool operator!=(const ThisType &other);

    Atom_ * operator*();

    // pre-increment
    ThisType &operator++();
    ThisType operator++(int);

    // pre-decrement
    ThisType &operator--();
    ThisType operator--(int);
  private:
    int _end,_pos;
    Mol_ * _mol;

    int _findNext(int from);  
    int _findPrev(int from);
  };
  
  //! Iterate over atoms matching a query. This is bidirectional.
  template <class Atom_, class Mol_>
  class QueryAtomIterator_ {
  public:
    typedef QueryAtomIterator_<Atom_,Mol_> ThisType;
    QueryAtomIterator_() : _mol(0),_qA(0) {}; 
    QueryAtomIterator_(Mol_ * mol,QueryAtom const *what);
    QueryAtomIterator_(Mol_ * mol,int pos);
    ~QueryAtomIterator_();
    QueryAtomIterator_(const ThisType &other);
    QueryAtomIterator_ &operator=(const ThisType &other);
    bool operator==(const ThisType &other);
    bool operator!=(const ThisType &other);

    Atom_ * operator*();

    // pre-increment
    ThisType &operator++();
    ThisType operator++(int);

    // pre-decrement
    ThisType &operator--();
    ThisType operator--(int);
  private:
    int _end,_pos;
    Mol_ * _mol;
    QueryAtom *_qA;

    int _findNext(int from);  
    int _findPrev(int from);
  };

  //! Iterate over atoms matching a query function. This is bidirectional.
  template <class Atom_, class Mol_>
  class MatchingAtomIterator_ {
  public:
    typedef MatchingAtomIterator_<Atom_,Mol_> ThisType;
    MatchingAtomIterator_() : _end(-1),_pos(-1),_mol(0),_qF(0) {}; 
    MatchingAtomIterator_(Mol_ * mol,bool (*fn)(Atom_ *));
    MatchingAtomIterator_(Mol_ * mol,int pos);
    ~MatchingAtomIterator_();
    MatchingAtomIterator_(const ThisType &other);
    MatchingAtomIterator_ &operator=(const ThisType &other);
    bool operator==(const ThisType &other);
    bool operator!=(const ThisType &other);

    Atom_ * operator*();

    // pre-increment
    ThisType &operator++();
    ThisType operator++(int);

    // pre-decrement
    ThisType &operator--();
    ThisType operator--(int);
  private:
    int _end,_pos;
    Mol_ * _mol;
    bool (*_qF)(Atom_ *);

    int _findNext(int from);  
    int _findPrev(int from);
  };

  

}  /* end o namespace */

#endif
