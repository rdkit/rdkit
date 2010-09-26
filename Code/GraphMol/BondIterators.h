//
//  Copyright (C) 2002-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file BondIterators.h

  \brief various tools for iterating over a molecule's Bonds

   <b>WARNING:</b> If you go changing the molecule underneath one of
   these iterators you will be sad...
*/  
#ifndef _RD_BOND_ITERATORS_H
#define _RD_BOND_ITERATORS_H

#include "ROMol.h"

namespace RDKit{

  //! \brief iterator for a molecule's bonds, currently BiDirectional,
  //! but it theoretically ought to be RandomAccess.
  class BondIterator_ {
    //FIX: I'm not pleased with the lack of internal testing code
    //  (PREs and the like) in here
  public:
    BondIterator_() : _mol(NULL) {}; 
    BondIterator_(ROMol *mol);
    BondIterator_(ROMol *mol,ROMol::EDGE_ITER pos);
    BondIterator_(const BondIterator_ &other);
    BondIterator_ &operator=(const BondIterator_ &other);
    bool operator==(const BondIterator_ &other);
    bool operator!=(const BondIterator_ &other);
    Bond *operator*();
    // pre-increment
    BondIterator_ &operator++();
    BondIterator_ operator++(int);
    // pre-decrement
    BondIterator_ &operator--();
    BondIterator_ operator--(int);
  private:
    ROMol::EDGE_ITER _beg,_end,_pos;
    ROMol *_mol;
  };
  //! \brief const iterator for a molecule's bonds, currently BiDirectional,
  //! but it theoretically ought to be RandomAccess.
  class ConstBondIterator_ {
  public:
    ConstBondIterator_() : _mol(NULL) {}; 
    ConstBondIterator_(ROMol const *mol);
    ConstBondIterator_(ROMol const *mol,ROMol::EDGE_ITER pos);
    ConstBondIterator_(const ConstBondIterator_ &other);
    ConstBondIterator_ &operator=(const ConstBondIterator_ &other);
    bool operator==(const ConstBondIterator_ &other);
    bool operator!=(const ConstBondIterator_ &other);
    Bond const *operator*();
    // pre-increment
    ConstBondIterator_ &operator++();
    ConstBondIterator_ operator++(int);
    // pre-decrement
    ConstBondIterator_ &operator--();
    ConstBondIterator_ operator--(int);
  private:
    ROMol::EDGE_ITER _beg,_end,_pos;
    ROMol const *_mol;
  };

}  


#endif
